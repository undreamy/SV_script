# ==============================
# 0. 加载包
# ==============================

library(data.table)
library(stringr)

# ==============================
# 1. 读取 bcftools 输出
# ==============================

gt_raw <- fread(
  "sample_gt_raw.txt",
  header = FALSE,
  sep = "\t",
  col.names = c("CHROM","POS","END","ID","REF","ALT")
)

# 检查样本列
n_samples <- ncol(gt_raw) - 6
if (n_samples <= 0) stop("未检测到样本列，请检查 bcftools 输出。")

cat("检测到样本数量:", n_samples, "\n")

# ==============================
# 2. 提取样本名
# ==============================

sample_names <- sub(
  "=.*",
  "",
  unlist(gt_raw[1, 7:ncol(gt_raw), with = FALSE])
)

# 重命名样本列
setnames(gt_raw, 7:ncol(gt_raw), sample_names)

# ==============================
# 3. 提取 GT 值
# ==============================

for (s in sample_names) {
  gt_raw[, (s) := sub(".*=", "", get(s))]
}

# ==============================
# 4. 定义 carrier 判断函数
# ==============================

is_carrier <- function(gt) {
  !is.na(gt) & !gt %in% c("./.", ".|.", "0/0", "0|0")
}

# ==============================
# 5. 生成 carrier 样本列表
# ==============================

gt_raw[, Samples := {
  carriers <- sample_names[sapply(.SD, is_carrier)]
  paste(carriers, collapse = ",")
}, .SDcols = sample_names]

# ==============================
# 6. 计算 carrier 数量
# ==============================

gt_raw[, Carrier_count :=
        ifelse(Samples == "", 0, lengths(strsplit(Samples, ",")))]

# ==============================
# 7. 计算 SV 频率
# ==============================

gt_raw[, Frequency := Carrier_count / n_samples]

# ==============================
# 8. 保存 SV 样本信息
# ==============================

sample_info <- gt_raw[, .(
  CHROM,
  POS,
  END,
  ID,
  REF,
  ALT,
  Samples,
  Carrier_count,
  Frequency
)]

fwrite(sample_info,
       "sample_info_per_sv.txt",
       sep = "\t",
       quote = FALSE)

cat("SV样本信息已保存: sample_info_per_sv.txt\n")

# ==============================
# 9. 读取 AGC SV 表
# ==============================

sv_final <- fread("AGC_SV_v3.txt")

# 如果 AGC 表中叫 START
if ("START" %in% names(sv_final)) {
  setnames(sv_final, "START", "POS")
}

# ==============================
# 10. 合并样本信息
# ==============================

if (all(!is.na(sv_final$ID)) & !any(duplicated(sv_final$ID))) {

  merged <- merge(
    sv_final,
    sample_info[, .(ID, Samples, Carrier_count, Frequency)],
    by = "ID",
    all.x = TRUE
  )

} else {

  merged <- merge(
    sv_final,
    sample_info[, .(CHROM, POS, END, Samples, Carrier_count, Frequency)],
    by = c("CHROM","POS","END"),
    all.x = TRUE
  )
}

# ==============================
# 11. 保存最终结果
# ==============================

fwrite(
  merged,
  "AGC_SV_with_samples.txt",
  sep = "\t",
  quote = FALSE
)

cat("最终结果已保存: AGC_SV_with_samples.txt\n")

# ==============================
# 12. 简单统计
# ==============================

cat("\nSV carrier 数量分布:\n")
print(table(merged$Carrier_count))

cat("\n平均SV频率:", mean(merged$Frequency, na.rm = TRUE), "\n")
