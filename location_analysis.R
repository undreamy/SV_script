library(data.table) 
library(stringr) 
# 读取数据 
sv_file <- "<您的 R筛选结果路径>/target_SV.txt" 
gff_file <- "<您的文件路径>snpEff/data/Zm-B73-NAM-5.0/genes.gff" 
sv <- fread(sv_file) 
gff <- fread(gff_file, skip = "#") 
# 提取基因信息 
gene_gff <- gff[V3 == "gene"] 
gene_gff[, GeneID := str_extract(V9, "Zm\\d+eb\\d+")] 
gene_gff <- gene_gff[!is.na(GeneID)] 
# 创建基因坐标表 
gene_coords <- gene_gff[, .( 
GENE_CHROM = V1, 
gene_start = V4, 
gene_end = V5, 
GeneID 
)] 
# 合并基因坐标 
sv2 <- merge(sv, gene_coords, by = "GeneID", all.x = TRUE) 
# 检查合并后的列名 
cat("合并后的列名：\n") 
print(names(sv2)) 
# 重命名重复的染色体列 
if ("CHROM.x" %in% names(sv2) && "CHROM.y" %in% names(sv2)) { 
setnames(sv2, "CHROM.x", "SV_CHROM") 
setnames(sv2, "CHROM.y", "GENE_CHROM") 
} else if ("CHROM" %in% names(sv2) && "CHROM.1" %in% names(sv2)) { 
setnames(sv2, "CHROM", "SV_CHROM") 
setnames(sv2, "CHROM.1", "GENE_CHROM") 
} else if ("CHROM" %in% names(sv2) && "GENE_CHROM" %in% names(sv2)) { 
setnames(sv2, "CHROM", "SV_CHROM") 
} 
# 进行overlap判断 
sv2[, overlap := 
(START <= gene_end) & 
(END >= gene_start) & 
(SV_CHROM == GENE_CHROM)] 
# 统计总体情况 
total_count <- nrow(sv2) 
overlap_count <- sum(sv2$overlap, na.rm = TRUE) 
non_overlap_count <- sum(!sv2$overlap, na.rm = TRUE) 
cat("\n===== 总体统计 =====\n") 
cat("总 SV记录数:", total_count, "\n") 
cat("与基因重叠的SV数:", overlap_count, "\n") 
cat("不与基因重叠的SV数:", non_overlap_count, "\n") 
cat("重叠比例:", round(overlap_count/total_count*100, 2), "%\n") 
cat("不重叠比例:", round(non_overlap_count/total_count*100, 2), "%\n") 
# 筛选出不重叠的SV记录（可能存在错误的记录） 
non_overlap_sv <- sv2[overlap == FALSE, ] 
cat("\n===== 不重叠的 SV记录详情 =====\n") 
if (nrow(non_overlap_sv) > 0) { 
cat("找到", nrow(non_overlap_sv), "条不与基因重叠的 SV记录\n\n") 
# 输出详细信息 
cat("不重叠 SV记录详情（前20条）：\n") 
print(non_overlap_sv[1:min(20, nrow(non_overlap_sv)),  
.(GeneID, START, END, SV_CHROM, GENE_CHROM,  
gene_start, gene_end, Effect, Impact, overlap)]) 
# 按染色体统计 
cat("\n 按染色体统计不重叠记录：\n") 
chrom_stats <- non_overlap_sv[, .(count = .N), by = SV_CHROM] 
print(chrom_stats[order(-count)]) 
    # 按影响类型统计 
    if ("Effect" %in% names(non_overlap_sv)) { 
        cat("\n按Effect类型统计不重叠记录：\n") 
        effect_stats <- non_overlap_sv[, .(count = .N), by = Effect] 
        print(effect_stats[order(-count)]) 
    } 
     
    if ("Impact" %in% names(non_overlap_sv)) { 
        cat("\n按Impact类型统计不重叠记录：\n") 
        impact_stats <- non_overlap_sv[, .(count = .N), by = Impact] 
        print(impact_stats[order(-count)]) 
    } 
     
    # 分析不重叠的可能原因 
    cat("\n===== 不重叠原因分析 =====\n") 
     
    # 1. 染色体不匹配 
    chrom_mismatch <- non_overlap_sv[SV_CHROM != GENE_CHROM, .N] 
    cat("1. 染色体不匹配的记录数:", chrom_mismatch,  
        sprintf("(%.1f%%)", chrom_mismatch/nrow(non_overlap_sv)*100), "\n") 
     
    if (chrom_mismatch > 0) { 
        cat("   示例：\n") 
        examples <- non_overlap_sv[SV_CHROM != GENE_CHROM,  
                                   .(GeneID, SV_CHROM, GENE_CHROM, START, END, 
gene_start, gene_end)][1:5] 
        print(examples) 
    } 
     
    # 2. 位置完全在基因之前 
    before_gene <- non_overlap_sv[END < gene_start & SV_CHROM == 
GENE_CHROM, .N] 
    cat("\n2. 完全在基因之前的记录数:", before_gene,  
        sprintf("(%.1f%%)", before_gene/nrow(non_overlap_sv)*100), "\n") 
     
    if (before_gene > 0) { 
        cat("   示例：\n") 
        examples <- non_overlap_sv[END < gene_start & SV_CHROM == 
GENE_CHROM, 
                                   .(GeneID, SV_CHROM, START, END, gene_start, gene_end,  
                                     distance_to_gene = gene_start - END)][1:5] 
        print(examples) 
    } 
     
    # 3. 位置完全在基因之后 
    after_gene <- non_overlap_sv[START > gene_end & SV_CHROM == 
GENE_CHROM, .N] 
    cat("\n3. 完全在基因之后的记录数:", after_gene,  
        sprintf("(%.1f%%)", after_gene/nrow(non_overlap_sv)*100), "\n") 
     
    if (after_gene > 0) { 
        cat("   示例：\n") 
        examples <- non_overlap_sv[START > gene_end & SV_CHROM == 
GENE_CHROM, 
                                   .(GeneID, SV_CHROM, START, END, gene_start, gene_end, 
                                     distance_to_gene = START - gene_end)][1:5] 
        print(examples) 
    } 
     
    # 4. 染色体匹配但位置不重叠 
    same_chrom_no_overlap <- non_overlap_sv[SV_CHROM == GENE_CHROM, .N] 
    cat("\n4. 染色体匹配但位置不重叠的记录数:", same_chrom_no_overlap,  
        sprintf("(%.1f%%)", same_chrom_no_overlap/nrow(non_overlap_sv)*100), "\n") 
     
    # 保存不重叠的记录到文件 
    output_file_non_overlap <- "<您的文件目录>/ non_overlap.txt" 
    fwrite(non_overlap_sv, output_file_non_overlap, sep = "\t") 
    cat("\n不重叠的SV记录已保存到：", output_file_non_overlap, "\n") 
     
    # 保存详细分析报告 
    output_file_analysis <- "<您的文件目录>/ non_overlap_analysis.txt" 
    sink(output_file_analysis) 
    cat("=== 不与基因重叠的SV记录分析报告 ===\n\n") 
    cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n") 
    cat("总记录数:", nrow(non_overlap_sv), "\n\n") 
     
    cat("按染色体统计：\n") 
    print(chrom_stats[order(-count)]) 
    cat("\n") 
     
    cat("按不重叠原因统计：\n") 
    cat("1. 染色体不匹配:", chrom_mismatch, "\n") 
    cat("2. 完全在基因之前:", before_gene, "\n") 
    cat("3. 完全在基因之后:", after_gene, "\n") 
    cat("4. 染色体匹配但位置不重叠:", same_chrom_no_overlap, "\n") 
    sink() 
    cat("详细分析报告已保存到：", output_file_analysis, "\n") 
     
} else { 
cat("没有找到不与基因重叠的SV记录，所有SV都与基因坐标匹配。\n") 
} 
# 保存所有记录（包含overlap标记） 
output_file_all <- "<您的文件目录>/sv_ with_overlap_marked.txt" 
fwrite(sv2, output_file_all, sep = "\t") 
cat("\n 所有记录（含overlap标记）已保存到：", output_file_all, "\n") 
cat("\n===== 分析完成 =====\n")
