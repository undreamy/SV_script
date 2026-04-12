library(data.table) 
library(stringr) 
# 1 文件路径 
anno_file  <- "<您的文件目录>/maize_long_table.txt" 
id_file    <- "<您的参考数据库的目的基因如AGC的id文件目录>/B73_id.txt" 
output_file<- "<您的输出目录>/target_SV.txt" 
# 2 读取数据 
dt_anno <- fread(anno_file) 
dt_ids  <- fread(id_file, header = FALSE) 
target_ids <- dt_ids$V1 
# 3 提取 GeneID 
dt_anno[, extracted_id := str_extract(GeneID, "Zm\\d+eb\\d+")] 
# 删除NA 
dt_anno <- dt_anno[!is.na(extracted_id)] 
# 4 AGC基因筛选 
dt_result <- dt_anno[extracted_id %in% target_ids] 
# 5 去掉 intergenic 
dt_result <- dt_result[tolower(Effect) != "intergenic_region"] 
# 6 删除临时列 
dt_result[, extracted_id := NULL] 
# 7 去重 
dt_result <- unique(dt_result) 
# 8 POS 转数值 
dt_result[, POS := as.numeric(POS)] 
# 9 输出 
fwrite(dt_result, output_file, sep="\t", quote=FALSE) 
cat("筛选完成！结果文件：", output_file, "\n")
