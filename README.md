# SV_script
基于玉米26个NAM品系的SV数据，经snpEff注释，python、R等处理。筛选得到AGC家族的SV数据，并给出overlap验证。接着进行表达关联分析。  

**输入数据**  
参考基因组的fa和gff3文件、待处理的SV数据（VCF）、目的基因列表（txt）、表达矩阵（TPM）  

**使用简介**  
1.snpEff的配置、建库和注释。教程详见（https://zhuanlan.zhihu.com/p/625865035）。  
2.调用snpSift提取注释信息（snpeff.txt）  
3.conver_v3.py转长表  
4.filter.R筛选目的基因  
5.location_analysis.R验证是否有重叠  
6.merge_sanmple.R进行样本名关联  
7.sv_express.py进行表达关联分析
