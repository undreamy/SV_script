# SV_script
基于玉米26个NAM品系的SV数据，经snpEff注释，python、R等处理。筛选得到AGC家族的SV数据，并给出overlap验证。接着进行表达关联分析。  

**依赖环境和包**  
Java Runtime Environment 8+，Python 3.7+（"sys""argparse""pandas""numpy"等），R 4.0+（"data.table""stringr""optparse"等），snpEff  

**输入数据**  
1.参考基因组的fa和gff3文件  
2.待处理的SV数据（VCF），需包含标准VCF字段（如 CHROM, POS, END, ID, REF, ALT）。 
3.目的基因列表和样本（品系）名称列表（txt）  
4.表达矩阵（TPM）  

**使用顺序与简介**  
1.snpEff的配置、建库和注释，调用snpSift提取注释信息（snpeff.txt）。教程详见（https://zhuanlan.zhihu.com/p/625865035）。  
2.conver_v3.py转长表  
3.filter.R筛选目的基因  
4.location_analysis.R验证是否有重叠  
5.merge_sample.R进行样本名关联  
6.sv_express.py进行表达关联分析

**结果示例**  
1.筛选得到的目的基因SV数据应如下图所示：  
<img width="1008" height="207" alt="局部截取_20260509_121009" src="https://github.com/user-attachments/assets/31a64921-fb7f-4b91-85e8-f6c0d0543d3b" />  

2.与样本关联后应如下图所示：  
<img width="1140" height="243" alt="局部截取_20260509_121303" src="https://github.com/user-attachments/assets/b0b65f37-a94a-43fb-82f4-0581d2b67fb8" />  

3.表达量关联分析结果如下图所示，仅展示全组织合并分析的csv文件截图，还应有显著基因的箱线图：  
<img width="800" height="223" alt="局部截取_20260509_121656" src="https://github.com/user-attachments/assets/bd88ac69-b646-4e30-82f4-22906c4af61d" />

