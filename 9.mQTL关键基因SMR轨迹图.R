m_e <- data.table::fread('../result/smr_result/smr_mqtl_eqtl.txt')
prob <- unique(m_e$Expo_ID)
write.table(prob,'../data/m_e_mQTL_prob.txt',row.names = F,col.names = F,quote = F)

# 画图文件输出
```{bash}
mkdir -p ../data/SMR_plot

export PATH=$PATH:/1.Tools/parallel-20240422/src

sort -u ../data/m_e_mQTL_prob.txt | parallel -j 7 "/software/SMR/smr-1.3.1-linux-x86_64/smr \
--bfile /mnt/bak/smr_data/g1000_eur/g1000_eur \
--gwas-summary ../data/finngen_R10_O15_PREECLAMPS.txt \
--beqtl-summary /public/SMR/mqtl/LBC_BSGS_meta/LBC_BSGS_meta_merge \
--out ../data/SMR_plot/cis-mQTLs_LBC_BSGS_meta_finngen_R10_O15_PREECLAMPS \
--plot \
--probe {} \
--probe-wind 500 \
--gene-list /public/SMR/glist-hg19.txt \
--thread-num 40"
```

prob <- read.table('../data/m_e_mQTL_prob.txt')
prob <- prob$V1
length(prob)

source("/2.Common_data/plot_SMR_refine.r") 
i <- prob[6]
# Read the data file in R:
for (i in i) {
  SMRData = ReadSMRData(paste0("../data/SMR_plot/plot/cis-mQTLs_LBC_BSGS_meta_finngen_R10_O15_PREECLAMPS.",i,".txt"))
  
  probe = data.table::fread('/mnt/bak/smr_data/mqtl/LBC_BSGS_meta/meta/GPL16304_LBC_BSGS.txt')
  head(probe$Closest_TSS_gene_name[match(SMRData$SMR$V1,probe$ID)])
  SMRData$SMR$V4 <- probe$Closest_TSS_gene_name[match(SMRData$SMR$V1,probe$ID)]
  # 替换 SMRData中的 p 值，以筛选标准为主，这里以p_SMR_multi为例
  GWAS_QTL <- data.table::fread('../result/smr_mQTLs_LBC_BSGS_meta_finngen_R10_O15_PREECLAMPS.msmr')
  SMRData$SMR$V8 <- GWAS_QTL$p_SMR_multi[match(SMRData$SMR$V1,GWAS_QTL$probeID)] #GWAS_QTL为QTL对应的SMR结果
  
  # 以m-e的分析结果为例，提取要显示基因对应的甲基化位点probeid.list
  eQTL_mQTL_pass_multi <- data.table::fread('../result/smr_result/smr_mqtl_eqtl.txt') #该数据需经过筛选条件，筛选之后的列表
  probeid.list <- eQTL_mQTL_pass_multi$Expo_ID
  
  
  jpeg(paste0('../figure/SMRLocusPlot_mQTLs_LBC_BSGS_meta_finngen_R10_O15_PREECLAMPS_',i,'.jpeg'), res = 300,units = "cm",
       width = 30, height = 20)
  SMRLocusPlot(data=SMRData, ## 原始参数，输入数据
               smr_thresh=0.05, ## 原始参数，p值阈值标准
               heidi_thresh=0.01, ## 原始参数，HEIDI检测阈值标准
               plotWindow=500,  ## 原始参数，展示窗口区域大小
               max_anno_probe=15, ## 原始参数，最多展示位点数
               smr_thresh_type = "p_SMR_multi", ## 更新参数，阈值类型标注
               smrindx = probeid.list) ## 更新参数，目标高亮位点
  dev.off()
  
  jpeg(paste0('../figure/SMREffectPlot_mQTLs_LBC_BSGS_meta_finngen_R10_O15_PREECLAMPS_',i,'.jpeg'), res = 300,units = "cm",
       width = 15, height = 15)
  SMREffectPlot(data=SMRData, 
                # trait_name=paste0("Effect Plot for ","Your_SYMBOL"), ## 主标题，可修改，不加标题也可以
                title1 = "mQTLs", ## 坐标轴标题，按数据类型修改
                title2 = "GWAS", ## 坐标轴标题，一般不修改
                Pval_title = "p_SMR_multi") ## 更新参数，判定标准类型
  dev.off()
  
}
dev.off()

