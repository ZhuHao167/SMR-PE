m_e <- data.table::fread('../result/smr_result/smr_mqtl_eqtl.txt')
prob <- unique(m_e$Outco_Gene)
write.table(prob,'../data/m_e_eQTL_prob.txt',row.names = F,col.names = F,quote = F)


# 画图文件输出
```{bash}
mkdir -p ../data/SMR_plot

export PATH=$PATH:/1.Tools/parallel-20240422/src

sort -u ../data/m_e_eQTL_prob.txt | parallel -j 7 "/software/SMR/smr-1.3.1-linux-x86_64/smr \
--bfile /mnt/bak/smr_data/g1000_eur/g1000_eur \
--gwas-summary ../data/finngen_R10_O15_PREECLAMPS.txt \
--beqtl-summary /mnt/bak/smr_data/eqtl/cis-eQTLs_eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense \
--out ../data/SMR_plot/eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS \
--plot \
--probe {} \
--probe-wind 500 \
--gene-list /public/SMR/glist-hg19.txt \
--thread-num 40"
```
cat ../data/m_e_eQTL_prob.txt 

prob <- read.table('../data/m_e_eQTL_prob.txt')
prob <- prob$V1
length(prob)
source("/2.Common_data/plot_SMR_refine.r") 
i <- prob[4]
for (i in i) {
  SMRData = ReadSMRData(paste0("../data/SMR_plot/plot/eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.",i,".txt"))
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  name <- bitr(SMRData$SMR$V4,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
  colnames(name) = c("Gene","Symbol")
  name = name[!duplicated(name$Gene),]
  smr=SMRData$SMR
  smr$ind = as.integer(rownames(smr))
  smr = merge(smr,name,by.x='V4',by.y='Gene',all.x = TRUE)
  smr = smr[order(smr$ind),]
  # smr$Symbol[which(is.na(smr$Symbol))] <- smr$V1
  SMRData$SMR$V4 <- smr$Symbol
  # 替换 SMRData中的 p 值，以筛选标准为主，这里以p_SMR_multi为例
  GWAS_QTL <- data.table::fread('../result/smr_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.msmr')
  SMRData$SMR$V8 <- GWAS_QTL$p_SMR_multi[match(SMRData$SMR$V1,GWAS_QTL$probeID)] #GWAS_QTL为QTL对应的SMR结果
  
  eQTL_mQTL_pass_multi <- data.table::fread('../result/smr_result/smr_mqtl_eqtl.txt') #该数据需经过筛选条件，筛选之后的列表
  probeid.list <- eQTL_mQTL_pass_multi$Outco_Gene
  
  jpeg(paste0('../figure/SMRLocusPlot_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS_',i,'.jpeg'),
       res = 300,units = "cm",
       width = 30, height = 20)
  SMRLocusPlot(data=SMRData, ## 原始参数，输入数据
               smr_thresh=0.05, ## 原始参数，p值阈值标准
               heidi_thresh=0.01, ## 原始参数，HEIDI检测阈值标准
               plotWindow=500,  ## 原始参数，展示窗口区域大小
               max_anno_probe=15, ## 原始参数，最多展示位点数
               smr_thresh_type = "p_SMR_multi", ## 更新参数，阈值类型标注
               smrindx = probeid.list## 更新参数，目标高亮位点
               ) 
  dev.off()
  
  jpeg(paste0('../figure/SMREffectPlot_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS_',i,'.jpeg'),
       res = 300,units = "cm",
       width = 15, height = 15)
  SMREffectPlot(data=SMRData, 
                # trait_name=paste0("Effect Plot for ","Your_SYMBOL"), ## 主标题，可修改，不加标题也可以
                title1 = "eQTLs", ## 坐标轴标题，按数据类型修改
                title2 = "GWAS", ## 坐标轴标题，一般不修改
                Pval_title = "p_SMR_multi") ## 更新参数，判定标准类型
  dev.off()
}
dev.off()


tmp <- data.table::fread('../result/smr_result/smr_eQTLs_GTEx_Uterus_finngen_R10_O15_PREECLAMPS.txt')
prob <- unique(tmp$probeID)
write.table(prob,'../data/tmp_eQTL_prob.txt',row.names = F,col.names = F,quote = F)


# 画图文件输出
```{bash}
mkdir -p ../data/SMR_plot

export PATH=$PATH:/1.Tools/parallel-20240422/src

sort -u ../data/tmp_eQTL_prob.txt | parallel -j 7 "/software/SMR/smr-1.3.1-linux-x86_64/smr \
--bfile /mnt/bak/smr_data/g1000_eur/g1000_eur \
--gwas-summary ../data/finngen_R10_O15_PREECLAMPS.txt \
--beqtl-summary /mnt/bak/smr_data/eqtl/cis-eQTLs_eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense \
--out ../data/SMR_plot/eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS \
--plot \
--probe {} \
--probe-wind 500 \
--gene-list /public/SMR/glist-hg19.txt \
--thread-num 40"
```
cat ../data/tmp_eQTL_prob.txt 

prob <- read.table('../data/tmp_eQTL_prob.txt')
prob <- prob$V1
length(prob)
source("/2.Common_data/plot_SMR_refine.r") 
# i <- prob[4]
for (i in prob) {
  SMRData = ReadSMRData(paste0("../data/SMR_plot/plot/eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.",i,".txt"))
  
  library(clusterProfiler)
  load('/public/biomaRt_ensembl_annotation/biomaRt_ensembl_annotation.Rdata')
  name <- bitr(SMRData$SMR$V4, fromType = "ensembl_gene_id",
               toType = c("hgnc_symbol"),
               OrgDb = human) #小鼠换成mouse
  colnames(name) = c("Gene","Symbol")
  name = name[!duplicated(name$Gene),]
  smr=SMRData$SMR
  smr$ind = as.integer(rownames(smr))
  smr = merge(smr,name,by.x='V4',by.y='Gene',all.x = TRUE)
  smr = smr[order(smr$ind),]
  # smr$Symbol[which(is.na(smr$Symbol))] <- smr$V1
  SMRData$SMR$V4 <- smr$Symbol
  # 替换 SMRData中的 p 值，以筛选标准为主，这里以p_SMR_multi为例
  GWAS_QTL <- data.table::fread('../result/smr_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.msmr')
  SMRData$SMR$V8 <- GWAS_QTL$p_SMR_multi[match(SMRData$SMR$V1,GWAS_QTL$probeID)] #GWAS_QTL为QTL对应的SMR结果
  
  #eQTL_mQTL_pass_multi <- data.table::fread('../result/smr_result/smr_mqtl_eqtl.txt') #该数据需经过筛选条件，筛选之后的列表
  #probeid.list <- eQTL_mQTL_pass_multi$Outco_Gene
  probeid.list <- prob
  
  jpeg(paste0('../figure/SMRLocusPlot_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS_',i,'.jpeg'),
       res = 300,units = "cm",
       width = 30, height = 20)
  SMRLocusPlot(data=SMRData, ## 原始参数，输入数据
               smr_thresh=0.05, ## 原始参数，p值阈值标准
               heidi_thresh=0.01, ## 原始参数，HEIDI检测阈值标准
               plotWindow=500,  ## 原始参数，展示窗口区域大小
               max_anno_probe=15, ## 原始参数，最多展示位点数
               smr_thresh_type = "p_SMR_multi", ## 更新参数，阈值类型标注
               smrindx = probeid.list## 更新参数，目标高亮位点
  ) 
  dev.off()
  
  
  tryCatch({
    jpeg(paste0('../figure/SMREffectPlot_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS_',i,'.jpeg'),
         res = 300,units = "cm", width = 15, height = 15)
    SMREffectPlot(data=SMRData, 
                  # trait_name=paste0("Effect Plot for ","Your_SYMBOL"), ## 主标题，可修改，不加标题也可以
                  title1 = "eQTLs", ## 坐标轴标题，按数据类型修改
                  title2 = "GWAS", ## 坐标轴标题，一般不修改
                  Pval_title = "p_SMR_multi") ## 更新参数，判定标准类型
    dev.off()
  }, error = function(e) {
    dev.off()
    cat("Error in SMREffectPlot:", conditionMessage(e), "\n")
  })
  
}

