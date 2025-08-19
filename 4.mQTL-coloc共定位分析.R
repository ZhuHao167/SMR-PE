# mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS----
library(data.table)
library(tidyverse)

dat1 <- fread('../result/smr_result/smr_mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS.txt')
finngen_R10_O15_PREECLAMPS <- fread('../rawdata/finngen_R10_O15_PREECLAMPS.gz')
head(finngen_R10_O15_PREECLAMPS)
finngen_R10_O15_PREECLAMPS$af_alt <- as.numeric(finngen_R10_O15_PREECLAMPS$af_alt)
table(finngen_R10_O15_PREECLAMPS$af_alt > 0.5)
finngen_R10_O15_PREECLAMPS$MAF = ifelse(finngen_R10_O15_PREECLAMPS$af_alt<0.5, finngen_R10_O15_PREECLAMPS$af_alt, 1 - finngen_R10_O15_PREECLAMPS$af_alt)
hist(finngen_R10_O15_PREECLAMPS$MAF)
table(finngen_R10_O15_PREECLAMPS$MAF == 0)
finngen_R10_O15_PREECLAMPS <- finngen_R10_O15_PREECLAMPS[finngen_R10_O15_PREECLAMPS$MAF > 0,]
finngen_R10_O15_PREECLAMPS <- dplyr::rename(finngen_R10_O15_PREECLAMPS,
                                            chromosome = `#chrom`,
                                            base_pair_location = pos,
                                            p_value = pval,
                                            variant_id = rsids ,
                                            MAF = MAF)

gwas1 = finngen_R10_O15_PREECLAMPS[order(finngen_R10_O15_PREECLAMPS$p_value),]
gwas1 = gwas1[!duplicated(gwas1$variant_id),]
head(gwas1)

gwas2 <- fread('../data/QTL_Coloc_Prepare/mQTLs_coloc.txt')
## 添加一下列名
colnames(gwas2) <- c("SNP","Chr","BP","A1","A2","Freq","Probe","Probe_Chr","Probe_bp","Gene","Orientation","b","SE","p")
head(gwas2)
gwas2 = gwas2[order(gwas2$p),]
# gwas2 = gwas2[!duplicated(gwas2$SNP),]
gwas2$MAF1 = ifelse(gwas2$Freq<0.5, gwas2$Freq, 1 - gwas2$Freq) 


output_table <- data.frame()

for (i in unique(dat1$probeID)) { # 循环对每个筛选到的结果进行共定位，eQTL对每个基因基因共定位分析，mQTL对每个甲基化位点进行共定位分析，pQTL对每个蛋白进行共定位分析。
  SNP <- dat1$topSNP[dat1$probeID == i]
  CHR <- dat1$ProbeChr[dat1$probeID == i]
  BP <- dat1$Probe_bp[dat1$probeID == i]
  gwas1_c <- gwas1[gwas1$chromosome == CHR,]
  head(gwas1_c)

  gwas2_c <- gwas2[gwas2$Chr == CHR,]
  # 以QTL数据中的关键probeID筛选上下游1000 kb的SNP
  gwas2_c <- gwas2_c[gwas2_c$BP >  BP - 1000000 &
                       gwas2_c$BP < BP + 1000000,]
  head(gwas2_c)
  gwas2_c <- gwas2_c[which(gwas2_c$Probe == i),] ## 确保共定位所用的QTL数据中，每个SNP的影响对象是目标探针
  dat_merge <- merge(gwas1_c,gwas2_c,
                     by.x = 'variant_id',by.y = 'SNP',
                     suffixes = c("_gwas1","_gwas2"))
  head(dat_merge)
  
  table(is.na(dat_merge$p_value))
  min(dat_merge$p_value)
  dat_merge$p_value[dat_merge$p_value==0] <- min(dat_merge$p_value[dat_merge$p_value!=0])
  
  table(is.na(dat_merge$p))
  min(dat_merge$p)
  dat_merge$p[dat_merge$p==0] <- min(dat_merge$p[dat_merge$p!=0])
  
  # 共定位分析
  library("coloc")
  ngwas = 7877+211957  #gwas样本量，根据GWAS队列作修改
  case = 7877 #得病数，根据GWAS队列作修改  
  neqtl = 1980 #eqtl样本量31684、mQTL样本量1980 、pqtl样本量有3个数据来源（参考以下），具体需根据实际分析使用数据来确定。
  
  #名称  数据来源  作者  队列人种  样本量
  #Decode_2550 PMID:34857953 Ferkingstad et al.  Icelanders  35,559
  #Fenland_4674_protein  PMID:34648354 Pietzner et al. European  10,708
  #UKB-PPP_2937  PMID:37794186 Sun et al.  UKB participants  54,219
  
  
  result <- coloc.abf(
    dataset1 = list(pvalues=dat_merge$p_value, #GWAS pvalue
                    snp=dat_merge$variant_id,
                    type="cc", N=ngwas, s=case/ngwas,
                    MAF = dat_merge$MAF),
    dataset2 = list(pvalues=dat_merge$p, #QTL pvalue
                    snp=dat_merge$variant_id,
                    type="quant", N=neqtl,
                    MAF = dat_merge$MAF1),
    p12 = 5e-05) # p12 default 1e-05，该设定会影响PPH4的值，选用5e-05得到的PPH4结果会优于默认值。
  
  tmp <- data.frame(ProbID = i, 
                    PPH4 = result$summary[6], # 森林图PPH4选用这个值
                    SNP.PP.H4 = max(result$results$SNP.PP.H4)) 
  
  
  output_table <- rbind(tmp, output_table) # 保存各个结果的PPH4值
  if (!file.exists('../result/coloc_result')) {
    dir.create('../result/coloc_result')
  } 
  if(result$summary[6] > 0.5){
    library(dplyr)
    
    need_result=result$results #输出参与共定位分析的全部SNP的共定位结果
    write.table(need_result,
                paste0("../result/coloc_result/coloc-cis-mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS-",SNP,".tsv"),
                col.names = T,row.names = F,sep="\t",quote = F)
    
    
    gwas = cbind(dat_merge$variant_id,dat_merge$p_value)
    colnames(gwas) = c("rsid","pval")
    eqtl = cbind(dat_merge$variant_id,dat_merge$p)
    colnames(eqtl) = c("rsid","pval")
    if (!file.exists('../data/locusdata')) {
      dir.create('../data/locusdata')
    } 
    write.table(gwas,
                paste0("../data/locusdata/gwas-cis-mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS-",i,".tsv"),
                col.names = T,row.names = F,sep="\t",quote = F)
    write.table(eqtl,
                paste0("../data/locusdata/qtl-cis-mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS-",i,".tsv"),
                col.names = T,row.names = F,sep="\t",quote = F)
    
    locuscompare(paste0("../data/locusdata/gwas-cis-mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS-",i,".tsv"),
                 paste0("../data/locusdata/qtl-cis-mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS-",i,".tsv"),
                 title1 = 'GWAS',title2 = 'mQTL',#title1参数表示可视化结果的横坐标标签，默认是eQTL，因输入数据1为GWAS数据，因此title1 = GWAS，eQTL和pQTL的共定位分析，记得修改title2参数，如title2 = “eQTL"
                 legend =T)  
    ggsave(paste0("../figure/locus-cis-mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS-",i,".jpeg"),
           width=10,height = 6)
  }
  
}


table(output_table$PPH4>0.5)
write.table(output_table,file = '../result/coloc_result/coloc-merge-mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS.txt',
            row.names = F,quote = F)


ls ../figure/locus-cis-mQTLs*|wc -l