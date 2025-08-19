library(data.table)
dat <- fread('../result/smr_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.msmr',
             header = T)
library(clusterProfiler)
library(org.Hs.eg.db)
gene_ENSEMBL<-bitr(dat$probeID,
                   fromType="ENSEMBL",
                   toType=c('ENTREZID','SYMBOL'),
                   OrgDb = "org.Hs.eg.db")
dat <- merge(dat,gene_ENSEMBL[,c('ENSEMBL','SYMBOL')],by.x = 'probeID',by.y = 'ENSEMBL')
qtl <- fread('../result/smr_result/smr_mqtl_eqtl.txt')

gene <- unique(qtl$SYMBOL)
gene <- gene[-2]

highlight_names <- qtl$Outco_Gene[qtl$SYMBOL %in% gene]
gwas <- dat[,c('SYMBOL','probeID','ProbeChr','Probe_bp','p_SMR_multi')]

source("/2.Common_data/refine_manh.R") 
jpeg("../figure/manhattan_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.jpeg", 
     res = 300,units = "cm", width = 30, height = 10)
## 绘图，直接使用 Manhattan_refine()即可，参数一般来说可不进行修改
Manhattan_refine(Background_data = gwas, ## 输入数据
                 highlight_names = highlight_names, ## 需要highlight的重点结果
                 thresholds_line_color = "#E58601", ## 阈值线颜色
                 points_color = "#B0AFA2",  ## 背景点颜色
                 highlight_color = "#046C9A",  ## 高亮点颜色
                 pointsize = 1.5, ## 背景点大小
                 highlight_size = 2, ## 高亮点大小
                 label_size = 3.5, ## 标签大小
                 force=0.2, ## 标签间斥力
                 pval_names = "p_SMR_multi", ## 阈值类型标注
                 thresholds = 0.05, ## 阈值标准
                 legend_position=c(0.9,0.8),
                 legend_labels = "SMR: cis-eQTL Plasma") ## 图例内容，注意修改不同QTL

dev.off()


library(data.table)
dat <- fread('../result/smr_mQTLs_LBC_BSGS_meta_finngen_R10_O15_PREECLAMPS.msmr',
             header = T)
probe = data.table::fread('/mnt/bak/smr_data/mqtl/LBC_BSGS_meta/meta/GPL16304_LBC_BSGS.txt')
names(probe)[5] <- 'SYMBOL'
dat <- merge(dat,probe[,c(1,5)],by.x = 'probeID',by.y = 'ID')

qtl <- fread('../result/smr_result/smr_mqtl_eqtl.txt')
dat1 <- fread('../result/smr_result/smr_mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS.txt')

gene <- unique(qtl$SYMBOL)
gene <- gene[-2]

highlight_names <- qtl$Expo_ID[qtl$SYMBOL %in% gene]
highlight_names <- highlight_names[-2]

gwas <- dat[,c('SYMBOL','probeID','ProbeChr','Probe_bp','p_SMR_multi')]

source("/2.Common_data/refine_manh.R") 
jpeg("../figure/manhattan_mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS.jpeg", 
     res = 300,units = "cm", width = 30, height = 15)
## 绘图，直接使用 Manhattan_refine()即可，参数一般来说可不进行修改
Manhattan_refine(Background_data = gwas, ## 输入数据
                 highlight_names = highlight_names, ## 需要highlight的重点结果
                 thresholds_line_color = "#E58601", ## 阈值线颜色
                 points_color = "#A5C2A3",  ## 背景点颜色
                 highlight_color = "#972D15",  ## 高亮点颜色
                 pointsize = 1.5, ## 背景点大小
                 highlight_size = 2, ## 高亮点大小
                 label_size = 3.5, ## 标签大小
                 force=0.2, ## 标签间斥力
                 pval_names = "p_SMR_multi", ## 阈值类型标注
                 thresholds = 0.05, ## 阈值标准
                 legend_position=c(0.9,0.8),
                 legend_labels = "SMR: cis-mQTL Plasma") ## 图例内容，注意修改不同QTL

dev.off()


library(data.table)
dat <- fread('../result/smr_pQTLs_Fenland_4674_protein_finngen_R10_O15_PREECLAMPS.msmr',
             header = T)
dat$SYMBOL <- dat$probeID
# eqtl <- fread('../result/smr_result/smr_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.txt')
# mqtl <- fread('../result/smr_result/smr_mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS.txt')
pqtl <- fread('../result/smr_result/smr_pQTLs_Fenland_4674_protein_finngen_R10_O15_PREECLAMPS.txt')

# gene <- intersect(intersect(eqtl$SYMBOL,mqtl$SYMBOL),pqtl$probeID)
gene <- pqtl$probeID
highlight_names <- pqtl$probeID[pqtl$probeID %in% gene]
gwas <- dat[,c('SYMBOL','probeID','ProbeChr','Probe_bp','p_SMR_multi')]

source("/2.Common_data/refine_manh.R") 
jpeg("../figure/manhattan_pQTLs_Fenland_4674_protein_finngen_R10_O15_PREECLAMPS.jpeg", 
     res = 300,units = "cm", width = 30, height = 10)
## 绘图，直接使用 Manhattan_refine()即可，参数一般来说可不进行修改
Manhattan_refine(Background_data = gwas, ## 输入数据
                 highlight_names = highlight_names, ## 需要highlight的重点结果
                 thresholds_line_color = "#E58601", ## 阈值线颜色
                 points_color = "#C6CDF7",  ## 背景点颜色
                 highlight_color = "#02401B",  ## 高亮点颜色
                 pointsize = 1.5, ## 背景点大小
                 highlight_size = 2, ## 高亮点大小
                 label_size = 3.5, ## 标签大小
                 force=0.2, ## 标签间斥力
                 pval_names = "p_SMR_multi", ## 阈值类型标注
                 thresholds = 0.05, ## 阈值标准
                 legend_position=c(0.9,0.8),
                 legend_labels = "SMR: cis-pQTL Plasma") ## 图例内容，注意修改不同QTL

dev.off()
