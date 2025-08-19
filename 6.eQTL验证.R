library(data.table)
eqtl <- fread('../result/smr_result/smr_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.txt')
# eQTLs_eQTLGen_GCST90301704----
library(data.table)
library(tidyverse)
dat <- read.table('../result/smr_eQTLs_eQTLGen_GCST90301704.msmr',
                  header = T)
head(dat)
gene <- fread('../data/genelist.txt')
table(dat$Gene %in% gene$ENSEMBL)
dat <- dat[dat$Gene %in% eqtl$probeID,]
dat <- merge(dat,gene[,c(1,3)],by.x = 'probeID',by.y = 'ENSEMBL')

dat$fdr_SMR <- p.adjust(dat$p_SMR,method = 'BH')

head(dat)
table(dat$p_HEIDI > 0.01 & dat$p_SMR_multi < 0.05 & dat$p_SMR < 0.05)
dat1 <- dat[dat$p_HEIDI > 0.01 & dat$p_SMR_multi < 0.05 & dat$p_SMR < 0.05,] 

table(is.na(dat1$p_HEIDI))
dat1 <- dat1[!is.na(dat1$p_HEIDI),]
sort(dat1$p_eQTL)

dat1$OR_SMR <- round(exp(dat1$b_SMR),3)
dat1$`95% CI_SMR` <- paste(round(exp(dat1$b_SMR - dat1$se_SMR * 1.96),3),
                           '-',
                           round(exp(dat1$b_SMR + dat1$se_SMR * 1.96),3))
dat1 <- arrange(dat1,SYMBOL)
dat1 <- dat1 %>%
  dplyr::select(probeID,SYMBOL,p_SMR,p_SMR_multi,fdr_SMR,OR_SMR,`95% CI_SMR`,everything())
write.table(dat1,'../result/smr_result/smr_eQTLs_eQTLGen_GCST90301704.txt',
            row.names = F, sep = '\t',quote = F)


# eQTLs_GTEx_Uterus_finngen_R10_O15_PREECLAMPS----
library(data.table)
library(tidyverse)
dat <- read.table('../result/smr_eQTLs_GTEx_Uterus_finngen_R10_O15_PREECLAMPS.msmr',
                  header = T)
head(dat)
gene <- fread('../data/genelist.txt')
table(dat$Gene %in% gene$SYMBOL)
dat <- dat[dat$Gene %in% eqtl$SYMBOL,]
# dat <- merge(dat,gene[,c(1,3)],by.x = 'probeID',by.y = 'ENSEMBL')

dat$fdr_SMR <- p.adjust(dat$p_SMR,method = 'BH')

head(dat)
table(dat$p_HEIDI > 0.01 & dat$p_SMR_multi < 0.05 & dat$p_SMR < 0.05)
dat1 <- dat[dat$p_HEIDI > 0.01 & dat$p_SMR_multi < 0.05 & dat$p_SMR < 0.05,]
# dat1 <- dat
table(is.na(dat1$p_HEIDI))
dat1 <- dat1[!is.na(dat1$p_HEIDI),]
sort(dat1$p_eQTL)

dat1$OR_SMR <- round(exp(dat1$b_SMR),3)
dat1$`95% CI_SMR` <- paste(round(exp(dat1$b_SMR - dat1$se_SMR * 1.96),3),
                           '-',
                           round(exp(dat1$b_SMR + dat1$se_SMR * 1.96),3))
dat1 <- arrange(dat1,Gene)
dat1 <- dat1 %>%
  dplyr::select(probeID,Gene,p_SMR,p_SMR_multi,fdr_SMR,OR_SMR,`95% CI_SMR`,everything())
write.table(dat1,'../result/smr_result/smr_eQTLs_GTEx_Uterus_finngen_R10_O15_PREECLAMPS.txt',
            row.names = F, sep = '\t',quote = F)
