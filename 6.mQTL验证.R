library(data.table)
mqtl <- fread('../result/smr_result/smr_mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS.txt')
# mQTLs_LBC_BSGS_meta_GCST90301704----
library(data.table)
library(tidyverse)
dat <- fread('../result/smr_mQTLs_LBC_BSGS_meta_GCST90301704.msmr')
table(dat$probeID %in% mqtl$probeID)
dat <- dat[dat$probeID %in% mqtl$probeID,]
head(dat)
dat <- merge(dat,mqtl[,c('probeID','SYMBOL')], by = 'probeID')
dat$fdr_SMR <- p.adjust(dat$p_SMR,method = 'BH')
head(dat)

table(dat$p_HEIDI > 0.01 & dat$p_SMR_multi < 0.05 & dat$p_SMR < 0.05)
dat1 <- dat[dat$p_HEIDI > 0.01 & dat$p_SMR_multi < 0.05 & dat$p_SMR < 0.05,]
head(dat1)
table(is.na(dat1$p_HEIDI))
dat1 <- dat1[!is.na(dat1$p_HEIDI),]

dat1$OR_SMR <- round(exp(dat1$b_SMR),3)
dat1$`95% CI_SMR` <- paste(round(exp(dat1$b_SMR - dat1$se_SMR * 1.96),3),
                           '-',
                           round(exp(dat1$b_SMR + dat1$se_SMR * 1.96),3))
dat1 <- arrange(dat1,SYMBOL)
dat1 <- dat1 %>%
  dplyr::select(probeID,SYMBOL,p_SMR,p_SMR_multi,fdr_SMR,OR_SMR,`95% CI_SMR`,everything())
if (!file.exists('../result/smr_result')) {
  dir.create('../result/smr_result')
} 
write.table(dat1,'../result/smr_result/smr_mQTLs_LBC_BSGS_meta_GCST90301704.txt',
            row.names = F, sep = '\t',quote = F)

