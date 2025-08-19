
# pQTLs_Fenland_4674_protein_GCST90301704----
library(data.table)
pqtl <- fread('../result/smr_result/smr_pQTLs_Fenland_4674_protein_finngen_R10_O15_PREECLAMPS.txt')
library(data.table)
library(tidyverse)
dat <- read.table('../result/smr_pQTLs_Fenland_4674_protein_GCST90301704.msmr',
                  header = T)
head(dat)
table(dat$probeID %in% pqtl$probeID)
dat <- dat[dat$probeID %in% pqtl$probeID,]

dat$fdr_SMR <- p.adjust(dat$p_SMR,method = 'BH')
head(dat)
table(dat$p_HEIDI > 0.01 & dat$p_SMR_multi < 0.05 & dat$p_SMR < 0.05)
dat1 <- dat[dat$p_HEIDI > 0.01 & dat$p_SMR_multi < 0.05 & dat$p_SMR < 0.05,]
# dat1 <- dat
head(dat1)
table(is.na(dat1$p_HEIDI))
dat1 <- dat1[!is.na(dat1$p_HEIDI),]
sort(dat1$p_eQTL)

dat1$OR_SMR <- round(exp(dat1$b_SMR),3)
dat1$`95% CI_SMR` <- paste(round(exp(dat1$b_SMR - dat1$se_SMR * 1.96),3),
                           '-',
                           round(exp(dat1$b_SMR + dat1$se_SMR * 1.96),3))
dat1 <- arrange(dat1,probeID)
dat1 <- dat1 %>%
  dplyr::select(probeID,p_SMR,p_SMR_multi,fdr_SMR,OR_SMR,`95% CI_SMR`,everything())
if (!file.exists('../result/smr_result')) {
  dir.create('../result/smr_result')
} 
write.table(dat1,'../result/smr_result/smr_pQTLs_Fenland_4674_protein_GCST90301704.txt',
            row.names = F, sep = '\t',quote = F)

