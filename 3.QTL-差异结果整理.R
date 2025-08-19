library(data.table)
library(tidyverse)
# eqtl---- 
eqtl <- read.table('../result/smr_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.msmr',
                   header = T)
eqtl$fdr_SMR <- p.adjust(eqtl$p_SMR,method = 'BH')
gene <- fread('../data/genelist.txt')
table(eqtl$Gene %in% gene$ENSEMBL)
eqtl <- merge(eqtl,gene[,c(1,3)],by.x = 'probeID',by.y = 'ENSEMBL')
table(eqtl$p_HEIDI > 0.01 & eqtl$p_SMR_multi < 0.05 & eqtl$p_SMR < 0.05)
eqtl1 <- eqtl[eqtl$p_HEIDI > 0.01 & eqtl$p_SMR_multi < 0.05 & eqtl$p_SMR < 0.05,]
table(is.na(eqtl1$p_HEIDI))
eqtl1 <- eqtl1[!is.na(eqtl1$p_HEIDI),]
# sort(eqtl1$p_eQTL)

eqtl1$OR_SMR <- round(exp(eqtl1$b_SMR),3)
eqtl1$`95% CI_SMR` <- paste(round(exp(eqtl1$b_SMR - eqtl1$se_SMR * 1.96),3),
                      '-',
                      round(exp(eqtl1$b_SMR + eqtl1$se_SMR * 1.96),3))
eqtl1 <- arrange(eqtl1,SYMBOL)
eqtl1 <- eqtl1 %>%
  dplyr::select(probeID,SYMBOL,p_SMR,p_SMR_multi,fdr_SMR,OR_SMR,`95% CI_SMR`,everything())
if (!file.exists('../result/smr_result')) {
  dir.create('../result/smr_result')
} 
write.table(eqtl1,'../result/smr_result/smr_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.txt',
            row.names = F, sep = '\t',quote = F)

# mqtl---- 
mqtl <- read.table('../result/smr_mQTLs_LBC_BSGS_meta_finngen_R10_O15_PREECLAMPS.msmr',
                   header = T)
mqtl$fdr_SMR <- p.adjust(mqtl$p_SMR,method = 'BH')
GPL16304 <- fread('/mnt/bak/smr_data/mqtl/LBC_BSGS_meta/meta/GPL16304_LBC_BSGS.txt')
head(GPL16304)
GPL16304 <- GPL16304[,c(1,5)]
names(GPL16304)[2] <- 'SYMBOL' 
mqtl <- merge(mqtl,GPL16304,by.x = 'probeID',by.y = 'ID')
head(mqtl)
table(mqtl$SYMBOL %in% gene$SYMBOL)
mqtl <- mqtl[mqtl$SYMBOL %in% gene$SYMBOL,]
table(mqtl$p_HEIDI > 0.01 & mqtl$p_SMR_multi < 0.05 & mqtl$p_SMR < 0.05)
mqtl1 <- mqtl[mqtl$p_HEIDI > 0.01 & mqtl$p_SMR_multi < 0.05 & mqtl$p_SMR < 0.05,]

table(is.na(mqtl1$p_HEIDI))
mqtl1 <- mqtl1[!is.na(mqtl1$p_HEIDI),]
# sort(mqtl1$p_eQTL)

mqtl1$OR_SMR <- round(exp(mqtl1$b_SMR),3)
mqtl1$`95% CI_SMR` <- paste(round(exp(mqtl1$b_SMR - mqtl1$se_SMR * 1.96),3),
                      '-',
                      round(exp(mqtl1$b_SMR + mqtl1$se_SMR * 1.96),3))
mqtl1 <- arrange(mqtl1,SYMBOL)
mqtl1 <- mqtl1 %>%
  dplyr::select(probeID,SYMBOL,p_SMR,p_SMR_multi,fdr_SMR,OR_SMR,`95% CI_SMR`,everything())

if (!file.exists('../result/smr_result')) {
  dir.create('../result/smr_result')
} 
write.table(mqtl1,'../result/smr_result/smr_mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS.txt',
            row.names = F, sep = '\t',quote = F)
# pqtl---- 
pqtl <- read.table('../result/smr_pQTLs_pqtl_p0.05_PMID_34857953_finngen_R10_O15_PREECLAMPS.msmr',
                   header = T)
pqtl$fdr_SMR <- p.adjust(pqtl$p_SMR,method = 'BH')
table(pqtl$probeID %in% gene$SYMBOL)
pqtl <- pqtl[pqtl$probeID %in% gene$SYMBOL,]
table(pqtl$p_HEIDI > 0.01 & pqtl$p_SMR_multi < 0.05 & pqtl$p_SMR < 0.05)
# table(pqtl$p_HEIDI > 0.01 & pqtl$p_SMR < 0.05)
pqtl1 <- pqtl[pqtl$p_HEIDI > 0.01 & pqtl$p_SMR_multi < 0.05 & pqtl$p_SMR < 0.05,]
head(pqtl1)
table(is.na(pqtl1$p_HEIDI))
pqtl1 <- pqtl1[!is.na(pqtl1$p_HEIDI),]
# sort(pqtl1$p_eQTL)

pqtl1$OR_SMR <- round(exp(pqtl1$b_SMR),3)
pqtl1$`95% CI_SMR` <- paste(round(exp(pqtl1$b_SMR - pqtl1$se_SMR * 1.96),3),
                      '-',
                      round(exp(pqtl1$b_SMR + pqtl1$se_SMR * 1.96),3))
pqtl1 <- arrange(pqtl1,probeID)
pqtl1 <- pqtl1 %>%
  dplyr::select(probeID,Gene,p_SMR,p_SMR_multi,fdr_SMR,OR_SMR,`95% CI_SMR`,everything())

if (!file.exists('../result/smr_result')) {
  dir.create('../result/smr_result')
} 
write.table(pqtl1,'../result/smr_result/smr_pQTLs_pqtl_p0.05_PMID_34857953_finngen_R10_O15_PREECLAMPS.txt',
            row.names = F, sep = '\t',quote = F)
# pqtl_Fenland_4674_protein---- 
pqtl_Fenland_4674_protein <- read.table('../result/smr_pQTLs_Fenland_4674_protein_finngen_R10_O15_PREECLAMPS.msmr',
                                        header = T)
pqtl_Fenland_4674_protein$fdr_SMR <- p.adjust(pqtl_Fenland_4674_protein$p_SMR,method = 'BH')
table(pqtl_Fenland_4674_protein$probeID %in% gene$SYMBOL)
pqtl_Fenland_4674_protein <- pqtl_Fenland_4674_protein[pqtl_Fenland_4674_protein$probeID %in% gene$SYMBOL,]
table(pqtl_Fenland_4674_protein$p_HEIDI > 0.01 & 
        pqtl_Fenland_4674_protein$p_SMR_multi < 0.05 & 
        pqtl_Fenland_4674_protein$p_SMR < 0.05)
# table(pqtl_Fenland_4674_protein$p_HEIDI > 0.01 & pqtl_Fenland_4674_protein$p_SMR < 0.05)
pqtl1_Fenland_4674_protein <- pqtl_Fenland_4674_protein[pqtl_Fenland_4674_protein$p_HEIDI > 0.01 & 
                                                          pqtl_Fenland_4674_protein$p_SMR_multi < 0.05 & 
                                                          pqtl_Fenland_4674_protein$p_SMR < 0.05,]
head(pqtl1_Fenland_4674_protein)
table(is.na(pqtl1_Fenland_4674_protein$p_HEIDI))
pqtl1_Fenland_4674_protein <- pqtl1_Fenland_4674_protein[!is.na(pqtl1_Fenland_4674_protein$p_HEIDI),]
# sort(pqtl1_Fenland_4674_protein$p_eQTL)

pqtl1_Fenland_4674_protein$OR_SMR <- round(exp(pqtl1_Fenland_4674_protein$b_SMR),3)
pqtl1_Fenland_4674_protein$`95% CI_SMR` <- paste(round(exp(pqtl1_Fenland_4674_protein$b_SMR - pqtl1_Fenland_4674_protein$se_SMR * 1.96),3),
                                           '-',
                                           round(exp(pqtl1_Fenland_4674_protein$b_SMR + pqtl1_Fenland_4674_protein$se_SMR * 1.96),3))
pqtl1_Fenland_4674_protein <- arrange(pqtl1_Fenland_4674_protein,probeID)
pqtl1_Fenland_4674_protein <- pqtl1_Fenland_4674_protein %>%
  dplyr::select(probeID,Gene,p_SMR,p_SMR_multi,fdr_SMR,OR_SMR,`95% CI_SMR`,everything())

if (!file.exists('../result/smr_result')) {
  dir.create('../result/smr_result')
} 
write.table(pqtl1_Fenland_4674_protein,'../result/smr_result/smr_pQTLs_Fenland_4674_protein_finngen_R10_O15_PREECLAMPS.txt',
            row.names = F, sep = '\t',quote = F)

# pqtl_UKB_PPP_2937---- 
pqtl_UKB_PPP_2937 <- read.table('../result/smr_pQTLs_UKB_PPP_2937_finngen_R10_O15_PREECLAMPS.msmr',
                                header = T)
pqtl_UKB_PPP_2937$fdr_SMR <- p.adjust(pqtl_UKB_PPP_2937$p_SMR,method = 'BH')
table(pqtl_UKB_PPP_2937$probeID %in% gene$SYMBOL)
pqtl_UKB_PPP_2937 <- pqtl_UKB_PPP_2937[pqtl_UKB_PPP_2937$probeID %in% gene$SYMBOL,]
table(pqtl_UKB_PPP_2937$p_HEIDI > 0.01 & 
        pqtl_UKB_PPP_2937$p_SMR_multi < 0.05 & 
        pqtl_UKB_PPP_2937$p_SMR < 0.05)
# table(pqtl_UKB_PPP_2937$p_HEIDI > 0.01 & pqtl_UKB_PPP_2937$p_SMR < 0.05)
pqtl1_UKB_PPP_2937 <- pqtl_UKB_PPP_2937[pqtl_UKB_PPP_2937$p_HEIDI > 0.01 & 
                                          pqtl_UKB_PPP_2937$p_SMR_multi < 0.05 & 
                                          pqtl_UKB_PPP_2937$p_SMR < 0.05,]
head(pqtl1_UKB_PPP_2937)
table(is.na(pqtl1_UKB_PPP_2937$p_HEIDI))
pqtl1_UKB_PPP_2937 <- pqtl1_UKB_PPP_2937[!is.na(pqtl1_UKB_PPP_2937$p_HEIDI),]
# sort(pqtl1_UKB_PPP_2937$p_eQTL)

pqtl1_UKB_PPP_2937$OR_SMR <- round(exp(pqtl1_UKB_PPP_2937$b_SMR),3)
pqtl1_UKB_PPP_2937$`95% CI_SMR` <- paste(round(exp(pqtl1_UKB_PPP_2937$b_SMR - pqtl1_UKB_PPP_2937$se_SMR * 1.96),3),
                                   '-',
                                   round(exp(pqtl1_UKB_PPP_2937$b_SMR + pqtl1_UKB_PPP_2937$se_SMR * 1.96),3))
pqtl1_UKB_PPP_2937 <- arrange(pqtl1_UKB_PPP_2937,probeID)
pqtl1_UKB_PPP_2937 <- pqtl1_UKB_PPP_2937 %>%
  dplyr::select(probeID,Gene,p_SMR,p_SMR_multi,fdr_SMR,OR_SMR,`95% CI_SMR`,everything())

if (!file.exists('../result/smr_result')) {
  dir.create('../result/smr_result')
} 
write.table(pqtl1_UKB_PPP_2937,'../result/smr_result/smr_pQTLs_UKB_PPP_2937_finngen_R10_O15_PREECLAMPS.txt',
            row.names = F, sep = '\t',quote = F)

intersect(pqtl1_Fenland_4674_protein$probeID,pqtl1$probeID)
intersect(intersect(eqtl1$SYMBOL,mqtl1$SYMBOL),pqtl1_Fenland_4674_protein$probeID)
intersect(intersect(eqtl1$SYMBOL,mqtl1$SYMBOL),pqtl1$probeID)
intersect(intersect(eqtl1$SYMBOL,mqtl1$SYMBOL),pqtl1_UKB_PPP_2937$probeID)
