
```{bash}
/mnt/bak/smr-1.3.1-linux-x86_64/smr \
--bfile /mnt/bak/smr_data/g1000_eur/g1000_eur \
--beqtl-summary /mnt/bak/smr_data/mqtl/LBC_BSGS_meta/meta/LBC_BSGS_meta_merge \
--beqtl-summary /mnt/bak/smr_data/eqtl/cis-eQTLs_eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense \
--extract-exposure-probe ../data/mQTL_prob.txt \
--extract-outcome-probe ../data/eQTL_prob.txt \
--out ../result/mqtl_eqtl \
--smr-multi
```
dat <- fread('../result/mqtl_eqtl.msmr')
head(dat)
gene <- fread('../data/genelist.txt')
table(dat$Outco_Gene %in% gene$ENSEMBL)
dat <- merge(dat,gene[,c(1,3)],by.x = 'Outco_Gene',by.y = 'ENSEMBL')
dat$fdr_SMR <- p.adjust(dat$p_SMR,method = 'BH')
dat$OR_SMR <- round(exp(dat$b_SMR),3)
dat$`95% CI_SMR` <- paste(round(exp(dat$b_SMR - dat$se_SMR * 1.96),3),
                     '-',
                     round(exp(dat$b_SMR + dat$se_SMR * 1.96),3))
dat <- arrange(dat,SYMBOL)
head(dat)
table(dat$p_HEIDI > 0.01 & dat$p_SMR_multi < 0.05 & dat$p_SMR < 0.05)
dat1 <- dat[dat$p_HEIDI > 0.01 & dat$p_SMR_multi < 0.05 & dat$p_SMR < 0.05,]
table(is.na(dat$p_HEIDI))
dat1 <- dat1[!is.na(dat$p_HEIDI),]
# sort(dat$p_dat)

dat1 <- dat1 %>%
  dplyr::select(Expo_ID,Outco_Gene,SYMBOL,p_SMR,p_SMR_multi,fdr_SMR,OR_SMR,`95% CI_SMR`,everything())

if (!file.exists('../result/smr_result')) {
  dir.create('../result/smr_result')
} 
write.table(dat1,'../result/smr_result/smr_mqtl_eqtl.txt',
            row.names = F, sep = '\t',quote = F)

write.csv(dat,'../result/sup_table/smr_mqtl_eqtl.csv',row.names = F,quote = F)



