# 下载GWAS数据 ----

```{bash}
mkdir rawdata data main result figure
cd rawdata
/1.Tools/axel-2.17.14
axel -a -n 50 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90301001-GCST90302000/GCST90301704/harmonised/GCST90301704.h.tsv.gz

ls /mnt/bak/finngen_r10/|grep "PREECLAMPS"

/1.Tools/axel-2.17.14
axel -a -n 50 https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_O15_PREECLAMPS.gz

```

# GWAS catalog数据整理----

library(data.table)
library(R.utils)
data <- fread('../rawdata/GCST90301704.h.tsv.gz')
names(data)

#提取需要的列
f <- data[, c("SNPID", "effect_allele", "other_allele",
              "effect_allele_frequency", "beta", "standard_error",
              "p_value")]
f$n <- NA
#更改为SMR分析的列名
colnames(f) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "n")
library(tidyverse)
table(grepl(',',f$SNP))
f <- f %>%
  filter(SNP != ""& !duplicated(SNP) & SNP != ".")
write.table(f, "../data/GCST90301704.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)


# 芬兰数据整理----

library(data.table)
library(R.utils)
finngen <- fread('../rawdata/finngen_R10_O15_PREECLAMPS.gz')
names(finngen)
#如果只有OR值，需要将OR转换为BETA，OR转换为b（b=log（OR））
# finngen$b <- log(finngen$OR)
#提取需要的列
f <- finngen[, c("rsids", "alt", "ref", "af_alt", "beta", "sebeta",
                 "pval")]
f$n <- NA
#更改为SMR分析的列名
colnames(f) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "n")
library(tidyverse)
table(grepl(',',f$SNP))
# f <- f[!grepl(',',f$SNP),]
f <- f %>% 
  separate_rows(SNP, sep = ',') %>% 
  filter(SNP != ""& !duplicated(SNP) & SNP != ".") %>%
  data.frame() 
write.table(f, "../data/finngen_R10_O15_PREECLAMPS.txt", sep = "\t", 
            row.names = FALSE, quote = FALSE)


