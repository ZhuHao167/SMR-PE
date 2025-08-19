# eQTL共定位数据准备----
dat1 <- fread('../result/smr_result/smr_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.txt')#SMR分析筛选之后的数据

# 基于SMR筛选的位点，准备QTL数据
prob <- unique(dat1$probeID)
write.table(prob,file = "../data/eQTL_prob.txt",sep = "\t",row.names = F,col.names = F,quote = F)

```{bash}
### 创建对应文件夹
mkdir ../data/QTL_Coloc_Prepare
mkdir ../data/QTL_Coloc_Prepare/Single

### 批量提取目标探针，并行运行加速完成，以eQTLs为例
export PATH=$PATH:/mnt/data/jiaozicong/1.Tools/parallel-20240422/src

sort -u ../data/eQTL_prob.txt | parallel -j 10 "cat /mnt/bak/smr_data/Coloc_QTLs_Data/cis-eQTLs_eQTLGen.txt | grep {} > ../data/QTL_Coloc_Prepare/Single/{}.txt"



### 合并上一步输出的所有探针数据
cat ../data/QTL_Coloc_Prepare/Single/*.txt > ../data/QTL_Coloc_Prepare/eQTLs_coloc.txt

### 将单独的探针数据清理干净
rm ../data/QTL_Coloc_Prepare/Single/*.txt

```
# mQTL共定位数据准备----
dat1 <- fread('../result/smr_result/smr_mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS.txt')#SMR分析筛选之后的数据

# 基于SMR筛选的位点，准备QTL数据
prob <- unique(dat1$probeID)
write.table(prob,file = "../data/mQTL_prob.txt",sep = "\t",row.names = F,col.names = F,quote = F)

```{bash}
### 创建对应文件夹
# mkdir ../data/QTL_Coloc_Prepare
# mkdir ../data/QTL_Coloc_Prepare/Single

### 批量提取目标探针，并行运行加速完成
export PATH=$PATH:/mnt/data/jiaozicong/1.Tools/parallel-20240422/src

sort -u ../data/mQTL_prob.txt | parallel -j 10 "cat /mnt/bak/smr_data/Coloc_QTLs_Data/LBC_BSGS_meta.txt | grep {} > ../data/QTL_Coloc_Prepare/Single/{}.txt"


### 合并上一步输出的所有探针数据
cat ../data/QTL_Coloc_Prepare/Single/*.txt > ../data/QTL_Coloc_Prepare/mQTLs_coloc.txt

### 将单独的探针数据清理干净
rm ../data/QTL_Coloc_Prepare/Single/*.txt

```