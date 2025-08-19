# eQTLs_eQTLGen----
library(data.table)
PPH4 <- fread('../result/coloc_result/coloc-merge-eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.txt')
dat1 <- fread('../result/smr_result/smr_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.txt')
dat1 <- merge(dat1,PPH4,by.x = 'probeID', by.y = 'ProbID')
dat1 <- arrange(dat1,probeID)
dat1$PPH4 <- round(dat1$PPH4,3)
dat1$PPH4 <- ifelse(dat1$PPH4 > 0.5, 
                    paste0(dat1$PPH4," *"),
                    dat1$PPH4)
library(stringr)
dat1$CI_SMR_lower <- str_split(dat1$`95% CI_SMR`,' - ',simplify = T)[,1] %>% as.numeric()
dat1$CI_SMR_upper <- str_split(dat1$`95% CI_SMR`,' - ',simplify = T)[,2] %>% as.numeric()

label <- cbind(c("probeID", dat1$probeID), 
               c("Gene", dat1$SYMBOL),
               c("OR (95%CI)", paste0(dat1$OR_SMR,'(',dat1$`95% CI_SMR`,')')),
               c("p_SMR", round(dat1$p_SMR,5)),
               # c("fdr_SMR", round(dat1$fdr_SMR,5)),
               c("p_SMR_multi", round(dat1$p_SMR_multi,5)),
               c("PPH4", dat1$PPH4))
library(forestplot)
jpeg("../figure/forestplot_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.jpeg", 
     res = 300,units = "cm",
     width = 25, height = 15)
forestplot(labeltext = label,#设置用于文本展示的列
           mean = c(NA,dat1$OR_SMR), #设置均值
           lower = c(NA,dat1$CI_SMR_lower), #设置均值的lowlimits限
           upper = c(NA,dat1$CI_SMR_upper), #设置均值的uplimits限
           zero = 1, #设置参照值，此处我们展示的是OR值，故参照值是1
           boxsize = 0.2, #设置点估计的方形大小
           # lineheight = unit(4,'mm'),#设置图形中的行距
           # colgap = unit(15,'mm'),#设置图形中的列间距
           lwd.zero = 2,#设置参考线的粗细
           lwd.ci = 2,#设置区间估计线的粗细
           col=fpColors(box='black',lines = 'black',zero = 'grey'),
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           #xlog = T,
           # xlab="The estimates",#设置x轴标签
           hrzl_lines = list("1" = gpar(lty=1,lwd=1,col="black"),
                             "2" = gpar(lty=1,lwd=1,col="black")),#在第一、二行上方添加黑色横线
           lwd.xaxis=2,#设置X轴线的粗细
           lty.ci = "solid",#线条的形状
           graph.pos = 3,#设置森林图的位置，此处设置为3，则出现在第三列列
           clip=c(0,3),#箭头
           xticks = c(0,1,2),#刻度线
           txt_gp=fpTxtGp(label=gpar(fontfamily="serif",cex=1),#设置字体和大小
                          ticks=gpar(fontfamily="serif",cex=1), 
                          xlab=gpar(fontfamily="serif",cex=0.8), 
                          title=gpar(fontfamily="serif",cex=1))
)|> 
  fp_set_zebra_style("#EFEFEF")
dev.off()

# mQTLs_LBC_BSGS_meta----
library(data.table)
PPH4 <- fread('../result/coloc_result/coloc-merge-mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS.txt')
PPH4 <- unique(PPH4)
dat1 <- fread('../result/smr_result/smr_mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS.txt')
dat1 <- merge(dat1,PPH4,by.x = 'probeID', by.y = 'ProbID')
eqtl <- fread('../result/smr_result/smr_eQTLs_eQTLGen_finngen_R10_O15_PREECLAMPS.txt')
# dat1 <- dat1 %>%
#   filter(PPH4 > 0.5)
dat1 <- dat1[dat1$SYMBOL %in% eqtl$SYMBOL,]
dat1 <- arrange(dat1,SYMBOL)
dat1$PPH4 <- round(dat1$PPH4,3)
dat1$PPH4 <- ifelse(dat1$PPH4 > 0.5, 
                    paste0(dat1$PPH4," *"),
                    dat1$PPH4)
library(stringr)
dat1$CI_SMR_lower <- str_split(dat1$`95% CI_SMR`,' - ',simplify = T)[,1] %>% as.numeric()
dat1$CI_SMR_upper <- str_split(dat1$`95% CI_SMR`,' - ',simplify = T)[,2] %>% as.numeric()

label <- cbind(c("probeID", dat1$probeID), 
               c("Gene", dat1$SYMBOL),
               c("OR (95%CI)", paste0(dat1$OR_SMR,'(',dat1$`95% CI_SMR`,')')),
               c("p_SMR", round(dat1$p_SMR,5)),
               # c("fdr_SMR", round(dat1$fdr_SMR,5)),
               c("p_SMR_multi", round(dat1$p_SMR_multi,5)),
               c("PPH4", dat1$PPH4))
library(forestplot)
jpeg("../figure/forestplot_mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS.jpeg", 
     res = 300,units = "cm",
     width = 25, height = 30)
jpeg("../figure/forestplot_mQTLs_LBC_BSGS_finngen_R10_O15_PREECLAMPS_1.jpeg", 
     res = 300,units = "cm",
     width = 25, height = 10)
forestplot(labeltext = label,#设置用于文本展示的列
           mean = c(NA,dat1$OR_SMR), #设置均值
           lower = c(NA,dat1$CI_SMR_lower), #设置均值的lowlimits限
           upper = c(NA,dat1$CI_SMR_upper), #设置均值的uplimits限
           zero = 1, #设置参照值，此处我们展示的是OR值，故参照值是1
           boxsize = 0.2, #设置点估计的方形大小
           # lineheight = unit(4,'mm'),#设置图形中的行距
           # colgap = unit(15,'mm'),#设置图形中的列间距
           lwd.zero = 2,#设置参考线的粗细
           lwd.ci = 2,#设置区间估计线的粗细
           col=fpColors(box='black',lines = 'black',zero = 'grey'),
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           #xlog = T,
           # xlab="The estimates",#设置x轴标签
           hrzl_lines = list("1" = gpar(lty=1,lwd=1,col="black"),
                             "2" = gpar(lty=1,lwd=1,col="black")),#在第一、二行上方添加黑色横线
           lwd.xaxis=2,#设置X轴线的粗细
           lty.ci = "solid",#线条的形状
           graph.pos = 3,#设置森林图的位置，此处设置为3，则出现在第三列列
           clip=c(0,3),#箭头
           xticks = c(0,1,2),#刻度线
           txt_gp=fpTxtGp(label=gpar(fontfamily="serif",cex=1),#设置字体和大小
                          ticks=gpar(fontfamily="serif",cex=1), 
                          xlab=gpar(fontfamily="serif",cex=0.8), 
                          title=gpar(fontfamily="serif",cex=1))
)|> 
  fp_set_zebra_style("#EFEFEF")
dev.off()

# pQTLs_Fenland_4674_protein----
library(data.table)
PPH4 <- fread('../result/coloc_result/coloc-merge-pQTLs_Fenland_4674_protein_finngen_R10_O15_PREECLAMPS.txt')
dat1 <- fread('../result/smr_result/smr_pQTLs_Fenland_4674_protein_finngen_R10_O15_PREECLAMPS.txt')
dat1 <- merge(dat1,PPH4,by.x = 'probeID', by.y = 'ProbID')
dat1 <- arrange(dat1,probeID)

dat1$PPH4 <- round(dat1$PPH4,3)
dat1$PPH4 <- ifelse(dat1$PPH4 > 0.5, 
                    paste0(dat1$PPH4," *"),
                    dat1$PPH4)
library(stringr)
dat1$CI_SMR_lower <- str_split(dat1$`95% CI_SMR`,' - ',simplify = T)[,1] %>% as.numeric()
dat1$CI_SMR_upper <- str_split(dat1$`95% CI_SMR`,' - ',simplify = T)[,2] %>% as.numeric()

label <- cbind(c("Gene", dat1$probeID), 
               c("Protein", dat1$Gene),
               c("OR (95%CI)", paste0(dat1$OR_SMR,'(',dat1$`95% CI_SMR`,')')),
               c("p_SMR", round(dat1$p_SMR,5)),
               # c("fdr_SMR", round(dat1$fdr_SMR,5)),
               c("p_SMR_multi", round(dat1$p_SMR_multi,5)),
               c("PPH4", dat1$PPH4))
library(forestplot)
jpeg("../figure/forestplot_pQTLs_Fenland_4674_protein_finngen_R10_O15_PREECLAMPS.jpeg", 
     res = 300,units = "cm",
     width = 25, height = 5)
forestplot(labeltext = label,#设置用于文本展示的列
           mean = c(NA,dat1$OR_SMR), #设置均值
           lower = c(NA,dat1$CI_SMR_lower), #设置均值的lowlimits限
           upper = c(NA,dat1$CI_SMR_upper), #设置均值的uplimits限
           zero = 1, #设置参照值，此处我们展示的是OR值，故参照值是1
           boxsize = 0.2, #设置点估计的方形大小
           # lineheight = unit(4,'mm'),#设置图形中的行距
           # colgap = unit(15,'mm'),#设置图形中的列间距
           lwd.zero = 2,#设置参考线的粗细
           lwd.ci = 2,#设置区间估计线的粗细
           col=fpColors(box='black',lines = 'black',zero = 'grey'),
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           #xlog = T,
           # xlab="The estimates",#设置x轴标签
           hrzl_lines = list("1" = gpar(lty=1,lwd=1,col="black"),
                             "2" = gpar(lty=1,lwd=1,col="black")),#在第一、二行上方添加黑色横线
           lwd.xaxis=2,#设置X轴线的粗细
           lty.ci = "solid",#线条的形状
           graph.pos = 3,#设置森林图的位置，此处设置为3，则出现在第三列列
           clip=c(0,3),#箭头
           xticks = c(0,1,2),#刻度线
           txt_gp=fpTxtGp(label=gpar(fontfamily="serif",cex=1),#设置字体和大小
                          ticks=gpar(fontfamily="serif",cex=1), 
                          xlab=gpar(fontfamily="serif",cex=0.8), 
                          title=gpar(fontfamily="serif",cex=1))
)|> 
  fp_set_zebra_style("#EFEFEF")
dev.off()
