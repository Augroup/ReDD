library(ggplot2)

alldata <- read.table("allsample.TE-Gene.TE_nonTE.editLevel.txt.4R",header=T, sep="\t")
pdf("allsample.TE-Gene.TE_nonTE.editLevel.pdf", w=7,h=4)

ggplot(alldata, aes(x=Sample, y=EditLevel,fill = Type)) + 
geom_violin(position="dodge", alpha=0.7, outlier.colour="transparent")+ 
#geom_boxplot(width=0.1, outlier.shape = NA)+
theme_classic()+
scale_y_continuous(limits = c(0, 1))+
scale_fill_manual(values=c('#A6A6A6','#ED7D31','#FFC000','#70AD47','#F90000','#0068B3'))+

theme(axis.text.x = element_text( color = "black", 
                           size=12, angle=0),
          axis.text.y = element_text( color = "black", 
                           size=12, angle=0))+
#scale_x_discrete(limits=c("NonTE", "DNA","LTR","LINE","SINE","Retroposon"))
scale_x_discrete(limits=c("H1-hESC","DE-H1","AFG-H1","PGC-H1","H9-hESC","DE-H9","AFG-H9"))
dev.off()