library(ggplot2)


pdf("allsample.editSite.gene_vs_te.pdf",w=6,h=2)
alldata <- read.table("allsample.editSite.gene_vs_te.4R",header=T, sep="\t")

# grouped boxplot
ggplot(alldata, aes(x=sample, y=editSite, fill=class))+ 
    geom_boxplot(outlier.shape = NA)+
    scale_fill_manual(values=c("grey","#971F15"))+
    scale_x_discrete(limits=c("H1-hESC","DE-H1","AFG-H1","PGC-H1","H9-hESC","DE-H9","AFG-H9"))+
    coord_cartesian(ylim = c(0, 30))+
    theme_classic()

dev.off()

pdf("allsample.editLevel.gene_vs_te.pdf",w=6,h=2)
alldata <- read.table("allsample.editSite.gene_vs_te.4R",header=T, sep="\t")

# grouped boxplot
ggplot(alldata, aes(x=sample, y=editLevel, fill=class))+ 
    geom_boxplot(outlier.shape = NA)+
    scale_fill_manual(values=c("grey","#971F15"))+
    scale_x_discrete(limits=c("H1-hESC","DE-H1","AFG-H1","PGC-H1","H9-hESC","DE-H9","AFG-H9"))+
    coord_cartesian(ylim = c(0, 0.5))+
    theme_classic()

dev.off()