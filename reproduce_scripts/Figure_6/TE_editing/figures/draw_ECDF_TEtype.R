library(ggplot2)

pdf("H1-hESC.TE_vs_gene.TEtype.ECDF.pdf", w=4,h=3)
alldata <- read.table("H1-hESC.TE_vs_gene.TEtype.ECDF.4R",header=T, sep="\t")
ggplot(alldata, aes(editLevel, color = class)) + stat_ecdf(geom = "line")+
scale_x_continuous(limits = c(0.1, 0.5))+
scale_color_manual(values=c('#ED7D31','#A6A6A6','#FFC000','#70AD47','#0068B3'))+
theme_classic()+
theme(axis.text.x = element_text( color = "black", 
                           size=14, angle=0),
          axis.text.y = element_text( color = "black", 
                           size=14, angle=0))

dev.off()

pdf("DE-H1.TE_vs_gene.TEtype.ECDF.pdf", w=4,h=3)
alldata <- read.table("DE-H1.TE_vs_gene.TEtype.ECDF.4R",header=T, sep="\t")
ggplot(alldata, aes(editLevel, color = class)) + stat_ecdf(geom = "line")+
scale_x_continuous(limits = c(0.1, 0.5))+
scale_color_manual(values=c('#ED7D31','#A6A6A6','#FFC000','#70AD47','#0068B3'))+
theme_classic()+
theme(axis.text.x = element_text( color = "black", 
                           size=14, angle=0),
          axis.text.y = element_text( color = "black", 
                           size=14, angle=0))

dev.off()

pdf("AFG-H1.TE_vs_gene.TEtype.ECDF.pdf", w=4,h=3)
alldata <- read.table("AFG-H1.TE_vs_gene.TEtype.ECDF.4R",header=T, sep="\t")
ggplot(alldata, aes(editLevel, color = class)) + stat_ecdf(geom = "line")+
scale_x_continuous(limits = c(0.1, 0.5))+
scale_color_manual(values=c('#ED7D31','#A6A6A6','#FFC000','#70AD47','#0068B3'))+
theme_classic()+
theme(axis.text.x = element_text( color = "black", 
                           size=14, angle=0),
          axis.text.y = element_text( color = "black", 
                           size=14, angle=0))

dev.off()


pdf("PGC-H1.TE_vs_gene.TEtype.ECDF.pdf", w=4,h=3)
alldata <- read.table("PGC-H1.TE_vs_gene.TEtype.ECDF.4R",header=T, sep="\t")
ggplot(alldata, aes(editLevel, color = class)) + stat_ecdf(geom = "line")+
scale_x_continuous(limits = c(0.1, 0.5))+
scale_color_manual(values=c('#ED7D31','#A6A6A6','#FFC000','#70AD47','#0068B3'))+
theme_classic()+
theme(axis.text.x = element_text( color = "black", 
                           size=14, angle=0),
          axis.text.y = element_text( color = "black", 
                           size=14, angle=0))

dev.off()




pdf("H9-hESC.TE_vs_gene.TEtype.ECDF.pdf", w=4,h=3)
alldata <- read.table("H9-hESC.TE_vs_gene.TEtype.ECDF.4R",header=T, sep="\t")
ggplot(alldata, aes(editLevel, color = class)) + stat_ecdf(geom = "line")+
scale_x_continuous(limits = c(0.1, 0.5))+
scale_color_manual(values=c('#ED7D31','#A6A6A6','#FFC000','#70AD47','#0068B3'))+
theme_classic()+
theme(axis.text.x = element_text( color = "black", 
                           size=14, angle=0),
          axis.text.y = element_text( color = "black", 
                           size=14, angle=0))

dev.off()

pdf("DE-H9.TE_vs_gene.TEtype.ECDF.pdf", w=4,h=3)
alldata <- read.table("DE-H9.TE_vs_gene.TEtype.ECDF.4R",header=T, sep="\t")
ggplot(alldata, aes(editLevel, color = class)) + stat_ecdf(geom = "line")+
scale_x_continuous(limits = c(0.1, 0.5))+
scale_color_manual(values=c('#ED7D31','#A6A6A6','#FFC000','#70AD47','#0068B3'))+
theme_classic()+
theme(axis.text.x = element_text( color = "black", 
                           size=14, angle=0),
          axis.text.y = element_text( color = "black", 
                           size=14, angle=0))

dev.off()

pdf("AFG-H9.TE_vs_gene.TEtype.ECDF.pdf", w=4,h=3)
alldata <- read.table("AFG-H9.TE_vs_gene.TEtype.ECDF.4R",header=T, sep="\t")
ggplot(alldata, aes(editLevel, color = class)) + stat_ecdf(geom = "line")+
scale_x_continuous(limits = c(0.1, 0.5))+
scale_color_manual(values=c('#ED7D31','#A6A6A6','#FFC000','#70AD47','#0068B3'))+
theme_classic()+
theme(axis.text.x = element_text( color = "black", 
                           size=14, angle=0),
          axis.text.y = element_text( color = "black", 
                           size=14, angle=0))

dev.off()

