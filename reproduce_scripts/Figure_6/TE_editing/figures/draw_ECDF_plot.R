library(ggplot2)

pdf("H1-hESC.Alu_type.ECDF.pdf", w=4,h=3)
alldata <- read.table("H1-hESC.Alu_type.Edit.cmp",header=T, sep="\t")
ggplot(alldata, aes(EditLevel, color = Type)) + stat_ecdf(geom = "line")+
scale_x_continuous(limits = c(0.1, 0.5))+
scale_color_manual(values = c("black","orange",  "blue"))+
theme_classic()+
theme(axis.text.x = element_text( color = "black", 
                           size=14, angle=0),
          axis.text.y = element_text( color = "black", 
                           size=14, angle=0))
dev.off()

pdf("H9-hESC.Alu_type.ECDF.pdf", w=4,h=3)
alldata <- read.table("H9-hESC.Alu_type.Edit.cmp",header=T, sep="\t")
ggplot(alldata, aes(EditLevel, color = Type)) + stat_ecdf(geom = "line")+
scale_x_continuous(limits = c(0.1, 0.5))+
scale_color_manual(values = c("black","orange",  "blue"))+
theme_classic()+
theme(axis.text.x = element_text( color = "black", 
                           size=14, angle=0),
          axis.text.y = element_text( color = "black", 
                           size=14, angle=0))
dev.off()


pdf("DE-H1.Alu_type.ECDF.pdf", w=4,h=3)
alldata <- read.table("DE-H1.Alu_type.Edit.cmp",header=T, sep="\t")
ggplot(alldata, aes(EditLevel, color = Type)) + stat_ecdf(geom = "line")+
scale_x_continuous(limits = c(0.1, 0.5))+
scale_color_manual(values = c("black","orange",  "blue"))+
theme_classic()+
theme(axis.text.x = element_text( color = "black", 
                           size=14, angle=0),
          axis.text.y = element_text( color = "black", 
                           size=14, angle=0))
dev.off()


pdf("DE-H9.Alu_type.ECDF.pdf", w=4,h=3)
alldata <- read.table("DE-H9.Alu_type.Edit.cmp",header=T, sep="\t")
ggplot(alldata, aes(EditLevel, color = Type)) + stat_ecdf(geom = "line")+
scale_x_continuous(limits = c(0.1, 0.5))+
scale_color_manual(values = c("black","orange",  "blue"))+
theme_classic()+
theme(axis.text.x = element_text( color = "black", 
                           size=14, angle=0),
          axis.text.y = element_text( color = "black", 
                           size=14, angle=0))
dev.off()




