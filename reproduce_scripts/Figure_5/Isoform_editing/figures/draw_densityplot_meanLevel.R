library(ggplot2)

pdf("H1.combined.isoform.meanLevel.txt.densityplot.pdf",w=5,h=4)
data <- read.table("H1.combined.isoform.meanLevel.txt", header = TRUE, sep = "\t")

p1 <-ggplot(data, aes(x=editlevel, color=sample)) +
geom_density(linewidth=1.2,alpha=0.5)+theme_classic()+
scale_x_continuous(limits = c(0, 1))+
#scale_fill_manual(values=c('#ED7D31'))+
scale_color_manual(values=c('#682B95','#2B66AD','#509BEA','#00B050'))

p1+theme(axis.text.x = element_text(face="bold", color = "black", 
                           size=18, angle=0),
          axis.text.y = element_text(face="bold", color = "black", 
                           size=18, angle=0))
dev.off()

summary(data$editlevel)


pdf("H9.combined.isoform.meanLevel.txt.densityplot.pdf",w=5,h=4)
data <- read.table("H9.combined.isoform.meanLevel.txt", header = TRUE, sep = "\t")

p1 <-ggplot(data, aes(x=editlevel, color=sample)) +
geom_density(linewidth=1.2,alpha=0.5)+theme_classic()+
scale_x_continuous(limits = c(0, 1))+
scale_y_continuous(limits = c(0, 12))+
#scale_fill_manual(values=c('#ED7D31'))+
scale_color_manual(values=c('#682B95','#2B66AD','#509BEA','#00B050'))

p1+theme(axis.text.x = element_text(face="bold", color = "black", 
                           size=20, angle=0),
          axis.text.y = element_text(face="bold", color = "black", 
                           size=20, angle=0))
dev.off()