library(ggplot2)


pdf("H1_3cell.editdiversity.densityplot.pdf",w=5,h=4)
data <- read.table("H1_3cell.editdiversity", header = TRUE, sep = "\t")

p1 <-ggplot(data, aes(x=cv, color=sample)) +
geom_density(linewidth=1.5,alpha=0.2)+theme_classic()+
scale_x_continuous(limits = c(0, 1))+
scale_color_manual(values=c('#682B95','#2B66AD','#509BEA','#00B050'))
#scale_fill_manual(values=c('#ED7D31'))+
#scale_color_manual(values=c('#ED7D31'))

p1+theme(axis.text.x = element_text(face="bold", color = "black", 
                           size=20, angle=0),
          axis.text.y = element_text(face="bold", color = "black", 
                           size=20, angle=0))
dev.off()


pdf("H9_3cell.editdiversity.densityplot.pdf",w=5,h=4)
data <- read.table("H9_3cell.editdiversity", header = TRUE, sep = "\t")

p1 <-ggplot(data, aes(x=cv, color=sample)) +
geom_density(linewidth=1.5,alpha=0.2)+theme_classic()+
scale_x_continuous(limits = c(0, 1))+
scale_color_manual(values=c('#682B95','#2B66AD','#509BEA','#00B050'))
#scale_fill_manual(values=c('#ED7D31'))+
#scale_color_manual(values=c('#ED7D31'))

p1+theme(axis.text.x = element_text(face="bold", color = "black", 
                           size=20, angle=0),
          axis.text.y = element_text(face="bold", color = "black", 
                           size=20, angle=0))
dev.off()