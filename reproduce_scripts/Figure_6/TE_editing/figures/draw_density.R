library(ggpubr)
library(dplyr) 
library(ggplot2)
theme_set(theme_pubclean())

file_names <- list ("H1-hESC","DE-H1","AFG-H1","PGC-H1","H9-hESC","DE-H9","AFG-H9")

for (file in file_names) {

data <- read.table(paste0(file,".DNA_vs_gene.edit.4R"), header = TRUE, sep = "\t")
p1<-ggplot(data, aes(x=editLevel, color=class)) + 
geom_density(linewidth=0.5,alpha=0.3)+theme_classic() +theme(legend.position = "none", axis.line = element_line(color="white"), axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
scale_color_manual(values=c('#ED7D31','#A6A6A6'))

data <- read.table(paste0(file,".LTR_vs_gene.edit.4R"), header = TRUE, sep = "\t")
p2<-ggplot(data, aes(x=editLevel, color=class)) +
geom_density(linewidth=0.5,alpha=0.3)+theme_classic() +theme(legend.position = "none", axis.line = element_line(color="white"),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
scale_color_manual(values=c('#A6A6A6','#70AD47'))


data <- read.table(paste0(file,".LINE_vs_gene.edit.4R"), header = TRUE, sep = "\t")
p3<-ggplot(data, aes(x=editLevel, color=class)) +
geom_density(linewidth=0.5,alpha=0.3)+theme_classic() +theme(legend.position = "none", axis.line = element_line(color="white"),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
scale_color_manual(values=c('#A6A6A6','#FFC000'))

data <- read.table(paste0(file,".SINE_vs_gene.edit.4R"), header = TRUE, sep = "\t")
p4<-ggplot(data, aes(x=editLevel, color=class)) +
geom_density(linewidth=0.5,alpha=0.3)+theme_classic() +theme(legend.position = "none", axis.line = element_line(color="white"),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
scale_color_manual(values=c('#A6A6A6','#0068B3'))

pdf(paste0(file,".TE_vs_gene.density.pdf"),w=8,h=2)
print(ggarrange(p4, p3, p2, p1, ncol = 4, nrow = 1))
dev.off()
}




file_names <- list ("H1-hESC","DE-H1","AFG-H1","PGC-H1","H9-hESC","DE-H9","AFG-H9")

for (file in file_names) {

data <- read.table(paste0(file,".DNA_vs_gene.edit.4R"), header = TRUE, sep = "\t")
p1<-ggplot(data, aes(x=editLevel, color=class)) + 
geom_density(linewidth=0.5,alpha=0.3)+theme_classic() +theme(legend.position = "none",  axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
scale_color_manual(values=c('#ED7D31','#A6A6A6'))

data <- read.table(paste0(file,".LTR_vs_gene.edit.4R"), header = TRUE, sep = "\t")
p2<-ggplot(data, aes(x=editLevel, color=class)) +
geom_density(linewidth=0.5,alpha=0.3)+theme_classic() +theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
scale_color_manual(values=c('#A6A6A6','#70AD47'))


data <- read.table(paste0(file,".LINE_vs_gene.edit.4R"), header = TRUE, sep = "\t")
p3<-ggplot(data, aes(x=editLevel, color=class)) +
geom_density(linewidth=0.5,alpha=0.3)+theme_classic() +theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
scale_color_manual(values=c('#A6A6A6','#FFC000'))

data <- read.table(paste0(file,".SINE_vs_gene.edit.4R"), header = TRUE, sep = "\t")
p4<-ggplot(data, aes(x=editLevel, color=class)) +
geom_density(linewidth=0.5,alpha=0.3)+theme_classic() +theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
scale_color_manual(values=c('#A6A6A6','#0068B3'))

pdf(paste0(file,".TE_vs_gene.density.withlable.pdf"),w=8,h=2)
print(ggarrange(p4, p3, p2, p1, ncol = 4, nrow = 1))
dev.off()
}