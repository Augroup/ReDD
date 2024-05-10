library(ggplot2)
#library(ggpubr)

pdf ("TEfamily.enrichment_analysis.pdf",h=8,w=8)
mydata <- read.table("enrichment_analysis.4R",header = T, sep = "\t")

p<- ggplot(mydata, aes(x = cell, y = family,
                      size = count,
                      color = pvalue))+
geom_point(alpha = 0.9)+theme_light()+
#scale_size_binned(name = "Fold Enrichment")+

scale_size_binned_area(
    limits = c(0, 4000),
    breaks = c(0, 100, 500, 1000, 4000),
  )+
#scale_size(breaks = waiver(),guide = "legend", name = "log2(Enrichment)")+
#scale_color_manual(values=c("#00afbb","#e7b800","#fc4e07"))+
#scale_colour_gradient(low="#91D1C2FF",high = "#E64B35FF", name="p-value")
scale_colour_gradient(low ="red", high="blue",space="Lab")

p+theme(axis.text.x = element_text(color = "black",
                           size=14, angle=0),
          axis.text.y = element_text(color = "black",  #face="bold",
                           size=14, angle=0))

dev.off()