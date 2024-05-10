library(ggplot2)
library(ggpubr)

pdf ("stemcell.top500.enrichment.pdf",h=12,w=12)
mydata <- read.table("stemcell.GO",header = T, sep = "\t")

p<- ggplot(mydata, aes(x = Type, y = Term,
                      size = Fold_Enrichment,
                      color = FDR_log))+
geom_point(alpha = 0.9)+theme_light()+
scale_size_binned(name = "Fold Enrichment")+
#scale_size(breaks = waiver(),guide = "legend", name = "log2(Enrichment)")+
#scale_color_manual(values=c("#00afbb","#e7b800","#fc4e07"))+
scale_colour_gradient(low="#91D1C2FF",high = "#E64B35FF", name="-log10(FDR)")

p+theme(axis.text.x = element_text(color = "black",
                           size=5, angle=0),
          axis.text.y = element_text(color = "black",  #face="bold",
                           size=14, angle=0))

dev.off()