library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(gplots)

paletteLength <- 50
myColor2 <- colorRampPalette(c("grey","white", "#ED7D31"))(paletteLength)

pdf("H1hESC_DE.editSwitch.Heatmap.pdf",w=4,h=4)
mydata <- read.table("H1hESC_DE.editSwitch.txt.4R", header=F, sep="\t",row.names=1)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, cellwidth = 20, cellheight = 15)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, border_color = "NA")
pheatmap(mydata, scale = "row", color = myColor2, cluster_rows=FALSE,cluster_cols=FALSE, border_color = "NA",show_rownames = F)
dev.off()

pdf("H1DE_PE.editSwitch.Heatmap.pdf",w=4,h=4)
mydata <- read.table("H1DE_PE.editSwitch.txt.4R", header=F, sep="\t",row.names=1)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, cellwidth = 20, cellheight = 15)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, border_color = "NA")
pheatmap(mydata, scale = "row", color = myColor2, cluster_rows=FALSE,cluster_cols=FALSE, border_color = "NA",show_rownames = F)
dev.off()

pdf("H9hESC_DE.editSwitch.Heatmap.pdf",w=4,h=4)
mydata <- read.table("H9hESC_DE.editSwitch.txt.4R", header=F, sep="\t",row.names=1)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, cellwidth = 20, cellheight = 15)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, border_color = "NA")
pheatmap(mydata, scale = "row", color = myColor2, cluster_rows=FALSE,cluster_cols=FALSE, border_color = "NA",show_rownames = F)
dev.off()

pdf("H9DE_PE.editSwitch.Heatmap.pdf",w=4,h=4)
mydata <- read.table("H9DE_PE.editSwitch.txt.4R", header=F, sep="\t",row.names=1)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, cellwidth = 20, cellheight = 15)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, border_color = "NA")
pheatmap(mydata, scale = "row", color = myColor2, cluster_rows=FALSE,cluster_cols=FALSE, border_color = "NA",show_rownames = F)
dev.off()







pdf("H1hESC_DE.editSwitch.isoform1.Heatmap.pdf",w=4,h=4)
mydata <- read.table("H1hESC_DE.editSwitch.txt.isoform1.4R", header=F, sep="\t",row.names=1)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, cellwidth = 20, cellheight = 15)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, border_color = "NA")
pheatmap(mydata, scale = "row", color = myColor2, cluster_rows=FALSE,cluster_cols=FALSE, border_color = "NA",show_rownames = F)
dev.off()

pdf("H1hESC_DE.editSwitch.isoform2.Heatmap.pdf",w=4,h=4)
mydata <- read.table("H1hESC_DE.editSwitch.txt.isoform2.4R", header=F, sep="\t",row.names=1)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, cellwidth = 20, cellheight = 15)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, border_color = "NA")
pheatmap(mydata, scale = "row", color = myColor2, cluster_rows=FALSE,cluster_cols=FALSE, border_color = "NA",show_rownames = F)
dev.off()




pdf("H1DE_PE.editSwitch.isoform1.Heatmap.pdf",w=4,h=4)
mydata <- read.table("H1DE_PE.editSwitch.txt.isoform1.4R", header=F, sep="\t",row.names=1)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, cellwidth = 20, cellheight = 15)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, border_color = "NA")
pheatmap(mydata, scale = "row", color = myColor2, cluster_rows=FALSE,cluster_cols=FALSE, border_color = "NA",show_rownames = F)
dev.off()

pdf("H1DE_PE.editSwitch.isoform2.Heatmap.pdf",w=4,h=4)
mydata <- read.table("H1DE_PE.editSwitch.txt.isoform2.4R", header=F, sep="\t",row.names=1)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, cellwidth = 20, cellheight = 15)
#pheatmap(mydata, scale = "row", color = myColor2, cluster_cols=FALSE, border_color = "NA")
pheatmap(mydata, scale = "row", color = myColor2, cluster_rows=FALSE,cluster_cols=FALSE, border_color = "NA",show_rownames = F)
dev.off()