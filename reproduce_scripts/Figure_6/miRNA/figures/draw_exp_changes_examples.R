
library(ggplot2)

data <- read.table("ENST00000378714.7.TPM", header = TRUE, sep = "\t")


data_summary <- aggregate(TPM ~ Stage, data,       # Create summary data
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[ , 1], data_summary$TPM)
data_summary


pdf("ENST00000378714.7.exp_hESC_DE.pdf")
ggplot(data, aes(TPM, Stage, fill = Stage)) +      # ggplot2 barplot with error bars
  coord_flip() +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.4)+
  scale_fill_manual(values=c('#C00000','#005493'))+
  theme_classic()
  #scale_x_discrete(limits=c("hESC","DE"))

dev.off()


data <- read.table("hsa-miR-1260b.tpm", header = TRUE, sep = "\t")

pdf("hsa-miR-1260b.exp_hESC_DE.pdf", w=3,h=6)
ggplot(data, aes( Stage,TPM, fill = Stage)) +  
geom_bar(stat="identity")+    # ggplot2 barplot with error bars
  scale_fill_manual(values=c('#C00000','#005493'))+
  scale_x_discrete(limits=c("hESC","DE"))+
  theme_classic()

dev.off()



data <- read.table("ENST00000360768.5.TPM", header = TRUE, sep = "\t")


data_summary <- aggregate(TPM ~ Stage, data,       # Create summary data
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[ , 1], data_summary$TPM)
data_summary


pdf("ENST00000360768.5.exp_hESC_DE.pdf")
ggplot(data, aes(TPM, Stage, fill = Stage)) +      # ggplot2 barplot with error bars
  coord_flip() +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.4)+
  scale_fill_manual(values=c('#C00000','#005493'))+
  theme_classic()
  #scale_x_discrete(limits=c("hESC","DE"))

dev.off()


data <- read.table("hsa-miR-106a-5p.tpm", header = TRUE, sep = "\t")

pdf("hsa-miR-106a-5p.exp_hESC_DE.pdf", w=3,h=6)
ggplot(data, aes( Stage,TPM, fill = Stage)) +  
geom_bar(stat="identity")+    # ggplot2 barplot with error bars
  scale_fill_manual(values=c('#C00000','#005493'))+
  scale_x_discrete(limits=c("hESC","DE"))+
  theme_classic()

dev.off()