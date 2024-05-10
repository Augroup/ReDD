library(ggpubr)
library(dplyr) 
library(ggplot2)
theme_set(theme_pubclean())

# List of your tables, replace 'your_table_list' with the actual list of your tables

file_names <- list("H1-hESC.TE_vs_gene.TEtype.4R", "DE-H1.TE_vs_gene.TEtype.4R", "AFG-H1.TE_vs_gene.TEtype.4R","PGC-H1.TE_vs_gene.TEtype.4R","H9-hESC.TE_vs_gene.TEtype.4R", "DE-H9.TE_vs_gene.TEtype.4R", "AFG-H9.TE_vs_gene.TEtype.4R")

sink("KS_test.results.out", append = TRUE)

for (file_name in file_names){

data<- read.table(file_name,header= T, sep="\t")
gene_data <- data[data$class == "gene", "editLevel"]
ltr_data <- data[data$class == "LTR", "editLevel"]
dna_data <- data[data$class == "DNA", "editLevel"]
line_data <- data[data$class == "LINE", "editLevel"]
sine_data <- data[data$class == "SINE", "editLevel"]

# Perform the Kolmogorov-Smirnov test
ltr_result <- ks.test(gene_data, ltr_data)
dna_result <- ks.test(gene_data, dna_data)
line_result <- ks.test(gene_data, line_data)
sine_result <- ks.test(gene_data, sine_data)

sum_gene <-summary(gene_data)
sum_ltr <-summary(ltr_data)
sum_line <-summary(line_data)
sum_sine <-summary(sine_data)
sum_dna <-summary(dna_data)


print(file_name)
print(dna_result)
print(ltr_result)
print(line_result)
print(sine_result)
cat("Gene:\n")
print (sum_gene)
cat("LTR:\n")
print (sum_ltr)
cat("LINE:\n")
print (sum_line)
cat("SINE:\n")
print (sum_sine)
cat("DNA:\n")
print (sum_dna)
}
sink()



for (file_name in file_names){

dataset<- read.table(file_name,header= T, sep="\t")
data <- rbind (dataset(class = "gene", editLevel=dataset$editLevel),
				data.frame(class=DNA, editLevel=dataset$editLevel))

pdf(paste (file_name,".dist.pdf"),w=6,h=6)
p1<-ggplot(data, aes(x=editLevel, color=class)) +
geom_density(linewidth=1.5,alpha=0.2)+theme_classic()
p1+theme(axis.text.x = element_text(face="bold", color = "black", 
                           size=14, angle=0),
          axis.text.y = element_text(face="bold", color = "black", 
                           size=14, angle=0))
print (p1)
dev.off()
}



