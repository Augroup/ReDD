#install.packages("UpSetR")

library(UpSetR)

pdf ("edit_isoform_cellline.upsetplot.pdf", w=14,h=7)
# Dataset
input <- c(

"H1-hESC" = 39691,
"H1-DE" = 16161,
"H1-PE" = 14914,
"H1-PGC" = 17083,
"H9-hESC" = 24472,
"H9-DE" = 32528,
"H9-PE" = 26781,
"H1-hESC&H1-DE" = 15246,
"H1-DE&H1-PE" = 11812,
"H1-PE&H1-PGC" = 12138,
"H1-hESC&H1-DE&H1-PE" = 11604,
"H1-PGC&H1-hESC&H1-DE&H1-PE" = 10437,
"H9-hESC&H9-DE" = 21812,
"H9-DE&H9-PE" = 23683,
"H9-hESC&H9-PE" = 20451,
"H9-hESC&H9-DE&H9-PE" = 19408,
"H1-PGC&H1-hESC&H1-DE&H1-PE&H9-hESC&H9-DE&H9-PE" = 10068

)


# Plot
upset(fromExpression(input), sets = c("H1-PGC", "H1-hESC", "H1-DE","H1-PE","H9-hESC", "H9-DE","H9-PE"), sets.bar.color = "#56B4E9",
      nintersects = 20, 
      nsets = 9, 
      order.by = "freq", keep.order=TRUE,
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 2.0, 
      point.size = 2.8, 
      line.size = 1
      )

dev.off()