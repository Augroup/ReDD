#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))


option_list <- list(
  make_option(c("-s", "--siteFile"), type="character", default=NULL, help="Sites file name", metavar="character"),
  make_option(c("-b", "--bedFile"), type="character", default=NULL, help="Bed file name", metavar="character"),
  make_option(c("-o", "--outDir"), type="character", default="output", help="output file path [default = %default]", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

######
# 01.load data
######

import_files <- function(file_path, type = c("sites","beds")){
  type <- match.arg(type)
  tab <- data.table::fread(input = file_path, sep="\t", header=FALSE)
  if(type =="sites"){
    colnames(tab) <- c('NGS_editing','read_id','transcript_id','transcript_loc','strand','editing_probability')
  }else if(type =="beds"){
    colnames(tab) <- c('transcript_id','transcript_loc','editing_cov','total_cov','editing_ratio','pred_ratio')
  }
  tab
}

######
# 02.pre-filter data
######

filter_editing_bed <- function(bed_tab,
                               edit_cov_cutoff = 1,
                               total_cov_cutoff = 20,
                               edit_ratio_cutoff = 0.01){
  bed_tab %>% filter(editing_cov >= edit_cov_cutoff & total_cov >= total_cov_cutoff & editing_ratio >= edit_ratio_cutoff) %>% arrange(transcript_id, transcript_loc)
}

filter_editing_site <- function(site_tab, bed_tab, norm_val = FALSE){
  if(norm_val){
    site_tab$editing_probability <- ifelse(site_tab$editing_probability>=0.5, 1, 0)
  }
  semi_join(x = site_tab, y = bed_tab, by = c('transcript_id', 'transcript_loc')) %>% arrange(transcript_id, transcript_loc, read_id)
}

######
# 04.get read~site matrix
######

get_read_site_mat <- function(site_tab, min_reads = 20, min_sites = 2, ncores=1){
  site_by_iso <- split(site_tab, site_tab$transcript_id)
  #--filtering isoforms without enough reads
  is_iso_keep <- purrr::map_lgl(.x = site_by_iso, .f = function(tab, min_read, min_site){
    length(unique(tab$read_id)) >= min_read & length(unique(tab$transcript_loc)) >= min_site
  }, min_read =  min_reads, min_site =  min_sites)
  #--generating matrix
  generate_mat <- function(tab){
    editing_tab <- data.table::dcast(data = tab, formula = read_id~transcript_loc, fill = NA, value.var = "editing_probability")
    editing_mat <- as.matrix(editing_tab[,-1])
    row.names(editing_mat) <- editing_tab$read_id
    Matrix::Matrix(editing_mat)
  }
  #---parallel calculation
  future::plan(multicore, workers = ncores)
  furrr::future_map(.x = site_by_iso[is_iso_keep],.f = generate_mat, .options = furrr::furrr_options())
}

######
# 05.filter read~site matrix
######
# 05.1 remove reads and sites with NAs.
rm_NAs <- function(mat,
                   max_na_site = 0.3,
                   max_na_read = 0.3){
  mat <- as.matrix(mat)
  #---cols
  NA_cols <- apply(X = mat, MARGIN = 2, FUN = function(x){length(x[is.na(x)])})/nrow(mat)
  
  #---rows
  NA_rows <- apply(X = mat, MARGIN = 1, FUN = function(x){length(x[is.na(x)])})/ncol(mat)
  
  if(all(NA_cols <=  max_na_site) & all(NA_rows <= max_na_read)){
    return(mat)
  }else{
    mat <- mat[, NA_cols <= max_na_site, drop=FALSE]
    mat <- mat[NA_rows <= max_na_read,, drop=FALSE]
    if(any(dim(mat)==0)){
      return(mat)
    }else{
      rm_NAs(mat = mat, max_na_site = max_na_site, max_na_read = max_na_read)
    }
  }
}

# 05.2 filter matrix 
mat_filter <- function(mat,
                       prob_cutoff = 0.5,
                       min_site_num = 1,
                       min_read_num = 2,
                       na.rm = TRUE){
  mat <- as.matrix(mat)
  site_num <- apply(X = mat, MARGIN = 1, FUN = function(x, prob_cutoff=0.5){
    x <- na.omit(object = x)
    length(x[x>=prob_cutoff])
  }, prob_cutoff=prob_cutoff)
  read_num <- apply(X = mat, MARGIN = 2, FUN = function(x, prob_cutoff=0.5){
    x <- na.omit(object = x)
    length(x[x>=prob_cutoff])
  }, prob_cutoff=prob_cutoff)
  if(all(site_num >= min_site_num) & all(read_num >= min_read_num)){
    return(mat)
  }else{
    row_keep <- site_num >= min_site_num
    mat <- mat[row_keep,, drop=FALSE]
    col_keep <- read_num>= min_read_num
    mat <- mat[,col_keep, drop=FALSE]
    mat_filter(mat = mat,prob_cutoff = prob_cutoff,
               min_site_num = min_site_num, min_read_num = min_read_num)
  }
}


#' Title
#'
#' @param mat The matrix being filtered.
#' @param editing_probablity_cutoff The cutoff value used to decide editing or not, default: 0.5.
#' @param min_editing_site_per_read  The minimum editing sites on each read. default: 1.
#' @param min_editing_read_per_site The minimum number of reads with editing sites. default: 2.
#' @param max_missing_site_rate The maximum ratio of reads with missing value. default: 0.3.
#' @param max_missing_read_rate The maximum ratio of sites with missing value. default: 0.3.
#' @param normalize_val Normalize the probability.
#' @param rm_outlier_read Remove reads with super editing.
#'
#' @return
#' @export
#'
#' @examples
read_site_qc <- function(mat,
                         editing_probablity_cutoff = 0.5,
                         min_editing_site_per_read = 1,
                         min_editing_read_per_site = 2,
                         max_missing_site_rate = 0.3,
                         max_missing_read_rate = 0.3,
                         rm_outlier_read = TRUE,
                         normalize_val = FALSE){
  mat <- as.matrix(mat)
  #---1. remove reads & sites with large proportion of missing value.
  mat <- rm_NAs(mat = mat, max_na_site = max_missing_site_rate, max_na_read = max_missing_read_rate)
  #---2. remove outlier
  if(rm_outlier_read){
    editing_site_per_read <- apply(X = mat, MARGIN = 1, function(x){
      x <- x[!is.na(x)]
      length(x[x>=0.5])
    })
    row_keep <- editing_site_per_read < quantile(x = editing_site_per_read, probs = 0.99)
    mat <- mat[row_keep,, drop=FALSE]
  }
  #---3. remove reads & sites with unwanted editing site number
  mat <- mat_filter(mat = mat, prob_cutoff = editing_probablity_cutoff, min_site_num =  min_editing_site_per_read,
                    min_read_num = min_editing_read_per_site)
  #---4. normalization
  if(normalize_val){
    mat <- ifelse(mat >= 0.5,1,0)
  }
  mat
}

######
# 06.visualize read~site matrix 
######

# 1. load data and processing

site_files <- import_files(file_path = opt$siteFile, type = 'sites')
bed_files <- import_files(file_path = opt$bedFile, type = 'beds')
#data.table::fwrite(x = bed_files ,file = opt$out,sep = '\t')
bed_files <- filter_editing_bed(bed_tab = bed_files)
site_files <- filter_editing_site(site_tab = site_files, bed_tab = bed_files)
site_mat_obj <- get_read_site_mat(site_tab = site_files)

site_mat_obj <- pbapply::pblapply(X = site_mat_obj, FUN = read_site_qc)
site_mat_obj <- site_mat_obj[sapply(X = site_mat_obj, FUN = nrow) > 1 & sapply(X = site_mat_obj, FUN = ncol) > 1]
dir.create(path = file.path(opt$outDir,"data"),recursive = T,showWarnings = FALSE)
saveRDS(object = site_mat_obj, file = file.path(opt$outDir,"data",'site_matrix_filtered.rds'))

# visualization
fig_out <- file.path(opt$outDir,"figures")
dir.create(path = file.path(fig_out,"heatmap"),recursive = T,showWarnings = FALSE)
dir.create(path = file.path(fig_out,"correlation"),recursive = T,showWarnings = FALSE)

isoform_id <- names(site_mat_obj)
visul_func <- pbapply::pblapply(X = isoform_id, FUN = function(iso, obj, fig_out){
  plot_mat <- as.matrix(obj[[iso]])
  clust_mat <- plot_mat
  clust_mat[is.na(clust_mat)] <- 0
  clust_res <- hclust(d = dist(clust_mat), method = 'ward.D2')
  p <- pheatmap::pheatmap(mat = plot_mat[clust_res$order,],cluster_rows = F, cluster_cols = F, show_rownames = F, 
                          clustering_method = 'ward.D2', fontsize = 6, silent=T)
  ggplot2::ggsave(filename = file.path(fig_out,"heatmap", paste(iso,'.svg',sep='')), plot = p, device = 'svg', width = 3,height = 3, units = 'in')
  # corr
  plot_cor <- clust_mat
  svg(file =  file.path(fig_out,"correlation",paste(iso,'.svg',sep='')))
  corrplot::corrplot(cor(plot_cor),method='color',type='upper',col = rev(COL2('RdBu', 200)),outline = FALSE,tl.cex = 0.5)
  dev.off()
}, obj = site_mat_obj,fig_out=fig_out)

save.image(file = file.path(opt$outDir,"data",'output.RData'))

