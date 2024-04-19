#!/usr/bin/env /programs/R-4.2.1-r9/bin/Rscript

def_colnames <- function(dataframe, name, suffixes){
  base_suffixes <- c("chr","start","end")
  dynamic_names <- sapply(suffixes, function(suffix) paste0(name, "_", suffix))
  new_colnames <- c(base_suffixes, dynamic_names)
  
  colnames(dataframe) = new_colnames
  return(dataframe)
}

DESeq_with_ctrl <- function(all, ctrl, outdir, model, num_rep_DNA, num_rep_RNA){
  all <- all[, -(1:3)]
  filter_full <- all[1:(nrow(all)-nrow(ctrl)),]
  filter_full <- which(rowSums(filter_full)>10)
  all_fil <- all[which(rowSums(all)>10),]
  rm(all)
  
  cts <- all_fil
  rownames(cts) <- rownames(all_fil)
  controlgene = seq(length(filter_full)+1, length=nrow(cts)-length(filter_full))
  
  coldata <- data.frame(colnames(cts))
  coldata$condition <- c(rep("DNA", num_rep_DNA),rep("RNA", num_rep_RNA))
  
  rownames(coldata) <- coldata$colnames.cts.
  coldata[,1] <- NULL
  coldata$condition <- factor(coldata$condition)
  
  ## reset cts row index
  #row.names(cts) <- 1:nrow(cts)
  
  all(rownames(coldata) %in% colnames(cts))
  all(rownames(coldata) == colnames(cts))
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)
  dds
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds <- estimateSizeFactors(dds, controlGenes=controlgene)#, type="poscounts")
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  resultsNames(dds)
  res <- results(dds, contrast=c("condition", "RNA", "DNA"))
  
  resLFC <- lfcShrink(dds, coef="condition_RNA_vs_DNA", type="apeglm")
  resLFC
  
  output <- resLFC[which((resLFC[1:length(filter_full),]$padj<0.05) == TRUE),]
  output <- resLFC[which((resLFC[1:length(filter_full),]$log2FoldChange>=1) == TRUE),]
  # print(output)
  
  # Define the directory path
  deseq_directory <- paste0(outdir, "/DESeq")
  
  # Check if the directory exists
  if (!file.exists(deseq_directory)) {
    # Create the directory
    dir.create(deseq_directory, recursive = TRUE)
  }
  
  write.table(output, file=paste0(outdir, "/DESeq/DE_results_", model, ".txt"),
              append = FALSE, sep = "\t", 
              row.names = TRUE, col.names = TRUE)
  message("DESeq results saved.")
}


save_enhancers <- function(deseq, outdir, design) {
  ### find the row index for full_f elements
  f_idx <- length(which(substring(rownames(deseq), 1, nchar(design)+2)== paste0(design,"_f")))

  ### save the elements for each subgroup
  full_f = substring(rownames(deseq)[1:f_idx], first=nchar(design)+3)
  full_r = substring(rownames(deseq)[f_idx+1:nrow(deseq)], first=nchar(design)+3)
  
  common = intersect(full_f, full_r)
  # print(common)
  
  ### need to find the original element
  ori_f <- read.table(paste0(outdir, "/", design, "/srt_", design, "_f_pairwise.bed"))
  ori_r <- read.table(paste0(outdir, "/", design, "/srt_", design, "_r_pairwise.bed"))
  
  rownames(ori_f) <- paste0("Full_f", 1:nrow(ori_f))
  rownames(ori_r) <- paste0("Full_r", 1:nrow(ori_r))
  
  names_f <- paste0("Full_f", common)
  names_r <- paste0("Full_r", common)
  
  # Filter dataframe based on the numbers in names
  filtered_rows_f <- ori_f[rownames(ori_f) %in% names_f, ]
  filtered_rows_r <- ori_r[rownames(ori_r) %in% names_r, ]
  
  write.table(filtered_rows_f, file=paste0(outdir, "/",design,"/srt_",design,"_f_e.bed"),
              append = FALSE, sep = "\t", 
              row.names = FALSE, col.names = FALSE)
  
  write.table(filtered_rows_r, file=paste0(outdir, "/",design,"/srt_",design,"_r_e.bed"),
              append = FALSE, sep = "\t", 
              row.names = FALSE, col.names = FALSE)
  
  if (nrow(filtered_rows_f)>=1) {
    message("Enhancer elements successfully identified and saved.")
  }
  else {
    message("No enhancer elements captured.")
  }
}

################################################################################
suppressPackageStartupMessages({
  library(argparse)
  library(DESeq2)
  library(dplyr)
});

# create parser object and add parser arguments
parser <- ArgumentParser()

parser$add_argument("-o", "--output", required=T, help="Output file directory")
parser$add_argument("--num_rep_DNA", type="integer", default=6, help="number of replicate")
parser$add_argument("--num_rep_RNA", type="integer", default=7, help="number of replicate")
parser$add_argument("--name", nargs="+", required=T, help="design of element, e.g. 5p; 3p; TSS_b")

args <- parser$parse_args()

################################################################################
if (!file.exists(paste0(args$output, "/full/srt_full_f_e.bed"))){
  message("Started searching for enhancers ... ")
  # Output full enhancer elements with orientation-independence
  full_f <- read.table(paste0(args$output,"/full/srt_full_f.bed"))
  full_r <- read.table(paste0(args$output,"/full/srt_full_r.bed"))
  
  ### read in control gene set
  control_gene <- read.table("/fs/cbsuhy01/storage/yz2676/data/STARR-seq/normalization/srt_deep_ATAC_exon_ctrl.bed")
  
  rownames(full_f) <- paste0("Full_f", 1:nrow(full_f))
  rownames(full_r) <- paste0("Full_r", 1:nrow(full_r))
  rownames(control_gene) <- paste0("ctrl", 1:nrow(control_gene))
  control_gene <- control_gene[, -(4)]
  colnames(control_gene) <- colnames(full_f)
  all <- bind_rows(full_f, full_r, control_gene)
  
  DESeq_with_ctrl(all, control_gene, args$output, "full", args$num_rep_DNA, args$num_rep_RNA)
  
  full <- read.table(paste0(args$output, "/DESeq/DE_results_full.txt"))
  save_enhancers(full, args$output)
}

# Also calculate activities for partial elements
for (n in args$name) {
  message(paste0("Calculating activities for ", n))
  partial_f <- read.table(paste0(args$output,"/", n, "/srt_", n, "_f_pairwise.bed"))
  partial_r <- read.table(paste0(args$output,"/", n, "/srt_", n, "_r_pairwise.bed"))

  ### read in control gene set
  control_gene <- read.table("/fs/cbsuhy01/storage/yz2676/data/STARR-seq/normalization/srt_deep_ATAC_exon_ctrl.bed")

  rownames(partial_f) <- paste0(n, "_f", 1:nrow(partial_f))
  rownames(partial_r) <- paste0(n, "_r", 1:nrow(partial_r))
  rownames(control_gene) <- paste0("ctrl", 1:nrow(control_gene))
  control_gene <- control_gene[, -(4)]
  colnames(control_gene) <- colnames(partial_f)
  all <- bind_rows(partial_f, partial_r, control_gene)

  DESeq_with_ctrl(all, control_gene, args$output, n, args$num_rep_DNA, args$num_rep_RNA)
  full <- read.table(paste0(args$output, "/DESeq/DE_results_", n, ".txt"))
  save_enhancers(full, args$output, n)
}