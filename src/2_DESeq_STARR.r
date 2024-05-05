#!/usr/bin/env /programs/R-4.2.1-r9/bin/Rscript

# Sample command to run the script: 
# /programs/R-4.2.1-r9/bin/Rscript 2_DESeq_STARR.r -o /path_to_file/EnhancerNet --name 5p, 3p
# /programs/R-4.2.1-r9/bin/Rscript 2_DESeq_STARR.r -o /path_to_file/new_data --name TSS_b TSS_p TSS_n

def_colnames <- function(dataframe, name, suffixes){
  base_suffixes <- c("chr","start","end")
  dynamic_names <- sapply(suffixes, function(suffix) paste0(name, "_", suffix))
  new_colnames <- c(base_suffixes, dynamic_names)
  
  colnames(dataframe) = new_colnames
  return(dataframe)
}


# Function to bind each pair with control_gene
bind_each_pair <- function(name_data) {
  # Binding the forward, reverse data frames
  bind_rows(name_data$forward, name_data$reverse)
}


DESeq_regular <- function(ctrl, num_rep_DNA, num_rep_RNA, output_file){
  ctrl <- ctrl[which(rowSums(ctrl)>10),]
  
  cts <- ctrl
  rownames(cts) <- rownames(ctrl)
  
  coldata <- data.frame(colnames(cts))
  coldata$condition <- c(rep("DNA", num_rep_DNA),rep("RNA", num_rep_RNA))
  
  rownames(coldata) <- coldata$colnames.cts.
  coldata[,1] <- NULL
  coldata$condition <- factor(coldata$condition)
  
  all(rownames(coldata) %in% colnames(cts))
  all(rownames(coldata) == colnames(cts))
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("condition", "RNA", "DNA"))
  # resLFC <- lfcShrink(dds, coef="condition_RNA_vs_DNA", type="apeglm")
  
  write.table(res, file=output_file,
              append = FALSE, sep = "\t",
              row.names = TRUE, col.names = TRUE)
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
  # dds
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds <- estimateSizeFactors(dds, controlGenes=controlgene)#, type="poscounts")
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  resultsNames(dds)
  res <- results(dds, contrast=c("condition", "RNA", "DNA"))
  
  resLFC <- lfcShrink(dds, coef="condition_RNA_vs_DNA", type="apeglm")
  # resLFC
  
  if (model == "full"){
    ## first option: LFC>=1 & padj<0.05
    output <- resLFC[1:length(filter_full), ][resLFC[1:length(filter_full), ]$padj < 0.05 & 
                                                resLFC[1:length(filter_full), ]$log2FoldChange >= 1, ]
    ## second option: z-score > 1.64 - 95th percentile for a one-sided Gaussian distribution
    # mean_neg_ctrl <- mean(resLFC[length(filter_full)+1:length(controlgene),]$log2FoldChange)
    # 
    # std_neg_ctrl = sd(resLFC[length(filter_full)+1:length(controlgene),]$log2FoldChange)
    # 
    # resLFC$"z-score" = (resLFC$log2FoldChange-mean_neg_ctrl)/std_neg_ctrl
    # logFC_cutoff = 1.64*std_neg_ctrl+mean_neg_ctrl
    # output <- resLFC[1:length(filter_full), ][resLFC[1:length(filter_full), ]$padj < 0.05 & 
    #                                             resLFC[1:length(filter_full), ]$log2FoldChange >= logFC_cutoff, ]
  }
  else {
    output <- resLFC[1:length(filter_full),]
  }
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
  ori_f <- read.table(paste0(outdir, "/", design, "/srt_", design, "_f.bed"))
  ori_r <- read.table(paste0(outdir, "/", design, "/srt_", design, "_r.bed"))
  
  rownames(ori_f) <- paste0(design, "_f", 1:nrow(ori_f))
  rownames(ori_r) <- paste0(design, "_r", 1:nrow(ori_r))
  
  names_f <- paste0(design, "_f", common)
  names_r <- paste0(design, "_r", common)
  
  # Filter dataframe based on the numbers in names
  filtered_rows_f <- ori_f[rownames(ori_f) %in% names_f, ]
  filtered_rows_r <- ori_r[rownames(ori_r) %in% names_r, ]
  
  # re-indexing the dataframes so that they can be added together
  rownames(filtered_rows_r) <- rownames(filtered_rows_f)
  
  result_values <- filtered_rows_f[,4:ncol(filtered_rows_f)] + filtered_rows_r[,4:ncol(filtered_rows_r)]
  result <- cbind(filtered_rows_f[,1:3], result_values)
  
  write.table(result, file=paste0(outdir, "/",design,"/srt_",design,"_e.bed"),
              append = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  
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
# Validate that the negative ctrl is normally distributed
if (!file.exists("/fs/cbsuhy01/storage/yz2676/data/STARR-seq/normalization/DESeq_result_ctrl.txt")){
  message("Testing negative controls ... ")
  ctrl <- read.table("/fs/cbsuhy01/storage/yz2676/data/STARR-seq/normalization/srt_deep_ATAC_exon_ctrl.bed")
  ctrl <- ctrl[, -(1:4)]
  output_file_path <- "/fs/cbsuhy01/storage/yz2676/data/STARR-seq/normalization/DESeq_result_ctrl.txt"
  DESeq_regular(ctrl, args$num_rep_DNA, args$num_rep_RNA, output_file_path)
  message("Negative ctrl results saved.")
}


# Output full enhancer elements with orientation-independence
if (!file.exists(paste0(args$output, "/full/srt_full_e.bed"))){
  message("Started searching for enhancers ... ")
  full_f <- read.table(paste0(args$output,"/full/srt_full_f.bed"))
  full_r <- read.table(paste0(args$output,"/full/srt_full_r.bed"))

  ### read in control gene set
  control_gene <- read.table("/fs/cbsuhy01/storage/yz2676/data/STARR-seq/normalization/srt_deep_ATAC_exon_ctrl.bed")

  rownames(full_f) <- paste0("full_f", 1:nrow(full_f))
  rownames(full_r) <- paste0("full_r", 1:nrow(full_r))
  rownames(control_gene) <- paste0("ctrl", 1:nrow(control_gene))
  control_gene <- control_gene[, -(4)]
  colnames(control_gene) <- colnames(full_f)
  all <- bind_rows(full_f, full_r, control_gene)

  DESeq_with_ctrl(all, control_gene, args$output, "full", args$num_rep_DNA, args$num_rep_RNA)
  
  full <- read.table(paste0(args$output, "/DESeq/DE_results_full.txt"))
  save_enhancers(full, args$output, "full")
}

# Also calculate activities for partial elements
partial_list <- list()
for (n in args$name) {
  message(paste0("Calculating activities for ", n))
  partial_f <- read.table(paste0(args$output,"/", n, "/srt_", n, "_f.bed"))
  partial_r <- read.table(paste0(args$output,"/", n, "/srt_", n, "_r.bed"))

  rownames(partial_f) <- paste0(n, "_f", 1:nrow(partial_f))
  rownames(partial_r) <- paste0(n, "_r", 1:nrow(partial_r))

  partial_list[[n]] <-list(forward=partial_f, reverse=partial_r)
}
control_gene <- read.table("/fs/cbsuhy01/storage/yz2676/data/STARR-seq/normalization/srt_deep_ATAC_exon_ctrl.bed")
rownames(control_gene) <- paste0("ctrl", 1:nrow(control_gene))
control_gene <- control_gene[, -(4)]
colnames(control_gene) <- colnames(partial_list[[n]]$forward)

# Combine full, partial, control as a single matrix
full <- read.table(paste0(args$output, "/full/srt_full_e.bed"))
rownames(full) <- paste0("full", 1:nrow(full))

all <- bind_rows(lapply(partial_list, bind_each_pair))
all <- bind_rows(full, all, control_gene)

DESeq_with_ctrl(all, control_gene, args$output, substr(n,1,nchar(n)-2), args$num_rep_DNA, args$num_rep_RNA)