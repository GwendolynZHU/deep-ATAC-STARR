#!/usr/bin/env /programs/R-4.2.1-r9/bin/Rscript

# Sample command to run the script: 
# /programs/R-4.2.1-r9/bin/Rscript 2_DESeq_STARR.r -o /path_to_file/EnhancerNet --name 5p, 3p
# /programs/R-4.2.1-r9/bin/Rscript 2_DESeq_STARR.r -o /path_to_file/new_data --name TSS_b TSS_p TSS_n

write_params <- function(args, file) {
  conn <- file(file, "w")
  on.exit(close(conn))
  
  lapply(names(args), function(arg) {
    write(paste(arg, args[[arg]]), conn)
  })
}


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


DESeq_regular <- function(ctrl, DNA_rep, RNA_rep, output_file){
  ctrl <- ctrl[which(rowSums(ctrl)>10),]
  
  cts <- ctrl
  rownames(cts) <- rownames(ctrl)
  
  coldata <- data.frame(colnames(cts))
  coldata$condition <- c(rep("DNA", length(DNA_rep)),rep("RNA", length(RNA_rep)))
  
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


DESeq_with_ctrl <- function(all, ctrl, outdir, model, DNA_rep, RNA_rep, cutoff, cpm_filter){
  all <- all[, -(1:3)]
  ori_DNA_columns <- paste0("DNA_", c(1,2,3,4,5,6))
  ori_RNA_columns <- paste0("RNA_", c(1,2,3,4,5,6,7))
  colnames(all) <- c(ori_DNA_columns, ori_RNA_columns)
  colnames(ctrl) <- c(ori_DNA_columns, ori_RNA_columns)
  

  ### add in flexibility to allow a subset of replicates
  DNA_columns <- paste0("DNA_", DNA_rep)
  RNA_columns <- paste0("RNA_", RNA_rep)
  all <- select(all, all_of(DNA_columns), all_of(RNA_columns))
  ctrl <- select(ctrl, all_of(DNA_columns), all_of(RNA_columns))
  full <- all[1:(nrow(all)-nrow(ctrl)),]
  
  cts <- all
  rownames(cts) <- rownames(all)
  
  coldata <- data.frame(colnames(cts))
  coldata$condition <- c(rep("DNA", length(DNA_rep)),rep("RNA", length(RNA_rep)))
  
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

  if (cpm_filter) { # at least 4 DNA libs >= cpm cutoff
    cpm = cpm(counts(dds), log=FALSE)
    lib_size = mean(colSums(counts(dds))[1:length(DNA_rep)])
    message(paste("mean_lib_size", lib_size))
    cpm_cutoff = cutoff/(lib_size/10^6)
    message(paste("CPM_filter_cutoff", cpm_cutoff))
    if(length(DNA_rep) == 1){
        keep = cpm[, DNA_columns] >= cpm_cutoff
    } else {
        keep = rowSums(cpm[, DNA_columns] >= cpm_cutoff) >= 4
    }
  } else { #raw counts filter
    keep <- rowSums(counts(dds)) >= cutoff
  }
  
  dds_filtered <- dds[keep,]
  control_genes <- rownames(dds_filtered) %in% rownames(ctrl)
  control_indexes <- which(control_genes)

  filtered_full <- rownames(dds_filtered) %in% rownames(full)
  filtered_full_indexes <- length(which(filtered_full))

  dds_filtered <- estimateSizeFactors(dds_filtered, controlGenes=control_indexes)#, type="poscounts")
  dds_filtered <- estimateDispersions(dds_filtered)
  dds_filtered <- nbinomWaldTest(dds_filtered)

  png(filename = paste0(outdir,"/DESeq/plots/", model, "_", cutoff, "_dispersion_plot.png"))
  plotDispEsts(dds_filtered)
  dev.off()

#   resultsNames(dds_filtered)
  res <- results(dds_filtered, contrast=c("condition", "RNA", "DNA"))
  
  resLFC <- lfcShrink(dds_filtered, coef="condition_RNA_vs_DNA", type="apeglm")
  # resLFC
  
  if (model == "full"){
    ## first option: LFC>=1 & padj<0.1
    output <- resLFC[1:filtered_full_indexes, ][resLFC[1:filtered_full_indexes, ]$padj < 0.1 & 
                                                resLFC[1:filtered_full_indexes, ]$log2FoldChange >= 1, ]
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
    output <- resLFC[1:filtered_full_indexes,]
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
  ori_f <- read.table(paste0(outdir, "/data/", design, "/srt_", design, "_f.bed"))
  ori_r <- read.table(paste0(outdir, "/data/", design, "/srt_", design, "_r.bed"))
  
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
  
  write.table(result, file=paste0(outdir, "/data/",design,"/srt_",design,"_e.bed"),
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
  library(edgeR)
});

# create parser object and add parser arguments
parser <- ArgumentParser()

parser$add_argument("-o", "--output", required=T, help="Output file directory")
parser$add_argument("--DNA_rep", type="integer", nargs="+", default=c(1,2,3,4,5,6), help="number of replicate")
parser$add_argument("--RNA_rep", type="integer", nargs="+", default=c(1,2,3,4,5,6,7), help="number of replicate")
parser$add_argument("-c", "--cutoff", type="integer", default=10)
parser$add_argument("--cpm_filter", type="logical", default=FALSE)
parser$add_argument("--name", nargs="+", required=T, help="design of element, e.g. 5p; 3p; TSS_b")

args <- parser$parse_args()

################################################################################
# Save parameters:
write_params(args, paste0(args$output, "/DESeq/params.txt"))

# Validate that the negative ctrl is normally distributed
if (!file.exists("/fs/cbsuhy01/storage/yz2676/data/STARR-seq/normalization/DESeq_result_ctrl.txt")){
  message("Testing negative controls ... ")
  ctrl <- read.table("/fs/cbsuhy01/storage/yz2676/data/STARR-seq/normalization/srt_deep_ATAC_exon_ctrl.bed")
  ctrl <- ctrl[, -(1:4)]
  output_file_path <- "/fs/cbsuhy01/storage/yz2676/data/STARR-seq/normalization/DESeq_result_ctrl.txt"
  DESeq_regular(ctrl, args$DNA_rep, args$RNA_rep, output_file_path)
  message("Negative ctrl results saved.")
}


# Output full enhancer elements with orientation-independence
if (!file.exists(paste0(args$output, "/data/full/srt_full_e.bed"))){
  message("Started searching for enhancers ... ")
  full_f <- read.table(paste0(args$output,"/data/full/srt_full_f.bed"))
  full_r <- read.table(paste0(args$output,"/data/full/srt_full_r.bed"))

  ### read in control gene set
  control_gene <- read.table("/fs/cbsuhy01/storage/yz2676/data/STARR-seq/normalization/srt_deep_ATAC_exon_ctrl.bed")

  rownames(full_f) <- paste0("full_f", 1:nrow(full_f))
  rownames(full_r) <- paste0("full_r", 1:nrow(full_r))
  rownames(control_gene) <- paste0("ctrl", 1:nrow(control_gene))
  control_gene <- control_gene[, -(4)]
  colnames(control_gene) <- colnames(full_f)
  all <- bind_rows(full_f, full_r, control_gene)

  DESeq_with_ctrl(all, control_gene, args$output, "full", args$DNA_rep, args$RNA_rep, args$cutoff, args$cpm_filter)
  
  full <- read.table(paste0(args$output, "/DESeq/DE_results_full.txt"))
#   save_enhancers(full, args$output, "full")
}

# Also calculate activities for partial elements
partial_list <- list()
for (n in args$name) {
  message(paste0("Calculating activities for ", n))
  partial <- read.table(paste0(args$output,"/data/", n, "/srt_", n, ".bed"))

  rownames(partial) <- paste0(n, 1:nrow(partial))

  partial_list[[n]] <- partial
}
control_gene <- read.table("/fs/cbsuhy01/storage/yz2676/data/STARR-seq/normalization/srt_deep_ATAC_exon_ctrl.bed")
rownames(control_gene) <- paste0("ctrl", 1:nrow(control_gene))
control_gene <- control_gene[, -(4)]
colnames(control_gene) <- colnames(partial_list[[n]])

# Combine full, partial, control as a single matrix
full <- read.table(paste0(args$output, "/data/full/srt_full_e.bed"))
rownames(full) <- paste0("full", 1:nrow(full))

all <- bind_rows(c(list(full), partial_list, list(control_gene)))

DESeq_with_ctrl(all, control_gene, args$output, substr(n,1,nchar(n)-2), args$DNA_rep, args$RNA_rep, args$cutoff, args$cpm_filter)