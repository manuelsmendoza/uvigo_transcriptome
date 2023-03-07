# SETUP -------------------------------------------------------------------------------------------
options(warn = -1)

# Command line parser accept short and long flag/options
suppressPackageStartupMessages(library("optparse", character.only = TRUE))

# Parse the options
opt_list   <- list(
  make_option(opt_str = c("-r", "--ref"), 
              type    = "character", 
              default = NULL, 
              help    = "Reference sequence (.fasta)"),
  make_option(opt_str = c("-a", "--ann"), 
              type    = "character", 
              default = NULL, 
              help    = "Reference annotation (.gtf)"),
  make_option(opt_str = c("-p", "--spl"), 
              type    = "character", 
              default = NULL, 
              help    = "Samples information"),
  make_option(opt_str = c("-o", "--out"), 
              type    = "character", 
              default = NULL, 
              help    = "Output directory"),
  make_option(opt_str = c("-t", "--tec"), 
              type    = "character", 
              default = NULL, 
              help    = "High throughput sequencing technique"),
  make_option(opt_str = c("-c", "--cpu"), 
              type    = "double", 
              default = 1, 
              help    = "Number of threads to use"))
opt_parser <- OptionParser(option_list = opt_list)
opt        <- parse_args(opt_parser)


if (is.null(opt$ref) || is.null(opt$spl) || is.null(opt$tec)) {
  stop(paste0("[", Sys.time(), "] [FAIL]: Information required missed"))
}

if (is.null(opt$out)) {
  opt$out <- getwd()
}

if (!opt$tec %in% c("tran", "geno")) {
  stop(paste0("[", Sys.time(), "] [FAIL]: Reference type wrong"))
}

# EXPRESSION QUANTIFICATION -----------------------------------------------------------------------
# Command line parser accept short and long flag/options
suppressPackageStartupMessages(library("Rsubread", character.only = TRUE))
suppressPackageStartupMessages(library("readr",    character.only = TRUE))
suppressPackageStartupMessages(library("stringr",  character.only = TRUE))

# Build reference index
ref_index   <- str_replace(string = opt$ref, pattern = ".fasta$", replacement = "")
buildindex(
  reference = opt$ref, 
  basename  = ref_index
)

# Align the reads to the reference
samples_info <- read_delim(file = opt$spl, delim = "\t", col_names = FALSE)
if (ncol(samples_info) == 3 ) {
  for (I in 1:nrow(samples_info)) {
    align(
      index                  = ref_index, 
      readfile1              = samples_info[I, 3], 
      type                   = "rna", 
      output_file            = file.path(opt$out, paste0(samples_info[I, 1], ".bam")), 
      nthreads               = as.numeric(opt$cpu), 
      sortReadsByCoordinates = TRUE, 
      useAnnotation          = TRUE, 
      annot.ext              = opt$ann, 
      isGTF                  = TRUE, 
      GTF.attrType           = "gene_name")
  }
} else if (ncol(samples_info) == 4) {
  for (I in 1:nrow(samples_info)) {
    align(
      index                  = ref_index, 
      readfile1              = samples_info[I, 3], 
      readfile2              = samples_info[I, 4], 
      type                   = "rna", 
      output_file            = file.path(opt$out, paste0(samples_info[I, 1], ".bam")), 
      nthreads               = as.numeric(opt$cpu), 
      sortReadsByCoordinates = TRUE, 
      useAnnotation          = TRUE, 
      annot.ext              = opt$ann, 
      isGTF                  = TRUE, 
      GTF.attrType           = "gene_name")
  }
}

# Quantify gene expression
paired_seq <- ncol(samples_info) == 4

gene_expression <- featureCounts(
  files               = list.files(path = opt$out, pattern = ".bam$", full.names = TRUE), 
  genome              = opt$ref , 
  annot.ext           = opt$ann, 
  isGTFAnnotationFile = TRUE, 
  GTF.attrType        = "gene_name", 
  useMetaFeatures     = TRUE, 
  isLongRead          = FALSE, 
  isPairedEnd         = paired_seq, 
  minMQS              = 20, 
  nthreads            = as.numeric(opt$cpu))
gene_counts <- gene_expression$counts
colnames(gene_counts) <- str_replace(string = colnames(gene_counts), pattern = ".bam$", replacement = "")

# Convert raw counts to TPM
gene_length <- gene_expression$annotation[, c("GeneID", "Length")]

count_tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
gene_tpms <- apply(X = gene_counts, MARGIN = 2, FUN = function(x){count_tpm(x, gene_length$Length)})


# Export gene expression matrix and remove tmp files
write_delim(x = gene_counts, file = file.path(opt$out, "gene_expression.cts.tsv"), delim = "\t", append = FALSE, col_names = FALSE)
write_delim(x = gene_tpms,   file = file.path(opt$out, "gene_expression.tpm.tsv"), delim = "\t", append = FALSE, col_names = FALSE)

tmp_files <- list.files(path = opt$out, pattern = ".bam", full.names = TRUE)
for (TF in tmp_files) {
  out <- file.remove(TF)
}
