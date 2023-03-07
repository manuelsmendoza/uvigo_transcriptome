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

# BUILD SEQUENCE INDEX ----------------------------------------------------------------------------
# Command line parser accept short and long flag/options
suppressPackageStartupMessages(library("Rsubread", character.only = TRUE))
suppressPackageStartupMessages(library("stringr",  character.only = TRUE))

# Build reference index
ref_index   <- str_replace(string = opt$ref, pattern = ".fasta$", replacement = "")
buildindex(
  reference = opt$ref, 
  basename  = ref_index
)

samples_info <- read_delim(file = opt$dpl, delim = "\t", col_names = FALSE)
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
