# SETUP -------------------------------------------------------------------------------------------
options(warn = -1)

# Command line parser accept short and long flag/options
suppressPackageStartupMessages(library("optparse", character.only = TRUE))

# Parse the options
opt_list   <- list(
  make_option(opt_str = c("-s", "--seq"), 
              type    = "character", 
              default = NULL, 
              help    = "Sequence to download"),
  make_option(opt_str = c("-o", "--out"), 
              type    = "character", 
              default = NULL, 
              help    = "Output file in FASTA fomat (.fasta)"),
  make_option(opt_str = c("-a", "--ann"), 
              type    = "logical", 
              default = FALSE, 
              help    = "Download the annotation of the sequence"),
  make_option(opt_str = c("-p", "--cds"),
              type    = "logical",
              default = FALSE, 
              help    = "Extract protein coding sequences"))
opt_parser <- OptionParser(option_list = opt_list)
opt        <- parse_args(opt_parser)

chr_names  <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY")
if (is.null(opt$seq) || !opt$seq %in% c("genome", chr_names)) {
  print_help(opt_parser)
  stop(paste0("[", Sys.time(), "] [FAIL]: Sequence name to downlaod is missed or not allowed"))
}

if (!opt$ann && opt$cds) {
  stop(paste0("[", Sys.time(), "] [FAIL]: Sequence annotation required to extract protein-coding sequences"))
}

if (is.null(opt$out)) {
  opt$out <- file.path(getwd(), paste0(opt$seq, ".fasta"))
}

# Define the chromosomes accession
hsa_accessions <- c(
  "NC_000001", "NC_000002", "NC_000003", "NC_000004", "NC_000005", "NC_000006",
  "NC_000007", "NC_000008", "NC_000009", "NC_000010", "NC_000011", "NC_000012", 
  "NC_000013", "NC_000014", "NC_000015", "NC_000016", "NC_000017", "NC_000018", 
  "NC_000019", "NC_000020", "NC_000021", "NC_000022", "NC_000023", "NC_000024")
names(hsa_accessions) <- chr_names

# DOWNLOAD THE SEQUENCE ---------------------------------------------------------------------------
# Manipulation of large biological sequences
suppressPackageStartupMessages(library("Biostrings", character.only = TRUE))
suppressPackageStartupMessages(library("stringr",    character.only = TRUE))
suppressPackageStartupMessages(library("readr",      character.only = TRUE))
suppressPackageStartupMessages(library("dplyr",      character.only = TRUE))

# Download and unzip the genome sequence
sequence_url <- "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz"
download.file(url = sequence_url, destfile = file.path(tempdir(), "hsa38.fasta.gz"), quiet = TRUE)
system(command = paste("gunzip", file.path(tempdir(), "hsa38.fasta.gz")))
all_sequence <- readDNAStringSet(filepath = file.path(tempdir(), "hsa38.fasta"))

# Select the sequences to export
if (opt$seq %in% chr_names) {
  message(paste0("[", Sys.time(), "] [INFO]: Downloading ", opt$seq, "sequence"))
  chr_sequence <- all_sequence[which(str_detect(string = names(all_sequence), pattern = hsa_accessions[opt$seq]))]
} else {
  message(paste0("[", Sys.time(), "] [INFO]: Downloading human genome GRCh38 sequence"))
  chr_sequence <- list()
  for (I in 1:length(hsa_accessions)) {
    chr_sequence[[I]] <- all_sequence[which(str_detect(string = names(all_sequence), pattern = hsa_accessions[I]))]
  }
  chr_sequence <- unlist(DNAStringSetList(chr_sequence))
}

# Export the sequence
writeXStringSet(x = chr_sequence, filepath = opt$out, append = FALSE)



# DOWNLOAD THE ANNOTATION -------------------------------------------------------------------------
if (opt$ann) {
  message(paste0("[", Sys.time(), "] [INFO]: Downloading sequence annotation"))
  # Download and unzip the genome annotation
  annotation_url <- "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz"
  download.file(url = annotation_url, destfile = file.path(tempdir(), "hsa38.gff.gz"), quiet = TRUE)
  annotation_tbl <- read_delim(file = file.path(tempdir(), "hsa38.gff.gz"), delim = "\t", col_names = FALSE, col_types = "ccciicccc", comment = "#")
  
  # Filter the annotation to extract only the chromosome of interest (if needed)
  if (opt$seq %in% chr_names) {
    annotation_tbl <- annotation_tbl %>%
      dplyr::filter(str_detect(X1, hsa_accessions[opt$seq]))
  }
  
  # Export the annotation
  annotation_gff <- str_replace(string = opt$out, pattern = ".fasta$", replacement = ".gff")
  annotation_gtf <- str_replace(string = opt$out, pattern = ".fasta$", replacement = ".gtf")
  write_delim(x = annotation_tbl, file = annotation_gff, delim = "\t", append = FALSE, col_names = FALSE)
  
  # Convert the annotation format from gff to gtf
  cmd <- paste("gffread -T", annotation_gff, ">", annotation_gtf)
  system(cmd)
}



# EXTRACT THE PROTEIN-CODING SEQUENCES ------------------------------------------------------------
message(paste0("[", Sys.time(), "] [INFO]: Extracting protein-coding sequences"))
# Extract protein-coding sequences (CDS)
if (opt$cds) {
  coding_sequences <- str_replace(string = opt$out, pattern = ".fasta$", replacement = ".CDS.fasta")
  cmd <- paste("gffread -g", opt$out, "-x", coding_sequences, annotation_gff)
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  out <- file.remove(paste0(opt$out, ".fai"))
}
