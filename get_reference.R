# SETUP -------------------------------------------------------------------------------------------
# Command line parser accept short and long flag/options
library("optparse", character.only = TRUE, verbose = FALSE) 

# Parse the options
opt_list   <- list(
  make_option(opt_str = c("-s", "--seq"), 
              type    = "character", 
              default = NULL, 
              help    = "Sequence to download"),
  make_option(opt_str = c("-o", "--out"), 
              type    = "character", 
              default = NULL, 
              help    = "Output file"),
  make_option(opt_str = c("-a", "--ann"), 
              type    = "logical", 
              default = TRUE, 
              help    = "Download the annotation of the sequence"),
  make_option(opt_str = c("-p", "--cds"),
              type    = "logical",
              default = TRUE, 
              help    = "Extract protein coding sequences"))
opt_parser <- OptionParser(option_list = opt_list)
opt        <- parse_args(opt_parser)

chr_names  <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY")
if (is.null(opt$seq) || !opt$seq %in% c("genome", chr_names)) {
  print_help(opt_parser)
  stop("Sequence name to downlaod is missed or not allowed")
}

if (is.null(opt$out)) {
  opt$out <- file.path(getwd(), paste0(opt$seq, ".fasta"))
}



# DOWNLOAD THE SEQUENCE ---------------------------------------------------------------------------
# R interface to search in the different NCBI databases
library(rentrez)

# Define the accession of the human chromosomes (autosomes, X and Y)
hsa_accessions <- c(
  "NC_000001", "NC_000002", "NC_000003", "NC_000004", "NC_000005", "NC_000006",
  "NC_000007", "NC_000008", "NC_000009", "NC_000010", "NC_000011", "NC_000012", 
  "NC_000013", "NC_000014", "NC_000015", "NC_000016", "NC_000017", "NC_000018", 
  "NC_000019", "NC_000020", "NC_000021", "NC_000022", "NC_000023", "NC_000024")
names(hsa_accessions) <- chr_names

# Download and unzip the genome sequence
sequence_url <- "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz"

download.file(url = sequence_url, destfile = file.path(tempdir(), "hsa38.fasta.gz"))
system(command = paste("gunzip", file.path(tempdir(), "hsa38.fasta.gz")))
all_sequence <- readDNAStringSet(filepath = file.path(tempdir(), "hsa38.fasta"))

if (opt$seq %in% chr_names) {
  chr_sequence <- all_sequence[which(str_detect(string = names(all_sequence), pattern = opt$seq))]
} else {
  chr_sequence <- list()
  for (I in 1:length(hsa_accessions)) {
    chr_sequence[[I]] <- all_sequence[which(str_detect(string = names(all_sequence), pattern = hsa_accessions[I]))]
  }
  chr_sequence <- unlist(DNAStringSetList(chr_sequence))
}

writeXStringSet(x = chr_sequence, filepath = opt$out, append = FALSE)


