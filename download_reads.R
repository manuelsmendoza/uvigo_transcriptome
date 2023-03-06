# SETUP -------------------------------------------------------------------------------------------
options(warn = -1)

# Command line parser accept short and long flag/options
suppressPackageStartupMessages(library("optparse", character.only = TRUE))

# Parse the options
opt_list   <- list(
  make_option(opt_str = c("-a", "--acc"), 
              type    = "character", 
              default = NULL, 
              help    = "Accession number"),
  make_option(opt_str = c("-o", "--out"), 
              type    = "character", 
              default = NULL, 
              help    = "Accession number"),
  make_option(opt_str = c("-c", "--cpu"), 
              type    = "double", 
              default = 1, 
              help    = "Accession number"))
opt_parser <- OptionParser(option_list = opt_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$acc)) {
  stop(paste0("[", Sys.time(), "] [FAIL]: Sample accession number required"))
}

if (is.null(opt$out)) {
  opt$out <- getwd()
}



# DOWNLOAD READS AND COMPRESS FILES ---------------------------------------------------------------
cmd <- paste("fasterq-dump --force", 
             "--outdir", opt$out,
             "--temp", opt$out, 
             "--threads", opt$cpu,
             opt$acc)
system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

sra_files <- list.files(path = opt$out, pattern = opt$acc, full.names = TRUE)

for (FILE in sra_files) {
  cmd <- paste("pigz --force --best --processes", opt$cpu, FILE)
}
