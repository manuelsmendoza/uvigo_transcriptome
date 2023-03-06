# SETUP -------------------------------------------------------------
# Command line parser accept short and long flag/options
library("optparse", character.only = TRUE, verbose = FALSE) 

# Parse the options
opt_list   <- list(
  make_option(opt_str = c("-s", "--seq"), 
              type    = "character", 
              default = "genome", 
              help    = "Sequence to download"),
  make_option(opt_str = c("-a", "--ann"), 
              type    = "logical", 
              default = TRUE, 
              help    = "Download the annotation of the sequence"),
  make_option(opt_str = c("-p", "--cds"),
              type    = "logical",
              default = TRUE, 
              help    = "Extract protein coding sequences")
  )

opt_parser <- OptionParser(option_list = opt_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$seq) || !opt$seq %in% c("genome", paste("chr", 1:22, sep = ""))) {
  print_help(opt_parser)
  stop("Sequence name to downlaod is missed or not allowed")
}

print(length(opt_list))
print(opt_list)
