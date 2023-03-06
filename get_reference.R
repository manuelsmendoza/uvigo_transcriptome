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
              type    = "logic", 
              default = TRUE, 
              help    = "Download the annotation of the sequence"),
  make_option(opt_str = c("-p", "--cds"),
              type    = "logic",
              default = TRUE, 
              help    = "Extract protein coding sequences")
  )

opt_parser <- OptionParser(option_list = opt_list)
opt        <- parse_args(opt_parser)


