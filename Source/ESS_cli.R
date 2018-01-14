### Command-line interface for ESS calculations

library(optparse)
source("ESS_input.r")

# Set defaults ------------------------------------------------------------

default.nthreads <- 1
default.verbose <- FALSE

# Parsing arguments -------------------------------------------------------

options <- list (
  
  make_option(
    opt_str = c("-t", "--threads"),
    dest    = "nthreads",
    type    = "integer",
    default = default.nthreads,
    help    = paste0("number of threads to use, defaults to ", default.nthreads),
    metavar = "4")
  
)

parser <- OptionParser(
  usage       = "Rscript %prog [options] input output",
  option_list = options,
  description = "\nan awesome R script"
)

cli <- parse_args(parser, positional_arguments = 2)

# Shortcuts ---------------------------------------------------------------

input    <- cli$args[1]
output   <- cli$args[2]
nthreads <- cli$options$nthreads

# Program -----------------------------------------------------------------

curdate <- format(Sys.Date(),"%d%b%Y")
ESS_input()
ESS_program()

foo(input, output, nthreads)