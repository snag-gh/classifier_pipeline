#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0) {
        stop("USAGE: ./3Dplot.R /path/to/plot.html sampleID /path/to/pipeline\n", call.=FALSE)
}

source(file.path(args[3], "scripts/classifier/plotting/tSNE.R"))

output_dir <- dirname(args[1])
output_file <- basename(args[1])

plotSNE(sampleID=args[2],output_file=output_file,output_dir=output_dir,pipeline=args[3])
