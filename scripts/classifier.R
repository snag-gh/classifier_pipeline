#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0) {
        stop("USAGE: ./classifier.R /path/to/samplesheet /path/to/output.pdf sampleID sentrix_id /path/to/pipeline version\n", call.=FALSE)
} 

.libPaths(c(.libPaths(), "/data/Compass/Tools/R_libraries/3.5/library"))
#cwd <- getwd()
source(file.path(args[5], "scripts/classifier/MNPreport.R"))
suppressPackageStartupMessages(library(mnp.v11b4))
suppressPackageStartupMessages(library(mnpqc))

path <- dirname(args[1])
sampleEpic <- args[4]
pathEpic <- paste(path,sampleEpic,sep='/')

RGsetEpic <- read.metharray(pathEpic, verbose=FALSE)

output_dir <- dirname(args[2])
output_file <- basename(args[2])

target <- read.csv(args[1], stringsAsFactors = FALSE, skip = 7)
sheet <- target[target$Sample_Name == args[3], ]


MNPreport(RGsetEpic, sample=args[4], sampleID=args[3], case=sheet$Surgical_Case, material=sheet$Material_Type, sex=sheet$Gender, output_file=output_file, output_dir=output_dir, pipeline=args[5], version=args[6])


#scripts/classifier.R /data/Compass/iScan_raw/202827620168/Sample_Sheet.csv report.pdf BSC-4897 202827620168_R01C01 `pwd`
