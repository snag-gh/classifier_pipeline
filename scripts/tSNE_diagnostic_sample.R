#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0) {
        stop("Please specify SampleName.\nUSAGE: ./tSNE_diagnostic_sample.R SampleID /path/to/samplesheet /path/to/pipeline /path/to/output\n", call.=FALSE)
}

.libPaths(c("/data/Compass/Tools/R_libraries/3.5/library", .libPaths()))
suppressPackageStartupMessages(library(mnp.v11b4))
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(Rtsne))
suppressPackageStartupMessages(library(RSpectra))
source(file.path(args[3], "scripts/classifier/RSpectra_pca.R"))

targets <- read.csv(args[2], stringsAsFactors=FALSE, skip=7)
sheet <- targets[targets$Sample_Name == args[1], ]
basedir <- dirname(args[2])
sheet$Basename <- file.path(basedir, paste(sheet$Sentrix_ID, sheet$Sentrix_Position, sep = "_"))
RGset <- read.metharray.exp(targets = sheet)
Mset <- MNPpreprocessIllumina(RGset)
Mset450 <- convertArray(Mset, outType = "IlluminaHumanMethylation450k", verbose = TRUE)
colnames(Mset450) <- sheet$Sample_Name

load(file.path(args[3], "scripts/classifier/RData/ba.coef.RData"))
load(file.path(args[3], "scripts/classifier/RData/betas.transposed.RData"))

filtered_probes <- rownames(tbetas)
keep <- rownames(Mset450) %in% filtered_probes
Mset450_filtered <- Mset450[keep,]
meths <- getMeth(Mset450_filtered)
unmeths <- getUnmeth(Mset450_filtered)

methy.coef <- as.data.frame(methy.coef)
keep <- rownames(methy.coef) %in% rownames(meths)
methy.coef.red <- methy.coef[keep, ]
meths <- meths[match(rownames(methy.coef.red), rownames(meths)),]
meth.b <- log2(meths + 1) + methy.coef.red$FFPE

unmethy.coef <- as.data.frame(unmethy.coef)
keep <- rownames(unmethy.coef) %in% rownames(unmeths)
unmethy.coef.red <- unmethy.coef[keep, ]
unmeths <- unmeths[match(rownames(unmethy.coef.red), rownames(unmeths)),]
unmeth.b <- log2(unmeths + 1) + unmethy.coef.red$FFPE

meth.b[meth.b < 0] <- 0
unmeth.b[unmeth.b <0] <- 0
meth.ba <- 2^meth.b
unmeth.ba <- 2^unmeth.b

betas.test <- meth.ba / (meth.ba + unmeth.ba +100)
betas.test <- as.data.frame(betas.test)

keep <- rownames(tbetas) %in% rownames(betas.test)
tbetas_keep <- tbetas[keep, ]
tbetas_keep$Sample <- betas.test[match(rownames(tbetas_keep), rownames(betas.test)), "betas.test"]
new_betas <- as.data.frame(t(tbetas_keep))

new_betas <- new_betas[,order(-apply(new_betas,2,sd))[1:32000]]
pca <- prcomp_svds(new_betas,k=94)
res <- Rtsne(pca$x,pca=F,max_iter=2500,theta=0,dim=3,verbose=T)
fname = file.path(args[4], paste0(args[1], "_tSNE.txt"))
write.table(res$Y, file = fname, row.names = FALSE, sep = "\t")

