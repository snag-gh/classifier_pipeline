#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0) {
	stop("\nUSAGE: ./methylArrayQC.R /path/to/arrayID /path/to/Sample_Sheet.csv /path/to/qcReport.pdf /path/to/suplementary_plots.pdf\n", call.=FALSE)
} 

suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICmanifest))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(ggplot2))

gen_plot <- function(data, controlType) {

my_colors <- c("forestgreen" ,"darkorange", "darkorchid1", "gold", "dodgerblue", "orangered", "darkblue", "cyan4", "pink2", "firebrick", "coral4", "darkorchid4")
x <- ggplot(data, aes(x=variable, y=value, color=ExtendedType)) + geom_point() + facet_grid(. ~ channel) + theme_bw() + theme(axis.text.x = element_text(angle = 90), panel.border = element_rect(size = 0.25, color = "darkgrey"), panel.spacing = unit(0, "lines"), legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = my_colors) + labs(title = controlType, y = "Intensity", x = NULL)
return(x)

}

plot_control <- function(allcontrols, controlType, r, g) {

  control <- allcontrols[allcontrols$Type==controlType, ]
  control_red <- r[control$Address, , drop = FALSE]
  control_red <- as.data.frame(control_red)
  control_red$channel <- rep("red", nrow(control_red))
  control_red$Address <- row.names(control_red)

  control_green <- g[control$Address, , drop = FALSE]
  control_green <- as.data.frame(control_green)
  control_green$channel <- rep("green", nrow(control_green))
  control_green$Address <- row.names(control_green)

  control_both <- rbind(control_red, control_green)
  control_both$ExtendedType <- control[match(control_both$Address, control$Address), "ExtendedType"]
  control_both$Color <- control[match(control_both$Address, control$Address), "Color"]

  control_melt <- melt(control_both, id=c("Address", "channel", "Color", "ExtendedType"))
  x <- gen_plot(control_melt, controlType)
  return(x)

}

#baseDir <- file.path("/data/Compass/iScan_result", args[1], fsep = "/")
#my_file <- list.files(path = baseDir, pattern = 'Sample_Sheet.*csv$')
baseDir <- args[1]
my_file <- args[2]

#targets <- read.csv(file.path(baseDir, my_file), stringsAsFactors = FALSE, skip = 7)
targets <- read.csv(my_file, stringsAsFactors = FALSE, skip = 7)
targets$Basename <- file.path(baseDir, paste(targets$Sentrix_ID, targets$Sentrix_Position, sep = "_"))

RGset <- read.metharray.exp(targets = targets)
pd <- pData(RGset)
colnames(RGset) <- targets$Sample_Name
Mset <- preprocessRaw(RGset)

#outfile = paste(args[1], "qcReport", "pdf", sep = ".")
#outpath = file.path("/data/Compass/Methylation/QCreports", args[1])

#directory creation can probably be skipped when put in snakemake pipeline.
#dir.create(outpath)

#pdf(file = file.path("/data/Compass/Methylation/QCreports", args[1], outfile, fsep = "/")) 
pdf(file = args[3])

controls <- getProbeInfo(RGset, type = "Control")
controls <- as.data.frame(controls)
r <- getRed(RGset)
g <- getGreen(RGset)

clist <- c("STAINING", "EXTENSION", "HYBRIDIZATION", "TARGET REMOVAL", "BISULFITE CONVERSION I", "BISULFITE CONVERSION II", "SPECIFICITY I", "SPECIFICITY II", "NON-POLYMORPHIC")
for (D in clist) {
  my_plot <- plot_control(controls, D, r, g)
  print(my_plot)
}

negative <- controls[controls$Type=="NEGATIVE", ]
negative_red <- r[negative$Address, , drop = FALSE]
negative_red <- as.matrix(negative_red)
means <- apply(negative_red, 2, mean)
means <- as.data.frame(means)
means$ExtendedType <- rep("Average", nrow(means))
means$Sample <- rownames(means)
colnames(means) <- c("value", "ExtendedType", "variable")
sds <- apply(negative_red, 2, sd)
sds <- as.data.frame(sds)
sds$ExtendedType <- rep("StdDev", nrow(sds))
sds$Sample <- rownames(sds)
colnames(sds) <- c("value", "ExtendedType", "variable")
neg_red_summ <- rbind(means, sds)
neg_red_summ$channel <- rep("red", nrow(neg_red_summ))

negative_green <- g[negative$Address, , drop = FALSE]
negative_green <- as.matrix(negative_green)
means <- apply(negative_green, 2, mean)
means <- as.data.frame(means)
means$ExtendedType <- rep("Average", nrow(means))
means$Sample <- rownames(means)
colnames(means) <- c("value", "ExtendedType", "variable")
sds <- apply(negative_green, 2, sd)
sds <- as.data.frame(sds)
sds$ExtendedType <- rep("StdDev", nrow(sds))
sds$Sample <- rownames(sds)
colnames(sds) <- c("value", "ExtendedType", "variable")
neg_green_summ <- rbind(means, sds)
neg_green_summ$channel <- rep("green", nrow(neg_green_summ))

neg_all <- rbind(neg_red_summ, neg_green_summ)
neg_all$variable <- factor(RGset$Sample_Name, levels = RGset$Sample_Name)

x <- gen_plot(neg_all, "NEGATIVE")
print(x)

my_plot <- plot_control(controls, "RESTORATION", r, g)
print(my_plot)

dev.off()

pdf(file = args[4])
#Detection p-value plot
detP <- detectionP(RGset)
pal <- brewer.pal(8, "Dark2")
#temp <- colMeans(detP)
#temp <- as.data.frame(temp)
#temp$x <- seq(1, nrow(temp))
#detPplot <- ggplot(temp, aes(x = x, y = temp))
#detPplot + geom_point(col = pal[factor(targets$Sample_Name)], size = 2) + theme_bw() + scale_x_continuous(breaks = seq(1, nrow(temp)), labels = rownames(temp)) + coord_cartesian(ylim = c(0,0.01)) + geom_hline(yintercept = 0.01, linetype = "dashed") + theme(axis.text.x = element_text(angle = 90), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5)) + labs(y = "Mean detection p-values", x = NULL, title = "Mean detection p-values")

#Density plot of beta values
densityPlot(Mset, sampGroups = pd$Sample_Name)

#Bean plot of beta values
par(mar=c(5,6,4,2))
densityBeanPlot(Mset, sampGroups = pd$Sample_Name)

#QC plot
qc <- getQC(Mset)
#plotQC(qc)
qc <- as.data.frame(qc)
qcplot <- ggplot(qc, aes(x = mMed, y = uMed, col = rownames(qc)))
qcplot + geom_point(size = 3, alpha = 0.8) + theme_bw() + coord_cartesian(xlim = c(8, 14), ylim = c(8,14)) + scale_color_manual(values = pal) + labs(x = "Meth median intensity (log2)", y = "Unmeth median intensity (log2)", title = "Median intensities in the Methylated and Unmethylated channels") + theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) + scale_x_continuous(breaks = c(9, 11, 13)) + scale_y_continuous(breaks = c(9, 11, 13))

#Detection p-value plot
failed <- detP > 0.01
#head(failed)
fprobes <- (colSums(failed)/nrow(failed))*100
barplot(fprobes, las=2, col = pal[factor(names(fprobes))], cex.names = 0.8, cex.axis = 0.6, ylim = c(0,50), ylab = "Failed percentage")
abline(h = 10, col = "Red", lty = 5)

dev.off()

