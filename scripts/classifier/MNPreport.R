#' Generate a MNP Methylation Profiling Report
#'
#' \code{MNPreport} generates a html report with classifier prediction, t-SNE clustering, copy number plot and MGMT prediction.
#' 
#' @param RGset an object of class RGset
#' @param sample the sample for which the scores are plotted
#' @param FFPE preparation protocol, will be estimated if not set
#' @param sampleID will be Sample\code{sample} if not set
#' @param tsne set true to get interactive plots and t-SNE clustering
#' @param output_file will be \code{sampleID}.html if not set
#' @param output_fir will be current directory if not set
#' @export
#' @import minfi
#' @import conumee
#' @import ggplot2

MNPreport <- function(RGset,sample=1,FFPE=NULL,sampleID=NULL,case=NULL,material=NULL,sex=NULL,tsne=FALSE,output_file=NULL,output_dir=getwd(),pipeline=NULL,version=NULL,...){
  if(is.null(sampleID)) sampleID <- paste0("Sample",sample)
  if(is.null(output_file)) output_file <- paste0(sampleID,".pdf")
#  if(tsne){
  report <- file.path(output_dir, "/NCIreport.Rmd")
  rmarkdown::render(report, output_dir=output_dir, output_file=output_file, intermediates_dir=output_dir, knit_root_dir=output_dir, params=list(sample=sample, sampleID=sampleID, case=case, material=material, sex=sex, output_dir=output_dir, version=version, pipeline=pipeline),...)
#  }else{
#  rmarkdown::render(system.file('report_2.Rmd',package = "mnp.v11b4"),output_dir=output_dir,
#                      output_file=output_file,...)  
#  }
}
