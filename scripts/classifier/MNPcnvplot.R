#' Generates copy number variation plot
#'
#' \code{MNPcnvplot} uses function implemented in the \code{conumee} package to generate copy number profile plots
#' 
#' @param Mset an object of class Mset
#' @param sex will be predicted if not set
#' @details See \code{\link[conumee]{CNV.segment}} for more details
#' @import minfi
#' @import conumee
#' @export
MNPcnvplot <- function(Mset,sex=NULL,fname=NULL,output_dir=NULL,...){
  path <- paste(path.package('mnp.v11b4'),'/ext/',sep='')
  chiptype <- annotation(Mset)[1]
  if(is.null(sex)){
    Rset <- ratioConvert(Mset, what = "both", keepCN = TRUE)
    # convert to genomic ratio set
    sex <- ifelse(MNPgetSex(Rset)$predictedSex=="M","Male","Female")
  }
  if(chiptype=="IlluminaHumanMethylationEPIC"){
    load(paste(path,"IlluminaArrayDBconumee_annotation_EPIC_B4.2017-06-07.RData",sep="/"))
    cndata <- CNV.load(Mset)
    #if(sex=="Male"){ 
    #  load(paste(path,"CNanalysis4_conumee_REF-M.vh20150715.RData",sep="/"))
    #  x <- CNV.fit(cndata, refM.data, annoXY)
    #}
    #if(sex=="Female"){
    #  load(paste(path,"CNanalysis4_conumee_REF-F.vh20150715.RData",sep="/"))
    #  x <- CNV.fit(cndata, refF.data, annoXY)
    #}
    load(paste(path,"CNanalysis5_conumee_REF.2017-02-10.RData",sep="/"))
    x <- CNV.fit(cndata, refEPIC.data, annoEPICxy)
    x <- CNV.bin(x)
    x <- CNV.detail(x)
    x <- CNV.segment(x)
    par(mar=c(1, 1, 1, 1))
    CNV.genomeplot(x, chrY=ifelse(sex=="Female", FALSE, TRUE),...)
    CNV.write(x, what = "segments", file = paste0(file.path(output_dir, fname), ".seg"))
    CNV.write(x, what = "detail", file = paste0(file.path(output_dir, fname), ".detail.txt"))
    CNV.write(x, what = "probes", file = paste0(file.path(output_dir, fname), ".CNVprobes.igv"))
  }else{
    load(paste(path,"CNanalysis4_conumee_ANNO.vh20150715.RData",sep="/"))
  #  load(paste(path,"CNanalysis4_conumee_ANNO.vh20150715.RData",sep="/"))
    cndata <- CNV.load(Mset)
    if(sex=="Male"){ 
      load(paste(path,"CNanalysis4_conumee_REF-M.vh20150715.RData",sep="/"))
      x <- CNV.fit(cndata, refM.data, annoXY)
    }
    if(sex=="Female"){
      load(paste(path,"CNanalysis4_conumee_REF-F.vh20150715.RData",sep="/"))
      x <- CNV.fit(cndata, refF.data, annoXY)
    }
    x <- CNV.bin(x)
    x <- CNV.detail(x)
    x <- CNV.segment(x)
    par(mar=c(1, 1, 1, 1))
    CNV.genomeplot(x, chrY=ifelse(sex=="Female", FALSE, TRUE),...)
  }
}
