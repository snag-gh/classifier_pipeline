codel_plot <- function(data) {
  library(ggplot2)
  library(reshape2)
  dtmelt <- melt(data, id.vars = "Sample", measure.vars = c("chr1p", "chr19q"))
  dtmelt$variable <- ifelse(dtmelt$variable == "chr1p", "1p", "19q")
  dtmelt2 <- melt(data, id.vars = "Sample", measure.vars = c("threshold_1p", "threshold_19q"))
  dtmelt2$variable <- ifelse(dtmelt2$variable == "threshold_1p", "1p", "19q")
  dtmelt$variable <- factor(dtmelt$variable, levels = c("1p", "19q"))
  dtmelt2$variable <- factor(dtmelt2$variable, levels = c("1p", "19q"))

  g <- ggplot(dtmelt, aes(x=value, y=Sample)) + geom_point(size = 5) + xlim(-1.0, 0.2) + xlab("Weighted Segmentation mean") + geom_vline(data=dtmelt2, aes(xintercept = value), color = "darkred", size = 1) + theme(strip.text.y = element_text(angle = 180), axis.title.y = element_blank(), axis.text.y = element_blank())

  g+facet_wrap(~ variable, ncol = 1, strip.position = "left") + scale_y_discrete("",breaks=c(0))
}
