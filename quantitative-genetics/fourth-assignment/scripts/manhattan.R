########################################################
##  Author: Isabelle Cleynen                          ##
##  Date last modified: 12/12/2017                    ##
##  create manhattan plots of outcome  				  ##
########################################################

# based on code from gettinggeneticsdone
# Stephen Turner
# http://StephenTurner.us/
# http://GettingGeneticsDone.blogspot.com/
# See license at http://gettinggeneticsdone.blogspot.com/p/copyright.html

# R code for making manhattan plots from plink output files.

setwd("/home/james/Documents/KULeuven/Bioinformatics/second_year/Evolutionary-and-Quantitative-Genetics/quantitative-genetics/fourth-assignment/data")

# manhattan plot using base graphics. 
# Copy the code below (lines 17 to 67) to R
manhattan <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), annotate=NULL, ...) {

    d=dataframe
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
    
    if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
    d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
    d$logp = -log10(d$P)
    d$pos=NA
    ticks=NULL
    lastbase=0
    colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
    if (ymax=="max") ymax<-ceiling(max(d$logp))
    if (ymax<8) ymax<-8
    
    numchroms=length(unique(d$CHR))
    if (numchroms==1) {
        d$pos=d$BP
        ticks=floor(length(d$pos))/2+1
    } else {
        for (i in unique(d$CHR)) {
          if (i==1) {
     d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
     } else {
     lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
     d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
     }
     ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
     }
    }
    
    if (numchroms==1) {
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
    } else {
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
        axis(1, at=ticks, lab=unique(d$CHR), ...)
        icol=1
        for (i in unique(d$CHR)) {
            with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
            icol=icol+1
     }
    }
    
    if (!is.null(annotate)) {
        d.annotate=d[which(d$SNP %in% annotate), ]
        with(d.annotate, points(pos, logp, col="green3", ...))
    }
    
    if (suggestiveline) abline(h=suggestiveline, col="blue")
    if (genomewideline) abline(h=genomewideline, col="red")
}


## to acually make the manhattan plot:
results <- read.table("gwa_clean.assoc.logistic", header=T)
head(subset(results, select=c(SNP, CHR, BP, P)))

# jpeg(filename="gwa_manhattan.jpg")
manhattan(results, colors=c("#008B8B","#A9A9A9"), pch=20)
# dev.off()

# make QQ plot
data <- read.table("gwa_clean.assoc.logistic", header=TRUE)
# jpeg("pvalue.qq.plot.jpg")
obs <- -log10(sort(data$P))
exp <- -log10(c(1:length(obs))/(length(obs) + 1))
plot(exp, obs, ylab="Observed (−logP)", xlab="Expected(−logP)", ylim=c(0,20), xlim=c(0,7))
lines(c(0,7), c(0,7), col=1, lwd=2)
# dev.off()