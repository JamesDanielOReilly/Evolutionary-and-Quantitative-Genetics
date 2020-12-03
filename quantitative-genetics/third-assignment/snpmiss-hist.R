## if working with notepad++, set the language to R if not automatically done so.

x<-read.table("GWA_sampleqc.lmiss",header=T)
ylabels=c("0","20K","40K","60K","80K","100K")
xlabels=c("0.0001","0.001","0.01","0.1","1")
#par(mfrow=c(1,1))
pdf("GWA_sampleqc.lmiss.pdf")
hist(log10(x$F_MISS),axes=F,xlim=c(-4,0),col="RED",ylab="Number of SNPs",xlab="Fraction of missing data",main="All SNPs",ylim=c(0,100000))
axis(side=2,labels=F)
mtext(ylabels,side=2,las=2, at=c(0,20000,40000,60000,80000,100000),line=1)
axis(side=1,labels=F)
mtext(xlabels,side=1,at=c(-4,-3,-2,-1,0),line=1)
abline(v=log10(0.05),lty=2)
dev.off()

# if you do not want to save this plot as a pdf, you can mute lines 7 and 14 of this code (using #)
# you can also choose to save it as a figure (png, tiff, jpeg...). Then adapt line 7 of the code accordingly.