# loading libraries
library(dplyr)

# set working directory
setwd("/home/james/Documents/KULeuven/Bioinformatics/second_year/Evolutionary-and-Quantitative-Genetics/quantitative-genetics/fourth-assignment/data")

#### PCA
dat <- read.table("gwa_hapmap.pruned.pca.eigenvec", h=F)
colnames(dat)[1] <- "FID" #change column header of first column to FID 
colnames(dat)[2] <- "IID" #change column header of second column to IID 
# add a column in the dat file indicating who is CEU, CHB, JPT, YRI; and who is from your own gwa dataset. This way we can colour the samples in the plot according to their etnicity/dataset. This information is present in the fam file of the merged dataset, in which the last column indicates whether the individual is a case (=='2'), control (=='1'), or which HapMap population (=='3','4','5',or '6'). 
pop <- read.table("gwa_hapmap.pruned.fam", header = F)
colnames(pop)[1] <- "FID"
colnames(pop)[2] <- "IID"
colnames(pop)[6] <- "group"
pops <- subset(pop, select=c(1,2,6))
head(pops)
# merge dat and pops
data <- merge(dat, pops, by="IID", all=FALSE)
names(data)
head(data)
dim(data)
str(data)
data$group <- as.factor(data$group)
cont <- which(data$group=="1")
case <- which(data$group=="2")
CEU <- which(data$group=="3") #note that this indicates that those individuals with '3' in phenotype column are the CEU (european population) individuals
CHB <- which(data$group=="4") 
JPT <- which(data$group=="5")
YRI <- which(data$group=="6")
#pdf("pca-ancestry-plot.pdf")
plot(0, 0, pch="", xlim=c(-0.1,0.05), ylim=c(-0.1,0.1), xlab="principal component 1", ylab="principal component 2")
points(data$V3[JPT], data$V4[JPT], pch=20, col="PURPLE")
points(data$V3[CHB], data$V4[CHB], pch=20, col="PURPLE")
points(data$V3[YRI], data$V4[YRI], pch=20, col="GREEN")
points(data$V3[CEU], data$V4[CEU], pch=20, col="RED")

par(cex=0.5)
points(data$V3[cont], data$V4[cont], pch="+", col="BLACK")
points(data$V3[case], data$V4[case], pch="+", col="BLACK")
# abline(h=0.072, col="gray32", lty=2)
abline(v=0.004, col='blue', lty=2)
legend("topright", 
       legend = c("JPT/CHB", "YRI", "CEU", "CASE/CONTROL"), 
       col = c("purple", "green", "red", "black"), 
       pch = c(19, 19, 19, 3), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.06, 0.06))

# View(data)
fails <- data %>% filter(group %in% c(1,2), V3 < 0.004) %>% select(IID, FID.x)
write.table(as.matrix(fails), file='fail-ancestry-QC.txt', col.names=F, row.names=F, quote=F)
# View(fails)
