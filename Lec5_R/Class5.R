#Bioinformatics Class 5
#Plots

x <- rnorm(1000,0)

summary(x)

#see this data as a graph
boxplot(x)


#Good old histogram
hist(x)

#Section 1 from lab sheet

baby <- read.table("bggn213_05_rstats/weight_chart.txt", header = TRUE)

plot(baby,type="b",pch=15, cex=1.5,lwd=2,ylim=c(2,10),xlab="Age (months)",ylab="Weight(kg)",main="Some title")

feat <- read.table("bggn213_05_rstats/feature_counts.txt", sep="\t", header = TRUE)

par(mar=c(5,11,4,2))
barplot(feat[,2],horiz = TRUE, ylab="Atitle", names.arg = feat[,1], main="Some title", las=1)
 
hist(c(rnorm(10000),rnorm(10000)+4),breaks=50) 

#Section 2 from lab sheet
#male <- read.table("bggn213_05_rstats/male_female_counts.txt", sep="\t", header = TRUE)
male <- read.delim("bggn213_05_rstats/male_female_counts.txt")

barplot(male[,2],horiz = TRUE, ylab="Atitle", names.arg = male[,1], main="Some title", las=1,col=rainbow(10) )
barplot(male[,2],horiz = TRUE, ylab="Atitle", names.arg = male[,1], main="Some title", las=1,col=c("blue2","red2") )


expression <- read.delim("bggn213_05_rstats/up_down_expression.txt")
palette(c("blue","grey","red"))
plot(expression$Condition1, expression$Condition2, col=expression$State)



methyl <- read.delim("bggn213_05_rstats/expression_methylation.txt")
map.colors <- colorRampPalette(c("grey","red"))(100)
map.colors
plot(methyl$promoter.meth,methyl$gene.meth, col=map.colors)

#Section 3
chrom <- read.delim("bggn213_05_rstats/chromosome_position_data.txt")
plot(chrom$Position, chrom$WT, type="l")





