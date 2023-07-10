library("methylKit")
library("genomation")
library("GenomicRanges")


samples = list( "DG-109", "DG-112", "DG-115", "DG-128", "DG-129", "DG-132", "DG-133", "DG-134", "DG-143", "DG-156", "DG-160", "DG-163", "NO-105","NG-94","NG-100","NG-102", "NG-106", "NG-108", "NG-114", "NG-120", "NG-127", "NG-135", "NG-139", "NG-147") 
N = 12

file.list <- list(
"SLX-21202.DG-109.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.DG-112.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.DG-115.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.DG-128.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.DG-129.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.DG-132.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.DG-133.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.DG-134.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.DG-143.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.DG-156.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.DG-160.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.DG-163.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.NO-105.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.NG-94.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.NG-100.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.NG-102.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.NG-106.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.NG-108.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.NG-114.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.NG-120.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.NG-127.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.NG-135.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz",
"SLX-21202.NG-139.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz", 
"SLX-21202.NG-147.s_2.r_1_trimmed_bismark_bt2.bismark.cov.gz")

samples 
samples[[1]]
# read the listed files into a methylRawList object making sure the other
# parameters are filled in correctly.
myobj <- methRead(file.list,
           sample.id=samples,
           pipeline = "bismarkCoverage",
           assembly="hg38",
           treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0),
           mincov = 10
           )

# check number of samples
# What type of data is stored here?

for (itr in 1:N)
        {	
	outname = paste("methylstats",samples[itr], sep="_")
	fig_name = paste(outname,"pdf", sep=".")
        print(fig_name)
        pdf(fig_name) 
        getMethylationStats(myobj[[itr]], plot=TRUE, both.strands=FALSE)
        dev.off()
        }
for (itr in 1:N)
        {
        outname = paste("coveragestats",samples[itr], sep="_")
        fig_name = paste(outname,"pdf", sep=".")
        print(fig_name)
        pdf(fig_name)
	getCoverageStats(myobj[[itr]], plot=TRUE, both.strands=FALSE)
	dev.off() 
        }
for (itr in 1:N)
        {
	getCoverageStats(myobj[[itr]], plot=FALSE, both.strands=FALSE)
        } 

myobj.filt <- filterByCoverage(myobj,
                      lo.count=10,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)
myobj.filt.norm <- normalizeCoverage(myobj.filt, method = "median")

meth <- unite(myobj.filt.norm, destrand=FALSE)

# get percent methylation matrix
pm=percMethylation(meth)

# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)

# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
pdf("CpGSD.pdf")
hist(sds, breaks = 100)
dev.off() 

# keep only CpG with standard deviations larger than 2%
meth <- meth[sds > 2]

# This leaves us with this number of CpG sites
nrow(meth)

pdf("correlations.pdf")
getCorrelation(meth,plot=TRUE)
dev.off() 

pdf("clustering.pdf") 
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
dev.off() 

pdf("PCA.pdf") 
PCASamples(meth)
dev.off()


myDiff <- calculateDiffMeth(meth,
                            overdispersion = "MN",
                            adjust="qvalue")
myDiff

write.csv(myDiff, file = "diff.csv", row.names = FALSE)

#head(myDiff[order(-abs(myDiff$meth.diff))], n = 20)

# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.9,type="hyper") #qvalue is high 
#
# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.9,type="hypo") #qvalue is high
#
#
# get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.9) #qvalue is high

diffMethPerChr <- diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.9, meth.cutoff=25) #qvalue is high 
write.csv(diffMethPerChr, file ="diffMethPerChr.csv", row.names =FALSE) 

gene.obj=genomation::readTranscriptFeatures("../test.bed",remove.unusual=TRUE,
                              up.flank=1000,down.flank=1000,unique.prom=TRUE) 

promoters=regionCounts(myobj,gene.obj$promoters)
exons = regionCounts(myobj, gene.obj$exon) 
introns = regionCounts(myobj, gene.obj$intron) 
#intergenic - regionCounts(myobj, gene.obj$intergenic) 
print("head promoters then rest") 
head(promoters) 
head(exons[[1]]) 
head(introns[[1]]) 
#head(intergenic[[1]]) 


annot <- genomation::annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj ) 

print("getting members") 
members <- getMembers(annot)
members 
print("getting tss") 
tss <- getAssociationWithTSS(annot)
tss 

annotations <- cbind(members, tss) 
print("Here is annotations") 
annotations 
write.csv(annotations, file ="annotations.csv", row.names=FALSE)
stats <- getTargetAnnotationStats(annot,percentage=TRUE,precedence=TRUE) 
print("head of stats") 
head(stats) 


pdf("annotations.pdf")
plotTargetAnnotation(annot,precedence=TRUE,
    main="differential methylation annotation")
dev.off()


cpg.obj=genomation::readFeatureFlank("../cpg_bed.txt",
                           feature.flank.name=c("CpGi","shores"))

# convert methylDiff object to GRanges and annotate
diffCpGann=genomation::annotateWithFeatureFlank(as(myDiff25p,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                         feature.name="CpGi",flank.name="shores")

pdf("CpGannot.pdf")
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),
       main="differential methylation annotation")
dev.off()

members <- getMembers(diffCpGann)
members 

tss <- getAssociationWithTSS(diffCpGann)
tss

