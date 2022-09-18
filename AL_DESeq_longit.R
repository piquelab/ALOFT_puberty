# this script does longitudinal analysis of puberty in ALOFT kids
# based on ./AL_DESeq_longit.R but model adapted for testing effect of age
# 3/5/2021 JR

require(ggplot2)
library(edgeR)
library(DESeq2)
library(qvalue)
library(dplyr)
library(reshape)
require(parallel)
require(BiocParallel)

myvar <- "cage1"
sex <- "Male"
# get command line arguments:
args = commandArgs(trailingOnly=TRUE)
if (length(args)>0){
      myvar <- args[1]
      sex <- args[2]
      }

# significance:
FDR <- 0.1

myvar
cat("##Processing ",myvar,"\t",sex,"\n")

cores            <- as.integer(Sys.getenv("NCPUS"))

timestamp()

## To run DESeq2 in parallel, using the
## BiocParallel library
register(MulticoreParam(cores))
ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

# load the covariate file:
cv <- read.table("/nfs/rprdata/ALOFT/AL1-6_ln/covariates/ALOFT_RNA-seq_sample_masterfile_corrected_merged_8-23-19_puberty_4-1-2021_updated.txt",sep="\t", stringsAsFactors=F, header=T, comment="", quote='"')

# susbset to the selected sex:
cv <- cv[cv$csex1_0==sex,]

# load the gene expression file:
data <- read.table("../../../counts/GC/gene_counts/AL_HTseq_gene_counts_only.txt",sep=" ", stringsAsFactors=F, header=T)
colnames(data) <- gsub("[.]","-",colnames(data))

# load annotation so you can remove X and Y chrs:
anno <- read.table("/nfs/rprdata/ALOFT/AL/counts/GC/gene_counts/HTSeq.anno.txt", sep="\t")
keep <- (rownames(data) %in% anno[anno$V1 %in% c(1:22),"V9"])
data <- data[keep,]

# keep only samples which have GE data (all do, though):
cv <- cv[cv$barcode %in% colnames(data),]

# keep only samples which have data for the variable:
cv <- cv[!is.na(cv[,myvar]),]

# keep only indidividuals with at least 2 timepoints:
summ <- summary(as.factor(cv$dbgap.ID), maxsum = 500)
cv <- cv[cv$dbgap.ID %in% names(summ)[summ>1],]

# keep only the first 2 timepoints for individuals with three timepoints:
cv <- data.frame(cv %>% group_by(dbgap.ID) %>% top_n(-2,sample.ID))

### longit height/puberty QC:
## all Males and Females age forward only

data <- data[,colnames(data) %in% cv$barcode]

# save the samples used in analysis:
write.table(cv[order(cv$sample.ID),c("sample.ID", "barcode")], paste0("logs/samples_", sex, "_",myvar,".txt"), sep="\t", col.names=F, row.names=F,quote=F)

##################################################################
## assign variables, load data, and load experiment information
topDirectory <- 'all-out_data'
outDir <- paste(topDirectory, "/",myvar, sep='')
system(paste("mkdir -p",outDir))

summary(as.numeric(cv[,myvar]))

sum(is.na(cv[,myvar]))

rownames(cv) <- cv$barcode

data <- data[,cv$barcode]

# drop lowly-expressed genes:
cpm <- cpm(data)
samples <- dim(data)[2]
table(rowSums(data==0)==samples) #Shows how many transcripts have 0 count across all samples
keep.exprs <- rowSums(cpm>0.25)>=(0.1*ncol(data)) #Only keep transcripts that have cpm>1 in at least 10% of samples
data <- data[keep.exprs,]

dim(data)
n <- dim(data)[2]

ddsFull <- DESeqDataSetFromMatrix(
    countData = round(data),
    colData = cv,
    design = as.formula(paste0(" ~
#as.factor(Plate) +
# as.numeric(genPC1) + as.numeric(genPC2) + as.numeric(genPC3) +
# as.factor(Wave) +
as.numeric(fraction_mapped) + as.numeric(percent_clean) + as.numeric(RIN) +
as.factor(dbgap.ID) + as.numeric(",myvar,")")))
## keep <- rowSums(counts(ddsFull)) > 0
## ddsFull <- ddsFull[keep,]
# colnames(ddsFull) <- cv$Filename

ddsFull<-DESeq(ddsFull,test="Wald")#,#parallel=TRUE,
               ## reduced =as.formula(paste0(" ~
#as.factor(Plate) +
# as.numeric(genPC1) + as.numeric(genPC2) + as.numeric(genPC3) +
#as.factor(Wave) +
## as.numeric(fraction_mapped) + as.numeric(percent_clean) + as.numeric(RIN)  + as.factor(dbgap.ID)")))

res <- results(ddsFull)#,parallel=TRUE) 

# some plots:
pdf(paste0("plots/QQplots/QQplot_", myvar, "_",sex,".pdf"))
library(qqman)
qq(res$pvalue)
dev.off()
pdf(paste0("./plots/histograms/histogram_",myvar,"_",sex,".pdf"))
hist(as.numeric(cv[,myvar]), xlab=NULL, main=paste0(myvar))
dev.off()

## cat("##20%FDR:\t",myvar,"\t",sum(res$padj<0.20 & abs(res$log2FoldChange)>0.25,na.rm=TRUE),"\n")
## cat("##15%FDR:\t",myvar,"\t",sum(res$padj<0.15 & abs(res$log2FoldChange)>0.25,na.rm=TRUE),"\n")
cat("##10%FDR:\t",myvar,sex,"\t",sum(res$padj<0.1,na.rm=TRUE),"\t",dim(cv)[1],"\n") #,"\t",sum(res$padj<0.1,na.rm=TRUE),"\n")
## cat("##5%FDR:\t",myvar,"\t",sum(res$padj<0.05 & abs(res$log2FoldChange)>0.25,na.rm=TRUE),"\n")
## cat("##1%FDR:\t",myvar,"\t",sum(res$padj<0.01 & abs(res$log2FoldChange)>0.25,na.rm=TRUE),"\n")

#Save results. 
write.table(res,file=gzfile(paste(outDir,"/",myvar,"_",sex,"_stats.txt.gz",sep="")),quote=F,sep="\t")

# volcano plot:
resp <- read.table(paste0("all-out_data/", myvar, "/", myvar,"_",sex, "_stats.txt.gz"), sep="\t", header=T, stringsAsFactors=F)
resp[!is.na(resp$padj) & resp$padj<=FDR,"significant"] <- "yes"
resp[!is.na(resp$padj) & resp$padj>FDR,"significant"] <- "no"
resp[is.na(resp$padj),"significant"] <- "no"
p <- ggplot(resp, aes(x=log2FoldChange, y=-log10(pvalue), col=significant)) +   geom_point() + theme_bw() + scale_color_manual(values=c("yes"="black", "no"="grey")) + ggtitle(paste0(myvar, " ", sex, ", longitudinal model"))
ggsave(paste0("plots/volcano/", myvar, "_", sex, ".png"),p)

# drop untested:
res <- res[!is.na(res$padj),]

DEensg <- rownames(res[res$padj<0.1,])# & res$abs(res$log2FoldChange)>0.25,])
DEensgup <- rownames(res[res$padj<0.1 & res$log2FoldChange>0,])# & res$log2FoldChange>0.25,])
DEensgdown <- rownames(res[res$padj<0.1 & res$log2FoldChange<0,])# & res$log2FoldChange<(-0.25),])
write.table(DEensg,file=paste0("./DEGs/",myvar,"_",sex,"_DEGs_10FDR.txt"),quote=F,sep="\t", row.names=F,col.names=F)
write.table(DEensgup,file=paste0("./DEGs/",myvar,"_",sex,"_up_DEGs_10FDR.txt"),quote=F,sep="\t", row.names=F,col.names=F)
write.table(DEensgdown,file=paste0("./DEGs/",myvar,"_",sex,"_down_DEGs_10FDR.txt"),quote=F,sep="\t", row.names=F,col.names=F)

# save all genes to use as background:
bckgrnd <- rownames(res[!is.na(res$padj),])
write.table(bckgrnd,file=paste0("./background_ensgs_",sex,"_",myvar,".txt"),quote=F,sep="\t", row.names=F,col.names=F)

sessionInfo()


### END 3/5/2021
