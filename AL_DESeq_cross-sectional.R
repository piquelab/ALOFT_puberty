# this script does longitudinal analysis of puberty in ALOFT kids
# based on ./AL_DESeq_cross-sectional_cage1.R but adapted for cross-sectional analysis
# 3/12/2021 JR

require(ggplot2)
library(edgeR)
library(DESeq2)
library(qvalue)
library(dplyr)
library(reshape)
require(parallel)
require(BiocParallel)
library(stringr)

myvar <- "cgpd"
## sex <- "Female"
# get command line arguments:
args = commandArgs(trailingOnly=TRUE)
if (length(args)>0){
      myvar <- args[1]
      sex <- args[2]
      }

myvar
cat("##Processing ",myvar )#,"\t",sex,"\n")

cores            <- as.integer(Sys.getenv("NCPUS"))

timestamp()

## To run DESeq2 in parallel, using the
## BiocParallel library
register(MulticoreParam(cores))
ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

# load the covariate file:
cv <- read.table("/nfs/rprdata/ALOFT/AL1-6_ln/covariates/ALOFT_RNA-seq_sample_masterfile_corrected_merged_8-23-19_puberty_4-1-2021.txt",sep="\t", stringsAsFactors=F, header=T, comment="", quote='"')

# add cohort variable:
cv[str_detect(cv$Wave,"[[:upper:]]"),"cohort"] <- "original"
cv[!str_detect(cv$Wave,"[[:upper:]]"),"cohort"] <- "refresher"

## #subset to the desired sex:
## cv <- filter(cv, csex1_0 == sex)

# keep only samples which have data for the variable:
cv <- cv[!is.na(cv[,myvar]),]

# keep only the first sample for each individual:
cv <- cv[!duplicated(cv$dbgap.ID),]

# relevel cgpd5:
if(myvar=="cgpd5"){
    cv[cv$cgpd5==2,"cgpd5"] <- 0
    }

# load the gene expression file:
data <- read.table("../../../counts/GC/gene_counts/AL_HTseq_gene_counts_only.txt",sep=" ", stringsAsFactors=F, header=T)
colnames(data) <- gsub("[.]","-",colnames(data))

# load annotation so you can remove X and Y chrs:
anno <- read.table("/nfs/rprdata/ALOFT/AL/counts/GC/gene_counts/HTSeq.anno.txt", sep="\t")
keep <- (rownames(data) %in% anno[anno$V1 %in% c(1:22),"V9"])
data <- data[keep,]

# keep only samples which have GE data (all do, though):
cv <- cv[cv$barcode %in% colnames(data),]

### longit height/puberty QC:
## all Males and Females age forward only

# save the values for QC:
write.table(cv[,c("sample.ID",myvar)], paste0("DESeq_values_", myvar,".txt"), sep="\t", col.names=TRUE, row.names=FALSE)

data <- data[,colnames(data) %in% cv$barcode]

# save the samples used in analysis:
write.table(cv[order(cv$sample.ID),c("sample.ID", "barcode")], paste0("logs/samples_",myvar,".txt"), sep="\t", col.names=F, row.names=F,quote=F)

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
as.factor(Plate) +
as.numeric(genPC1) + as.numeric(genPC2) + as.numeric(genPC3) +
# as.factor(Wave) +
as.numeric(fraction_mapped) + as.numeric(percent_clean) + as.numeric(RIN) +
as.numeric(",myvar,")")))
## keep <- rowSums(counts(ddsFull)) > 0
## ddsFull <- ddsFull[keep,]
# colnames(ddsFull) <- cv$Filename

ddsFull<-DESeq(ddsFull,test="Wald")#,#parallel=TRUE,
               ## reduced =as.formula(paste0(" ~
#as.factor(Plate) +
# as.numeric(genPC1) + as.numeric(genPC2) + as.numeric(genPC3) +
#as.factor(Wave) +
## as.numeric(fraction_mapped) + as.numeric(percent_clean) + as.numeric(RIN)+ as.factor(cohort)  + as.factor(dbgap.ID)")))

res <- results(ddsFull)#,parallel=TRUE) 

# some plots:
pdf(paste0("plots/QQplots/QQplot_", myvar, ".pdf"))
library(qqman)
qq(res$pvalue)
dev.off()
pdf(paste0("./plots/histograms/histogram_",myvar,".pdf"))
hist(as.numeric(cv[,myvar]), xlab=NULL, main=paste0(myvar))
dev.off()

## cat("##20%FDR:\t",myvar,"\t",sum(res$padj<0.20 & abs(res$log2FoldChange)>0.25,na.rm=TRUE),"\n")
## cat("##15%FDR:\t",myvar,"\t",sum(res$padj<0.15 & abs(res$log2FoldChange)>0.25,na.rm=TRUE),"\n")
cat("##10%FDR:\t",myvar,"\t",sum(res$padj<0.1,na.rm=TRUE),"\t",dim(cv)[1],"\n") #,"\t",sum(res$padj<0.1,na.rm=TRUE),"\n")
## cat("##5%FDR:\t",myvar,"\t",sum(res$padj<0.05 & abs(res$log2FoldChange)>0.25,na.rm=TRUE),"\n")
## cat("##1%FDR:\t",myvar,"\t",sum(res$padj<0.01 & abs(res$log2FoldChange)>0.25,na.rm=TRUE),"\n")

#Save results. 
write.table(res,file=gzfile(paste(outDir,"/",myvar,"_stats.txt.gz",sep="")),quote=F,sep="\t")

# drop untested:
res <- res[!is.na(res$padj),]

DEensg <- rownames(res[res$padj<0.1,])# & res$abs(res$log2FoldChange)>0.25,])
DEensgup <- rownames(res[res$padj<0.1 & res$log2FoldChange>0,])# & res$log2FoldChange>0.25,])
DEensgdown <- rownames(res[res$padj<0.1 & res$log2FoldChange<0,])# & res$log2FoldChange<(-0.25),])
write.table(DEensg,file=paste0("./DEGs/",myvar,"_DEGs_10FDR.txt"),quote=F,sep="\t", row.names=F,col.names=F)
write.table(DEensgup,file=paste0("./DEGs/",myvar,"_up_DEGs_10FDR.txt"),quote=F,sep="\t", row.names=F,col.names=F)
write.table(DEensgdown,file=paste0("./DEGs/",myvar,"_down_DEGs_10FDR.txt"),quote=F,sep="\t", row.names=F,col.names=F)

# save all genes to use as background:
bckgrnd <- rownames(res[!is.na(res$padj),])
write.table(bckgrnd,file=paste0("./background_ensgs_",myvar,".txt"),quote=F,sep="\t", row.names=F,col.names=F)

sessionInfo()


### END 3/12/2021
