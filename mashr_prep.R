# this script fits mash on all the preprocessed SCAIP FastQTL dispersion eQTL mapping results and generates the input for mashr (pvalues,SEs)
# based on mashr_cross-sectional_Wald_height-weight-lab-uncorr_final/mashr_prep.R
# 4/2/2021 JR

library(ashr)
library(data.table)

var <- 'cage1'
args = commandArgs(trailingOnly=TRUE)

if (length(args)>0){
      dataset <- args[1]
    }

dir <- "../cross-sectional_Wald_height-weight-lab-uncorr_bigger-sample/"

for(sex in c("Female","Male")){
# read in, grab effect estimate and its standard error
m <- fread(paste0(dir, "all-out_data/", var,"/",var, "_", sex,"_stats.txt.gz")) 
colnames(m)[1] <- "ENSG"
    
# calculate lfsr:
ashr=ash(m$log2FoldChange,m$lfcSE)
m$lfsr <- ashr$result$lfsr
    
# save the effect size estimate and the SE in separate files:
lfc <- m[,c("ENSG","log2FoldChange")]
colnames(lfc)[2] <- paste0(var, "_", sex)
write.table(lfc, paste0("input/", var, "_", sex, "_lfc.txt"), sep="\t", col.names=T, row.names=F, quote=F)
SE <- m[,c("ENSG","lfcSE")]
colnames(SE)[2] <- paste0(var, "_", sex)
write.table(SE, paste0("input/", var, "_", sex, "_SE.txt"), sep="\t", col.names=T, row.names=F, quote=F)
# also save pvalues to select "strong" signals for estimating covariances step:
pvals <- m[,c("ENSG","pvalue")]
colnames(pvals)[2] <- paste0(var, "_", sex)
# save lfsr to use instead of pvalues
write.table(pvals, paste0("input/", var, "_", sex, "_pvalue.txt"), sep="\t", col.names=T, row.names=F, quote=F)
lfsr <- m[,c("ENSG","lfsr")]
colnames(lfsr)[2] <- paste0(var, "_", sex)
write.table(lfsr, paste0("input/", var, "_", sex, "_lfsr.txt"), sep="\t", col.names=T, row.names=F, quote=F)

# gzip all:
system(paste0("for file in input/", var, "_", sex,"*txt ; do gzip -f $file > $file.gz; done"))
}

### END 4/2/2021 JR
