# this script fits mash on all the preprocessed SCAIP FastQTL eQTL mapping results and saves the model to be run on chunked data
# based on mashr_cross-sectional_Wald_height-weight-lab-uncorr/mashr-fit-model.R
# 4/2/2021 JR

library(ashr)
library(mashr)
library(data.table)
library(dplyr)
library(doParallel)
cores <- as.integer(Sys.getenv("SLURM_STEP_TASKS_PER_NODE"))
registerDoParallel(cores = cores)
library(RhpcBLASctl)
blas_set_num_threads(12)

## # done once and for all 4/2/2021:
## # 1. read in all the data and convert into a mashr object:
## lfc_files <-list.files(path=paste0("input/"),pattern=".*_lfc.txt",full.name=T)
## datasets <- gsub("_lfc.txt","",list.files(path=paste0("input/"),pattern=".*_lfc.txt"))
## SE_files <- list.files(path=paste0("input/"),pattern=".*_SE.txt", full.name=T)
## p_files <- list.files(path=paste0("input/"),pattern=".*_pvalue.txt", full.name=T)
## lfsr_files <- list.files(path=paste0("input/"),pattern=".*lfsr.txt",full.name=T)
## # read in all the files:
## ldf_lfc <- lapply(lfc_files, read.table, sep="\t", header=T)
## ldf_SE <- lapply(SE_files, read.table, sep="\t", header=T)
## ldf_p <- lapply(p_files, read.table, sep="\t", header=T)
## ldf_lfsr <- lapply(lfsr_files, read.table, sep="\t", header=T)
## # merge them into one data frame:
## # reduced_lfc <- Reduce(merge(ldf_lfc)) # WHY won't this work??!
## lfcs <- Reduce(function(x,y) merge(x = x, y = y, all=TRUE), ldf_lfc)
## # remove the duplicated row for now:
## rownames(lfcs) <- lfcs[,1]
## lfcs <- lfcs[,-1]
## SEs <- Reduce(function(x,y) merge(x = x, y = y, all=TRUE), ldf_SE)
## # remove the duplicated row for now:
## rownames(SEs) <- SEs[,1]
## SEs <- SEs[,-1]
## pvalues <- Reduce(function(x,y) merge(x = x, y = y, all=TRUE), ldf_p)
## # remove the duplicated row for now:
## rownames(pvalues) <- pvalues[,1]
## pvalues <- pvalues[,-1]
## lfsrs <- Reduce(function(x,y) merge(x = x, y = y, all=TRUE), ldf_lfsr)
## # remove the duplicated row for now:
## rownames(lfsrs) <- lfsrs[,1]
## lfsrs <- lfsrs[,-1]
## # save it for fututre use since this step takes forever:
## save(lfcs, SEs, pvalues, lfsrs,file=paste0("./mashr_input.Rd"))

load(paste0("./mashr_input.Rd"))
# some rows have some NAs, but I think it doesn't work with NAs
# so let's remove all NAs as well as Inf from SEs:
SEs <- SEs[!rowSums(is.finite(as.matrix(SEs)))<ncol(SEs),]
# remove the same rows from pvalue and lfc matrixes:
lfcs <- lfcs[rownames(SEs),]
pvalues <- pvalues[rownames(SEs),]
# check:
sum(rowSums(is.na(as.matrix(lfcs))>0))
sum(rowSums(is.na(as.matrix(pvalues))>0))
stopifnot(identical(rownames(lfcs),rownames(SEs)))
stopifnot(identical(rownames(lfcs),rownames(pvalues)))

# save the rownames:
write.table(rownames(pvalues),paste0("rownames_mash.txt"),row.names=F,col.names=F,quote=F)

data = mash_set_data(as.matrix(lfcs), as.matrix(SEs))
# calculate Vhat:
Vhat = estimate_null_correlation_simple(data)
data = mash_set_data(as.matrix(lfcs), as.matrix(SEs),V=Vhat)

# 2. Set up the DATA-DRIVEN covariance matrices # NEXT TRY USING BOTH OR JUST CANONICAL
# 2.1. select strong signals
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05)
# 2.2. Obtain initial data-driven covariance matrices
# this step estimates the covariances of the observed data
U.pca = cov_pca(data,2,subset=strong)
print(names(U.pca))
# canonical matrix:
U.c = cov_canonical(data) 
# 2.3. Apply Extreme Deconvolution
# this step estimates the covariances of the actual underlying effects
U.ed = cov_ed(data, U.pca, subset=strong)

# save all objects created up to now:
save(data,U.c,U.ed,Vhat, file="mashr-fit-model_objects.Rd")

# Step 3: fit the model
# this fits a mixture model to the data, estimating the mixture proportions
## Sys.time()
## m.ed = mash(data, U.ed,outputlevel=1)

Sys.time()
m   = mash(data, c(U.c,U.ed))
Sys.time()

save(data,m,Vhat,file=paste0("mash-model-fit.Rd"))

lfsr=data.frame(get_lfsr(m))

summary(as.factor(rowSums(lfsr<0.1)))
colSums(lfsr<0.1)

# save the lfsr results:
write.table(get_lfsr(m),paste0("output/lfsr/all.txt"),sep="\t",quote=F, row.names=T,col.names=T)
# save the posterior means results:
write.table(get_pm(m),paste0("output/posterior_mean/all.txt"),sep="\t",quote=F, row.names=T,col.names=T)
# save the posterior SDs results:
write.table(get_psd(m),paste0("output/posterior_SD/all.txt"),sep="\t",quote=F, row.names=T,col.names=T)
write.table(get_lfdr(m),paste0("output/lfdr/all.txt"),sep="\t",quote=F, row.names=T,col.names=T)
write.table(get_np(m),paste0("output/NegativeProb/all.txt"),sep="\t",quote=F, row.names=T,col.names=T)

### END 3/12/2021
