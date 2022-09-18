# this file uses lm to test for genotype by puberty interaction after qunatile-normalizing the trait
# 1/28/2021 JR

library(dplyr)
library(qvalue)

args = commandArgs(trailingOnly=TRUE)

sex <- "Female"
trait <- "cgpd1"
k=0

if (length(args)>0){
    sex <- args[1]
    trait <- args[2]
    }

FDR <- 0.1

# run through all the traits and sexes:
## for(sex in c("Females","Males")){
## if(sex=="Females"){
## traits=c("cgpd1","cgpd3","cgpd4","cgpd5","cgpd5a","cgpd5b")
##     }esle{
## traits=c("cbpd1","cbpd3","cbpd4","cbpd5")
##         }
##     for(trait in traits){

# load the needed files:
# 1. "bed" file with gene expression:
GE <- read.table(paste0("/nfs/rprdata/ALOFT/AL1-6_ln/counts/GC/gene_counts/residuals/for_plotting/residuals_cross-sectional_cage1_", sex,"_age-weight-height-retained.txt"), sep="\t", stringsAsFactors=FALSE)
colnames(GE) <- gsub("[.]","-", colnames(GE))
# 2. txt file with dosages:
dosbed <- read.table("../../permutations/analysis/AL_PC0_signif-dosages.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(dosbed) <- dosbed[,3]
colnames(dosbed)[4:length(colnames(dosbed))] <- gsub("[.]","-", colnames(dosbed)[4:length(colnames(dosbed))])
dos <- dosbed[,-c(1:3)]

# 3. txt file with list of testable pairs:
pairs <- read.table("../../permutations/analysis/PC0_significant_topeeQTL_pairs.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# 4. covariate file:
cv1 <- read.table("/nfs/rprdata/ALOFT/AL1-6_ln/covariates/ALOFT_RNA-seq_sample_masterfile_corrected_merged_8-23-19_puberty_4-1-2021.txt", sep="\t", header=T, stringsAsFactors=F, quote='"',comment="")
# keep only samples selected for eQTL analysis and order by dbgap:
cv2 <- cv1[cv1$expanded_eQTL=="TRUE",]
cv2 <- cv2[order(cv2$dbgap.ID),]
rownames(cv2) <- cv2$barcode

# subset everything to the given sex and individuals with measured trait:
cv <- filter(cv2, csex1_0==sex)
cv <- cv[!is.na(cv[,trait]),]
# subset GE:
GE <- GE[,cv$barcode]
# subset pairs to genes present in GE:
pairs <- pairs[pairs$pid %in% rownames(GE),]
# subset the dos:
dos <- dos[,cv$dbgap.ID]

# quantile-normalize the trait:
trait_norm <- qqnorm(rank(cv[,trait], ties.method = "random"), plot = F)$x

# subset GE df to only tested genes and sort the same way as in pairs:
expression <- GE[pairs$pid,]

# duplicate dosages that are eQTLs for several genes (and order as pairs):
dosages <- dos[pairs$sid,]

# transpose expression and dosages (because lm takes columns):
expression <- t(expression)
dosages <- t(dosages)

# make a data frame that will take the results:
pairs$Intercept_pval <- NA
pairs$dosage_pval <- NA
pairs$trait_pval <- NA
pairs$interaction_pval <- NA
pairs$Intercept_effect <- NA
pairs$dosage_effect <- NA
pairs$trait_effect <- NA
pairs$interaction_effect <- NA
pairs$Intercept_SE <- NA
pairs$dosage_SE <- NA
pairs$trait_SE <- NA
pairs$interaction_SE <- NA

# loop with the best number of PCs - 0:
 for (i in 1:ncol(expression)){
      model <- lm(expression[,i]~dosages[,i]*trait_norm)
      pvalues <- summary(model)$coefficients[,4]
      pairs[i,3] <- pvalues[1]
      pairs[i,4] <- pvalues[2]
      pairs[i,5] <- pvalues[3]
      pairs[i,6] <- pvalues[4]
      effects <- summary(model)$coefficients[,1]
      pairs[i,7] <- effects[1]
      pairs[i,8] <- effects[2]
      pairs[i,9] <- effects[3]
      pairs[i,10] <- effects[4]
      ses  <- summary(model)$coefficients[,2]
      pairs[i,11] <- ses[1]
      pairs[i,12] <- ses[2]
      pairs[i,13] <- ses[3]
      pairs[i,14] <- ses[4]
      }

# multiple test correct:
pairs$Intercept_padj <- p.adjust(pairs$Intercept_pval)
pairs$dosage_padj <- p.adjust(pairs$dosage_pval)
pairs$trait_padj <- p.adjust(pairs$trait_pval)
pairs$interaction_padj <- p.adjust(pairs$interaction_pval)

# save the results:
write.table(pairs, file=paste0("./GxE_results/GxE_", sex, "_",trait, "_0PCs.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# report number of significant interactions:
signif <- sum(pairs$interaction_padj<FDR)
tab <- t(c(sex, trait, signif, nrow(cv)))

write.table(tab, file="./signif_interactions_list_uncorrected.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

## make a qqplot and p-value histogram:
pdf(paste0("./plots/QQplots/GxE_pvalues_QQplot_",sex, "_",trait,"_0PCs.pdf"))
library(qqman)
qq(pairs$interaction_pval)
dev.off()
pdf(paste0("./plots/pval_histograms/GxE_pvalues_histogram_",sex, "_", trait, "_0PCs.pdf"))
hist(pairs$interaction_pval)
dev.off()
## }
##     }

### END 1/28/2021
