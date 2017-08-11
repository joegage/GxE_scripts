## Joe Gage
## 22 June 2015

## Format SNPs.HYBRIDS.60 as plink format

load("./files/SNPs.HYBRIDS.60.rda")

hyb <- TFI2
rm(TFI2); gc()

head(rownames(hyb))

# Get rid of monomorphic SNPs
n_alleles <- apply(hyb, 2, function(x) length(unique(x[!is.na(x)])))
mono <- n_alleles < 2

hyb <- hyb[, !mono]

missing <- apply(hyb, 2, function(x) sum(is.na(x))/length(x)) ##
hyb <- hyb[,missing < 0.8] ##


# Make GM data frame
names <- colnames(hyb)
chr <- gsub("S([0-9]{1,2})_[0-9]*", "\\1", names)
pos <- gsub("S[0-9]{1,2}_([0-9]*)", "\\1", names)

GM <- data.frame(Name=names, Chromosome=chr, Position=pos)

# Remove any SNPs on chromosome "0"
chr0 <- chr == 0
hyb <- hyb[,!chr0]
GM <- GM[!chr0, ]

## dim(hyb)
missing <- missing[missing < 0.8]
missing <- missing[chr != 0]
summary(missing)

# Make and impute GD
GD <- hyb
impMean <- function(x){
    u <- mean(x, na.rm=TRUE)
    x[is.na(x)] <- u
    x
}
GD <- apply(GD, 2, impMean)
GD <- as.data.frame(GD)
taxa <- tolower(gsub("[^[:alnum:]]", "", rownames(GD)))
GD <- cbind(taxa, GD)


save(list=c("GD", "GM"), file="hybrid_geno_forGAPIT.rda")

## Write out SNP list for null distribution of distance to genes
writeSNPs <- GM
writeSNPs$alleles <- NA
writeSNPs <- writeSNPs[,c(1,4,2,3)]
colnames(writeSNPs) <- c("rs", "alleles", "chr", "pos")
write.table(writeSNPs, "all_gwas_snps.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
