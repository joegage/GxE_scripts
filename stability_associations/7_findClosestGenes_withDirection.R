library(data.table)
library(gridExtra)
library(ggplot2)
library(cowplot)

nSNPs <- 50
window <- 0
version <- 2

window <- window/1000
files <- list.files(path="./", pattern=paste0(".*top",nSNPs,"_positions_",window,"kbWindow_additive.txt"), full.names=TRUE)
files <- c(files,
           list.files(path="./perse_gwas/", pattern=paste0(".*top",nSNPs,"_positions_",window,"kbWindow_additive.txt"), full.names=TRUE))

if(version == 2){
    # AGPv2 from:
    # ftp://ftp.gramene.org/pub/gramene/maizesequence.org/release-5b/filtered-set/ZmB73_5b_FGS.gff.gz
    annotation <- fread("./ZmB73_5b.60_FGS.gff", stringsAsFactors=FALSE, header=FALSE)
} else if(version == 3){
    # AGPv3 from:
    # ftp://ftp.gramene.org/pub/gramene/maizesequence.org/release-5b+/zea_mays.protein_coding.gff
    # system("grep -vF '#' zea_mays.protein_coding.gff > zea_mays.protein_coding.txt") # Remove metadata lines
    annotation <- fread("./zea_mays.protein_coding.txt", stringsAsFactors=FALSE, header=FALSE, sep="\t")
}

annotation <- as.data.frame(annotation)
annotation <- annotation[annotation[,3] == "gene", ]
annotation <- annotation[annotation[,1] %in% 1:10, ]

annotation[,1] <- as.integer(annotation[,1])
annotation[,4] <- as.integer(annotation[,4])
annotation[,5] <- as.integer(annotation[,5])


find_closest_gene <- function(snp){
    snp <- suppressWarnings(as.numeric(as.vector(snp)))
    chr <- snp[3]
    pos <- snp[4]

    subset <- annotation[annotation[,1] == chr, ]
    #print(dim(subset))

    intra <- any((pos >= subset[,4]) & (pos <= subset[,5]))
    if(intra){
        dist <- 0
    } else{
        distStart <- pos - subset[,4]
        distEnd <- pos - subset[,5]

        minStart <- min(abs(distStart))
        minStartInd <- which.min(abs(distStart))
        minEnd <- min(abs(distEnd))
        minEndInd <- which.min(abs(distEnd))
                
        dist <- min(minStart, minEnd)

        if(pos < subset[minStartInd, 4]){
            dist <- dist * -1
            if(subset[minStartInd, 7] == "-"){
                dist <- dist * -1
            }
        } else{
            if(subset[minEndInd, 7] == "-"){
                dist <- dist * -1
            }
        }
    }

    return(dist)
}

## Get null distribution
## allSNPs <- fread("./all_gwas_snps.txt", header=TRUE, stringsAsFactors=FALSE)
## allSNPs <- as.data.frame(allSNPs)
## null <- apply(allSNPs, 1, find_closest_gene)
## save("null", file=paste0("null_distribution_distances_AGPv", version, ".RData"))
load(paste0("null_distribution_distances_AGPv", version, ".RData"))

for(file in files){
    print(file)
    snps <- read.table(file, header=TRUE, stringsAsFactors=FALSE)
    distances <- apply(snps, 1, function(x) find_closest_gene(x))
    assign(file, distances)
}

all_dist <- NULL
mse_dist <- NULL
slope_dist <- NULL
perse_dist <- NULL
for(file in files){
    all_dist <- c(all_dist, get(file))
    if(grepl("/mse_", file)){mse_dist <- c(mse_dist, get(file))}
    if(grepl("/slope_", file)){slope_dist <- c(slope_dist, get(file))}
    if(grepl("/PerSe_", file)){perse_dist <- c(perse_dist, get(file))}
}

proportions <- function(x){
    c(sum(x < -5000),
      sum(x >= -5000 & x < 0),
      sum(x == 0),
      sum(x > 0 & x <=5000),
      sum(x > 5000)) / length(x)
}

nullP <- proportions(null)
slopeP <- proportions(slope_dist)
mseP <- proportions(mse_dist)
perseP <- proportions(perse_dist)

props <- rbind(nullP, slopeP, mseP, perseP); props

counts <- function(x){
    c(sum(x < -5000),
      sum(x >= -5000 & x < 0),
      sum(x == 0),
      sum(x > 0 & x <=5000),
      sum(x > 5000))
}

nullC <- counts(null)
slopeC <- counts(slope_dist)
mseC <- counts(mse_dist)
perseC <- counts(perse_dist)

cts <- rbind(nullC, slopeC, mseC, perseC);cts


p.slope <- apply(cts, 2, function(x)
    binom.test(x["slopeC"], n=sum(cts["slopeC",]), p=x["nullC"]/sum(cts["nullC",]))$p.value
    )
ci.slope <- apply(cts, 2, function(x)
    binom.test(x["slopeC"], n=sum(cts["slopeC",]), p=x["nullC"]/sum(cts["nullC",]))$conf.int
    )
p.mse <- apply(cts, 2, function(x)
    binom.test(x["mseC"], n=sum(cts["mseC",]), p=x["nullC"]/sum(cts["nullC",]))$p.value
    )
ci.mse <- apply(cts, 2, function(x)
    binom.test(x["mseC"], n=sum(cts["mseC",]), p=x["nullC"]/sum(cts["nullC",]))$conf.int
    )
p.perse <- apply(cts, 2, function(x)
    binom.test(x["perseC"], n=sum(cts["perseC",]), p=x["nullC"]/sum(cts["nullC",]))$p.value
    )
ci.perse <- apply(cts, 2, function(x)
    binom.test(x["perseC"], n=sum(cts["perseC",]), p=x["nullC"]/sum(cts["nullC",]))$conf.int
    )

## Bonferroni adjusted p threshold would be:
.05/15

## Visualize results as bar plot
errorBars=FALSE

propsPlot <- melt(props)
propsPlot <- propsPlot[order(propsPlot$Var1), ]
propsPlot$pval <- c(rep(NA, 5), p.slope, p.mse, p.perse)
propsPlot$pval <- format(round(propsPlot$pval, 3), nsmall=2)
## Ugly hardcoding of changing pval display if it rounds to zero:
if(round(p.slope[2], 3) == 0){
    propsPlot$pval[7] <- "<0.001"
}
##
propsPlot$CI.low <- c(rep(NA, 5), ci.slope[1,], ci.mse[1,], ci.perse[1,])
propsPlot$CI.high <- c(rep(NA, 5), ci.slope[2,], ci.mse[2,], ci.perse[2,])
propsPlot$n <- c(rep(NA, 5), cts["slopeC",], cts["mseC",], cts["perseC",])
propsPlot$labs <- paste0(propsPlot$pval, "\nn=", propsPlot$n)
propsPlot$labs[1:5] <- NA
propsPlot$Var1 <- factor(propsPlot$Var1, levels(propsPlot$Var1)[c(1,4,2,3)])

if(errorBars==TRUE){
bars <- ggplot(propsPlot, aes(x=Var2, y=value, fill=Var1)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(name=NULL, labels=c("Upstream\nIntergenic", "Upstream\nGene-proximal", "Genic", "Downstream\nGene-proximal", "Downstream\nIntergenic"), breaks=1:5) +
    scale_fill_manual(name="", labels=c("All SNPs", "Per Se", "Slope","MSE"), values=c("gray", "seagreen3", "orange", "skyblue"))+
    geom_text(aes(label=labs), vjust=-4.25, position=position_dodge(0.9)) +
    geom_errorbar(aes(ymax=CI.high, ymin=CI.low), position=position_dodge(0.9), width=.1) +
scale_y_continuous(name="Proportion", limits=c(0,.65))
} else{
bars <- ggplot(propsPlot, aes(x=Var2, y=value, fill=Var1)) +
    geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(name=NULL, labels=c("Upstream\nIntergenic", "Upstream\nGene-proximal", "Genic", "Downstream\nGene-proximal", "Downstream\nIntergenic"), breaks=1:5) +
    scale_fill_manual(name="", labels=c("All SNPs", "Per Se", "Slope","MSE"), values=c("gray", "seagreen3", "orange", "skyblue"))+
    geom_text(aes(label=labs), vjust=-.75, position=position_dodge(0.9)) +
    scale_y_continuous(name="Proportion", limits=c(0,.65))    
}
bars

name <- ifelse(errorBars==TRUE, "regressionParameters_barPlots.png", "regressionParameters_barPlots_noError.png")
ggsave(name, bars, width=12, height=6)


## Figures for talks:
## figProps <- propsPlot[!grepl("perse", propsPlot$Var1), ]
## figProps$pval <- as.numeric(figProps$pval) * 15
## figProps$pval[figProps$pval > .05] <- NA

## null <- figProps
## null$value[!grepl("null", null$Var1)] <- 0
## null$pval[!grepl("null", null$Var1)] <- NA

## slope <- figProps
## slope$value[!grepl("null|slope", null$Var1)] <- 0
## slope$pval[!grepl("null|slope", null$Var1)] <- NA

## mse <- figProps;

## bars1 <- ggplot(null, aes(x=Var2, y=value, fill=Var1)) +
##     geom_bar(stat="identity", position="dodge") +
##     scale_x_continuous(name=NULL, labels=c("Upstream\nIntergenic", "Upstream\nGene-proximal", "Genic", "Downstream\nGene-proximal", "Downstream\nIntergenic"), breaks=1:5) +
##     scale_fill_manual(name="", labels=c("All SNPs", "Slope","MSE"), values=c("gray", "orange", "skyblue"))+
##     scale_y_continuous(name="Proportion", limits=c(0,.65)) +
##     geom_text(aes(label=pval), size=5, vjust=-.75, position=position_dodge(0.9)) +
##     theme(axis.text.x = element_text(size=16),
##           axis.text.y = element_text(size=16),
##           axis.title.y = element_text(size=16))
## ggsave("bars_forTalk_1.png", bars1, width=12, height=6)
## bars2 <- ggplot(slope, aes(x=Var2, y=value, fill=Var1)) +
##     geom_bar(stat="identity", position="dodge") +
##     scale_x_continuous(name=NULL, labels=c("Upstream\nIntergenic", "Upstream\nGene-proximal", "Genic", "Downstream\nGene-proximal", "Downstream\nIntergenic"), breaks=1:5) +
##     scale_fill_manual(name="", labels=c("All SNPs", "Slope","MSE"), values=c("gray", "orange", "skyblue"))+
##     scale_y_continuous(name="Proportion", limits=c(0,.65)) +
##     geom_text(aes(label=pval), size=5, vjust=-.75, position=position_dodge(0.9)) +
##     theme(axis.text.x = element_text(size=16),
##           axis.text.y = element_text(size=16),
##           axis.title.y = element_text(size=16))
## ggsave("bars_forTalk_2.png", bars2, width=12, height=6)
## bars3 <- ggplot(mse, aes(x=Var2, y=value, fill=Var1)) +
##     geom_bar(stat="identity", position="dodge") +
##     scale_x_continuous(name=NULL, labels=c("Upstream\nIntergenic", "Upstream\nGene-proximal", "Genic", "Downstream\nGene-proximal", "Downstream\nIntergenic"), breaks=1:5) +
##     scale_fill_manual(name="", labels=c("All SNPs", "Slope","MSE"), values=c("gray", "orange", "skyblue"))+
##     scale_y_continuous(name="Proportion", limits=c(0,.65)) +
##     geom_text(aes(label=pval), size=5, vjust=-.75, position=position_dodge(0.9)) +
##     theme(axis.text.x = element_text(size=16),
##           axis.text.y = element_text(size=16),
##           axis.title.y = element_text(size=16))
## ggsave("bars_forTalk_3.png", bars3, width=12, height=6)
