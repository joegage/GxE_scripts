getTopLoci <- function(dat, nLoci, bufferSize, LD, ld_thresh){
### dat is a data frame of all MLM results
### nLoci is the number of top loci to keep
### bufferSize is the minimum distance between selected loci
### LD_file is a previously read in output table from plink 1.9's r^2 analysis

### For Testing:
   # library(data.table)
   # dat <- "./map_hybrids_slope_hybrid_geno_g2f_chr0_1_2_3_4_5_6_7_8_9_10_523360-150107228_+_slope_hybrids_forGWAS_lowerNames_forTASSEL_+_hybrid_geno_g2f_PCs_stats.txt"
   # dat <- as.data.frame(fread(dat))
   # dat <- dat[dat$Trait == "pltht", ]
   # nLoci <- 100
   # bufferSize <- 0
   # LD_file <- "./hybrid_geno_g2f.ld"

    dat <- dat[order(dat$p), ]
    #topLoci <- dat[1:nLoci, c("Chr", "Pos", "add_p")]
                                       
    #placeholder <- nLoci + 1
    #old <- 1
    #new <- 2

    ## Iterate through SNPs from most to least significant,
    ## and remove any less significant SNPs that have are within
    ## 5k AND have r^2 > ld_thresh with that SNP.
    ## Have to screen both 'SNP_A' and 'SNP_B' columns for the
    ## significant target SNP, because PLINK's r^2 output does not
    ## have all pairwise combinations of SNPs
    for(snp in 1:nLoci){
        ld_subA <- ld[ld$SNP_A == dat[snp, "Marker"], ]
        ld_subB <- ld[ld$SNP_B == dat[snp, "Marker"], ]

        rm_A <- ld_subA[ld_subA$R2 > ld_thresh, "SNP_B"]
        rm_B <- ld_subB[ld_subB$R2 > ld_thresh, "SNP_A"]

        rm_A <- dat$Marker %in% rm_A
        dat <- dat[!rm_A, ]
        rm_B <- dat$Marker %in% rm_B
        dat <- dat[!rm_B, ]

    }

    topLoci <- dat[1:nLoci, c("Chr", "Pos", "p")]

    ## This section was for when we were using windows to buffer against
    ## picking SNPs that were in high LD. Pretty deprecated now.
    
    ## while(any(old != new)){
    ##     old <- topLoci
        
    ##     for(locus in 1:nLoci){
    ##         rownames(topLoci) <- 1:nLoci
    ##         remove <- which(as.numeric(rownames(topLoci)) > locus &
    ##                         topLoci$Chr == topLoci$Chr[locus] &
    ##                         topLoci$Pos > (topLoci$Pos[locus]-bufferSize) &
    ##                         topLoci$Pos < (topLoci$Pos[locus]+bufferSize))
    ##         if(!length(remove) == 0){
    ##             topLoci <- topLoci[-remove, ]
                
    ##             add <- nLoci - nrow(topLoci)

    ##             toAdd <-dat[placeholder:(placeholder + add -1), c("Chr", "Pos", "add_p")]
                
    ##             topLoci <- rbind(topLoci, toAdd)

    ##             placeholder <- placeholder + add

    ##         }
    ##     }

    ##     new <- topLoci

    ## }
    topLoci$rs <- paste(topLoci$Chr, topLoci$Pos, sep="_")
    topLoci$alleles <- NA
    topLoci[,c("rs", "alleles", "Chr", "Pos", "p")]
}
