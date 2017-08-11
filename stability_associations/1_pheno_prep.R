## Joe Gage
## 16 June 2015

## Prepare file containing slopes, intercepts, and MSEs of regressing
## hybrids against environmental check means

library(ggplot2)
library(gridExtra)

ggtheme <- theme(panel.background=element_rect(fill = "white", colour="white"),
                 panel.grid=element_blank(),
                 legend.key=element_rect(colour="black"),
                 #legend.position=c(.8,.5),
                 axis.text=element_text(size = 12, colour="black"),
                 axis.title=element_text(size=14, face="bold"),
                 #legend.title=element_text(size=16),
                 #legend.text=element_text(size=14),
                 plot.title=element_text(size=14, face="bold")
                 )

dat <- read.csv("./files/hybrid_data_clean_060915.csv", header=T, stringsAsFactors=F)
common <- read.table("./files/common_hybrids.txt", header=F, stringsAsFactors=F)
colnames(common) <- "hybrid"

traitCols <- c("Experiment", "Rep", "Pedigree", "Pollen.DAP..days.", "Silk.DAP..days.", "Plant.height..in.", "Ear.height..in.", "Grain.yield..bu.A.")
set <- dat$Trial
dat <- dat[,traitCols]
colnames(dat) <- c("experiment", "rep", "pedigree", "pollendap", "silkdap", "pltht", "earht", "yld_bu_ac")

dat$experiment <- gsub("([A-Z]*)H([0-9])", "\\1\\2", dat$experiment)
# Merge IA1a,b,c and NY1,2 because they were grown in the sampe place
dat$experiment[which(dat$experiment %in% c("IA1a", "IA1b", "IA1c"))] <- "IA1"
dat$experiment[dat$experiment %in% c("NY1", "NY2")] <- "NY1"

## Make geno by location figure
source("Hybrids_x_Locations_script.R")
## Specify common hybrids
nExp <- length(unique(dat$experiment))

common$nLoc <- NA

for(i in 1:length(common$hybrid)){
    hybridDat <- dat[which(dat$pedigree == common[i,"hybrid"]), ]
    common[i,"nLoc"] <- length(unique(hybridDat$experiment))
}

common <- common[which(common$nLoc >= 20), ]

## Pull out phenotypes for the checks and aggregate them
checkDat <- dat[which(dat$pedigree %in% common$hybrid), ]

means_by_check <- aggregate(checkDat, by=list(checkDat$experiment, checkDat$pedigree), FUN=mean, na.rm=T)
means_by_check$experiment <- means_by_check$Group.1
means_by_check$pedigree <- means_by_check$Group.2
means_by_check <- means_by_check[,-(1:2)]

means_by_exp <- aggregate(checkDat, by=list(checkDat$experiment), FUN=mean, na.rm=T)
means_by_exp$experiment <- means_by_exp$Group.1
means_by_exp <- means_by_exp[,-1]

traits <- colnames(dat)[4:ncol(dat)]
names <- c("Pollen DAP", "Silk DAP", "Plant Height", "Ear Height", "Yield (bu/ac)")
index <- 1

for(trait in traits) {

    dir.create(paste0("./", trait))

    traitName <- names[index]

    expMeans <- means_by_exp[order(means_by_exp[ ,trait]), ]
    expMeans$rank <- 1:nrow(expMeans)
    expMeans <- expMeans[complete.cases(expMeans[,trait]), ]

    checkVar <- aggregate(checkDat, by=list(checkDat$pedigree), FUN=var, na.rm=T)
    checkVar$pedigree <- checkVar$Group.1
    checkVar <- checkVar[,-1]
    checkVar <- checkVar[order(checkVar[ , trait]), ]
    checkVar$color <- c("blue", "blue", rep("gray", (nrow(checkVar)-4)), "red", "red")

    checkMeans <- means_by_check
    checkMeans <- checkMeans[complete.cases(checkMeans[, trait]), ]
    checkMeans$color <- checkVar$color[match(checkMeans$pedigree, checkVar$pedigree)]

    # Match environment rank back into means of data by environment and hybrid
    checkMeans$rank <- expMeans$rank[match(checkMeans$experiment, expMeans$experiment)]
    checkMeans$experiment <- reorder(checkMeans$experiment, checkMeans$rank)

    ## Plot checks across locations

    pdf(paste0("./", trait, "/checks_across_envs_", trait, "_boxplot.pdf"), width=8, height=4)
    print(ggplot(checkMeans, aes_string(x="experiment", y=trait)) +
          geom_boxplot(outlier.shape=NA) +
          geom_point(position=position_jitter(width=.0)) +
          labs(x="Experiment", y=traitName) + ggtheme
          )
    dev.off()

    ## Plot check means against overall means at each location
    overallExpMeans <- aggregate(dat, by=list(dat$experiment), FUN=mean, na.rm=TRUE)
    matchMeans <- merge(overallExpMeans[,c("Group.1",trait)], expMeans[,c("experiment", trait)], by.x="Group.1", by.y="experiment", all=TRUE)
    pdf(paste0("./", trait, "/check_vs_overall_location_means_",trait,".pdf"), width=6, height=6)
    plot(matchMeans[,2], matchMeans[,3],
         pch=NA,
         xlab="Overall Means", ylab="Check Means", main=trait)
    text(matchMeans[,2], matchMeans[,3], matchMeans[,1])
    abline(a=0, b=1, lty=2, col="red")
    dev.off()


    pdf(paste0("./",trait,"/checks_across_envs_",trait,"_variance.pdf"), width=15, height=6)
    print(ggplot(checkMeans, aes_string(x = "experiment", y = trait, group="pedigree")) +
          geom_point(size = 4, colour=checkMeans$color) +
          geom_line(, colour=checkMeans$color) +
          ggtheme + labs(y = traitName, x=""))
    dev.off()

    pdf(paste0("./",trait,"/checks_across_envs_",trait,"_points.pdf"), width=15, height=6)
    print(ggplot(checkMeans, aes_string(x = "experiment", y = trait, color="pedigree", group="pedigree")) +
          geom_point(size = 4) +
          ggtheme + labs(y = traitName, x="", colour="Pedigree"))
    dev.off()

    pdf(paste0("./",trait,"/checks_across_envs_",trait,"_lines.pdf"), width=15, height=6)
    print(ggplot(checkMeans, aes_string(x = "experiment", y = trait, color="pedigree", group="pedigree")) +
          geom_point(size = 4) +
          geom_line() +
          ggtheme + labs(y = traitName, x="", colour="Pedigree"))
    dev.off()

    ## Hybrid regression on location means

    hybrids <- unique(dat$pedigree)

    reg <- matrix(nrow=length(hybrids), ncol= 6, dimnames=list(hybrids, c("beta0", "beta1", "MSE", "nobs", "nenv", "Rsq")))

    for(hybrid in hybrids){

        hybDat <- dat[dat$pedigree == hybrid, c("experiment", trait)]

        if(all(is.na(hybDat[,trait]))){
            reg[hybrid, ] <- NA
            next
        }
        hybDat <- hybDat[!is.na(hybDat[,trait]), ]

        nenv <- length(unique(hybDat$experiment))

        hyb_loc_means <- aggregate(hybDat, by=list(hybDat$experiment), FUN=mean, na.rm=T)
        hyb_var <- var(hyb_loc_means[,trait], na.rm=T)

        regression_dat <- merge(expMeans[,c("experiment", trait)],
                                hybDat,
                                by="experiment",
                                all=T)
        if(all(is.na(regression_dat[,3]))){
            reg[hybrid, "beta0"] <- NA
            reg[hybrid, "beta1"] <- NA
            reg[hybrid, "MSE"] <- NA
        } else{
            trait_regression <- lm(regression_dat[,3] ~ regression_dat[,2])
            reg[hybrid, "beta0"] <- trait_regression$coefficients[1]
            reg[hybrid, "beta1"] <- trait_regression$coefficients[2] -1 ## Subtract one for deviation from checks
            reg[hybrid, "MSE"] <- sum(trait_regression$residuals^2) / trait_regression$df.residual
            reg[hybrid, "Rsq"] <- summary(trait_regression)$r.squared
        }


        reg[hybrid, "nobs"] <- nrow(hybDat)
        reg[hybrid, "nenv"] <- nenv
    }

    write.table(reg, paste0("./", trait, "/regression_", trait, ".txt"), sep="\t", col.names=NA, row.names=T)

    pdf(paste0("./", trait, "/stability_histograms_",trait,".pdf"), width=6, height=6)
    par(mfrow=c(3,2))
    hist(reg[,1], breaks=20, col="cadetblue", main="Intercept", ylab="Count", xlab="")
    hist(reg[,2], breaks=20, col="cadetblue", main="Slope", ylab="Count", xlab="")
    hist(reg[,3], breaks=20, col="cadetblue", main="MSE", ylab="Count", xlab="")
    hist(reg[,6], breaks=20, col="cadetblue", main="R squared", ylab="Count", xlab="")
    hist(reg[,4], breaks=20, col="cadetblue", main="Number of Obs", ylab="Count", xlab="")
    hist(reg[,5], breaks=20, col="cadetblue", main="Number of Environments", ylab="Count", xlab="")
    dev.off()

    pdf(paste0("./", trait, "/stability_regressions_", trait, ".pdf"), width=8, height=4)
    spread <- range(dat[,trait], na.rm=T)
    plot(1, xlim=c((min(spread)*.9),(max(spread)*1.1)), ylim=c((min(spread)*.9),(max(spread)*1.1)),
         xlab="Environmental Mean",
         ylab=paste0("Hybrid Phenotype"))
    for(i in 1:nrow(reg)){
        if(!is.na(reg[i,2])){
            abline(a=reg[i,1], b=(reg[i,2]+1)) ## For check slope correction
        }}
    abline(a=0, b=1, col="red")
    dev.off()

    index <- index + 1
}
