## Joe Gage
## 30 June 2015

## Analyze results from running stability regression parameters of
## GxE hybrid data through TASSEL's MLM

library(data.table)
library(ggplot2)
library(cowplot)

nLoci <- 50
window <- 0
ld_thresh <- 0.5

LD_file <- "./hybrid_geno_g2f.ld"
ld <- as.data.frame(fread(LD_file))

dir <- "./vanRaden/defGroups/"
files <- list.files(path=dir, pattern="GAPIT.*GWAS\\.Results\\.csv", full.names=T, recursive=TRUE)

params <- c("mse", "slope")

for(param in params){
    trait_param <- gsub(".*/GAPIT\\.\\.(.*)\\.GWAS\\.Results\\.csv", "\\1", files)
    
    read <- grepl(param, trait_param)
    all <- NULL
    for(file in files[read]){
        tmp <- as.data.frame(fread(file, stringsAsFactors=FALSE, sep=","))
        tmp$Trait <- gsub(".*/GAPIT\\.\\.([a-z]*)_.*", "\\1", file)
        all <- rbind(all, tmp)
    }

    colnames(all) <- c("Marker", "Chr", "Pos", "p", "maf", "nobs", "Rsq.without", "Rsq.with","FDR.adj.p", "Trait")

    traits <- unique(all$Trait)

    properNames <- c("Ear Height", "Plant Height", "Days to Anthesis", "Days to Silk", "Grain Yield", "MSE", "Slope")
    names(properNames) <- c(traits, params)

     all$cumulative <- all$Pos
     for(chr in 1:10){
         max <- max(all[all$Chr == chr, "Pos"])
         if(chr<10){
             all[all$Chr > chr, "cumulative"] <- all[all$Chr > chr, "cumulative"] + max
         }
     }

     borders <- rep(0,11)
     for(chr in 1:10){
         borders[chr+1] <- max(all$cumulative[all$Chr == chr], na.rm=T)
     }

     labelpos <- rep(1,10)
     for(i in 1:10){
         pos <- mean(all$cumulative[all$Chr == i], na.rm=T)
         labelpos[i] <- pos
     }

    for(trait in traits){
        subset <- all$Trait == trait
        dat <- all[subset, ]
       
        ## Assign colors by chromosome for plot
        dat$color <- dat$Chr %% 2
      
        ##ggplot Manhattan
        title <- paste0("GWAS of ", properNames[trait], " ", properNames[param], "\n", dat$nobs[1], " Individuals, ", nrow(dat), " Markers")
        man <- ggplot(dat, aes(x=cumulative, y=-log10(p), color=factor(color))) +
            geom_point() +
            scale_color_manual(values=c("blue4", "blue")) +
            scale_x_continuous("Chromosome", breaks=labelpos, labels=1:10) +
            labs(title=title) +
            guides(color=FALSE) +
            theme_classic()
        assign(paste0(trait, "_", param, "_man"), man)

        ## QQ plot
        pval <- dat$p
        pval <- pval[!is.na(pval)]
        pval <- pval[pval > 0 & pval <=1]
        pval <- sort(pval, decreasing=F)
        quantiles <- (1:length(pval))/(length(pval)+1)

        qq <- ggplot(data.frame(Observed=pval, Expected=quantiles), aes(x=-log10(Expected), y=-log10(Observed))) +
            geom_point() +
            geom_abline(slope=1, intercept=0, color="red") +
            labs(x="Expected", y="Observed") +
            theme_classic()

        assign(paste0(trait, "_", param, "_qq"), qq)

        source("./getTopLoci.R")
        
        top <- getTopLoci(dat, nLoci, window, ld, ld_thresh)
        write.table(top, paste0(param,"_", trait, "_top",nLoci,"_positions_",(window/1000),"kbWindow_additive.txt"), col.names=T, row.names=F, quote=F, sep="\t")             
    }

}

## Create file for each parameter with all manhattan and qq plots
for(param in params){
    g <- ggdraw()
    unit <- 1/length(traits)
    t <- 1 - unit
    for(trait in traits){
        p <- 1
        man <- get(paste0(trait, "_", param, "_man"))
        qq <- get(paste0(trait, "_", param, "_qq"))
        draw_plot(man, )  
        g <- g +
            draw_plot(man, 0, t, 4/5, unit) +
            draw_plot(qq, 4/5, t, 1/5, unit)
        t <- t - unit

    }

    save_plot(paste0("all_GWAS_", param, ".png"), g, ncol=2, nrow=1, base_height=12, base_width=8)
}
