library(data.table)

print("Reading Temperate file...")
temp = fread("window_nd_temp.txt", header=TRUE, stringsAsFactors=FALSE)
print("Reading Tropical file...")
trop = fread("window_nd_trop.txt", header=TRUE, stringsAsFactors=FALSE)

temp = as.data.frame(temp)
trop = as.data.frame(trop)

temp = temp[,c(1,2,4,8)]
trop = trop[,c(1,2,4,8)]

temp$Pairwise = exp(temp$Pairwise)
trop$Pairwise = exp(trop$Pairwise)

print("Aggregating...")
winTemp = aggregate(temp[,3], by=list(window=temp$window), FUN=sum)
winTrop = aggregate(trop[,3], by=list(window=trop$window), FUN=sum)


high = read.table("highFst_windows.txt", stringsAsFactors=FALSE)
high$Fst = "High Fst"
low = read.table("lowFst_windows.txt", stringsAsFactors=FALSE)
low$Fst = "Low Fst"

windows = rbind(high, low)
colnames(windows) = c("window", "Fst")

print("Merging....")
winTemp = merge(winTemp, windows, by="window", all=TRUE)
winTrop = merge(winTrop, windows, by="window", all=TRUE)

winTemp$Population = "Temperate"
#winTemp$x[is.na(winTemp$x)] = 0
winTrop$Population = "Tropical"
#winTrop$x[is.na(winTrop$x)] = 0

winThetas = rbind(winTemp, winTrop)

winThetas$start = as.numeric(gsub(".*:(.*)-.*", "\\1", winThetas$window))
winThetas$end = as.numeric(gsub(".*:.*-(.*)", "\\1", winThetas$window))
winThetas$size = winThetas$end - winThetas$start

winThetas$Pi = winThetas$x / winThetas$size

rm(temp, trop, winTemp, winTrop); gc()
save.image(file="plot_thetas.rda")

library(ggplot2)

print("Making plots...")

funSum <- function(x){
    data.frame(ymin = median(x) - sd(x),
               y = median(x),
               ymax = median(x) + sd(x))
}

ggplot(winThetas, aes(x=Fst, y=Pi, color=Population)) +
    stat_summary(fun.data=funSum, na.rm=TRUE)

NDplot = ggplot(winThetas, aes(x=Fst, y=Pi, fill=Population)) +
    geom_boxplot() +
    ylab("Per-site Nucleotide Diversity") +
    scale_fill_manual(values=c("dodgerblue3", "green")) +
    theme_classic()

ggsave("nucleotide_diversity.png", NDplot)

library(reshape2)

winThetas = dcast(winThetas, window + Fst ~ Population)

highNDplot = ggplot(winThetas[winThetas$Fst == "High Fst", ], aes(Temperate, Tropical)) +
    geom_point() +
    coord_fixed(ratio=1, xlim=c(0, .03), ylim=c(0,.03)) +
    theme_classic()

ggsave("highFst_nucleotide_diversity_compare.png", highNDplot)
