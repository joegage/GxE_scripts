fst <- read.csv("~/Dropbox/papers/GxE/analysis/var_explained_2/TEMP.ALL.INF.csv", stringsAsFactors=FALSE)
fst <- fst[!is.na(fst$Fst), ]

high <- fst[fst$Fst > .5, ]
low <- fst[fst$Fst < .15, ]
all <- rbind(high, low)
all <- all[order(all$Chr, all$Pos_i), ]

highWindows <- unique(paste0(high$Chr,":", high$Pos_i, "-", high$Pos_f))
lowWindows <- unique(paste0(low$Chr,":", low$Pos_i, "-", low$Pos_f))
allWindows <- unique(paste0(all$Chr,":", all$Pos_i, "-", all$Pos_f))

write.table(highWindows, "highFst_windows.txt", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(lowWindows, "lowFst_windows.txt", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(allWindows, "allFst_windows.txt", quote=FALSE, col.names=FALSE, row.names=FALSE)
