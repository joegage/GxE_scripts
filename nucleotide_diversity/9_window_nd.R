library(data.table)

print("Loading temperate file...")
temp = fread("temperate_nd.thetas", header=TRUE, stringsAsFactors=FALSE)
temp = as.data.frame(temp)
print("Loading tropical file...")
trop = fread("tropical_nd.thetas", header=TRUE, stringsAsFactors=FALSE)
trop = as.data.frame(trop)

##if( ! all(temp$Chromo == trop$Chromo & temp$Pos == trop$Pos)){
##	stop("Temperate and Tropical Thetas are not at the same positions")
##}
##print("All positions match")

#colnames(temp) = paste0(colnames(temp), ".temp")
#colnames(trop) = paste0(colnames(trop), ".trop")

high = read.table("highFst_windows.txt", stringsAsFactors=FALSE)
low = read.table("lowFst_windows.txt", stringsAsFactors=FALSE)

high$chromo = as.numeric(gsub("(.*):.*-.*", "\\1", high[,1]))
high$start = as.numeric(gsub(".*:(.*)-.*", "\\1", high[,1]))
high$end = as.numeric(gsub(".*:.*-(.*)", "\\1", high[,1]))

low$chromo = as.numeric(gsub("(.*):.*-.*", "\\1", low[,1]))
low$start = as.numeric(gsub(".*:(.*)-.*", "\\1", low[,1]))
low$end = as.numeric(gsub(".*:.*-(.*)", "\\1", low[,1]))

windows = rbind(high, low)
windows = windows[order(windows$chromo, windows$start), ]

trop$window = NA
temp$window = NA

library(doMC)

registerDoMC(cores=32)

winTheta = foreach(w=1:nrow(windows)) %dopar% {

	if(w %% 100 == 0){
		print(paste("Window", w))
	}
	
	chrW = windows[w, "chromo"]
	startW = windows[w, "start"]
	endW = windows[w, "end"]
	
	insideTemp = (temp$Chromo == chrW) &
			 (temp$Pos >= startW) & 
			 (temp$Pos <= endW)
	insideTrop = (trop$Chromo == chrW) &
			 (trop$Pos >= startW) & 
			 (trop$Pos <= endW)

	trop$window[insideTrop] = windows[w,1]
	temp$window[insideTemp] = windows[w,1]
	
	outTemp = temp[insideTemp, ]
	outTrop = trop[insideTrop, ]
	outName=paste0("window_nd/window_",chrW, ":", startW,"-",endW)
	
	write.table(outTemp, paste0(outName, "_temp.txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)	
	write.table(outTrop, paste0(outName, "_trop.txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)	

}

print("Done")

#write.table(winTheta, "thetas_persite_withWindows.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)