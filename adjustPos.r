############################
### adjustPos.r --- collapse chromosome & position into one
### by Tina
### created on 07.28.08
############################
adjustPos <- function(data) {

	thalStarts <- read.csv("dmel_sizes.txt",sep="\t", header=T, as.is=T);
	rownames(thalStarts) <- as.vector(thalStarts$chr)
#	thalStarts <- thalStarts[1:6,];

	lyrStarts  <- read.csv("contig_sizes.txt",sep="\t", header=T, as.is=T);
	rownames(lyrStarts)  <- as.vector(lyrStarts$contig)

#   ### reorder lyrStarts relative to thalStarts based on data
   data$chr_1 <- factor(data$chr_1, levels=unique(thalStarts$chr))
   lyrChrORDER <- unique(data[order(data$chr_1, data$start_1),]$chr_2)

#   lyrChrORDER <- c()
#   for (chr in thalStarts$chr) {
#      chrData <- data[data$chr_1 == chr,]
#      lyrChrORDER <- union(lyrChrORDER, unique(chrData[order(chrData$start_1),]$chr_2))
#   }

   lyrStarts <- lyrStarts[union(lyrChrORDER, lyrStarts$contig),]
   lyrStarts$plotStart <- cumsum(c(0,lyrStarts[-nrow(lyrStarts),]$length))

#   tackOn_1 <- as.vector(unlist(lapply(data$chr_1, function(x) { thalStarts[as.vector(x),]$plotStart})))
#   tackOn_2 <- as.vector(unlist(lapply(data$chr_2, function(x) { lyrStarts[as.vector(x),]$plotStart})))
   tackOn_1 <- thalStarts[as.vector(data$chr_1),]$plotStart
   tackOn_2 <- lyrStarts[as.vector(data$chr_2),]$plotStart

	data$pos_start_1	<- data[,"start_1"] + tackOn_1;
	data$pos_end_1		<- data[,"end_1"]   + tackOn_1;
	data$pos_start_2	<- data[,"start_2"] + tackOn_2;
	data$pos_end_2		<- data[,"end_2"]   + tackOn_2;

	rm(thalStarts,lyrStarts,tackOn_1,tackOn_2);
	return(data);
}
