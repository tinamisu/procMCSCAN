############################
### dotplot.r 
### by Tina
### created on 07.20.08
############################
readData <- 1;
if (readData == 1) {
	### READ IN DATA FROM MCSCAN
	########################################################
	#thalStarts <- read.csv("dmel_sizes.txt",sep="\t",header=T,as.is=T);
	#thalStarts <- thalStarts[1:6,];
   #lyrStarts <- read.csv("contig_sizes.txt",sep="\t",header=T,as.is=T);

	thal2lyr <-read.csv("contigs2dmel.aligns.postProcess.csv",header=T,as.is=T)#,nrows=5);
	thal2lyr_bounds <- read.csv("contigs2dmel.aligns.postProcess.boundaries.csv",header=T,as.is=T); 

#   ### REORDER THE SORTING by chromosome order as specified in thalStarts
#   thal2lyr$chr_1 <- factor(thal2lyr$chr_1, levels=unique(thalStarts$chr))
#   thal2lyr <- thal2lyr[order(thal2lyr$chr_1, thal2lyr$start_1),]

#	thal2lyr_file <- "contigs2dmel.blocks.csv";
#	print_order <- c("block_id","num_in_block",#"crossesAnother",
#		"chr_1","id_1_start","id_1_end","pos_start_1","pos_end_1","start_1","end_1","plot_start_pos_1","plot_end_pos_1",
#		"chr_2","id_2_start","id_2_end","pos_start_2","pos_end_2","start_2","end_2","plot_start_pos_2","plot_end_pos_2");
	
	print(paste("thal2ly",dim(thal2lyr)[1]));
	print(paste("thal2ly_bounds",dim(thal2lyr_bounds)[1]));
	
	### ADJUST POSITIONS --- order positions on contigs relative to their hit on dmel
	########################################################
	source("procMCSCAN/adjustPos.r");
	thal2lyr <- adjustPos(thal2lyr);
	print(".....adjusted thal2lyr");
	
#	source("adjustPosBounds.r");
	thal2lyr_bounds <- adjustPos(thal2lyr_bounds);
	#thal2lyr_bounds <- adjustPosBounds(thal2lyr_bounds,"thal","lyr");
	print(".....adjusted thal2lyr_bounds");
}   
	
if(1) {
	### PLOT MCSCAN RESULTS
	########################################################
	print("making plotDotplot.r");
	source("procMCSCAN/plotDotPlot.r");
	plotDotPlot(thal2lyr,thal2lyr_bounds,"dmel","contigs");

# 	source("plotDotPlot_global.r");
# 	plotDotPlot_global(thal2lyr,thal2lyr_bounds,"thal","lyr");
# 	
# 	
# #	### WRITE OUT RESULTS
# #	########################################################
# #	### NOW CHECK IF ANY BLOCKS INTERSECT OR CONTAINED WITHIN ONE ANOTHER
# #	#source("crossingBlocks.r");
# #	#thal2lyr_bounds <- crossingBlocks(thal2lyr_bounds);
# #	##thal2thal_bounds <- crossingBlocks(thal2thal_bounds);
# #	##lyr2lyr_bounds <- crossingBlocks(lyr2lyr_bounds);
# #	
# #	write.table(thal2lyr_bounds[,print_order],file=thal2lyr_file,append=F,col.names=T,row.name=F,quote=F,sep=",");
# #	#write.table(thal2thal_bounds[,print_order],file="thal2thal.blocks.csv",append=F,col.names=T,row.name=F,quote=F,sep=",");
# #	#write.table(lyr2lyr_bounds[,print_order]  ,file="lyr2lyr.blocks.csv",	append=F,col.names=T,row.name=F,quote=F,sep=",");
}


if (0) {
   #############################
   #### TESTING IF BOUNDARIES ARE SELECTED CORRECTLY
   plotCheckBlock <- function(check_block,thal2lyr,thal2lyr_bounds,check_block2=check_block) {
   	block2check<-thal2lyr[thal2lyr$block_id==check_block,]
   	block2check2<-thal2lyr[thal2lyr$block_id==check_block2,]
   	block2check_bounds<-thal2lyr_bounds[thal2lyr_bounds$block_id==check_block,]
   	block2check_bounds2<-thal2lyr_bounds[thal2lyr_bounds$block_id==check_block2,]
   
   	xrange <- range(c(block2check$pos_start_1,block2check$pos_end_1,block2check2$pos_start_1,block2check2$pos_end_1));
   	yrange <- range(c(block2check$pos_start_2,block2check$pos_end_2,block2check2$pos_start_2,block2check2$pos_end_2));
   	plot(c(block2check$pos_start_1,block2check$pos_end_1),c(block2check$pos_start_2,block2check$pos_end_2),cex=1.68,
   		xlim=c(xrange[1]-10^6,xrange[2]+10^6),
   		ylim=c(yrange[1]-10^6,yrange[2]+10^6),
   		main=paste("blue block",check_block,"red block",check_block2),
   		xlab=paste(range(c(block2check$pos_start_1,block2check$pos_end_1))),
   		ylab=paste(range(c(block2check$pos_start_2,block2check$pos_end_2))))
   	segments(block2check_bounds$pos_start_1,block2check_bounds$pos_start_2,block2check_bounds$pos_end_1,block2check_bounds$pos_end_2,lwd=20,lty=1,col="blue")
   	segments(block2check_bounds2$pos_start_1,block2check_bounds2$pos_start_2,block2check_bounds2$pos_end_1,block2check_bounds2$pos_end_2,lwd=20,lty=1,col="red")
   	segments(block2check$pos_start_1,block2check$pos_start_2,block2check$pos_end_1,block2check$pos_end_2,col="gray28",lwd=3)
   	points(c(thal2lyr$pos_start_1,thal2lyr$pos_end_1),c(thal2lyr$pos_start_2,thal2lyr$pos_end_2),pch=16,cex=.8,col="gray38");
   	segments(thal2lyr$pos_start_1,thal2lyr$pos_start_2,thal2lyr$pos_end_1,thal2lyr$pos_end_2,lwd=3,col="gray68");
   }
   
   
   source("breakupBlocks.r");
   check_block <- 135; ### block_to_break
   check_block2 <- 138; ### breakup_block
   plotCheckBlock(check_block,thal2lyr,thal2lyr_bounds,check_block2);
   breakupBlocks(thal2lyr,check_block,check_block2);
}
