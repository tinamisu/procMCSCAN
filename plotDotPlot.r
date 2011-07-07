############################
### dotplot.r 
### by Tina
### created on 07.20.08
############################
plotDotPlot <- function(data,bounds,species1,species2) {
#   contigSize <- read.csv("contig_sizes.txt",header=T,sep="\t")
#   contigSize$plotStart <- cumsum(c(0,contigSize[-nrow(contigSize),]$length))
#   write.table(file="contig_sizes.txt",contigSize,col.names=T,row.names=F,quote=F,sep="\t")
#   dmelSize <- read.csv("dmel_sizes.txt",header=T,sep="\t")
#   dmelSize$plotStart <- cumsum(c(0,dmelSize[-nrow(dmelSize),]$length))
#   write.table(file="dmel_sizes.txt",dmelSize,col.names=T,row.names=F,quote=F,sep="\t")

	thalStarts <- read.csv("dmel_sizes.txt",header=T,sep="\t",as.is=T);
	thalStarts <- thalStarts[1:6,];
   thalStarts$end <- thalStarts$plotStart + thalStarts$length;

	### put all together onto one plot
	blocks <- unique(as.vector(data$block_id));
	colors <- c("red","chartreuse","blue","darkorange","purple","gold","steelblue","orchid","darkolive"); ### COLORS FOR BLOCKS

	### START PLOTTING
	############################
	pdf(file=paste("matchBlocks.",species1,"2",species2,".all.pdf",sep=""),width=15,height=15);
	par(bg="transparent",col.main="black",col.lab="black",cex.lab=.8,cex.axis=.5,mar=c(3.5,3.5,1.5,3),xpd=F,lend=1,xaxs="i",yaxs="i");

	x_plotRange <- c(0,max(thalStarts$end));
	y_plotRange <- c(0,max(c(data$pos_start_2,data$pos_end_2)));
	plot(c(0,0),col="transparent",xlab="",ylab="",xlim=x_plotRange,ylim=y_plotRange,axes=F);
   box(col="gray68");
	mtext(paste(species1,"pos"),side=1,line=2,font=1);
	mtext(paste(species2,"pos"),side=2,line=2,font=1);
	mtext(paste("num blocks",length(blocks)),side=3,line=0.1,font=2,cex=1);
#	if (species1 == "dmel") { abline(v=thalStarts$end,lwd=2,col="black"); } else if (species1 == "dsim") { abline(v=lyrStarts$end,lwd=2,col="black"); }
#	if (species2 == "dmel") { abline(h=thalStarts$end,lwd=2,col="black"); } else if (species2 == "dsim") { abline(h=lyrStarts$end,lwd=2,col="black"); }
   abline(v=thalStarts$end,lwd=2,col="black");

   rect(bounds$pos_start_1,y_plotRange[1]-5000,bounds$pos_end_1,y_plotRange[1],border="transparent",col="yellow",xpd=T)
   mtext(thalStarts$chr,side=1,line=.5,at=thalStarts$plotStart,font=2,col="black");

   segments(bounds$pos_start_1,bounds$pos_start_2,bounds$pos_end_1,bounds$pos_end_2,col=c(1:nrow(bounds)),xpd=T,lwd=3)
#   segments(data$pos_start_1,data$pos_start_2,data$pos_end_1,data$pos_end_2,xpd=T,lwd=1)
	dev.off();

if (0) {
	############################
	pdf(file=paste("matchBlocks.",species1,"2",species2,".byChr.pdf",sep=""),width=5,height=3.8);
	par(bg="transparent",col.main="black",col.lab="black",cex.lab=.8,cex.axis=.5,mar=c(3.5,3.5,1,5),xpd=F,lend=1,xaxs="i",yaxs="i");

	for (chr1 in thalStarts$chr) {
#	for (chr1 in sort(unique(data$chr_1))) {
		for (chr2 in sort(unique(data$chr_2))) {
			subdata  <- data[data$chr_1==chr1 & data$chr_2==chr2,];
         subdata1not2 <- data[data$chr_1==chr1 & data$chr_2!=chr2,]
         subdata2not1 <- data[data$chr_2==chr2 & data$chr_1!=chr1,]

			blocks <- unique(as.vector(subdata$block_id));

			if(length(blocks)) {
				x_plotRange <- range(c(subdata$pos_start_1,subdata$pos_end_1));
				y_plotRange <- range(c(subdata$pos_start_2,subdata$pos_end_2));

				plot(c(0,0),col="transparent",xlab="",ylab="",xlim=x_plotRange,ylim=y_plotRange,axes=F);
            axis(1)
            axis(2)
				mtext(paste(species1," ",chr1," (",abs(x_plotRange[1]-x_plotRange[2]),": ",x_plotRange[1],"-",x_plotRange[2],")",sep=""),side=1,line=2,font=1,cex=.68);
				mtext(paste(chr2," (",abs(y_plotRange[1]-y_plotRange[2]),": ",y_plotRange[1],"-",y_plotRange[2],")",sep="",collapse="\n"),side=2,line=2,font=1,cex=.38);
				mtext(paste("num blocks",length(blocks)),side=3,line=0.05,font=2,cex=.8);

            ### plot external blocks in gray
            if (nrow(subdata1not2)>0) {
               rect(subdata1not2$pos_start_1,y_plotRange[1],subdata1not2$pos_end_1,y_plotRange[2],col=rgb(t(col2rgb("gray28")/255),alpha=.5),bor="transparent");
            }

            if (nrow(subdata2not1)>0) {
               rect(x_plotRange[1],subdata2not1$pos_start_2,x_plotRange[2],subdata2not1$pos_end_2,col=rgb(t(col2rgb("gray28")/255),alpha=.5),bor="transparent");
            }

				counter <- 1;
				for (b in blocks) {
					block_bounds <- bounds[bounds$block_id==b,];
					block_data <- subdata[subdata$block_id==b,];
					col4points <- colors[counter %% length(colors)];
					col4lines <- paste(colors[counter %% length(colors)],"4",sep="");

					first_entry <- block_data[1,];
					last_entry <- block_data[dim(block_data)[1],];

#					points(block_data[,c("pos_start_1","pos_start_2")],pch="|",cex=.8,col=col4points);
					segments(block_data$pos_start_1,block_data$pos_start_2,block_data$pos_end_1,block_data$pos_end_2,lwd=5,col=col4points);

					x_seg_pos <- y_plotRange[1];
					y_seg_pos <- x_plotRange[1];

					points(c(first_entry[,"pos_start_1"],last_entry[,"pos_start_2"]),rep(x_seg_pos,2),pch="|",cex=.6,col=col4points,xpd=T);
					points(rep(y_seg_pos,2),c(first_entry[,"pos_start_2"],last_entry[,"pos_end_2"]),pch="-",cex=.8,col=col4points,xpd=T);

					block_bounds <- bounds[bounds$block_id==b,];
					segments(block_bounds[,"pos_start_1"],block_bounds[,"pos_start_2"],block_bounds[,"pos_end_1"],block_bounds[,"pos_end_2"],
						lwd=2,xpd=T,col=col4lines); #col="black",
					segments(first_entry[,"pos_start_1"],x_seg_pos,last_entry[,"pos_end_1"],x_seg_pos,lwd=3,col=col4lines,xpd=T);
					segments(y_seg_pos,first_entry[,"pos_start_2"],y_seg_pos,last_entry[,"pos_end_2"],lwd=3,col=col4lines,xpd=T);
					
					label_entry <- first_entry;
					if (first_entry[,"pos_start_2"]>last_entry[,"pos_end_2"]) { label_entry <- last_entry; }
					text(max(x_plotRange),label_entry[,c("pos_start_2"),],labels=paste(" blk ",b," (",dim(block_data)[1],")",sep=""),
						adj=c(0,.5),offset=2,font=1,col=col4lines,xpd=T,cex=.8);

					counter <- counter+1;
				}
			}
		}
	}
	dev.off();
}

}
