#usage $0 minValue maxValue chip.peaks control.peaks outputFile.png output.txt
args <- commandArgs(TRUE)
minValue <- type.convert(args[2])
maxValue <- type.convert(args[3])
dataTable <-suppressWarnings(read.table(args[4], header=FALSE));
chip<-data.frame(dataTable)
dataTable <-suppressWarnings(read.table(args[5], header=FALSE));
control<-data.frame(dataTable)
x <-c((minValue-0.5):(maxValue+0.5))
breaks <- c(0,x,1000)
controlHist <-hist(control$V6, breaks=breaks, right=FALSE, plot=FALSE )
chipHist <-hist(chip$V6, breaks=breaks, right=FALSE, plot=FALSE )

ifPDF <- 0
if (length(args)>=8) {
	ifPDF=args[8]
}
if (ifPDF==1) {
	pdf(file = args[6], width = 8, height = 7, pointsize = 20, bg = "white")
} else {
	png(filename = args[6], width = 580, height = 580, units = "px", pointsize = 20, bg = "white", res = NA)
}
suppressWarnings(plot(controlHist$mids,controlHist$counts,xlab = "Peak height",xlim=c(minValue-0.5,maxValue+0.5), ylab = "Peak count",pch=17, col = colors()[328], log = "y"))
suppressWarnings(points(chipHist$mids,chipHist$counts,xlab = "peak height",xlim=c(minValue-0.5,maxValue+0.5),ylab = "peak count",pch=15, col = colors()[131], log = "y"))
legend(maxValue*0.7,y = max(chipHist$counts)*0.7, bty="n",c("ChIP","Control"), cex=1, col = c(colors()[131],colors()[328]), lty = c(-1, -1), pch = c(15, 17))

dev.off()

sink(args[7], append=FALSE, split=FALSE)
cat (paste("peak height","# peaks in ChIP","# peaks in Control","#control/chip","Cumulative #control/chip (FDR)","\n",sep='\t'))
for (xval in c(minValue:maxValue)) {
  for (i in (1:length(chipHist$mids))) {     
     if (xval==chipHist$mids[i]) {
      ychip <- chipHist$counts[i]
     }
  }
  for (i in (1:length(controlHist$mids))) {     
     if (xval==controlHist$mids[i]) {
      ycontrol <- controlHist$counts[i]
     }
  }
	ychipCum = 1;
	ycontrolCum = 1	;
  for (i in (1:length(chipHist$mids))) {     
     if (xval<=chipHist$mids[i]) {
      ychipCum <- ychipCum+chipHist$counts[i]
     }
  }
  for (i in (1:length(controlHist$mids))) {     
     if (xval<=controlHist$mids[i]) {
      ycontrolCum  <- ycontrolCum+controlHist$counts[i]
     }
  }
  cumRatio= ycontrolCum/ychipCum
 if(cumRatio>1){cumRatio = 1}
  cat (paste(xval,ychip,ycontrol,ycontrol/ychip,cumRatio,"\n",sep='\t'))
}
