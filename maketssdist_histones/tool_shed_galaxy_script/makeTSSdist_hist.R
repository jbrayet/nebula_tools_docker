#usage $0 STEP RIGHT chipPeaks outputFile.png output.txt [controlPeaks] [1 for pdf]
args <- commandArgs()
print (args)
myStep <- type.convert(args[4])
maxValue <- type.convert(args[5])

dataTable <-read.table(file=paste(args[6],".genes.ClosestPeakDist", sep=""), header=TRUE);
chip.genes.ClosestPeakDist<-data.frame(dataTable)
ifReg <- 0
if (length(unique(chip.genes.ClosestPeakDist$Reg))>1) {
 ifReg <- 1
}
ifControl <- 0


ifPDF <- 0
if (length(args)>=10) {
	ifPDF=args[10]
}
if (length(args)==9 & args[9]==1) {
	ifPDF=1
}

library(Hmisc)

if (length(args)>=9 & args[9]!=1 & args[9]!=0) {
  dataTable <-read.table(file=paste(args[9],".genes.ClosestPeakDist", sep=""), header=TRUE);
  control.genes.ClosestPeakDist<-data.frame(dataTable)
  ifControl <- 1
}
if (ifReg & ifControl) {
  if (ifPDF==1) {
	pdf(file = args[7], width = 19, height = 8, pointsize = 20, bg = "white")
  } else {
  	png(filename = args[7], width = 1440 , height = 680, units = "px", pointsize = 20, bg = "white", res = NA)
	plot(1:10)
  }
  op <- par(mfrow = c(2,3))
} else {
  if (ifPDF==1) {
	pdf(file = args[7], width = 10, height = 13, pointsize = 20, bg = "white")
  } else {
  	png(filename = args[7], width = 680, height = 880, units = "px", pointsize = 20, bg = "white", res = NA)
	plot(1:10)
  }
 # plot(1:10)
  op <- par(mfrow = c(2,1))
}
myColor <- 1
myColor[1] <- colors()[131]
myColor[2] <- "darkolivegreen3"
myColor[3] <- "azure4"
myColor[4] <- "royalblue3"
myColor[5] <- colors()[17]

myColorControl <- 1

myColorControl[1] <- colors()[24]
myColorControl[2] <- colors()[278]
myColorControl[3] <- colors()[305]
myColorControl[4] <- colors()[219]
myColorControl[5] <- colors()[343]



#for cumulative:
dist_real_f <- chip.genes.ClosestPeakDist

if (ifControl) {
 dist_control_f <- control.genes.ClosestPeakDist

}
step <- myStep
lim <- maxValue
x <- 0
count <- 1
countL <-1
n.types <- 1
myLevels <- 0
countTotalCont <- 0
countTotal <-0
countLCont <- 0
cumTotalCont <- 0


if (ifReg) {
	 n.types <- length(levels(chip.genes.ClosestPeakDist$Reg))
	 myLevels <- levels(chip.genes.ClosestPeakDist$Reg)
	 cum = matrix( 0, nrow=lim/step +1,  ncol=n.types, byrow = TRUE) 
	 for (i in c(1:n.types)) {
	 	t <- which ((dist_real_f$Reg==myLevels[i]))
	 	countL[i] <- length(t)
         }	
	count <-1 
	for (i in seq(length=lim/step +1, from=0, by=step)) {
		for (t in c(1:n.types)) {
			tt <- which ((dist_real_f$Reg==myLevels[t])&(dist_real_f$Dist<=i)&(dist_real_f$Dist>=-i))
			cum[count,t] <- length(tt)	 		
         	}		
		x[count] <- i
		count <- count + 1
	}
	ymax <- max(cum[,1]/countL[1], na.rm=TRUE)
	for (i in c(2:n.types)) {
	 	ymax <- max(ymax,max(cum[,i]/countL[i], na.rm=TRUE))
         }
	myLocCol <- myColor[2]	

	par(mar=c(5.1, 7.1, 4.1, 2.1)) 

	 plot (x,cum[,1]/countL[1] ,col = myColor[2],type="l", main="",xlab="",ylab="", lwd = 2, xlim = c(0, lim),xaxt="n" , ylim=c(0,ymax))
	 for (i in c(2:n.types)) {
	 	colorr <- i+1
		myLocCol <- c(myLocCol,myColor[colorr])
		lines (x,cum[,i]/countL[i] ,col = myColor[colorr],type="l", lwd = 2)
		print (myColor[colorr])
         }	
	
	gradi <- 1000
	if (lim>10000) {
		gradi <- 10000
	}
	if (lim<3000) {
		gradi <- 500
	}
	 axisx <- seq(length=lim/gradi+1, from=0, by=gradi)
  	 axisxlab <- paste(axisx/1000,"", sep = "")
   	 axis(1, at=axisx,labels=axisxlab , las=1, cex.axis=1)
  	 ymax <- max(cum[,i]/countL[i], na.rm=TRUE)

	  minor.tick(nx=5,tick.ratio=0.5)

	 legend(lim*0.45, ymax*0.45, myLevels,  cex=1, lwd = 2, bty = "n", col = myLocCol, lty = c(1),  pt.bg= c(myLocCol) , merge = TRUE)

	title( main="",xlab="Distance from TSS (Kb)",ylab="Proportion of genes with a peak\nat a given distance (cumulative)")
	
	 if (ifControl) {
                count <-1
 	 	n.types <- length(levels(control.genes.ClosestPeakDist$Reg))
	 	myLevels <- levels(control.genes.ClosestPeakDist$Reg)
		 cumCont = matrix( 0, nrow=lim/step +1,  ncol=n.types, byrow = TRUE) 
		 for (i in c(1:n.types)) {
		 	t <- which ((dist_control_f$Reg==myLevels[i]))
		 	countLCont[i] <- length(t)
		 }	 
		for (i in seq(length=lim/step +1, from=0, by=step)) {
			for (t in c(1:n.types)) {
				tt <- which ((dist_control_f$Reg==myLevels[t])&(dist_control_f$Dist<=i)&(dist_control_f$Dist>=-i))
				cumCont[count,t] <- length(tt)	 		
		 	}		
			x[count] <- i
			count <- count + 1
		}
		ymax <- max(cumCont[,1]/countLCont[1], na.rm=TRUE)
		for (i in c(2:n.types)) {
		 	ymax <- max(ymax,max(cumCont[,i]/countLCont[i], na.rm=TRUE))
		 }
		myLocColCntr <- myColorControl[2]	
		plot (x,cumCont[,1]/countLCont[1] ,col = myLocColCntr[1],type="l", main="",xlab="",ylab="", lwd = 2, xlim = c(0, lim),xaxt="n" , ylim=c(0,ymax))
		 for (i in c(2:n.types)) {
		 	colorr <- i+1
			myLocColCntr <- c(myLocColCntr,myColorControl[colorr])
			lines (x,cumCont[,i]/countLCont[i] ,col = myColorControl[colorr],type="l", lwd = 2)
		 }	
	   	   if (lim>10000) {
				gradi <- 10000
		   }
		   if (lim<3000) {
				gradi <- 500
		   }
		   axisx <- seq(length=lim/gradi+1, from=0, by=gradi)
	  	 axisxlab <- paste(axisx/1000, sep = "")
	   	 axis(1, at=axisx,labels=axisxlab , las=1, cex.axis=1)
  		minor.tick(nx=5,tick.ratio=0.5)
		legend(lim*0.45, ymax*0.45, myLevels,  cex=1 , lwd = 2, bty = "n", col = myLocColCntr, lty = c(1),  pt.bg= c(myLocCol) , merge = TRUE)
		title( main="",xlab="Distance from TSS (Kb)",ylab="Proportion of genes with a peak\nat a given distance (cumulative)")
		#real_vs_control_cumulative:
		count <-1
		countTotal <- length(dist_real_f$Reg)
	   	cumTotal  <- 0
		for (i in seq(length=lim/step +1, from=0, by=step)) {	
			t <- which ((dist_real_f$Dist<=i)&(dist_real_f$Dist>=-i))
			cumTotal[count] <- length(t)
			x[count] <- i
			count <- count + 1
		}
		   plot (x,cumTotal/countTotal ,col = myColor[1],type="l", main="",xlab="",ylab="", lwd = 2, xlim = c(0, lim),xaxt="n",ylim = c(0, max(cumTotal/countTotal,na.rm=TRUE)))
		   gradi <- 1000
		   if (lim>10000) {
				gradi <- 10000
		    }
		    if (lim<3000) {
				gradi <- 500
		    }
		    axisx <- seq(length=lim/gradi+1, from=0, by=gradi)
		   axisxlab <- paste(axisx/1000, sep = "")
		   axis(1, at=axisx,labels=axisxlab , las=1, cex.axis=1)
		   ymax <- max(cumTotal/countTotal, na.rm=TRUE)
    		 minor.tick(nx=5,tick.ratio=0.5)
		countTotalCont <- length(dist_control_f$Reg)
	   	cumTotalCont  <- 0
		count <- 1
		for (i in seq(length=lim/step +1, from=0, by=step)) {	
			t <- which ((dist_control_f$Dist<=i)&(dist_control_f$Dist>=-i))
			cumTotalCont[count] <- length(t)
			x[count] <- i
			count <- count + 1
		}
		lines (x,cumTotalCont/countTotalCont ,col = colors()[328],type="l", lwd = 2)
		legend(lim*0.45, ymax*0.45, c("ChIP","Control"),  cex=1 , lwd = 2, bty = "n", col = c(myColor[1], colors()[328]), lty = c(1),  pt.bg= c(myColor[1], colors()[328]) , merge = TRUE)
		title( main="",xlab="Distance from TSS (Kb)",ylab="Proportion of genes with a peak\nat a given distance (cumulative)")
	}
} else {
   countTotal <- length(dist_real_f$Reg)
   cumTotal  <- 0
   count <-1
   gradi <- 1000
   if (lim>10000) {
	gradi <- 10000
   }
   if (lim<3000) {
	gradi <- 500
   }

   for (i in seq(length=lim/step +1, from=0, by=step)) {	
	t <- which ((dist_real_f$Dist<=i)&(dist_real_f$Dist>=-i))
	cumTotal[count] <- length(t)
	x[count] <- i
	count <- count + 1
   }
   par(mar=c(5.1, 7.1, 4.1, 2.1)) 
  
   plot (x,cumTotal/countTotal ,col = myColor[1],type="l", main="",xlab="",ylab="", lwd = 2, xlim = c(0, lim),ylim = c(0, max(cumTotal/countTotal,na.rm=TRUE)),xaxt="n" )
   axisx <- seq(length=lim/gradi+1, from=0, by=gradi)
   axisxlab <- paste(axisx/1000, sep = "")
   axis(1, at=axisx,labels=axisxlab , las=1, cex.axis=1)
title( main="",xlab="Distance from TSS (Kb)",ylab="Proportion of genes with a peak\nat a given distance (cumulative)")
   ymax <- max(cumTotal/countTotal, na.rm=TRUE)

   if (ifControl) {
	countTotalCont <- length(dist_control_f$Reg)
   	cumTotalCont  <- 0
	count <- 1
	for (i in seq(length=lim/step +1, from=0, by=step)) {	
		t <- which ((dist_control_f$Dist<=i)&(dist_control_f$Dist>=-i))
		cumTotalCont[count] <- length(t)
		x[count] <- i
		count <- count + 1
	}
	lines (x,cumTotalCont/countTotalCont ,col = colors()[328],type="l", lwd = 2)
        legend(lim*0.45, ymax*0.45, c("ChIP","Control"),  cex=1 , lwd = 2, bty = "n", col = c(myColor[1], colors()[328]), lty = c(1),  pt.bg= c(myColor[1], colors()[328]) , merge = TRUE)
   } else {
        legend(lim*0.45, ymax*0.45, c("ChIP"),  cex=1 , lwd = 2, bty = "n", col = c(myColor[1]), lty = c(1),  pt.bg= c(myColor[1]) , merge = TRUE)
   }
}
sink(args[8], append=FALSE, split=FALSE)
if (ifReg) {
	if (ifControl) {
		cat (paste("Dist_TSS","% genes w/ a peak in ChIP","% genes w/ a peak in control",sep='\t'))
		cat("\t")
 		for (i in c(1:n.types)) {
		 	cat(paste("% ", myLevels[i]," genes w/ a peak in ChIP", sep=""))
                	cat("\t")
		}	

 		for (i in c(1:n.types)) {
		 	cat(paste("% ", myLevels[i]," genes w/ a peak in Control", sep=""))
			cat("\t")
		}
                cat("\n")
		for (i in c(1:length(x))) {
			cat(paste(x[i],cumTotal[i]/countTotal,cumTotalCont[i]/countTotalCont,sep="\t"))
                	cat("\t")
	 		for (t in c(1:n.types)) {
			 	cat(paste(cum[i,t]/countL[t],"\t", sep=""))
			}
	 		for (t in c(1:n.types)) {
			 	cat(paste(cumCont[i,t]/countLCont[t],"\t", sep=""))
			}	 		
                	cat("\n")
		}
	}else {
            cat (paste("Dist_TSS","\t",sep=''))
 	    for (i in c(1:n.types)) {
		 cat(paste("% ", myLevels[i]," genes w/ a peak in ChIP", "\t", sep=""))
	    }
 	    cat("\n")
	    for (i in c(1:length(x))) {
			cat(paste(x[i],"\t",sep=""))	 		
	 		for (t in c(1:n.types)) {
			 	cat(paste(cum[i,t]/countL[t],"\t", sep=""))
			}
                	cat("\n")
	    }
	}
} else {
	if (ifControl) {
		cat (paste("Dist_TSS","% genes w/ a peak in ChIP","% genes w/ a peak in control",sep='\t'))
                cat("\n")
		for (i in c(1:length(x))) {
			cat(paste(x[i],cumTotal[i]/countTotal,cumTotalCont[i]/countTotalCont,sep="\t"))
                	cat("\n")
		}
	}else {
		cat (paste("Dist_TSS","% genes w/ a peak in ChIP",sep='\t'))
                cat("\n")
		for (i in c(1:length(x))) {
			cat(paste(x[i],cumTotal[i]/countTotal,sep="\t"))
                	cat("\n")
		}

	}
}
#stop sinking:
sink() 

#around TSS:
lim <- maxValue
step <- myStep
xlabs <-seq(from = -lim, to = lim, by = step) 
chip.genes <- read.table(file=paste(args[6],".genes", sep=""), header=TRUE) ;
dist_real_f <- chip.genes
if (ifControl) {
   control.genes <- read.table(file=paste(args[9],".genes", sep=""), header=TRUE) ;   
   dist_control_f<-data.frame(control.genes)
}
if (ifReg) {
	#n.types <- length(levels(chip.genes.ClosestPeakDist$Reg))
	 #myLevels <- levels(dist_real_f$Reg)

	t<- which (dist_real_f$Reg==myLevels[1])
	values_real <-dist_real_f$Dist[t]

	values_Start <-dist_real_f$DistS[t]
	values_End <-dist_real_f$DistE[t]
  
	density <-NULL
	for (i in seq(from=-lim,to=lim,by=step)) {
	  density <- c(density,length(which(values_End>=i&values_Start<=i)))   
	}
	
	density<-density/length(unique(dist_real_f$GeneCoord[t]))
	  
  ymax <- max(density, na.rm=TRUE)
  
	for (i in c(2:n.types)) {
		t<- which (dist_real_f$Reg==myLevels[i])
		values_Start <-dist_real_f$DistS[t]
		values_End <-dist_real_f$DistE[t]
		density <-NULL
		for (i in seq(from=-lim,to=lim,by=step)) {
		  density <- c(density,length(which(values_End>=i&values_Start<=i)))   
		}
		density<-density/length(unique(dist_real_f$GeneCoord[t]))		
    ymax <- max(ymax,max(density, na.rm=TRUE))
  }


	t<- which (dist_real_f$Reg==myLevels[1])
	values_Start <-dist_real_f$DistS[t]
	values_End <-dist_real_f$DistE[t]
	density <-NULL
	for (i in seq(from=-lim,to=lim,by=step)) {
	  density <- c(density,length(which(values_End>=i&values_Start<=i)))   
	}
	density<-density/length(unique(dist_real_f$GeneCoord[t]))		
	
  plot (xlabs,density,col = myLocCol[1],type="l", main="",xlab="",ylab="", lwd = 2, xlim = c(-lim, lim),ylim = c(0, ymax), xaxt="n" ) 
	
	title( main="",xlab="Distance from TSS (Kb)",ylab="Proportion of genes with a peak\ncovering this area (density)")

	for (i in c(2:n.types)) {
		t<- which (dist_real_f$Reg==myLevels[i])
		values_Start <-dist_real_f$DistS[t]
		values_End <-dist_real_f$DistE[t]
		density <-NULL
		for (j in seq(from=-lim,to=lim,by=step)) {
		  density <- c(density,length(which(values_End>=j&values_Start<=j)))   
		}
		density<-density/length(unique(dist_real_f$GeneCoord[t]))		
  		lines (xlabs,density,col = myLocCol[i],type="l", lwd = 2)
  }
	legend(lim*0.1, ymax*0.9, myLevels,  cex=1 , lwd = 2, bty = "n", col = myLocCol, lty = c(1),  pt.bg= c(myLocCol) , merge = TRUE)
		
	   gradi <- 1000
	   if (lim>10000) {
		gradi <- 10000
	   }
	   if (lim<3000) {
		gradi <- 500
	   }

	   axisx <- seq(length=2*lim/gradi+1, from=-lim, by=gradi)
      	   axisxlab <- paste(axisx/1000, sep = "")
	   axis(1, at=axisx,labels=axisxlab , las=1, cex.axis=1)


  	#minor.tick(nx=10,tick.ratio=0.5)
       if (ifControl) {
		t<- which (dist_control_f$Reg==myLevels[1])
    
		values_Start <-dist_control_f$DistS[t]
		values_End <-dist_control_f$DistE[t]    
		density <-NULL
		for (i in seq(from=-lim,to=lim,by=step)) {
		  density <- c(density,length(which(values_End>=i&values_Start<=i)))   
		}
		density<-density/length(unique(dist_control_f$GeneCoord[t]))		
    ymax <- max(density, na.rm=TRUE)
		for (i in c(2:n.types)) {
			t<- which (dist_control_f$Reg==myLevels[i])
			values_Start <-dist_control_f$DistS[t]
			values_End <-dist_control_f$DistE[t]    
			density <-NULL
			for (j in seq(from=-lim,to=lim,by=step)) {
			  density <- c(density,length(which(values_End>=j&values_Start<=j)))   
			}
			density<-density/length(unique(dist_control_f$GeneCoord[t]))		
			ymax <- max(ymax,max(density, na.rm=TRUE))
		}
		t<- which (dist_control_f$Reg==myLevels[1])
		values_Start <-dist_control_f$DistS[t]
		values_End <-dist_control_f$DistE[t]    
		density <-NULL
		for (i in seq(from=-lim,to=lim,by=step)) {
		  density <- c(density,length(which(values_End>=i&values_Start<=i)))   
		}
		density<-density/length(unique(dist_control_f$GeneCoord[t]))	
		plot (xlabs,density,col = myLocColCntr[1],type="l", main="",xlab="",ylab="", lwd = 2, xlim = c(-lim, lim),ylim = c(0, ymax),xaxt="n" ) 
		title( main="",xlab="Distance from TSS (Kb)",ylab="Proportion of genes with a peak\ncovering this area (density)")		
		for (i in c(2:n.types)) {
			t<- which (dist_control_f$Reg==myLevels[i])
			values_Start <-dist_control_f$DistS[t]
			values_End <-dist_control_f$DistE[t]    
			density <-NULL
			for (j in seq(from=-lim,to=lim,by=step)) {
			  density <- c(density,length(which(values_End>=j&values_Start<=j)))   
			}
			density<-density/length(unique(dist_control_f$GeneCoord[t]))	
	  		lines (xlabs,density,col = myLocColCntr[i],type="l", lwd = 2)
		}

		gradi <- 1000
		   if (lim>10000) {
			gradi <- 10000
		   }
		   if (lim<3000) {
			gradi <- 500
		   }

		   axisx <- seq(length=2*lim/gradi+1, from=-lim, by=gradi)
	      	   axisxlab <- paste(axisx/1000, sep = "")
		   axis(1, at=axisx,labels=axisxlab , las=1, cex.axis=1)

		legend(lim*0.1, ymax*0.9, myLevels,  cex=1 , lwd = 2, bty = "n", col = myLocColCntr, lty = c(1),  pt.bg= c(myLocCol) , merge = TRUE)
		
	  	# minor.tick(nx=10,tick.ratio=0.5)
		#control vs real
    
		values_Start <-dist_real_f$DistS
		values_End <-dist_real_f$DistE   
		density <-NULL
		for (i in seq(from=-lim,to=lim,by=step)) {
		  density <- c(density,length(which(values_End>=i&values_Start<=i)))   
		}
		density<-density/length(unique(dist_real_f$GeneCoord))	
    
    
		 plot (xlabs,density,col = myColor[1],type="l", main="",xlab="",ylab="", lwd = 2, xlim = c(-lim, lim),xaxt="n",ylim=c(0,max(density,na.rm=TRUE))) 
		title( main="",xlab="Distance from TSS (Kb)",ylab="Proportion of genes with a peak\ncovering this area (density)")		
		 ymax <- max(density, na.rm=TRUE)
    
    
		values_Start <-dist_control_f$DistS
		values_End <-dist_control_f$DistE   
		density <-NULL
		for (i in seq(from=-lim,to=lim,by=step)) {
		  density <- c(density,length(which(values_End>=i&values_Start<=i)))   
		}
		density<-density/length(unique(dist_control_f$GeneCoord))    
    
		  lines (xlabs,density,col = colors()[328],type="l", lwd = 2)
		  legend(lim*0.2, ymax*0.9, c("ChIP","Control"),  cex=1 , lwd = 2, bty = "n", col = c(myColor[1], colors()[328]), lty = c(1),  pt.bg= c(myColor[1], colors()[328]) , merge = TRUE)		 
		
		   gradi <- 1000
		   if (lim>10000) {
			gradi <- 10000
		   }
		   if (lim<3000) {
			gradi <- 500
		   }

		   axisx <- seq(length=2*lim/gradi+1, from=-lim, by=gradi)
	      	   axisxlab <- paste(axisx/1000, sep = "")
		   axis(1, at=axisx,labels=axisxlab , las=1, cex.axis=1)



		 # minor.tick(nx=10,tick.ratio=0.5)
       }
} else {
 #values_real <-dist_real_f$Dist
 values_Start <-dist_real_f$DistS
 values_End <-dist_real_f$DistE
 density <-NULL
 xlabs <-seq(from = -lim, to = lim, by = step) 
 for (i in seq(from=-lim,to=lim,by=step)) {
   density <- c(density,length(which(values_End>=i&values_Start<=i)))   
 }
  
 density<-density/length(unique(dist_real_f$GeneCoord))
 

 plot (xlabs,density,col = myColor[1],type="l", main="",xlab="",ylab="", lwd = 2,xaxt="n",ylim=c(0,max(density,na.rm=TRUE))) 
		title( main="",xlab="Distance from TSS (Kb)",ylab="Proportion of genes with a peak\ncovering this area (density)")		
 ymax <- max(density, na.rm=TRUE)
 if (ifControl) {
  
  values_Start <-dist_control_f$DistS
  values_End <-dist_control_f$DistE
  density <-NULL
  for (i in seq(from=-lim,to=lim,by=step)) {
    density <- c(density,length(which(values_End>=i&values_Start<=i)))   
  }
  density<-density/length(unique(dist_control_f$GeneCoord))
  
  
  lines (xlabs,density_c,col = colors()[328],type="l", lwd = 2)
  legend(lim*0.2, ymax*0.9, c("ChIP","Control"),  cex=1 , lwd = 2, bty = "n", col = c(myColor[1], colors()[328]), lty = c(1),  pt.bg= c(myColor[1], colors()[328]) , merge = TRUE)
 } else {
  legend(lim*0.2, ymax*0.9, c("ChIP"),  cex=1 , lwd = 2, bty = "n", col = c(myColor[1]), lty = c(1),  pt.bg= c(myColor[1]) , merge = TRUE)
 }

 gradi <- 1000
 if (lim>10000) {
    gradi <- 10000
  }
  if (lim<3000) {
	gradi <- 500
  }

  axisx <- seq(length=2*lim/gradi+1, from=-lim, by=gradi)
  axisxlab <- paste(axisx/1000, sep = "")
  axis(1, at=axisx,labels=axisxlab , las=1, cex.axis=1)


 # minor.tick(nx=10,tick.ratio=0.5)
}
dev.off()
q();
cat (paste("peak height","# peaks in ChIP","# peaks in Control","#control/chip","\n",sep='\t'))
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
  cat (paste(xval,ychip,ycontrol,ycontrol/ychip,"\n",sep='\t'))
}
