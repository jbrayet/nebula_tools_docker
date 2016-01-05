args <- commandArgs(TRUE)
input <- args[2]
pngFile <- args[3]
dataTable <-read.table(file=input, header=TRUE);
chip.data<-data.frame(dataTable)
ifReg <- 0
if (length(unique(chip.data$Reg))>1) {
 ifReg <- 1
}

ifPDF <- 0
bootstrap <- 1

if (length(args)>=6) {
	ifPDF=args[6]
	bootstrap=args[7]
}
if (length(args)==5 & args[5]==1) {
	ifPDF=1
}

ifControl <- 0
if (length(args)>=5 & args[5]!=1 & args[5]!=0) {
  dataTable <-read.table(file=args[5], header=TRUE);
  control.data<-data.frame(dataTable)
  ifControl <- 1


}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
      if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
      stop("vectors must be same length")
      arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}



logFile <- args[4]
sink(logFile, append=FALSE, split=FALSE)

if (ifReg & ifControl) {

  if (ifPDF==1) {
       pdf(file = pngFile, width = 14, height = 13, pointsize = 20, bg = "white")
  } else {
       png(filename = pngFile, width =  1140, height =  840, units = "px", pointsize = 20, bg = "white", res = NA)
       plot(1:10)
  }
  op <- par(mfrow = c(3,1))
} else {
	if (ifReg | ifControl) {

	  if (ifPDF==1) {
		pdf(file = pngFile, width = 20, height = 17, pointsize = 20, bg = "white")
	  } else {
		png(filename = pngFile, width = 1580, height = 1230, units = "px", pointsize = 20, bg = "white", res = NA)
		plot(1:10)
	  }
	  op <- par(mfrow = c(2,1))
	} else {
	  if (ifPDF==1) {
		pdf(file = pngFile, width = 22, height = 8, pointsize = 20, bg = "white")
	  } else {
		png(filename = pngFile, width = 1580, height = 530, units = "px", pointsize = 20, bg = "white", res = NA)
		plot(1:10)
	  }
	  op <- par(mfrow = c(1,1))
	}
}
myColor <- 1
myColor[1] <- colors()[131]
myColor[2] <- colors()[59]
myColor[3] <- colors()[76]
myColor[4] <- colors()[88]
myColor[5] <- colors()[17]
myColor[6] <- colors()[565]
myColor[7] <- colors()[454]
myColor[8] <- colors()[401]
myColor[9] <- colors()[99]
myColorControl <- 1
myColorControl[1] <- colors()[24]
myColorControl[2] <- colors()[278]
myColorControl[3] <- colors()[305]
myColorControl[4] <- colors()[219]
myColorControl[5] <- colors()[343]
myColorControl[6] <- colors()[245]
myLevels <- 0
nn <- colnames(chip.data)
ifk27 <- 0
if (length(colnames(chip.data))==41) {
	cc <- c(1:41)
	colnames(chip.data) <- cc
}

if (length(colnames(chip.data))==45) {
	cc <- c(1:45)
	colnames(chip.data) <- cc
	ifk27 <- 1
}

if (!ifk27) {

	countTotal <- length(unique(chip.data$"1"))
	tt <- which(chip.data$"16">0)
	countTotalEnh <- length(unique(chip.data$"1"[tt]))
	tt <- which(chip.data$"10">0)
	countTotalProm <- length(unique(chip.data$"1"[tt]))
	tt <- which(chip.data$"12">0)
	countTotalImDown <- length(unique(chip.data$"1"[tt]))
	tt <- which(chip.data$"18">0)
	countTotalIntra <- length(unique(chip.data$"1"[tt]))
	tt <- which(chip.data$"20">0)
	countTotal5kbD <- length(unique(chip.data$"1"[tt]))
	tt <- which(chip.data$"28">0)
	countTotal1IntronI <- length(unique(chip.data$"1"[tt]))
	tt <- which(chip.data$"32">0)
	countTotalExonsI <- length(unique(chip.data$"1"[tt]))
	tt <- which(chip.data$"36">0)
	countTotalIntronsI <- length(unique(chip.data$"1"[tt]))
	tt <- which(chip.data$"40">0)
	countTotalJonctionsI <- length(unique(chip.data$"1"[tt]))

	names <- c("Enh.","Prom.","Imm.Down.","Intrag.","GeneDown.","F.Intron","Exons","2,3,etc.Introns","E.I.Junctions")
	yChIPTotal <- c (countTotalEnh/countTotal,countTotalProm/countTotal, countTotalImDown/countTotal,countTotalIntra/countTotal,countTotal5kbD/countTotal,countTotal1IntronI/countTotal,countTotalExonsI/countTotal,countTotalIntronsI/countTotal,countTotalJonctionsI/countTotal)

} else {
	countTotal <- length(unique(chip.data$"1"))
	tt <- which(chip.data$"42">0)
	countTSS <- length(unique(chip.data$"1"[tt]))
	tt <- which(chip.data$"44">0)
	countGene <- length(unique(chip.data$"1"[tt]))	
	names <- c("gene TSS","gene Body")
	yChIPTotal <- c (countTSS/countTotal,countGene/countTotal)
}

if(!ifReg&!ifControl) {
        par(mar=c(5.1, 7.1, 4.1, 2.1)) 
	barplot(yChIPTotal,xlab="",beside=TRUE, col=c(myColor), names.arg=c(names),ylab="Proportion of genes with a peak")
	cat ("Proportion of genes with a peak in a given genomic region:\n")
	cat (paste(c(names),sep='\t'))
	cat("\n")
	cat (paste(c(yChIPTotal),sep='\t'))
	cat("\n\n")
}
if (ifControl) {
	if (bootstrap>1) {
		yControlTotalMatrix <- NULL	
	}
	for (fileNumber in 1:bootstrap) {

		if (fileNumber>=2) {
			dataTable <-read.table(file=paste(args[5],fileNumber,sep=""), header=TRUE);
 			control.data<-data.frame(dataTable)
		}

                colnames(control.data) <- cc

		if (!ifk27) {

			countTotalCntr <- length(unique(control.data$"1"))
			tt <- which(control.data$"16">0)
			countTotalEnhCntr <- length(unique(control.data$"1"[tt]))
			tt <- which(control.data$"10">0)
			countTotalPromCntr <- length(unique(control.data$"1"[tt]))
			tt <- which(control.data$"12">0)
			countTotalImDownCntr <- length(unique(control.data$"1"[tt]))
			tt <- which(control.data$"18">0)
			countTotalIntraCntr <- length(unique(control.data$"1"[tt]))
			tt <- which(control.data$"20">0)
			countTotal5kbDCntr <- length(unique(control.data$"1"[tt]))
			tt <- which(control.data$"28">0)
			countTotal1IntronICntr <- length(unique(control.data$"1"[tt]))
			tt <- which(control.data$"32">0)
			countTotalExonsICntr <- length(unique(control.data$"1"[tt]))
			tt <- which(control.data$"36">0)
			countTotalIntronsICntr <- length(unique(control.data$"1"[tt]))
			tt <- which(control.data$"40">0)
			countTotalJonctionsICntr <- length(unique(control.data$"1"[tt]))
			yControlTotal <- c (countTotalEnhCntr/countTotalCntr,countTotalPromCntr/countTotalCntr, countTotalImDownCntr/countTotalCntr,countTotalIntraCntr/countTotalCntr,countTotal5kbDCntr/countTotalCntr,countTotal1IntronICntr/countTotalCntr,countTotalExonsICntr/countTotalCntr,countTotalIntronsICntr/countTotalCntr,countTotalJonctionsICntr/countTotalCntr)	

		} else {
			countTotalCntr <- length(unique(control.data$"1"))
			tt <- which(control.data$"42">0)
			countTSSCntr <- length(unique(control.data$"1"[tt]))
			tt <- which(control.data$"44">0)
			countGeneCntr <- length(unique(control.data$"1"[tt]))
			yControlTotal <- c (countTSSCntr/countTotalCntr,countGeneCntr/countTotalCntr)	

		}
		if (bootstrap>1) {
			yControlTotalMatrix <- rbind(yControlTotalMatrix,yControlTotal)	
		}
	}
	if (bootstrap>1) {
		yControlTotal <- colMeans(yControlTotalMatrix)
		Nrows <- nrow(yControlTotalMatrix)
		colVars <-  Nrows/(Nrows-1) * (colMeans(yControlTotalMatrix*yControlTotalMatrix)-colMeans(yControlTotalMatrix)^2)
	} 
	if (!ifReg) {
	 	cum = matrix( 0, nrow=2,  ncol=length(names), byrow = TRUE) 
		for (i in c(1:length(names))) {	
		  cum[1,i] <- yChIPTotal[i]
		  cum[2,i] <- yControlTotal[i]
		}
		if (bootstrap>1) {
			wiskers <- matrix(c(colVars-colVars,sqrt(colVars)),2,length(names),byrow=TRUE)*1.96
		} 
		par(mar=c(5.1, 7.1, 4.1, 2.1)) 
		barx <- barplot(cum,xlab="",beside=TRUE, col=c(myColor[6],myColor[5]),  legend = c("ChIP","Control"), names.arg=c(names),ylab="Proportion of genes with a peak")
		cat ("Proportion of genes with a peak from the ChIP dataset in a given genomic region:\n")
		cat (paste(c(names),sep='\t'))
		cat("\n")
		cat (paste(c(yChIPTotal),sep='\t'))
		cat("\n\n")
		cat ("Proportion of genes with a peak from the Control dataset in a given genomic region:\n")
		cat (paste(c(names),sep='\t'))
		cat("\n")
		cat (paste(c(yControlTotal),sep='\t'))
		cat("\n\n")
		if (bootstrap>1) {
			error.bar(barx,cum,wiskers)
			cat ("Standard deviation for the Control dataset in a given genomic region:\n")
			cat (paste(c(names),sep='\t'))
			cat("\n")
			cat (paste(c(sqrt(colVars)),sep='\t'))
			cat("\n\n")
		} 
		enrich <- NULL
		for (i in c(1:length(names))) {	
		  enrich[i] <- yChIPTotal[i]/yControlTotal[i];	  
		}
		barplot(enrich-1,xlab="",beside=TRUE, col=c(myColor), names.arg=c(names),ylab="Enrichment in comparison\nwith the control dataset",  yaxt="n")
		minX <- min(enrich-1)
		maxX <- max(enrich-1)
		x = seq(length=11, from=round(minX*10)/10, by=round((maxX-minX)*10)/100)
		axis(2, at=x,labels=x+1, las=2)

		cat ("Enrichment of genomic regions, ChIP peaks vs Control Peaks:\n")
		cat (paste(c(names),sep='\t'))
		cat("\n")
		cat (paste(c(yChIPTotal/yControlTotal),sep='\t'))
		cat("\n\n")
		if (bootstrap>1) {
			z <- (yChIPTotal-yControlTotal)/sqrt(colVars)
			pvalues <- 2*pnorm(-abs(z))
			cat ("Two-side P-values for each genomic category:\n")
			cat (paste(c(names),sep='\t'))
			cat("\n")
			cat (paste(c(yChIPTotal/yControlTotal),sep='\t'))
			cat("\n\n")
		}
	}
}
if (ifReg) { 	
	 n.types <- length(levels(chip.data$"6"))
	 myLevels <- levels(chip.data$"6")
	 nlev <- length(names)

	if (ifControl) {
		cum = matrix( 0, nrow=length(myLevels)+1,  ncol=nlev, byrow = TRUE) 
         	cumEnrichTotal = matrix( 0, nrow=length(myLevels),  ncol=nlev, byrow = TRUE) 
         	cumEnrichControl = matrix( 0, nrow=length(myLevels),  ncol=nlev, byrow = TRUE) 
	}else  {
		cum = matrix( 0, nrow=length(myLevels),  ncol=nlev, byrow = TRUE) 
         	cumEnrichTotal = matrix( 0, nrow=length(myLevels),  ncol=nlev, byrow = TRUE) 
	}
	colReg <-NULL
	for (r in c(1:length(myLevels))) {
		      tt <- which(chip.data$"6"==myLevels[r])
		      subset.data <- (chip.data[tt,])
		      countTotalSubset <- length(unique(subset.data$"1"))

			if (!ifk27) {

				tt <- which(subset.data$"16">0)
				countTotalEnhSubset <- length(unique(subset.data$"1"[tt]))
				tt <- which(subset.data$"10">0)
				countTotalPromSubset <- length(unique(subset.data$"1"[tt]))
				tt <- which(subset.data$"12">0)
				countTotalImDownSubset <- length(unique(subset.data$"1"[tt]))
				tt <- which(subset.data$"18">0)
				countTotalIntraSubset <- length(unique(subset.data$"1"[tt]))
				tt <- which(subset.data$"20">0)
				countTotal5kbDSubset <- length(unique(subset.data$"1"[tt]))
				tt <- which(subset.data$"28">0)
				countTotal1IntronISubset <- length(unique(subset.data$"1"[tt]))
				tt <- which(subset.data$"32">0)
				countTotalExonsISubset <- length(unique(subset.data$"1"[tt]))
				tt <- which(subset.data$"36">0)
				countTotalIntronsISubset <- length(unique(subset.data$"1"[tt]))
				tt <- which(subset.data$"40">0)
				countTotalJonctionsISubset <- length(unique(subset.data$"1"[tt]))
			      ySubsetTotal <- c (countTotalEnhSubset/countTotalSubset,countTotalPromSubset/countTotalSubset, countTotalImDownSubset/countTotalSubset,countTotalIntraSubset/countTotalSubset,countTotal5kbDSubset/countTotalSubset,countTotal1IntronISubset/countTotalSubset,countTotalExonsISubset/countTotalSubset,countTotalIntronsISubset/countTotalSubset,countTotalJonctionsISubset/countTotalSubset)
			} else {
				tt <- which(subset.data$"42">0)
				countTotalTSSSubset <- length(unique(subset.data$"1"[tt]))
				tt <- which(subset.data$"44">0)
				countTotalGeneSubset <- length(unique(subset.data$"1"[tt]))
				ySubsetTotal <- c (countTotalTSSSubset/countTotalSubset,countTotalGeneSubset/countTotalSubset)
			}	      
			for (i in c(1:nlev)) {	
			cum[r,i] <- ySubsetTotal[i]  
			cumEnrichTotal[r,i] <- ySubsetTotal[i]/yChIPTotal[i]   
			if (ifControl) {
				cumEnrichControl[r,i] <- ySubsetTotal[i]/yControlTotal[i]   
			}      
		      }
		      colReg[r]<-myColor[r+3]
	}
	if (ifControl) {
	 	for (i in c(1:nlev)) {	
			cum[4,i] <- yControlTotal[i]                
		}		
	}
	
	cat ("Proportion of genes with a peak from the ChIP dataset in a given genomic region:\n")
	for (r in c(1:length(myLevels))) {
		cat (paste(myLevels[r],":\n",sep=""))
		cat (paste(c(names),sep='\t'))
		cat("\n")
		cat (paste(c(cum[r,] ),sep='\t'))
		cat("\n")
	}
	cat("\n")
	par(mar=c(5.1, 7.1, 4.1, 2.1)) 	
	if (ifControl) {
		barx <- barplot(cum,xlab="",beside=TRUE, col=c(colReg, colors()[334]),  legend = c(myLevels,"Control"), names.arg=c(names),ylab="Proportion of genes with a peak")
		cat ("Proportion of genes with a peak from the Control dataset in a given genomic region:\n")
		cat (paste(c(names),sep='\t'))	
		cat("\n")
		cat (paste(c(yControlTotal),sep='\t'))
		cat("\n\n")				
		if (bootstrap>1) {
			wiskers <- cum-cum
			wiskers[nrow(wiskers),] <- sqrt(colVars)*1.96
			error.bar(barx,cum,wiskers)	
			cat ("Standard deviation for the Control dataset in a given genomic region:\n")
			cat (paste(c(names),sep='\t'))
			cat("\n")
			cat (paste(c(sqrt(colVars)),sep='\t'))
			cat("\n\n")	
		} 
	} else {
		barplot(cum,xlab="",beside=TRUE, col=c(colReg),  legend = c(myLevels), names.arg=c(names),ylab="Proportion of genes with a peak")
	}		
	barplot(cumEnrichTotal-1,xlab="",beside=TRUE, col=c(colReg), names.arg=c(names),ylab="Enrichment in comparison\nwith the total gene set",  yaxt="n")
	minX <- min(cumEnrichTotal-1)
	maxX <- max(cumEnrichTotal-1)
	x = seq(length=11, from=round(minX*10)/10, by=round((maxX-minX)*10)/100)
	axis(2, at=x,labels=x+1, las=2)

	cat ("Enrichment of genomic regions, Transcriptional categories vs All Genes:\n")
	for (r in c(1:length(myLevels))) {
			cat (paste(myLevels[r],":\n",sep=""))
			cat (paste(c(names),sep='\t'))	
			cat("\n")
			cat (paste(c(cumEnrichTotal[r,]),sep='\t'))			
			cat("\n")
	}	
	cat("\n")

	if (ifControl) {
		barplot(cumEnrichControl-1,xlab="",beside=TRUE, col=c(colReg), names.arg=c(names),ylab="Enrichment in comparison\nwith control",  yaxt="n")
		minX <- min(cumEnrichControl-1)
		maxX <- max(cumEnrichControl-1)
		x = seq(length=11, from=round(minX*10)/10, by=round((maxX-minX)*10)/100)
		axis(2, at=x,labels=x+1, las=2)
		cat ("Enrichment of genomic regions, ChIP peaks vs Control Peaks:\n")
		for (r in c(1:length(myLevels))) {
			cat (paste(myLevels[r],":\n",sep=""))
			cat (paste(c(names),sep='\t'))	
			cat("\n")
			cat (paste(c(cumEnrichControl[r,]),sep='\t'))			
			cat("\n")
		}	
		cat("\n")
		if (bootstrap>1) {
			cat ("Two-side P-values for each genomic category:\n")
			for (r in c(1:length(myLevels))) {
				z <- (cum[r,]-yControlTotal)/sqrt(colVars)
				pvalues <- 2*pnorm(-abs(z))
				cat (paste(myLevels[r],":\n",sep=""))
				cat (paste(c(names),sep='\t'))	
				cat("\n")
				cat (paste(c(pvalues),sep='\t'))			
				cat("\n")
			}	
		}
	}	
} 
sink() #stop sinking :)
dev.off()

