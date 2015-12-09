#usage $0 STEP RIGHT chipPeaks outputFile.png output.txt [controlPeaks]
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
if (length(args)>=5) {
	ifPDF=args[5]
}
if (length(args)==4 & args[4]==1) {
	ifPDF=1
}

ifControl <- 0
if (length(args)>=4 & args[4]!=1 & args[4]!=0) {
  dataTable <-read.table(file=args[4], header=TRUE);
  control.data<-data.frame(dataTable)
  ifControl <- 1
}
if (ifReg & ifControl) {

  if (ifPDF==1) {
	pdf(file = pngFile, width = 14, height = 13, pointsize = 20, bg = "white")
  } else {
       png(filename = pngFile, type="cairo" , width =  1440, height =  1040, units = "px", pointsize = 20, bg = "white", res = NA)
	plot(1:10)
  }
  op <- par(mfrow = c(3,2))
} else {

  if (ifPDF==1) {
	pdf(file = pngFile, width = 22, height = 8, pointsize = 20, bg = "white")
  } else {
        png(filename = pngFile, type="cairo" , width = 1580, height = 530, units = "px", pointsize = 20, bg = "white", res = NA)
	plot(1:10)
  }
  op <- par(mfrow = c(1,2))
}
myColor <- 1
myColor[1] <- colors()[131]
myColor[2] <- colors()[59]
myColor[3] <- colors()[76]
myColor[4] <- colors()[88]
myColor[5] <- colors()[17]
myColor[6] <- colors()[565]
myColorControl <- 1
myColorControl[1] <- colors()[24]
myColorControl[2] <- colors()[278]
myColorControl[3] <- colors()[305]
myColorControl[4] <- colors()[219]
myColorControl[5] <- colors()[343]
myColorControl[6] <- colors()[245]
myLevels <- 0
if (ifReg) {
 	if (ifControl) {
	 #control vs real:

		 countTotal <- length(chip.data$Reg)
		  totalChIP <- summary(chip.data$Type)/countTotal
		  tt <- which(chip.data$Type=="intragenic")	
		  subset.chip <- chip.data[tt,]
		  countIntra <- length(subset.chip$Reg)
		  intraChip<- summary(subset.chip$TypeIntra)/countTotal
		  nlev <- length(levels(chip.data$Type))
		  countTotalCont <- length(control.data$Reg)
		      totalContr <- summary(control.data$Type)/countTotalCont
		      tt <- which(control.data$Type=="intragenic")	
		      subset.control <- control.data[tt,]
		      countIntraCont <- length(subset.control$Reg)
		      intraControl<- summary(subset.control$TypeIntra)/countTotalCont
		      cum = matrix( 0, nrow=2,  ncol=nlev, byrow = TRUE) 
		      for (i in c(1:nlev)) {	
			cum[1,i] <- totalChIP[i]
			cum[2,i] <- totalContr[i]
		      }

			labels<-c("GeneDown.", "Enh.", "Imm.Down.", "Interg.", "Intrag.", "Prom.")
			if (length(labels)==length(levels(chip.data$Type))) {
				barplot(cum,xlab="",beside=TRUE, col=c(myColor[1],colors()[328]), names.arg=labels,ylab="Proportion of peaks")
			} else {
		      		barplot(cum,xlab="",beside=TRUE, col=c(myColor[1],colors()[328]), names.arg=c(levels(chip.data$Type)),ylab="Proportion of peaks")
			}

    			position <- 'topleft'
     			inset <- c(0.1, 0)
			legend(position, c("ChIP","Control"), bty="n",fill=c(myColor[1],colors()[328]), inset=inset)

		      nlev <- length(levels(subset.chip$TypeIntra))
		      cum = matrix( 0, nrow=2,  ncol=nlev, byrow = TRUE) 
		      for (i in c(1:nlev)) {	
			cum[1,i] <- intraChip[i]
			cum[2,i] <- intraControl[i]
		      }
		      barplot(cum,xlab="",beside=TRUE, col=c(myColor[1],colors()[328]), names.arg=c(levels(subset.chip$TypeIntra)),ylab="Proportion of peaks")

			position <- 'topleft'
     			inset <- c(0.1, 0)
			legend(position, c("ChIP","Control"), bty="n",fill=c(myColor[1],colors()[328]), inset=inset)
	}
	 n.types <- length(levels(chip.data$Reg))
	 myLevels <- levels(chip.data$Reg)
	 nlev <- length(levels(chip.data$Type))
         cum = matrix( 0, nrow=length(myLevels),  ncol=nlev, byrow = TRUE) 
         countTotal <- length(chip.data$Reg)
	 colReg <-NULL
         for (r in c(1:length(myLevels))) {
	      tt <- which(chip.data$Reg==myLevels[r])
              totalChIP <- summary(chip.data$Type[tt])/countTotal
	      for (i in c(1:nlev)) {	
		cum[r,i] <- totalChIP[i]                
	      }
	      colReg[r]<-myColor[r+3]
         }

	labels<-c("GeneDown.", "Enh.", "Imm.Down.", "Interg.", "Intrag.", "Prom.")
	if (length(labels)==length(levels(chip.data$Type))) {
		#barplot(cum,xlab="",beside=TRUE, col=c(myColor[1],myColor[5]), names.arg=labels,ylab="Proportion of peaks")
		barplot(cum,xlab="",beside=TRUE, col=c(colReg), names.arg=labels,ylab="Proportion of peaks")
	} else {
		barplot(cum,xlab="",beside=TRUE, col=c(colReg), names.arg=c(levels(chip.data$Type)),ylab="Proportion of peaks")
	}

    	position <- 'topleft'
     	inset <- c(0.1, 0)
	legend(position, c(myLevels), bty="n",fill=c(colReg), inset=inset)


	 nlev <- length(levels(chip.data$TypeIntra))
         cum = matrix( 0, nrow=length(myLevels),  ncol=nlev, byrow = TRUE) 	 
         for (r in c(1:length(myLevels))) {
	      tt <- which(chip.data$Reg==myLevels[r]&chip.data$Type=="intragenic")
              totalChIP <- summary(chip.data$TypeIntra[tt])/countTotal
	      for (i in c(1:nlev)) {	
		cum[r,i] <- totalChIP[i]                
	      }	      
         }
         barplot(cum,xlab="",beside=TRUE, col=c(colReg), names.arg=c(levels(chip.data$TypeIntra)),ylab="Proportion of peaks")
	position <- 'topleft'
     	inset <- c(0.1, 0)
	legend(position, c(myLevels), bty="n",fill=c(colReg), inset=inset)

	 
	 if (ifControl) {
		nlev <- length(levels(control.data$Type))
		cum = matrix( 0, nrow=length(myLevels),  ncol=nlev, byrow = TRUE) 
		 countTotal <- length(control.data$Reg)
		 colReg <-NULL
		 for (r in c(1:length(myLevels))) {
		      tt <- which(control.data$Reg==myLevels[r])
		      totalcontrol <- summary(control.data$Type[tt])/countTotal
		      for (i in c(1:nlev)) {	
			cum[r,i] <- totalcontrol[i]                
		      }
		      colReg[r]<-myColorControl[r+3]
		 }
		labels<-c("GeneDown.", "Enh.", "Imm.Down.", "Interg.", "Intrag.", "Prom.")
		if (length(labels)==length(levels(chip.data$Type))) {
			barplot(cum,xlab="",beside=TRUE, col=c(colReg), names.arg=labels,ylab="Proportion of peaks")
		} else {
		 	barplot(cum,xlab="",beside=TRUE, col=c(colReg), names.arg=c(levels(control.data$Type)),ylab="Proportion of peaks")
		}
		position <- 'topleft'
	     	inset <- c(0.1, 0)
		legend(position, c(myLevels), bty="n",fill=c(colReg), inset=inset)

		 nlev <- length(levels(control.data$TypeIntra))
		 cum = matrix( 0, nrow=length(myLevels),  ncol=nlev, byrow = TRUE) 	 
		 for (r in c(1:length(myLevels))) {
		      tt <- which(control.data$Reg==myLevels[r]&control.data$Type=="intragenic")
		      totalcontrol <- summary(control.data$TypeIntra[tt])/countTotal
		      for (i in c(1:nlev)) {	
			cum[r,i] <- totalcontrol[i]                
		      }		      
		 }
		 barplot(cum,xlab="",beside=TRUE, col=c(colReg),  names.arg=c(levels(control.data$TypeIntra)),ylab="Proportion of peaks")	   
		position <- 'topleft'
	     	inset <- c(0.1, 0)
		legend(position, c(myLevels), bty="n",fill=c(colReg), inset=inset)
	 }
} else {
  countTotal <- length(chip.data$Reg)
  totalChIP <- summary(chip.data$Type)/countTotal
  tt <- which(chip.data$Type=="intragenic")	
  subset.chip <- chip.data[tt,]
  countIntra <- length(subset.chip$Reg)
  intraChip<- summary(subset.chip$TypeIntra)/countTotal
  nlev <- length(levels(chip.data$Type))
 if (ifControl) {
      countTotalCont <- length(control.data$Reg)
      totalContr <- summary(control.data$Type)/countTotalCont
      tt <- which(control.data$Type=="intragenic")	
      subset.control <- control.data[tt,]
      countIntraCont <- length(subset.control$Reg)
      intraControl<- summary(subset.control$TypeIntra)/countTotalCont
      cum = matrix( 0, nrow=2,  ncol=nlev, byrow = TRUE) 
      for (i in c(1:nlev)) {	
	cum[1,i] <- totalChIP[i]
	cum[2,i] <- totalContr[i]
      }

      labels<-c("GeneDown.", "Enh.", "Imm.Down.", "Interg.", "Intrag.", "Prom.")
      if (length(labels)==length(levels(chip.data$Type))) {
		#barplot(cum,xlab="",beside=TRUE, col=c(myColor[1],myColor[5]), names.arg=labels,ylab="Proportion of peaks")
	 barplot(cum,xlab="",beside=TRUE, col=c(myColor[1],colors()[328]), names.arg=labels,ylab="Proportion of peaks")
       } else {
     	 barplot(cum,xlab="",beside=TRUE, col=c(myColor[1],colors()[328]), names.arg=c(levels(chip.data$Type)),ylab="Proportion of peaks")
       }

    	position <- 'topleft'
     	inset <- c(0.1, 0)
	legend(position,c("ChIP","Control"), bty="n",fill=c(myColor[1],colors()[328]), inset=inset)

      nlev <- length(levels(subset.chip$TypeIntra))
      cum = matrix( 0, nrow=2,  ncol=nlev, byrow = TRUE) 
      for (i in c(1:nlev)) {	
	cum[1,i] <- intraChip[i]
	cum[2,i] <- intraControl[i]
      }
      barplot(cum,xlab="",beside=TRUE, col=c(myColor[1],colors()[328]), names.arg=c(levels(subset.chip$TypeIntra)),ylab="Proportion of peaks")	
	position <- 'topleft'
     	inset <- c(0.1, 0)
	legend(position,c("ChIP","Control"), bty="n",fill=c(myColor[1],colors()[328]), inset=inset)

   } else {
	labels<-c("GeneDown.", "Enh.", "Imm.Down.", "Interg.", "Intrag.", "Prom.")
	if (length(labels)==length(levels(chip.data$Type))) {
		barplot(totalChIP,xlab="", col=myColor, names.arg=labels,ylab="Proportion of peaks")
	} else {
      		barplot(totalChIP,xlab="", col=myColor,ylab="Proportion of peaks")
	}
        barplot(intraChip,xlab="", col=myColor,ylab="Proportion of peaks")
   }
}
dev.off()
