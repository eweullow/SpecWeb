

#Author: Andrew Sila
#Email: asila@cgiar.org
#Purpose: Function for reading Bruker OPUS files from MPA, HTS-xt and Alpha
#Date: 30/March/2014
#Version: Alpha (Testing)

#To use the function you are required to change the file paths at the bottom of this file: 
#1. opus.folder=Folder with raw spectra to be converted
#2. converted.folder=Folder where converted tables will be kept
#3. run the last line to get conversions

# Begin by checking the required packages are installed. --------------------------------
is.installed <- function(anypkg){
	
  is.element(anypkg, installed.packages()[,1])
  
}

required.packages <- c("hexView","caret","httr","shiny","shinyFiles","jsonlite","DT","shinythemes","readr", "shinyjs", "prospectr", "soil.spec", "xlsx", "reshape2", "ggplot2", "plotly", "BalancedSampling")

installp <- which(!is.installed(required.packages) == TRUE)

#Install missing packages

if (length (installp) > 0){
	
  install.packages(required.packages[installp])

}

require(shinythemes)

require(hexView) #Needs readRaw function from this library

read.opus <- function(file.name, sp=NULL, codes=c("ZFF","RES","SNM","DAT","LWN","FXV","LXV","NPT","MXY","MNY","END","TIM"), 

  plot.spectra=FALSE, print.progress=FALSE, speclib="ICRAF")

{
	
	spec.opts <- new.env(hash=TRUE)
	
	spec.env <- function(MIR = c(390, 7500), NIRS = c(3900, 12500), NIRP = c(4000,10000), VISNIR1 = c(420, 960), VISNIR2 = c(1020, 1770), VISNIR3 = c(1830, 2480), icraf.htsxt = c(3578, 7497.964, 599.76), icraf.alpha = c(2542, 3998.12872, 399.387991),icraf.alpha.znse = c(1714, 3996.4810, 499.8151), icraf.mpa = c(2307, 12493.2, 3598.69), CO2.band = c(2350.8,2379.8), signif.digit=5, attributes = c("ORCCNS", "PHIHOX", "ALUM3S", "ECAM3S", "EXKM3S", "EMGM3S", "ENAM3S", "EXB", "NITCNS", "SNDLDF"), mdnames = c("Instrument_name", "Instrument_URL", "Laboratory_name", "Laboratory_contact", "Laboratory_URL", "Material_class", "Wavenumber_conversion", "Wavenlength_unit", "Location_error"), show.env = FALSE){

    pl.lst <- list(

     MIR = MIR,

     NIRS = NIRS,

     NIRP = NIRP,

     VISNIR1 = VISNIR1,

     VISNIR2 = VISNIR2,

     VISNIR3 = VISNIR3,

     icraf.htsxt = icraf.htsxt,

     icraf.alpha = icraf.alpha,

     icraf.mpa = icraf.mpa,

     icraf.alpha.znse = icraf.alpha.znse,

     CO2.band = CO2.band,

     signif.digit = signif.digit,

     attributes = attributes,

     mdnames = mdnames

     )

    x <- lapply(names(pl.lst), function(x){ assign(x, pl.lst[[x]], envir=spec.opts) })

    if(show.env){

     return(pl.lst)

   }

 }

 spec.env()

 if(!(speclib=="ICRAF"|speclib=="New")){ stop("'speclib' must be one of the following: 'ICRAF' or 'New'") }

 if(file.exists(file.name)){

    	## Read metainfo

      try( pa <- hexView::readRaw(file.name, offset = 0, nbytes = file.info(file.name)$size, human = "char", size = 1, endian = "little"), silent=TRUE )

      if(!class(.Last.value)[1]=="try-error"){

        pr <- pa$fileRaw

	   		## Get source of instrument

       ins <- grepRaw("INS", pr, all=TRUE)

       ins <- readRaw(file.name, offset = ins[length(ins)]+7, nbytes = 3, human = "char", size = 1, endian = "little")

       ins <- blockString(ins)

		 	## Get source of infrared to know if NIR or MIR

     src <- grepRaw("SRC", pr, all=TRUE)

     src <- readRaw(file.name, offset = src[length(src)]+4, nbytes = 3, human = "char", size = 1, endian = "little")

     src <- blockString(src)

     instr.range <- tolower(paste(ins, src, sep="-"))

  		  	## Get Beam Splitter

          bms <- grepRaw("BMS", pr, all=TRUE)

          bms <- readRaw(file.name, offset = bms[length(bms)]+4, nbytes = 4, human = "char", size = 1, endian = "little")

          bms <- blockString(bms)

         		## Wavenumbers for MIR spectra from Tensor are assigned prefix "m", MIR spectra from Alpha prefixed "a"; for NIR MPA "n"

        		if(instr.range=="ten-off"){ instr.range="ten-mir"} ## AS: Old ten-mir written as tensor-27

            pref <- ifelse(instr.range=="ten-mir", "m", ifelse(instr.range=="alp-mir", "a", ifelse(instr.range=="mpa-nir", "n", "X")))

            if(speclib=="ICRAF"){ 

             if(instr.range=="ten-mir"){


               wb <- rev(seq(get("icraf.htsxt", spec.opts)[3], get("icraf.htsxt", spec.opts)[2], (get("icraf.htsxt", spec.opts)[2]-get("icraf.htsxt", spec.opts)[3])/(get("icraf.htsxt", spec.opts)[1]-1)))

             }

             if(instr.range=="alp-mir"){

              if(bms!="ZnSe"){

                wb <- rev(seq(get("icraf.alpha", spec.opts)[3], get("icraf.alpha", spec.opts)[2], (get("icraf.alpha", spec.opts)[2]-get("icraf.alpha", spec.opts)[3])/(get("icraf.alpha", spec.opts)[1]-1)))

              }

              if(bms=="ZnSe"){

                wb <- rev(seq(get("icraf.alpha.znse", spec.opts)[3], get("icraf.alpha.znse", spec.opts)[2], (get("icraf.alpha.znse", spec.opts)[2]-get("icraf.alpha.znse", spec.opts)[3])/													(get("icraf.alpha.znse", spec.opts)[1]-1)))}

              }


              if(instr.range=="mpa-nir"){

                wb <- rev(seq(get("icraf.mpa", spec.opts)[3], get("icraf.mpa", spec.opts)[2], (get("icraf.mpa", spec.opts)[2]-get("icraf.mpa", spec.opts)[3])/(get("icraf.mpa", spec.opts)[1]-1)))

              }

            }

            if(!(instr.range=="ten-mir"|instr.range=="alp-mir"|instr.range=="mpa-nir"|instr.range=="ten-off"|bms=="ZnSe")){ stop("Unknown file format. See '?read.opus' for more info.") }  

            if(speclib=="New"){

             if(instr.range=="ten-mir"){ pref="m" }

             if(instr.range=="alp-mir"){ pref="a"}

             if(instr.range=="mpa-nir"){ pref="n"}

             if(bms=="ZnSe"){ pref="a"}

           }
			## Get positions where the following parameters are found in the file  

     z <- grepRaw(codes[1],pr,all=TRUE)[1]+5

     re <- grepRaw(codes[2],pr,all=TRUE)[1]+5

     snm <- grepRaw(codes[3],pr,all=TRUE)[1]+7

     dat <- grepRaw(codes[4],pr,all=TRUE)[1]+7

     lwn <- grepRaw(codes[5],pr,all=TRUE)[1]+7

     fx <- grepRaw(codes[6],pr,all=TRUE)[3]+7

     lx <- grepRaw(codes[7],pr,all=TRUE)[3]+7

     npt0 <- grepRaw(codes[8],pr,all=TRUE)[2]+3

     npt1 <- grepRaw(codes[8],pr,all=TRUE)[3]+7

     mxy <- grepRaw(codes[9],pr,all=TRUE)[1]+7 

     mny <- grepRaw(codes[10],pr,all=TRUE)[3]+7 

     end <- grepRaw(codes[11],pr,all=TRUE)+11

     tim <- grepRaw(codes[12],pr,all=TRUE)+11

			## calculate end and start of each block:

      offs <- sapply(5:10, function(x){end[x]})

      byts <- diff(offs)

      ZFF <- readRaw(file.name, offset=z, nbytes=4, human="int", size=2)[[5]][1]

      RES <- readRaw(file.name, offset=re, nbytes=4, human="int", size=2)[[5]][1]

      snm.lab.material <- blockString(readRaw(file.name, offset = snm, nbytes = 22, human = "char", size = 1, endian = "little"))

      if(!nzchar(snm.lab.material)){

        SSN <- ""

        Material <- ""

        warning("Product name not found inside OPUS file...")

      }

      else {

        if(!length(grep(snm.lab.material, pattern=";"))==0){

          snm.lab.material <- as.vector(strsplit(snm.lab.material,";"))[[1]]

          SSN <- paste0(snm.lab.material[2], snm.lab.material[1])

          Material <- snm.lab.material[3]
        } 
        else {

         if(!length(grep(snm.lab.material, pattern="_"))==0){

          SSN <- sub("_", "", snm.lab.material)

          Material <- ""          
        }

        else {

         if(!length(snm.lab.material)==0)

         { 
         	
           SSN <- snm.lab.material

           Material <- ""          

         }

       }

     }   

   }

		## Set three SSN first three characters to lower
		
		SSN <- paste0(tolower(substr(SSN,1,3)), substr(SSN,4,20))

   Scandate <- blockString(readRaw(file.name, offset = dat, nbytes = 10, human = "char", size = 1, endian = "little"))

   Scantime <- blockString(readRaw(file.name, offset = tim[2]-4, nbytes = 8, human = "char", size = 1, endian = "little"))

   Scandate <- paste(Scandate,Scantime)

   LWN <- readRaw(file.name, offset=lwn, nbytes=8, human="real", size=8)[[5]][1]

  	 	## Combine the above parameters

      spectrum.meta <- c(SSN, Material, Scandate, ZFF, RES, LWN)

		## Get number of data points for each spectra data block   		

   NPT0 <- readRaw(file.name, offset=npt0, nbytes=12, human="int", size=4)[[5]][2]

   NPT1 <- readRaw(file.name, offset=npt1, nbytes=4, human="int", size=4)[[5]][1]

  	 	fxv <- readRaw(file.name, offset=fx, nbytes=16, human="real", size=8)[[5]][1] ## fxv:	Frequency of first point
  	 	
  	 	lxv <- readRaw(file.name, offset=lx, nbytes=16, human="real", size=8)[[5]][1] ## lxv:	Frequency of last point
  	 	
  	 	Wavenumbers <- rev(seq(lxv, fxv, (fxv-lxv)/(NPT1-1)))
  	 	
  	 	## Read all through all the data blocks inside the OPUS file:
  	 	
  	 	nbytes1 <- NPT0*4 ## initial parameters
  	 	
  	 	smxa <- c()
  	 	
  	 	smna <- c()
  	 	
  	 	nbytes.f <- NPT1*4
  	 	
  	 	#Start reading available data blocks

      mxsa<-NULL

      for (f in 1:3){

        offs.f<-offs[f]

        try(opus.p <- readRaw(file.name,width=NULL,offset=offs.f-4,nbytes=nbytes.f,human="real",size=4,endian="little"), silent = TRUE)


        if(!class(.Last.value)[1]=="try-error"){

         spectra <- opus.p[[5]]

         mxs<-max(spectra)

         mxsa<-c(mxsa,mxs)

       }

     }

		#Determine valid data block to read
		
   fs<-which(mxsa==max(mxsa))

   if (!length(fs)<1){

    offs.f<-offs[fs]

    try(opus.p <- readRaw(file.name,width=NULL,offset=offs.f-4,nbytes=nbytes.f,human="real",size=4,endian="little"), silent=TRUE)

    if(!class(.Last.value)[1]=="try-error"){

     spectra <- opus.p[[5]]

   }

 }
			## Make compatible to ICRAF spectra:

			if(speclib=="ICRAF"){

       spectra <- spline(Wavenumbers, spectra, xout=wb, method="natural")$y

       Wavenumbers <- wb

     }

			## Specify if graphics showing spectra being converted is displayed
			
			if(plot.spectra==TRUE){

       plot(Wavenumbers, spectra, ylab="Absorabance", xlab=expression("Wavenumbers cm"^-1), type="l")   

       mtext(paste("File source: ", getwd(),file.name,sep="/"), side=3,line=2,cex=1)
     }
			## Print progress of conversion
			
			if(print.progress==TRUE){

       message(paste("Converting ", file.name, " file", sep=""))

     }


     samples <- data.frame(SAMPLEID=spectrum.meta[1], Material=spectrum.meta[2], Zero.Filing=spectrum.meta[4], Resolution=spectrum.meta[5], LWN=spectrum.meta[6], DateTime=as.POSIXct(spectrum.meta[3], format="%d/%m/%Y %H:%M:%S "))

     ab <- data.frame(as.list(spectra, signif.digit))

     spec<-c(spectrum.meta,round(spectra,5))

     names(spec)<-c("SAMPLEID","Material","Datetime","Zero.filling","Resolution","LWN",paste0(pref,round(Wavenumbers,1)))

   }

   out<-spec

   } else {

     warning(paste("File",file.name,"does not exist"))

   }
 }


 library(shiny)
 library(shinyFiles)
 library(xlsx)
 library(readr)

######################################## Load data ##########################################################
#Read Mua data
mua.avg<-readRDS("./data/Mua_std_averaged_spectra.rds")

########################################################################

#Read White sand data
wtsd.avg<-readRDS("./data/WhiteSand_std_averaged_spectra.rds")

###...Kenndard stone
n = 1
study <- 'afsis'


#### Calibration script for PLS and Rf methods
# Script for fitting calibration models for infrared spectroscopy data. ----------------------------------------------------------------
# Author: Andrew Sila , April 2016. --------------------------------

# Begin by checking the required packages are installed. --------------------------------
is.installed <- function(anypkg){
	
  is.element(anypkg, installed.packages()[,1])
  
}

required.packages <- c("hexView","caret","readr","mlbench","pROC","rpart","caretEnsemble","soil.spec","FitAR",

  "doParallel","ggplot2","dplyr","downloader","prospectr","randomForest","gridExtra")

installp <- which(!is.installed(required.packages) == TRUE)

#Install missing packages

if (length (installp) > 0){
	
  install.packages(required.packages[installp])
}

# Load packages.

suppressMessages(library(caret))
suppressMessages(library(readr))
suppressMessages(library(mlbench))
suppressMessages(library(pROC))
suppressMessages(library(rpart))
suppressMessages(library(caretEnsemble))
suppressMessages(library(soil.spec))
suppressMessages(library(FitAR))
suppressMessages(library(doParallel))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(downloader))
suppressMessages(library(prospectr))
suppressMessages(library(randomForest))
suppressMessages(library(gridExtra))

# Run in parallel to speed up the validation of the model fitting process.

registerDoParallel()

getDoParWorkers()

calibrate <- function(wd,infrared.data,reference.data,hout,method = c("RF","PLS")){
	
	if(method == "RF"){
		
  # Start by setting working directory. 
  
  setwd(wd)
  
  mir <- infrared.data
  
  ref <- reference.data
  
  hout <- hout
  
  # Exclude metadata variables.
  
  mir1 <- as.matrix(mir[,-1])
  
  wave <- as.numeric(substr(colnames(mir1),2,19))
  
  #prefx <- substr(colnames(mir1),1,1)[900]

  colnames(mir1) <- wave
  
  
  # First derivative.
  
  de1 <- trans(mir1,tr = "derivative",order = 1,gap = 23)
  
  der1 <- rev(as.data.frame(de1$trans))
  
  #colnames(der1) <- paste0(prefx,wave)
  
  # Save derivative spectra.
  der1.ssn <- as.data.frame(cbind(as.vector(mir[,1]),der1))
  
  colnames(der1.ssn) <- c("SSN",colnames(der1))
  
  write.table(der1.ssn,file = "First derivative.csv",sep = ",",row.names = FALSE)
  
  der1.ssn<-as.data.frame(read_csv("First derivative.csv"))
  
  # Merge with first derivative preprocessed spectra.
  
  ref.mir <- merge(ref,der1.ssn,by.x = "SSN",by.y = "SSN")
  
  rc <- colnames(ref)
  
  # Which columns contains reference data?
  
  ref <- ref.mir[,rc]
  
  # Extract spectral predictors
  
  mirp <- colnames(der1.ssn)[-1]
  
  spectra <- ref.mir[,mirp]
  
  #Create two new subfolders within the current working using:
  
  b <- getwd()
  
  if(!file.exists("Models")){dir.create("Models")}
  
  if(!file.exists("calibration_plots")){dir.create("calibration_plots")}
  
  # Fit calibration models for the training set and use the testing set to validate
  
  # the models.
  
  set.seed(67523)
  
  testing <- which(ref.mir$SSN%in%hout$SSN)
  
  #with hout
  
  #Loop for calibration of all soil properties in the reference set starts here
  
  msummary <- NULL
  
  hd <- colnames(ref)[-1]
  
  for (q in 1:length(hd)){
  	
  	#Get variable names from reference table
  	
  	refq <- which(colnames(ref)%in%hd[q])
  	
  	ref.q <- ref[,refq]
  	
  	pms.a <- NULL
  	
  	pred.all <- NULL
  	
  	cal <- na.omit(cbind(as.vector(ref.q),spectra)[-testing,])
  	
  	val <- na.omit(cbind(as.vector(ref.q),spectra)[testing,])
  	
  	colnames(cal) <- c(colnames(ref)[refq],colnames(spectra))
  	
  	colnames(val) <- colnames(cal)
  	
  	trainX <- cal[,-1]
  	
  	trainY <- cal[,1]
  	
  	rf.m <- randomForest(data = cal, x = cal[,-1], y = cal[,1],

      ntree = 200, importance = TRUE, na.action = na.omit)

  	predicted<-predict(rf.m )
  	
  	measured<-trainY
  	
  	# compute RMSE and R-squared values for the calibration set

  	training.parameters <- round(postResample(predicted,measured),2)
  	
  	RSQ <- training.parameters[2]
  	
  	RMSE <- training.parameters[1]
  	
  	predicted.test <- predict(rf.m,val[,-1])
  	
  	measured.test <- val[,1]
  	
  	# compute RMSE and R-squared values for the validation set
  	
  	testing.parameters <- round(postResample(predicted.test,measured.test),2)
  	
  	RSP <- testing.parameters[2]
  	
  	RMSEP <- testing.parameters[1]
  	
  	model.summary <- c(hd[q],training.parameters,testing.parameters)
  	
  	msummary <- rbind(msummary,model.summary)
  	
  	# Save obtained model
  	
  	saveRDS(rf.m,file = paste0(b,"/","models/",hd[q],".rds"))
  	
  	pm <- as.data.frame(cbind(measured,predicted))
  	
  	# Create scatter plot for the predicted versus the measured
  	
   p <- ggplot(pm, aes(x = measured,y = predicted)) + 

   geom_point(col = "black",size = 3,alpha = 0.5) + 

   ggtitle(paste0("Calibration for ",hd[q])) + 

   xlab("Measured") + 

   ylab("Predicted")

   p <- p + stat_smooth(method = lm, se = FALSE, color = 'black',alpha = 0.1)

   p <- p + theme(plot.title = element_text(lineheight = 3, face = "bold", color = "black", size = 20))

	p <- p + theme(text = element_text(size = 20)) # this will change all text size

	p <- p + annotate('text', label = paste('R^2 == ',RSQ), parse = TRUE,Inf, 

   -Inf,hjust = 2.5, vjust = -7.8) + annotate('text', label = paste('RMSE == ',RMSE),

   parse = TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)

	# Centre title
	p <- p + theme(plot.title = element_text(hjust  = 0.5))
	
	p <- p + xlim(range(pm)) + ylim(range(pm))
	
	# Validation data
	
	pmp <- as.data.frame(cbind(measured.test,predicted.test))
	
	p2 <- ggplot(pmp, aes(x = measured.test,y = predicted.test))  + 
	
	geom_point(col = "brown",size = 3,alpha = 0.5)  + 
	
	ggtitle(paste0("Validation for ",hd[q]))  + 
	
	xlab("Measured")  + 
	
	ylab("Predicted")
	
  p2 <- p2 + stat_smooth(method = lm, se = FALSE, color = 'brown',alpha = 0.1)

  p2 <- p2 + theme(plot.title = element_text(lineheight = 3, face = "bold",

   color = "black", size = 20))

    # this will change all text size 

    p2 <- p2 + theme(text = element_text(size = 20)) 
    
    p2 <- p2 + annotate('text', label = paste('R^2 == ',RSP),

      parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8) + 
    
    annotate('text', label = paste('RMSE == ',RMSEP), 

      parse = TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)

	# Centre title
 p2 <- p2 + theme(plot.title = element_text(hjust  = 0.5))

 p2 <- p2 + xlim(range(pmp)) + ylim(range(pmp))

    # Save produced scatter plots inro png files
    
    png(file = paste0(b,"/Calibration_plots/",hd[q],".png"),height = 400,width = 800)

    grid.arrange(p,p2,nrow = 1)

    dev.off()

  }
  
  	# Save the model summary of the calibration and validation sets.

  	colnames(msummary)<-c("Soil properties","OOB RMSEC","OOB Rsquared", "Holdout

     RMSEP","Holdout Rsquared")
  	
  	write.table(msummary,file = "Model_Summary.csv",sep = ",",row.names = FALSE)
  	
  	# Combine validation and calibration sets to get full model.

    setwd(b)

    if(!file.exists("Full_Models")){
      dir.create("Full_Models")}

      if(!file.exists("Full_calibration_plots")){

       dir.create("Full_calibration_plots")}

       msummary <- NULL

  	hd <- colnames(ref[,-1])#Exclude SSN 
  	
  	all.predicted <- NULL
  	
  	for (q in 1:length(hd)){

     refq <- which(colnames(ref)%in%hd[q])

     ref.q <- ref[,refq]

     pms.a <- NULL

     pred.all <- NULL

     cal <- cbind(as.vector(ref.q),spectra)

     colnames(cal) <- c(colnames(ref)[refq],colnames(spectra))

     p <- which(is.na(der1.ssn[,1]) ==  TRUE)

     ifelse(length(p)>0,ssn<-der1.ssn[-p,1],ssn<-der1.ssn[,1])

     ifelse(length(p)>0,der1.ssn<-der1.ssn[-p,],der1.ssn<-der1.ssn)

  	#Select training and testing sets
  	
  	cal <- na.omit(cal)
  	
  	trainX <- cal[, -1]
  	
  	trainY <- cal[,1]
  	
  	# Begin model
  	
  	rf.m <- randomForest(data = cal, x = cal[,-1], y = cal[,1], ntree = 200, 

     importance = TRUE, na.action = na.omit)
  	
  	predi <- predict(rf.m)
  	
  	y <- trainY
  	
  	training.parameters <- c(hd[q],round(postResample(predi,y),3))
  	
  	msummary <- rbind(msummary,training.parameters)
  	
  	# Save the model.
  	
  	saveRDS(rf.m,file = paste0(b,"/","Full_Models/",hd[q],".rds"))
  	
  	# Make prediction using the fitted model and all the spectra in the 
  	
  	# training data
  	
  	predicted <- predict(rf.m )
  	
  	measured <- trainY
  	
  	# Computes RMSE and R-squared values for the calibration set
  	
  	training.parameters <- round(postResample(predicted,measured),3)
  	
  	RSQ <- training.parameters[2]
  	
  	RMSE <- training.parameters[1]
  	
  	predicted.test <- predict(rf.m,val[,-1])
  	
  	measured.test <- val[,1]
  	
  	testing.parameters <- round(postResample(predicted.test,measured.test),2)
  	
  	#computes RMSE and R-squared values for the validation set
  	
  	RSP <- testing.parameters[2]
  	
  	RMSEP <- testing.parameters[1]
  	
  	model.summary <- c(hd[q],training.parameters,testing.parameters)
  	
  	msummary <- rbind(msummary,model.summary)
  	
  	pm <- as.data.frame(cbind(y,predicted))
  	
  	colnames(pm) <- c("measured","predicted")
  	
  	# Create scatter plot for the predicted versus the measured
  	
  	p1 <- ggplot(pm, aes(x = measured,y = predicted)) + 

   geom_point(col = "brown",size = 3,alpha = 0.2) + 
   
   ggtitle(paste0("Calibration for ",hd[q])) + 
   
   xlab("Measured") + 
   
   ylab("Predicted")
   
   p1 <- p1 + stat_smooth(method = lm, se = FALSE, color = 'brown',alpha = 0.1)
   
   p1 <- p1 + theme(plot.title = element_text(lineheight = 3,

     face = "bold", color = "black", size = 20))
   
   # this will change all text size 

   p1 <- p1 + theme(text = element_text(size = 20))
   
   p1 <- p1 + annotate('text', label = paste('R^2 == ',RSQ),

     parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8) + 
   
   annotate('text', label = paste('RMSE == ',RMSE), 

     parse = TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)

  # Centre title
  p1 <- p1 + theme(plot.title = element_text(hjust  = 0.5))
  
  p1 <- p1 + xlim(range(pm)) + ylim(range(pm))

  # Save the created models 
  
  par(mfrow = c(1,1))
  
  ggsave(file = paste0(b,"/","Full_calibration_plots/",hd[q],".png"),

    height = 6, width = 6, units = "in",p1)
  
  # Get prediction from the full model
  
  predicted.pq <- predict(rf.m,der1.ssn[,-1])

  # Combine all predicted values from the full models
  
  all.predicted <- cbind(all.predicted,predicted.pq)
  
}

  #Add sample ids (SSN) to the combined predicted values.

  all.predicted.SSN <- cbind(as.vector(ssn),all.predicted)

  colnames(all.predicted.SSN) <- c("SSN",hd)

  colnames(msummary) <- c("Soil properties","RMSEC","Rsquared")

  #Save full models' summary and predicted values

  write.table(msummary, file = "Full models summary.csv",

   sep = ",",row.names = FALSE)

  write.table(all.predicted.SSN, file = "All predictions.csv",

   sep = ",", row.names = FALSE)

}

  # PLS Regression method
  # ----------------------------------------------------------------
  if(method ==  "PLS"){
  	
  	setwd(wd)
  	
  	mir <- infrared.data
  	
  	ref <- reference.data
  	
  	#Exclude metadata variables

  	mir1 <- as.matrix(mir[,-1])
  	
  	wave <- as.numeric(substr(colnames(mir1),2,19))
  	
  	colnames(mir1) <- wave
  	
  	#First derivative
  	
  	de1 <- trans(mir1,tr = "derivative",order = 1,gap = 23)
  	
  	der1 <- rev(as.data.frame(de1$trans))
  	
  	colnames(der1) <- paste0("m",wave)
  	
  	#Save derivative spectra
  	
  	der1.ssn <- as.data.frame(cbind(as.vector(mir[,1]),der1))
  	
  	colnames(der1.ssn) <- c("SSN",colnames(der1))
  	
  	write.table(der1.ssn,file = "First derivative.csv",

     sep = ",",row.names = FALSE)
  	
  	#merge with first derivative preprocessed spectra
  	
  	ref.mir <- merge(ref,der1.ssn,by.x = "SSN",by.y = "SSN")
  	
  	rc <- colnames(ref)
  	
  	#which columns contains reference data?
  	
  	ref<-ref.mir[,rc]
  	
  	#Extract spectral predictors
  	
  	mirp<-colnames(der1.ssn)[-1]
  	
  	spectra<-ref.mir[,mirp]
  	
  	#Create two new subfolders within the current working using:
  	
  	b<-getwd()
  	
  	if(!file.exists("Models")){dir.create("Models")}
  	
  	if(!file.exists("calibration_plots")){dir.create("calibration_plots")}
  	
  	# Fit calibration models for the training set and
  	
  	# use the testing set to validate the models

  	set.seed(67523)
  	
  	testing <- which(ref.mir$SSN%in%hout$SSN) #with hout
  	
  	#Use Kennard_Stone.
  	
  	# This is an optional step just to show distribution of spectra in a PCA space.
  	
  	sel <- kenStone(spectra,k = round(0.33*nrow(spectra)),pc = .99)
  	
  	# To view selected samples, remove "#" below two lines to plot
  	
  	# plot(sel$pc[,1:2],xlab = 'PC1',ylab = 'PC2')
  	
  	# points(sel$pc[sel$model,1:2],pch = 19,col = 2)
  	
  	# points selected for calibration
  	
  	#Loop for calibration of all soil properties in the reference set starts here
  	
   msummary <- NULL

   hd <- colnames(ref)[-1]

   for (q in 1:length(hd)){

    refq <- which(colnames(ref)%in%hd[q])

    ref.q <- ref[,refq]

    pms.a <- NULL

    pred.all <- NULL

    cal <- cbind(as.vector(ref.q),spectra)[-testing,]

    val <- cbind(as.vector(ref.q),spectra)[testing,]

    colnames(cal) <- c(colnames(ref)[refq],colnames(spectra))

    colnames(val) <- colnames(cal)

    cal <- na.omit(cal)

    val <- na.omit(val)

    trainX <- cal[, -1]

    set.seed(100)

    colnames(cal) <- c("trainY", colnames(trainX))

    cal[,"trainY"] <- log(cal[,"trainY"])

    indx <- createFolds(cal[,"trainY"], returnTrain = TRUE)

    ctrl <- trainControl(method = "cv", index = indx)

    rf.m <- train(trainY~., method = "pls", data = cal,trControl =

      ctrl,tuneGrid = expand.grid(ncomp = 1:10),metric = "RMSE",preProc = 

      c("center", "scale"))

		# Get final model to compute coefficient for variation explained
		
		predi <- exp(predict(rf.m,rf.m$trainingData))
		
		y <- exp(cal[,"trainY"])
		
		#computes RMSE and R-squared values for the calibration set

		training.parameters <- round(postResample(predi,y),2)
		
		RSQ <- training.parameters[2]
		
		RMSE <- training.parameters[1]
		
		# Predict qth soil property of the holdoutset using
		
		# the MIR data and compare with the actual measurement
		
		predi.test <- exp(predict(rf.m,val[,-1]))
		
		y.test <- val[,1]
		
		#Get PCs used
		
		PCs <- rf.m$finalModel$ncomp
		
		#computes RMSE and R-squared values for the validation set

		testing.parameters <- round(postResample(predi.test,y.test),2)
		
		RSP <- testing.parameters[2]
		
		RMSEP <- testing.parameters[1]
		
		model.summary <- c(hd[q],PCs,training.parameters,testing.parameters)
		
		msummary <- rbind(msummary,model.summary)
		
		saveRDS(rf.m,file = paste0(b,"/","models/",hd[q],".rds"))
		
		pm <- as.data.frame(cbind(y,predi))
		
		colnames(pm) <- c("measured","predicted")
		
		# Create scatter plot for the predicted versus the measured - training data set
		
		p <- ggplot(pm, aes(x = measured,y = predicted)) + 
		
		geom_point(col = "black",size = 3,alpha = 0.5) + 
		
		ggtitle(paste0("Calibration for ",hd[q])) + 
		
		xlab("Measured") + 
		
		ylab("Predicted")
		
		p <- p + stat_smooth(method = lm, se = FALSE, color = 'black',alpha = 0.1)
		
		p <- p + theme(plot.title = element_text(lineheight = 3, face = "bold",

      color = "black", size = 20))
		
		 # this will change all text size 

     p <- p + theme(text = element_text(size = 20))

     p <- p + annotate('text', label = paste('R^2 == ',RSQ),

      parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8)  + 

     annotate('text', label = paste('RMSE == ',RMSE), 

      parse = TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)

		# Centre title
		
   p <- p + theme(plot.title = element_text(hjust  = 0.5))

   p <- p + xlim(range(pm)) + ylim(range(pm))

		#Validation data
		
		pmp <- as.data.frame(cbind(y.test,predi.test))
		
		colnames(pmp)<-c("measured.test","predicted.test")
		
		# Create scatter plot for the predicted versus the measured
		
		# the validation set
		
		p2 <- ggplot(pmp, aes(x = measured.test,y = predicted.test)) + 
		
		geom_point(col = "brown",size = 3,alpha = 0.5) + 
		
		ggtitle(paste0("Validation for ",hd[q])) + 
		
		xlab("Measured") + 
		
		ylab("Predicted")
		
		p2 <- p2 + stat_smooth(method = lm, se = FALSE, color = 'brown',

      alpha = 0.1)
		
		p2 <- p2 + theme(plot.title = element_text(lineheight = 3,

      face = "bold", color = "black", size = 20))
		
		# this will change all text size 
		
		p2 <- p2 + theme(text = element_text(size = 20))

		p2 <- p2 + annotate('text', label = paste('R^2 == ',RSP),

      parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8)  +
		
		annotate('text', label = paste('RMSE == ',RMSEP),

      parse = TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)
		
		# Centre title
   p2 <- p2 + theme(plot.title = element_text(hjust  = 0.5))

   p2 <- p2 + xlim(range(pmp)) + ylim(range(pmp))

		# Save calibration and validation plots
    
		png(file = paste0(b,"/Calibration_plots/",hd[q],".png"),

      height = 400,width = 800)
		
		grid.arrange(p,p2,nrow = 1)
		
		dev.off()
		
  }

  colnames(msummary) <- c("Soil properties","PCs","LOOCV RMSEC",

    "LOOCV Rsquared", "Holdout RMSEP","Holdout Rsquared")

  write.table(msummary,file = "Model_Summary.csv",sep = ",",row.names = FALSE)

      # All Samples
      
      b<-getwd()
      
      if(!file.exists("Full_Models")){dir.create("Full_Models")}
      
      if(!file.exists("Full_calibration_plots")){dir.create("Full_calibration_plots")}
      
      # Begin calibration 
      
      msummary<-NULL
      
      hd<-colnames(ref[,-1])#Exclude SSN 
      
      all.predicted<-NULL
      
      for (q in 1:length(hd)) {
      	
      	refq<-which(colnames(ref)%in%hd[q])
      	
      	ref.q<-ref[,refq]
      	
      	cal<-cbind(as.vector(ref.q),spectra)
      	
      	cal<-na.omit(cal)
      	
      	trainX <-cal[, -1]
      	
      	colnames (cal) <- c("trainY",colnames(trainX))
      	
      	cal[,"trainY"] <-log(cal[,"trainY"])
      	
      	#colnames(cal)<-c(colnames(ref)[refq],colnames(spectra))
      	
      	p<-which(is.na(der1.ssn[,1]) == TRUE)
      	
      	ifelse(length(p)>0,ssn<-der1.ssn[-p,1],ssn <- der1.ssn[,1])
      	
      	ifelse(length(p)>0,der1.ssn<-der1.ssn[-p,],der1.ssn<-der1.ssn)
      	
      	#Select training and testing sets
      	
      	set.seed(100)
      	
      	indx <- createFolds(cal[,"trainY"], returnTrain = TRUE)
      	
      	ctrl <- trainControl(method = "cv", index = indx)
      	
      	rf.m <- train(trainY~., method = "pls", data = cal,

         trControl = ctrl,tuneGrid = expand.grid(ncomp = 1:10),

         metric = "RMSE",preProc = c("center", "scale"))
      	
      	#Save the model
      	
      	saveRDS(rf.m,file = paste0(b,"/","Full_Models/",hd[q],".rds"))
      	
      	#Get final model to compute coefficient for variation explained
      	
      	predi <- exp(predict(rf.m,rf.m$trainingData))
      	
      	y <- exp(cal[,1])
      	
      	#Get PCs used
      	
      	PCs <- rf.m$finalModel$ncomp
      	
      	training.parameters <- c(hd[q],PCs,round(postResample(predi,y),3))
      	
      	RSQ <- round(as.numeric(training.parameters[4]),2)
      	
      	RMSE <- round(as.numeric(training.parameters[3]),2)
      	
      	msummary <- rbind(msummary,training.parameters)
      	
      	#Training
      	
      	pm <- as.data.frame(cbind(y,predi))
      	
      	colnames(pm) <-c ("measured","predicted")
      	
      	png(file = paste0(b,"/","Full_calibration_plots/",hd[q],".png"),

         height =  600,width = 600)
      	
      	p1 <- ggplot(pm, aes(x = measured,y = predicted)) + 
      	
      	geom_point(col = "brown",size = 5,alpha = 0.5) + 
      	
      	ggtitle(paste0("Calibration for ",hd[q])) + 
      	
      	xlab("Measured") + 
      	
      	ylab("Predicted")
      	
      	p1 <- p1 + stat_smooth(method = lm, se = FALSE, color = 'brown',

         alpha = 0.1,size = 1.5) + 
      	
      	theme(plot.title = element_text(lineheight = 3, 

         face = "bold", color = "black", size = 20)) + 
      	
      	# this will change all text size 

      	theme(text = element_text(size = 20)) +
      	
      	annotate('text', label = paste('R^2 == ',RSQ),

         parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8,size = 5) +
      	
      	annotate('text',label = paste('RMSE == ',RMSE), 

         parse = TRUE,Inf, -Inf,hjust = 2.0, vjust = -6.8,size = 5) + 
      	
      	annotate('text', label = paste('PCs == ',PCs), 

         parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -3.9,size = 5)
      	
      	# Centre title
      	
      	p1 <- p1 + theme(plot.title = element_text(hjust  = 0.5))
      	
      	# Create scatter plot for the predicted versus the measured 
      	
      	# the combined dataset
      	
      	p1 <- p1 + xlim(range(pm)) + ylim(range(pm))

      	ggsave(file = paste0(b,"/","Full_calibration_plots/",hd[q],".png"),

         height = 6, width = 6, units = "in", p1)
      	
      	prediction.f <- round(exp(predict(rf.m,der1.ssn[,-1])),2)
      	
      	all.predicted <- cbind(all.predicted,prediction.f)
      	
      }

      	#Combine the predicted values together
      	
      	all.predicted.SSN <- cbind(as.vector(ssn),all.predicted)
      	
      	colnames(all.predicted.SSN) <- c("SSN",hd)
      	
      	colnames(msummary)<-c("Soil properties","PCs","RMSEC","Rsquared")
      	
      	#Save full model summaries
      	
      	write.table(msummary, file = "Full models summary.csv",sep = ",",

         row.names = FALSE)
      	
      	#Save the linked file
      	
      	write.table(all.predicted.SSN, file = "All predictions.csv",

         sep = ",", row.names = FALSE)

      }

    }
   ################## end of PLS and Rf calibration method #########################