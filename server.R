library(shiny)
library(shinyjs)
library(prospectr)
library(shinyFiles)
library(DT)
library(caret)
library(httr)
#library(jsonlite)
library(xlsx)
library(soil.spec)
library(reshape2)
library(ggplot2)
library(dplyr)
library(plotly)
library(BalancedSampling)


server  <- shinyServer(function(input, output,session) {

	observeEvent(list(input$file,input$diropus,input$alphadt), {
 
   	volumes.opus <- getVolumes()
  
  	folderInput.opus  <-  reactive( {
  	
    	shinyDirChoose(input, 'diropus', roots  =  volumes.opus, session  =  session)
    
    	return(parseDirPath(volumes.opus, input$diropus))
  
  	}

  	)
########################################################################################	  	
 # Simulate work being done for 1 seconds

  Sys.sleep(1)

# Hide the loading message when the rest of the server function has executed

  hide(id = "loading-content", anim = TRUE, animType = "fade")  

########################################################################################	

 	output$fileopus  =  renderDataTable({
 
      
        lsf  <-  as.vector(list.files(path  =  folderInput.opus(), pattern  =  "*.0", full.names  =  FALSE))

    	spectra.a  <-  NULL
   
      
  		for (k in lsf){
    
    		spectra <-  read.opus(paste0(folderInput.opus(),"/",k))
      
    		spectra.a  <-  rbind(spectra.a,spectra)
    		
  			}
    
    	spectra  <- as.data.frame(spectra.a[,-c(2:7)])
    	
   
  		volumes.op<-getVolumes()
  		
  		shinyFileSave(input, "saveop",roots = volumes.op, session = session)
  
  		fileinfo  <- parseSavePath(volumes.op, input$saveop)
    	
    	spectra2<-spectra
  		
  		if (nrow(fileinfo) > 0) {
  			
      		write.csv(  spectra2, as.character(fileinfo$datapath),row.names = FALSE)
      		
      		}
      		
      	spectra
      	
      	}
      	
      	)
      	
      	output$raw0  <-  renderPlot({
      		
      		if (nrow(spectra)<1)
      		
      		return(NULL)
      		
      		wavenumbers  <-  as.numeric(substr(colnames(spectra.a[,-c(1:6)]),2,19))
      		
      		spectra  <- as.data.frame(spectra.a[,-c(1:7)])
      		
      		wavenumbers  <-  as.numeric(substr(colnames(spectra),2,19))
      		
      		plot(wavenumbers, spectra[1,],main  =  "Raw spectra", ylim  =  range(spectra), xlim  =  rev (range(wavenumbers)))
      		
      		for ( i in 2:nrow(spectra)){
      			
      			lines (wavenumbers, spectra[i,], col  =  i)}
      			
      			}
      			
      			)
      			
###########################################################
# Tab 2 Exploratory Analysis
###########################################################

### 2.1 Visualize converted OPUS files - raw spectra

options(shiny.maxRequestSize=150*1024^2)

volumes = getVolumes()

  output$view.raw <- renderTable({
  
    infile0 <- input$file0

    if (is.null(infile0))

    return(NULL)
    
    all.spec <- read.csv(infile0$datapath)
  
  wavenumbers<-round(as.numeric(substr(colnames(all.spec[,-1]),2,19)),1)
  
  colnames(all.spec)<-c("SSN", wavenumbers)
  
  all.spec[1:10,-1]

  }
  
  )
  
    output$raw.plot<-  renderPlot({
    
    infile0 <- input$file0

    if (is.null(infile0))
      return(NULL)
    
    all.spec<-read.csv(infile0$datapath)
    
  	wavenumbers<-round(as.numeric(substr(colnames(all.spec[,-1]),2,19)),1)
  	
  	colnames(all.spec)<-c("SSN", wavenumbers)

  	spec <- melt(all.spec, id.vars = c("SSN"))

  	p0 <- ggplot(data = spec, aes(x = as.numeric(as.vector(variable)),y = value,group = SSN)) +
    
    geom_line(size = 0.14,colour = "red", alpha = 0.6) +
    
    ggtitle("Raw MIR spectra") +
    
    xlab(expression("Wavenumbers cm"^-1)) +
    
    ylab("Absorbance") +
    
    theme_bw() +
  	
  	theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=32)) +
    
    theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22)) +
    
    theme(
        plot.background = element_blank()
        
        ,panel.grid.major = element_blank()
        
        ,panel.grid.minor = element_blank()
    )
    
  p0 <- p0 + theme(plot.title = element_text(hjust = 0.5))

  p0

    }
    
    )


# PCA start here???? To revisit like Uhuru!

options(shiny.maxRequestSize=150*1024^2)

volumes = getVolumes()

  output$pca.raw <- renderTable({
  
    inFile <- input$file1

    if (is.null(inFile))

    return(NULL)
    
    all.spec<-read.csv(inFile$datapath)
  
  wavenumbers<-round(as.numeric(substr(colnames(all.spec[,-1]),2,19)),1)
  
  colnames(all.spec)<-c("SSN", wavenumbers)
  
  all.spec[1:10,-1]

  }
  
  )
  
    output$pca.plot1<-  renderPlot({
    
    inFile <- input$file1

    if (is.null(inFile))
      return(NULL)
    
    all.spec<-read.csv(inFile$datapath)
    
  	wavenumbers<-round(as.numeric(substr(colnames(all.spec[,-1]),2,19)),1)
  	
  	colnames(all.spec)<-c("SSN", wavenumbers)

  	spec <- melt(all.spec, id.vars = c("SSN"))

  	p <- ggplot(data = spec, aes(x = as.numeric(as.vector(variable)),y = value,group = SSN)) +
    
    geom_line(size = 0.14,colour = "red", alpha = 0.6) +
    
    ggtitle("Raw MIR spectra") +
    
    xlab(expression("Wavenumbers cm"^-1)) +
    
    ylab("Absorbance") +
    
    theme_bw() +
  	
  	theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=32)) +
    
    theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22)) +
    
    theme(
      
        plot.background = element_blank()
        
        ,panel.grid.major = element_blank()
        
        ,panel.grid.minor = element_blank()
    )
    
  p <- p + theme(plot.title = element_text(hjust =0.5))

  p

    }
    
    )
    
# 2.2 Preprocess using first derivative and visualize the preprocessed spectra

	output$pca.plot2<-  renderPlot({
    
     inFile <- input$file1
     
     if (is.null(inFile))
     
     return(NULL)
    
    all.spec<-read.csv(inFile$datapath)
    
    wavenumbers<-round(as.numeric(substr(colnames(all.spec[,-1]),2,19)),1)
    
    colnames(all.spec)<-c("SSN", wavenumbers)
    
    spectra <- as.matrix(all.spec[,-1])
    
    der <- as.data.frame(trans(spectra)$trans)
    
    all.spec [,-1] <- rev(der)
    
    spec <- melt(all.spec, id.vars = c("SSN"))
    
    p2 <- ggplot(data = spec, aes(x = as.numeric(as.vector(variable)),y = value,group = SSN)) +
    
    geom_line(size = 0.14,colour = "red", alpha = 0.6) +
    
    ggtitle("Preprocessed MIR spectra") +
    
    xlab(expression("Wavenumbers cm"^-1)) +
    
    ylab("Absorbance") +

    theme_bw() +
    
    theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=32)) +
    
    theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22)) +
    
    theme(
    
        plot.background = element_blank()
    
        ,panel.grid.major = element_blank()
    
        ,panel.grid.minor = element_blank()
    )
    
	p2 <- p2 + theme(plot.title = element_text(hjust =0.5))
	
	p2

    }

    	)




#####################################################################################################



  methodnames<-c('Savtizky_Golay (1st Der)','Standard Normal Variate (SNV)','Multiplicative Scatter Correction (MSC)')
  output$preprospecsel <- renderUI({checkboxGroupInput("specsel", label = "Select a preprocessing Method", choices=methodnames)})



prepromethod <- eventReactive(input$specsel,{prepromethodinput<-as.vector(input$specsel)})


###################################################################################################################









# 2.3 Do a PCA

	output$pca.plot3<-  renderPlot({
    
	inFile <- input$file1

	if (is.null(inFile))
	
	return(NULL)
    
    all.spec<-read.csv(inFile$datapath)
    
    wavenumbers<-round(as.numeric(substr(colnames(all.spec[,-1]),2,19)),1)
    
    colnames(all.spec)<-c("SSN", wavenumbers)
    
    spectra <- as.matrix(all.spec[,-1])
    
    der <- as.data.frame(trans(spectra)$trans)
    
    pc  <- as.data.frame(prcomp(der)$x)[,1:10]

	
	# 2.3.1 Interactive outlier selection with mouse selection

	plot (pc[,1], pc[,2], pch = 16, col = "blue", xlab = "PC1", ylab = "PC2")
	
	}
	
	)
 
		# output (x, y) of the clicked point
		output$clicked <- renderPrint({
			
        c("input$clickPoint$x" = input$clickPoint$x, "input$clickPoint$y" = input$clickPoint$y)
        
        }
        
        )

		output$clicked <- renderPrint({
		
		cat("Click:\n")
		
		str(input$plot_click)
		
		}
		
		)

 # Selected outliers interactively; revisit
 
	output$outliers <- renderTable({
  
    inFile <- input$file1

    if (is.null(inFile))
    
    return(NULL)
    
    selc<-sample(all.spec[,1],0.1*nrow(all.spec[,1]))
    
    selc
   
  	}
  	
	)
  
	#Save  sel
	
	observe({
		
		volumesraw <-c(home="~/")
		
		shinyFileSave(input, "saveraw",roots=volumes, session=session)
		
		fileinfo <-parseSavePath(volumesraw, input$saveraw)
		
		pdf(file = paste0(fileinfo,".pdf"))
		
		dev.off()
		
		volumesprocessed <-c(home="~/")
		
		shinyFileSave(input, "saveprocessed",roots=volumes, session=session)
		
		fileinfo <-parseSavePath(volumesprocessed, input$saveprocessed)
		
		}
		
		) 
###########################################################      			
# Tab 3 Data quality
###########################################################

     
## 3.1 Internal standards QC
        
## read aggregated spectrums 
      
      avgd.std = reactive({switch(input$alphadt,"Mua" = mua.avg,"Whitesand" = wtsd.avg)})

      volumes.raw <- getVolumes()
      	
      shinyFileChoose(input,'file', session = session,roots = volumes.raw)
      	
      inFile  <-  parseFilePaths(roots = volumes.raw, input$file)
      	
      path <- as.character(inFile$datapath)
      	
      output$table1 <- renderText({path})
      	
      output$contents  <-  renderPlot({
      		
      avgd.stddt <- avgd.std()[,-1]
      		
      wvnbs <- as.numeric(substr(colnames(avgd.stddt),2,19))
      		
      if (is.null(inFile))
      		
      return(NULL)
      		
      spectra <- read.opus(path[1])
      		
      wavenumbers  <-  as.numeric(substr(names(spectra[-c(1:6)]),2,19))
      		
      z <- length(wavenumbers)
      		
      wavenumbers <- (wavenumbers[-z])
      		
      spectra  <- as.numeric(spectra[-c(1:6)])
      		
      spectra <- spectra[-z]
      		
      plot(wvnbs,avgd.stddt, type =  "l", col  =  "red", lwd  =  0.8,xlim = rev(range(wvnbs)),ylim = range(avgd.stddt,spectra),xlab = expression(paste("Wavenumbers cm^1")),ylab = "Absorbance",main =  paste0 (" Overlay of " ,input$alphadt, " reference with new spectrum(s)"))
         
      lines(wavenumbers,spectra,col = "black",lwd  =  0.6)
         
      hqi  <-  cor(avgd.stddt , t(spectra))
         
      legend("bottomright",c("Master","Test"),bty  =  "n", col  =  c("red", "black"),lty  =  1,lwd  =  0.5)
         
       }
         
       )
         
       }
         
       )
       
 ## 3.2 Replicate spectrum(s)
 
######################################
# Tabpanel 4.Tools
######################################

## Under the tools  tab add wrapper functions for:
 
### 4.1. QR_Codes  
  
	output$qrcodes <-  renderDataTable({
		
		if(n>0)
		
			{
			qrcodes<-c(1:input$n)
			
			for (i in 1:input$n){
				
				qrcodes[i]<-paste0(c(tolower(substr(input$study,1,3)),sample(letters,3), sample(c(0:9),2),sample(letters,3),sample(c(1:8),2)),collapse="",sep="")
				
				}

			}
  
	qrcodes<-matrix(qrcodes,ncol=1)
	
	colnames(qrcodes)<-"Codes"
	
	data<-as.data.frame(qrcodes)
  
	volumes.qr <- getVolumes()
	
	shinyFileSave(input, "saveqr",roots = volumes.qr, session = session)
	
	fileinfo<-parseSavePath(volumes.qr, input$saveqr)
	
	data2<-data
	
    if (nrow(fileinfo) > 0) {
    	
    	write.table(data2, as.character(fileinfo$datapath),row.names = FALSE)
    	
    	}
    	
    data  

}

)


### 2. Kennard starts here

	options(shiny.maxRequestSize=150*1024^2)
	
	volumes.ks=getVolumes()

	output$contents2 <- renderTable({
  
    inFile <- input$filek

    if (is.null(inFile))
    
    return(NULL)
    
    all.spec<-read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
    
    wavenumbers<-round(as.numeric(substr(colnames(all.spec[,-1]),2,19)),1)
    
    colnames(all.spec)<-c("SSN", wavenumbers)
    
    set.seed(19802)
    
    b <- round(0.1*nrow(all.spec))
         
    b <- sample(nrow(all.spec), b)
    
    all.spec[b,]#Exclude SSN
    
    }
    
    )
    
    output$plot1 <- renderPlot({
    	
    	 inFile <- input$filek
    	 
    	 if (is.null(inFile))
    	 
    	 return(NULL)
    	 
    	 all.spec<-read.csv(inFile$datapath, header=input$header, sep=input$sep, 
         quote=input$quote)
         
         wavenumbers<-round(as.numeric(substr(colnames(all.spec[,-1]),2,19)),1)
         
         colnames(all.spec)<-c("SSN", wavenumbers)
         
         plot(wavenumbers,all.spec[1,-1],type="l",lwd= 0.45, col="purple", ylim=range(all.spec[,-1]),xlab=expression("Wavenumbers cm"^-1),ylab="Absorbance")
         
         # Select 10% of the spectra randomly for visualization
         
         set.seed(19802)
         
         b <- round(0.1*nrow(all.spec))
         
         b <- sample(nrow(all.spec), b)
         
         allspec <- all.spec[b,]
         
         for (k in 1:nrow(allspec)){
         	
         	lines(wavenumbers,allspec[k,-1], col="purple", lwd = 0.45)
         	
         	}
         	
         	}

         	)
    
  # Use Kennard_Stone algorithm for selection 
  
	output$plot2<-  renderPlot({
  	
  	inFile <- input$filek

    if (is.null(inFile))
      
      return(NULL)
    
    all.spec<-read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
    
    wavenumbers<-round(as.numeric(substr(colnames(all.spec[,-1]),2,19)),1)
    
    colnames(all.spec)<-c("SSN", wavenumbers)
    
    spect<-all.spec[,-1]

#library(prospectr)
#library(soil.spec)
	
	spectra <- trans(spect)$trans
	
	sel <- kenStone(spectra,k=round(input$perc*nrow(spectra)),pc=.99)

	#View selected samples
	
	plot(sel$pc[,1:2],xlab='PC1',ylab='PC2')
	
	points(sel$pc[sel$model,1:2],pch=19,col=2) # points selected for calibration

	#Selected reference by Kennard_Stone
	
	output$selected <- renderTable({
  
    inFile <- input$filek

    if (is.null(inFile))
    
    return(NULL)
    
    
    selc<-cbind(1:length(sel$model),as.vector(all.spec[sel$model,1]))#SSN for selected samples
    
    colnames(selc)<-c("No","Reference SSN")
    
    selc<-as.data.frame(selc)

	shinyFileSave(input, "saveks", roots=volumes.ks, session=session)
 
	fileinfo <- parseSavePath(volumes.ks , input$saveks)
    
    selc2<-selc

    if (nrow(fileinfo) > 0) {
      
      write.table(selc2, as.character(fileinfo$datapath),row.names = FALSE, sep = ",")
      
      }

	selc
	
	}
	
	)
	
	}
	
	)
	
### 3. Fertilizer screening here
# ...

######################################	
# Tabpanel 5.Calibration
######################################
 
observeEvent(list(input$runaction), {

      #code here

predfiledata <- reactive({

  infile.pred <- input$ircalfile

  if (is.null(infile.pred)) {

  return(NULL)

   }

read.csv(infile.pred$datapath)

}

)

filerefdata <- reactive({

  infile.ref <- input$refcalfile

  if (is.null(infile.ref)) {

  return(NULL)

  }

  read.csv(infile.ref$datapath)

}

)

#shinyDirChoose(input, 'dir', roots = getVolumes())

dir  <-  reactive({
  	
    shinyDirChoose(input, 'dir', roots  =  getVolumes(), session  =  session)
    
    return(parseDirPath(getVolumes(), input$dir))
  
}
  
)
  
 
infiledir<-dir()
 
methodInput <- reactive({input$model}) 

actiondata <- eventReactive(input$runaction,{

  infile <- input$runaction

  if (is.null(infile)) {

  return(NULL)

  }
    
  output$ircalfiletable <- renderDataTable({

  predfiledata()[1:6,1:6]

	}

	)
	
output$refcalfiletable <- renderDataTable({

  filerefdata()[1:6,]

	}

	)
		
 #Calibration using PLS or RF methods begins here.

	#Select training and testing set using Balance sampling

set.seed(12345);

	if(nrow(input$refcalfile)>0){
    
	#Get reference variable with highest variance

	sdx <- apply(na.omit(filerefdata()[,-1]),2,sd)

	mdx <- apply(na.omit(filerefdata()[,-1]),2,mean)

	dsp <- as.vector(which(sdx/mdx==max(sdx/mdx)))+1

	n = round(0.7*nrow(filerefdata())); # sample size

	p = rep(n/nrow(filerefdata()),nrow(filerefdata())); # inclusion probabilities

	s = cube(p,as.matrix(filerefdata()[,dsp])); # select sample

	pkfr <- as.vector(c(1:nrow(filerefdata())))

	test <- which(!pkfr%in%s)

	hout<-filerefdata()[test,]

	}

	calibrate(infiledir,predfiledata(),filerefdata(),hout,method = methodInput())
     
  }

  )

actiondata()

}

)
  
###########################################################
# Tabpanel 6. Prediction
###########################################################


observeEvent(list(input$predaction,input$comparedata,input$dirmodels), {    

predfiledata <- reactive({

  infile.pred <- input$predfile

  if (is.null(infile.pred)) {

  return(NULL)

  }

  read.csv(infile.pred$datapath)

}

)


dirmodel  <-  reactive({
    
  shinyDirChoose(input, 'dirmodels', roots  =  getVolumes(), session  =  session)
    
  return(parseDirPath(getVolumes(), input$dirmodels))
  
  }
  
  )
  

infiledirmodel<-dirmodel()
  
comparepredactiondata <- eventReactive(input$comparedata,{

infile <- input$comparedata

if (is.null(infile)) {

return(NULL)

}
    
modelpaths<-list.files(path=infiledirmodel, pattern="*.rds", full.names  =  FALSE)

  varn <- gsub(".rds", "",modelpaths)

  spectradata<-readRDS(paste0(infiledirmodel,"/",varn[[1]],".rds"))$trainingData

	spectradata <- prcomp(as.data.frame(spectradata[,-1]))$x[,1:2]	
    
  output$predfile <- renderPlot({

  raw0 <- as.matrix(predfiledata()[,-1])

  colnames(raw0) <-as.numeric(substr(colnames(raw0),2,19))
     
  der <- trans(raw0,tr = "derivative",order = 1,gap = 23)
  
  der <- rev(as.data.frame(der$trans))

  pc.pred <- prcomp(der)$x[,1:2]

  pc.pred <- as.data.frame(pc.pred)

	spectradata <- as.data.frame(spectradata)

	pc.pred$set <- rep("Prediction",nrow(pc.pred))

  spectradata$set <- rep("Calibration", nrow(spectradata))

  pc.data <-rbind(spectradata,pc.pred)
        
  g <- ggplot(pc.data, aes(x = PC1, y = PC2))+
      
	geom_point(size = 3.0,alpha= 0.6, aes(col = set)) +

	ggtitle("PCA scores spectra") +

	xlab("PC1") +

	ylab("PC2") +

  theme_bw() +

  theme(
        
  plot.background = element_blank()

  ,panel.grid.major = element_blank()

  ,panel.grid.minor = element_blank()

  )

  g <- g + scale_color_manual(breaks = c("Calibration", "Prediction"),values=c("blue", "red"))
    
  g <- g + theme(plot.title = element_text(lineheight = 3, face = "bold", color = "black", size = 17))
		
	g <- g + theme(text = element_text(size = 12)) # this will change all text size +
	
	theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=26)) +

  theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22)) +

  theme(axis.text=element_text(size=12)) 
   
	g + theme(plot.title = element_text(hjust = 0.5)) 
	
  }

  )
      
  }

  )

comparepredactiondata()


predactiondata <- eventReactive(input$predaction,{

  infile <- input$predaction

  if (is.null(infile)) {

  return(NULL)

  }
  
  observeEvent(input$dirmodels,{

  modelpaths1<-list.files(path=infiledirmodel, pattern="*.rds", full.names  =  FALSE, ignore.case = TRUE)

  varn1 <- gsub(".rds", "",modelpaths1, ignore.case = TRUE)

  output$make_box <- renderUI({

  checkboxGroupInput("checkbox", label = "Select parameter(s) to predict", choices=varn1)

  }

  )

}

)


modelsetInput <- eventReactive(input$checkbox,{

  modelinput<-as.vector(input$checkbox)
     
  }

  ) 

output$predfiletable <- renderDataTable({

  raw0 <- as.matrix(predfiledata()[,-1]) 

  colnames(raw0) <-as.numeric(substr(colnames(raw0),2,19))

  der <- trans(raw0,tr = "derivative",order = 1,gap = 23)

  der <- rev(as.data.frame(der$trans))

	colnames(der) <- paste0("m",colnames(der))

	modata<-reactive({

  mod<-NULL

sapply(1:length(input$checkbox), function(i) mod <- cbind(mod,round(exp(predict(readRDS(paste0(infiledirmodel,"/",modelsetInput()[i],".rds")),der)),2)))

  }

  )
 
  modf <- modata()

  SSN <- predfiledata()[,1]

  modf <- cbind(as.vector(SSN),modf)

  colnames(modf) <- c("SSN", input$checkbox)

  modf

}

)

}

)

predactiondata()

}

)
	
}

)

  

