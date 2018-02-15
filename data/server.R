library(shinyjs)

library(shiny)

library(shinyFiles)

library(DT)

library(caret)

library(jsonlite)

server  <- shinyServer(function(input, output,session) {

observeEvent(list(input$file,input$alphadt), {

  avgd.std = reactive({switch(input$alphadt,"Mua" = mua.avg,"Whitesand" = wtsd.avg)})# Read aggregated spectrums 
  
  volumes <- getVolumes()
  
  folderInput  <-  reactive({
  	
    shinyDirChoose(input, 'directory', roots  =  volumes, session  =  session)
    
    return(parseDirPath(volumes, input$directory))
  
  }
  
  )
  	
 output$directorypath  =  renderPrint({  folderInput()  })

	
 output$file  =  renderDataTable({
 
      
        lsf  <-  as.vector(list.files(path  =  folderInput(), pattern  =  "*.0", full.names  =  FALSE))

    spectra.a  <-  NULL
      
	for (k in 1:100){
    
    spectra <-  read.opus(paste0(folderInput(),"/",lsf[k]))
      
      spectra.a  <-  rbind(spectra.a,spectra)
	}
    
    spectra  <- as.data.frame(spectra.a[,-c(2:7)])
  	
  	spectra
   
   }
   
   )
    
   observe({
	
	volumes  <- c(home = "~/")
	
	shinyFileSave(input, "save",roots = volumes, session = session)
	
	fileinfo  <- parseSavePath(volumes, input$save)

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

  volumes <- getVolumes()

  shinyFileChoose(input,'file', session = session,roots = volumes)

    inFile  <-  parseFilePaths(roots = volumes, input$file)
    
    path <- as.character(inFile$datapath)
    
    output$table1 <- renderText({
    	path }
    	)

   output$contents  <-  renderPlot(
   {
    
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
  	

    plot(wvnbs,avgd.stddt, type =  "l", col  =  "red", lwd  =  0.8,xlim = rev(range(wvnbs)),ylim = range(avgd.stddt,spectra),
         xlab = expression(paste("Wavenumbers cm^1")),ylab = "Absorbance",main =  paste0 (" Overlay of " ,input$alphadt, " reference with new spectrum(s)"))
         
    lines(wavenumbers,spectra,col = "black",lwd  =  0.6)
    
    hqi  <-  cor(avgd.stddt , t(spectra))
    
    legend("bottomright",c("Master","Test"),bty  =  "n", col  =  c("red", "black"),lty  =  1,lwd  =  0.5)
    
    
   	}
   	)
   
   	   
   ## ... QR_Codes   	
  output$qrcodes <-  renderDataTable({
		if(n>0)
		
	{
	qrcodes<-c(1:input$n)
	for (i in 1:input$n){
	qrcodes[i]<-paste0(c(tolower(substr(input$study,1,3)),sample(letters,3), sample(c(0:9),2),sample(letters,3),sample(c(1:8),2)),collapse="",sep="")}
	}
	
qrcodes<-matrix(qrcodes,ncol=1)
colnames(qrcodes)<-"Codes"
data<-as.data.frame(qrcodes)
  
    volumes <- getVolumes()
      
  folderInput  <-  reactive({
  	
    shinyDirChoose(input, 'directory', roots  =  volumes, session  =  session)
    
    return(parseDirPath(volumes, input$directory))
  
  }
  )
      volumes <- getVolumes()

    shinyFileSave(input, "save", roots=volumes, session=session)
    fileinfo <- parseSavePath(volumes, input$save)
   data2<-data
    if (nrow(fileinfo) > 0) {
    	
    	nsf <- paste0(as.character(fileinfo$datapath),"/png_files")
   	dir.create(nsf)
	#dir.create("png_files")
      write.table(data2, as.character(fileinfo$datapath),row.names = FALSE)
      
    }
  
data

	}
	

)

#Kennard_Stone

options(shiny.maxRequestSize=150*1024^2)
volumes = getVolumes()

  output$contents <- renderTable({
  
    inFile <- input$file1

    if (is.null(inFile))
    return(NULL)
    
    all.spec<-read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
	wavenumbers<-round(as.numeric(substr(colnames(all.spec[,-1]),2,19)),1)
	colnames(all.spec)<-c("SSN", wavenumbers)
	all.spec[1:2,-1]#Exclude SSN
	}
	)
	
  output$plot1<-	renderPlot({
  	 inFile <- input$file1

    if (is.null(inFile))
      return(NULL)
    
    all.spec<-read.csv(inFile$datapath, header=input$header, sep=input$sep, 
				 quote=input$quote)
	wavenumbers<-round(as.numeric(substr(colnames(all.spec[,-1]),2,19)),1)
	colnames(all.spec)<-c("SSN", wavenumbers)


  	plot(wavenumbers,all.spec[1,-1],type="l",lwd=2, col="brown", ylim=range(all.spec[,-1]))
  	
  	for (k in 4:9){
  		lines(wavenumbers,all.spec[k,-1], col=k)
  	}
  	}
  	)
  	
  # Use Kennard_Stone
  output$plot2<-	renderPlot({
  	 inFile <- input$file1

    if (is.null(inFile))
      return(NULL)
    
    all.spec<-read.csv(inFile$datapath, header=input$header, sep=input$sep, 
				 quote=input$quote)
	wavenumbers<-round(as.numeric(substr(colnames(all.spec[,-1]),2,19)),1)
	colnames(all.spec)<-c("SSN", wavenumbers)

spectra<-all.spec[,-1]
library(prospectr)
sel <- kenStone(spectra,k=round(input$perc*nrow(spectra)),pc=.99)
#View selected samples
plot(sel$pc[,1:2],xlab='PC1',ylab='PC2')
points(sel$pc[sel$model,1:2],pch=19,col=2) # points selected for calibration

  #Selected reference by Kennard_Stone
output$selected <- renderTable({
  
    inFile <- input$file1

    if (is.null(inFile))
    return(NULL)
    
    #all.spec<-read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
    selc<-cbind(1:length(sel$model),as.vector(all.spec[sel$model,1]))#SSN for selected samples
    colnames(selc)<-c("No","Reference SSN")
    selc<-as.data.frame(selc)
   	selc
   
	}
	
)	
 # Save  sel
observe({
	volumes <-c(home="~/")
	shinyFileSave(input, "save",roots=volumes, session=session)
	fileinfo <-parseSavePath(volumes, input$save)
	
}
)
  
  }
  )


  }
   	)
	}
	)



