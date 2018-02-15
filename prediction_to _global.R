
library(soil.spec)
library(tidyr)

new.spectra <- read.csv ("/Volumes/BIG CASS/test_data/pred_mir.csv")

#Preprocess the new spectra using first derivative

wavenumbers <- as.numeric(substr(colnames(new.spectra[,-1]),2,19))

colnames(new.spectra) <- c("SSN", wavenumbers)

der.new <- trans(new.spectra[,-1],tr = "derivative",order = 1,gap = 23)$trans

colnames(der.new) <-  paste0("m",colnames(der.new))


# Selects folder with the models 
model.f <- list.files("/Volumes/BIG CASS/Models", pattern = ".rds")

model.folder.path <- "/Volumes/BIG CASS/Models"

# Get the models in the folder
varn <- gsub(".rds", "",model.f)


#Extract IR data used to create models to overlay with the new spectra
mod <- readRDS(paste0(model.folder.path,"/",varn[1],".rds"))

#Extract data used for training

calib_data <- mod$trainingData[,-1]

cb <- rev(calib_data)

all <- rbind(cb, der.new)

pcs <- prcomp(all)$x[,1:10]

pcs <- as.data.frame(pcs)

pcs$set <- 0

pcs$set[1:nrow(cb)] <- "training"

pcs$set[-(1:nrow(cb))] <- "prediction"

imp <- summary(pc)$importance

p.calpred <- ggplot(pcs, aes(x = PC1, y = PC2, colour = set)) + 

geom_point(size = 1.4) + 

ggtitle("PCA scores spectra for calibration and preiction") +

xlab(paste0("PC1 explains ", round(imp[2,1],3)*100, " % total variation")) +

ylab(paste0("PC2 explains ", round(imp[2,2],3)*100, " % total variation")) +

    theme_bw() +
    
    theme(
        plot.background = element_blank()
        
        ,panel.grid.major = element_blank()
        
        ,panel.grid.minor = element_blank()
    )
    
p.calpred + theme(plot.title = element_text(hjust = 0.5))

#Selected by clicking the checkbox

varn.s <- varn[c(1,3)]

# If the PCA overlay shows the a good overlap, then proceed to predict

predic.all <- NULL

for (i in 1:length(varn.s)){

mod <- readRDS(paste0(model.folder.path,"/",varn.s[i],".rds"))

pr.i <- round(exp(predict(mod, der.new)),2)

predic.all <- cbind(predic.all, pr.i)
 }
predic.all <- cbind(as.vector(new.spectra[,1]),predic.all)
colnames(predic.all) <- c(colnames(new.spectra)[1],varn.s)

# Save predicted variables into a txt/csv file.

