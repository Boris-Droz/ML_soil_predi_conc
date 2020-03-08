####################################################################################################################################
####################################################################################################################################
###					
###  MODEL PREDICTION		    ###   
###                     	  ###					
####################################################################################################################################
####################################################################################################################################
## HISTORIC:
## --------
## v7.4 September 2016 - Boris DROZ & Gerrad Jones, ETHZ & EAWAG -- Projet Se around the world
## v8.0 Mars 2017 - Droz --- implementation of more techic
## v9.0 June 2017 - Droz --- ranger and Xboost package included

####################################################################################################################################
## DESCRIPTION
################      ####################################
#                   --- CHECK PARAMETER UNTIL LINE 150 ---
##                   #####################################

##    Tested under R 3.4 
#########################

## PREDICT MODEL ON A LARGE SCALE
##
# Input: - models calibrate
#       - spatial data variable (raster)
#       !!! RASTER OF ALL SCENARIO SHOULD HAVE SIMILAR NAME !!!!!
#       !!! Input parameter should be most similar as model calibration !!!!
# Output: 
#=======
# Need to prepare folders in order to have same folder than input
#==============================================================
# Structure: Output/Model  -> current
#                          -> future
#                          -> -.... (and other scenario)
#
####################################################################################################################################
## Set library path
#.libPaths("C:/Users/bdroz/Documents/R/win-library/3.4")

## library list
###############

# spatial analysis
library(raster)
library(gdalUtils)

# stat tool
library(ade4)

# Modelling technic
library(gam)
library(nnet)
library(neuralnet) 
library(NeuralNetTools)# importance and sensitivity analysis of nnet
library(glmnet)
library(randomForest)
library(AICcmodavg)
library(glmulti)
library(RSNNS)
library(elmNN)
library(gbm)
#library(adabag) # bagging method
#library(caret)
library(ipred) # bagging method
library(ranger)
library(xgboost)

# other
library(rms)  
library(foreign)
library(ggplot2)
library(gplots)
library(plyr)
library(Hmisc)
library(corrplot)
library(MASS)
library(AICcmodavg)
library(stringr)
library(boot)
library(pracma)
library(party)
library(h2o)

## library parrallel core
library(snowfall)
library(parallel)

####################################################################################################################################
## SCRIPT PARAMETER  --> NEED TO GO TROUH AND MAKE THE APPROPRIATE MODIFICATION !!!
#####################

# Set folder input folder content the different folder with raster for prediction
##                  if mask should be on the same input too
##           output folder content the model calibration
# 
## =========
inpath <-"D:/copper/input/pred"
output <- "D:/copper/output/20180330_model_v18_1000"

outpath <- output ## DON'T MODIFIED ##

## set your work space ## DON' MODIFIE ##
setwd(inpath) ## DON'T MODIFIED ##

# MASK OF THE PROJECTION --> raster NA or 1 to define area of prediction
##########################
mask.proj <- "YES"

mask.r <- "maskvine_wgs84_250m.tif"

## TRANSFORMATION 
## ==============  
transf.y <- alist(x=,(10^x )) # dep. variable write the back transformation

# Transformation list same number as predictive variable 
## # alist(x=,log10(abs(x))) , alist(x=,sqrt(abs(x)))) # EXEMPLE OF TRANSF.
#==============================================================

transf <- list( alist(x=,log10(abs(x))) , alist(x=,(x)) , alist(x=,(x)),  alist(x=,(1/x)), 
			 alist(x=,(x)) , alist(x=,(x)) , 
			 alist(x=,log10(abs(x))), alist(x=,log10(abs(x))),  alist(x=,(x)), alist(x=,(x)) )

length(transf)

## NUMBER OF MODEL BUILDING
# ========================
n.rep <- 1000 #should be similar or less as the calibration

## -- CHOOSE MODELLING TECHN:
## YES or NO
#############################
# machine learning method
NNET <- "NO"
RSNNS <-"NO"
ELM <- "NO"

RF <- "NO"
CSRF <- "NO"

BAG <- "YES"
EGB <- "NO"
GBM <- "NO"

# linear additive
GLM <- "NO"
GLM_sw <- "NO"
GAM <- "NO"

##### .....
####################################################################################################################################
####################################################################################################################################
#################################################################################################################################################
# FUNCTION check and produced subDir folder
###########################################
#February 2017

creat.subDir <- function (mainDir,subDir)
    {
      if (file.exists(subDir)){
        setwd(file.path(mainDir, subDir))
      } else {
        dir.create(file.path(mainDir, subDir))
        setwd(file.path(mainDir, subDir))
      }
    }
####################################################################################################################################
####################################################################################################################################

check.res <- function (r,res=res,ext=ext) # set raster, resolution and extent
          {
            empty.raster <- raster(res=res, ext=ext )
            
            if (extent(r) == ext & signif(res(r)[1],digits=7) == res) 
            { r <-r }else{
              r <- resample(r,empty.raster,method="ngb")
            }
            return(r)
          }

####################################################################################################################################
####################################################################################################################################

predict.mat.csrf <- function (r,model, cut)      
          {
            ######## FUNCTION DIVIDE A RASTER ######################
            # cut raster in egal part in function of the nb of cut
            # diffine  a function FUN to applied into a in.raster
            # put the data into an out raster !!!! should be eware of 
            ## resolution in = out otherwise dowscaling is performed!!!
            
            dx <- res(r[[1]])[1]*cut
            dy <- res(r[[1]])[2]*cut
            
            start.x <- extent(r)[1]
            start.y <- extent(r)[3]
            
            nb.loop <- 1
            
            print(round((extent(r)[2]-extent(r)[1])/dx) * round((extent(r)[4]-extent(r)[3])/dy))
            
            for (n in 1: ( round((extent(r)[2]-extent(r)[1])/dx) ) ) # dx
            {
              for (m in 1: ( round((extent(r)[4]-extent(r)[3])/dy) ) )# dy
              {
                cat("> LOOP ## NUMBER ##  ",nb.loop,"...", "\n",append = FALSE)
                
                ext <- c( (start.x + dx* (n-1)), (start.x + dx * n), # define the extent
                          (start.y+ dy* (m-1)), (start.y+ dy* m) )
                
                temp.raster <- crop(r, ext) # crop the raster into small one
                
                if ( length( temp.raster [!is.na(temp.raster[[1]] )] ) >= 1 ) # if all = NA skip the prediction
                    {
                      r.empty <- raster(res=res(temp.raster[[1]]), crs= crs(temp.raster[[1]]),
                                        ext= extent(temp.raster[[1]]))
                      
                      mat <- rasterToPoints(temp.raster)  # extract all values of the raster with point coordinate
                      
                      mat <- na.omit(mat) # so check again if empty row
                      
                      pred <- predict(model, mat[,3:ncol(mat)])$predictions # predict the point
                      
                      pred <- rasterize(mat[,1:2],r.empty,pred) # reinject into a raster
                      
                    }else { pred <- temp.raster }
                  
                  if (n == 1 & m ==1) {out.raster <- pred
                      }else {out.raster <- merge(out.raster,pred)}  
                
                nb.loop <- nb.loop + 1
                
              }  
            }
            
            return(out.raster)
        }
##################################################################################################################################

predict.mat.egb <- function (r,model, cut)      
          {
            ######## FUNCTION DIVIDE A RASTER ######################
            # cut raster in egal part in function of the nb of cut
            # diffine  a function FUN to applied into a in.raster
            # put the data into an out raster !!!! should be eware of 
            ## resolution in = out otherwise dowscaling is performed!!!
            
            dx <- res(r[[1]])[1]*cut
            dy <- res(r[[1]])[2]*cut
            
            start.x <- extent(r)[1]
            start.y <- extent(r)[3]
            
            nb.loop <- 1
            
            print(round((extent(r)[2]-extent(r)[1])/dx) * round((extent(r)[4]-extent(r)[3])/dy))
            
            for (n in 1: ( round((extent(r)[2]-extent(r)[1])/dx) ) ) # dx
            {
              for (m in 1: ( round((extent(r)[4]-extent(r)[3])/dy) ) )# dy
              {
                cat("> LOOP ## NUMBER ##  ",nb.loop,"...", "\n",append = FALSE)
                
                ext <- c( (start.x + dx* (n-1)), (start.x + dx * n), # define the extent
                          (start.y+ dy* (m-1)), (start.y+ dy* m) )
                
                temp.raster <- crop(r, ext) # crop the raster into small one
                
                if ( length( temp.raster [!is.na(temp.raster[[1]] )] ) >= 1 ) # if all = NA skip the prediction
                    {
                      r.empty <- raster(res=res(temp.raster[[1]]), crs= crs(temp.raster[[1]]),
                                        ext= extent(temp.raster[[1]]))
                      
                      mat <- rasterToPoints(temp.raster)  # extract all values of the raster with point coordinate
                      
                      mat <- na.omit(mat) # so check again if empty row
                      
                      pred <- predict(model, mat[,3:ncol(mat)]) # predict the point
                      
                      pred <- rasterize(mat[,1:2],r.empty,pred) # reinject into a raster
                      
                    }else { pred <- temp.raster }
                  
                  if (n == 1 & m ==1) {out.raster <- pred
                      }else {out.raster <- merge(out.raster,pred)}  
                
                nb.loop <- nb.loop + 1
                
              }  
            }
            
            return(out.raster)
        }



###################################################################################################################################
####################################################################################################################################
## SCRIPT START HERE
####################################################################################################################################
##############################################
# (1) Loading  datasets #
##############################################

ptm <- proc.time()# ignite timer

#beginCluster(detectCores()) # ACTIVATE THIS MULTI CORE CALCULATION 

# folder list of diff projection
proj <- list.dirs(path = inpath, full.names = FALSE, recursive = FALSE)
fn <- list.dirs(path = inpath, full.names = TRUE, recursive = FALSE)

# folder list of element
elem <- list.dirs(path = output, full.names = FALSE, recursive = FALSE)
fn.out <- list.dirs(path = output, full.names = TRUE, recursive = FALSE)

## start loop for diff element.
##############################
for (j in 1:length(fn.out))
    { 
    cat("> ",elem[j]," projection is running now...", "\n",append = FALSE)
    cat("##############################################", "\n",append = FALSE)
    
    output <- fn.out[j]
  
      # reload the rescaling parameter
      par.resc <- read.table(paste (output,"/rescal_para.txt",sep=""), header= TRUE )
      
      setwd(output)
      
      ## start loop for diff proj.
      ###########################
      for (d in 1:length(proj))
      { 
        cat("> ",proj[d]," projection is running now...", "\n",append = FALSE)
        cat("##############################################", "\n",append = FALSE)

        # STACK AND MODIF RASTER
        ########################
        
        creat.subDir (output,paste(proj[d],"_projection",sep="") )
        
        outpath <- paste(output,"/",proj[d],"_projection",sep="")
      
        ## data file (batch all file in input)
        #fn <- paste(file.path(getwd()),"/",proj[d],sep="")
        fns <- list.files(fn[d],pattern=".tif$",full.names = TRUE)
        
        r <- stack(fns) # stack all predictor
	      name.pred <- names(r)
        
        for (z in 1:length(fns))
          {
            ## APPLIED DATA MODIFICATION
            f <- as.function(transf[[z]]) # define function
            r[[z]] <- calc(r[[z]], fun=f) # applied modification
            names(r) <- name.pred # keep name!!

            ## RESCALE THE PREDICTOR 
            pred.mean <- par.resc[1,names(r[[z]]) ]
            pred.sd <- par.resc[2,names(r[[z]]) ]
            
            r[[z]] <- (r[[z]]-pred.mean)/pred.sd # function rescal
          } 
        # mask the area of prediction
       if  (mask.proj == "NO") { 
          }else{ r <- r * raster(paste(inpath,"/",mask.r,sep=""))
			      names(r) <- name.pred }
	
        writeRaster(r, filename = paste(output,"/",proj[d],"_projection/stack.rast.grd",sep=""), 
                         overwrite=TRUE)

        file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 
        
      ###############################################################################################################
        if (NNET == "NO") {}else{
           
          for (i in 1:n.rep)    
          {
            cat("> NNET MODEL run ",i,"predict  ...", "\n",append = FALSE)
            
            mod.tech <- "NNET"
            
            creat.subDir (paste(output,"/",proj[d],"_projection",sep=""),mod.tech)
          
            outpath <- paste(output,"/",proj[d],"_projection/",mod.tech,sep="")

	        	r <- stack( paste(output,"/",proj[d],"_projection/stack.rast.grd",sep="") )
        
            load(paste(output,"/",mod.tech,"/",mod.tech,"_model_",i,sep=''))
            
            pred <- predict(r, net.tr, type='raw', na.rm = TRUE)
            
            # unscale the data
            pred.mean <- par.resc[1,1]
            pred.sd <- par.resc[2,1]
            
            pred <- (pred * pred.sd) + pred.mean 
            
            # back transformed
            f <- as.function(transf.y) # define function
            pred <- calc(pred, fun=f)
            
            writeRaster(pred, filename = paste(outpath,"/",mod.tech,"_model_run_",i,".tif",sep=""), 
                        datatype="FLT8S", overwrite=TRUE)
            
            file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER    
            
            }# end loop model replicate
          } # end of if model 
        
      ###############################################################################################################
        if (RSNNS == "NO") {}else{
          
          for (i in 1:n.rep)    
          {
            cat("> RSNNS MODEL run ",i,"predict  ...", "\n",append = FALSE)
            
            mod.tech <- "RSNNS"
            
            creat.subDir (paste(output,"/",proj[d],"_projection",sep=""),mod.tech)
            
            outpath <- paste(output,"/",proj[d],"_projection/",mod.tech,sep="")
            
            r <- stack( paste(output,"/",proj[d],"_projection/stack.rast.grd",sep="") )
            
            load(paste(output,"/",mod.tech,"/",mod.tech,"_model_",i,sep=''))
            
            pred <- predict(r, rsnns.tr) 
            
            # unscale the data
            pred.mean <- par.resc[1,1]
            pred.sd <- par.resc[2,1]
            
            pred <- (pred * pred.sd) + pred.mean 
            
            # back transformed
            f <- as.function(transf.y) # define function
            pred <- calc(pred, fun=f)
            
            writeRaster(pred, filename = paste(outpath,"/",mod.tech,"_model_run_",i,".tif",sep=""), 
                        datatype="FLT8S", overwrite=TRUE)
            
            file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER    
            
          }# end loop model replicate
        } # end of if model 
        
      ############################################################################################################### 
         if (ELM == "NO") {}else{
          
          for (i in 1:n.rep)    
          {
            cat("> ELM MODEL run ",i,"predict  ...", "\n",append = FALSE)
            
            mod.tech <- "ELM"
            
            creat.subDir (paste(output,"/",proj[d],"_projection",sep=""),mod.tech)
            
            outpath <- paste(output,"/",proj[d],"_projection/",mod.tech,sep="")
            
            r <- stack( paste(output,"/",proj[d],"_projection/stack.rast.grd",sep="") )
            
            load(paste(output,"/",mod.tech,"/",mod.tech,"_model_",i,sep=''))
            
            pred <- predict(r, elm.tr, type="raw") 
            
            # unscale the data
            pred.mean <- par.resc[1,1]
            pred.sd <- par.resc[2,1]
            
            pred <- (pred * pred.sd) + pred.mean 
            
            # back transformed
            f <- as.function(transf.y) # define function
            pred <- calc(pred, fun=f)
            
            writeRaster(pred, filename = paste(outpath,"/",mod.tech,"_model_run_",i,".tif",sep=""), 
                        datatype="FLT8S", overwrite=TRUE)
            
            file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER    
            
          }# end loop model replicate
        } # end of if model  
      ############################################################################################################### 
        if (RF == "NO") {}else{
          
          for (i in 1:n.rep)    
          {
            cat("> RF MODEL run ",i,"predict  ...", "\n",append = FALSE)
            
            mod.tech <- "RF"
            
            creat.subDir (paste(output,"/",proj[d],"_projection",sep=""),mod.tech)
            
            outpath <- paste(output,"/",proj[d],"_projection/",mod.tech,sep="")
            
            r <- stack( paste(output,"/",proj[d],"_projection/stack.rast.grd",sep="") )
            
            load(paste(output,"/",mod.tech,"/",mod.tech,"_model_",i,sep=''))
            
            pred <- predict(r, rf ) 
            
            # unscale the data
            pred.mean <- par.resc[1,1]
            pred.sd <- par.resc[2,1]
            
            pred <- (pred * pred.sd) + pred.mean 
            
            # back transformed
            f <- as.function(transf.y) # define function
            pred <- calc(pred, fun=f)
            
            writeRaster(pred, filename = paste(outpath,"/",mod.tech,"_model_run_",i,".tif",sep=""), 
                        datatype="FLT8S", overwrite=TRUE)   
            
            file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER    
            
          }# end loop model replicate
        } # end of if model  
        
      ############################################################################################################### 
        if (CSRF == "NO") {}else{
          
          for (i in 1:n.rep)    ## actually don't work for large raster on perso computer
          {
            cat("> CSRF MODEL run ",i,"predict  ...", "\n",append = FALSE)
            
            mod.tech <- "CSRF"
            
            creat.subDir (paste(output,"/",proj[d],"_projection",sep=""),mod.tech)
            
            outpath <- paste(output,"/",proj[d],"_projection/",mod.tech,sep="")
            
            r <- stack( paste(output,"/",proj[d],"_projection/stack.rast.grd",sep="") )
            
            load(paste(output,"/",mod.tech,"/",mod.tech,"_model_",i,sep=''))
            
		pred <- predict.mat.csrf(r,model=csrf, cut=5000) 

		pred <- check.res(pred,res=res(r),ext=extent(r))
            
            # unscale the data
            pred.mean <- par.resc[1,1]
            pred.sd <- par.resc[2,1]
            
            pred <- (pred * pred.sd) + pred.mean 
            
            # back transformed
            f <- as.function(transf.y) # define function
            pred <- calc(pred, fun=f)
            
            writeRaster(pred, filename = paste(outpath,"/",mod.tech,"_model_run_",i,".tif",sep=""), 
                        datatype="FLT8S", overwrite=TRUE)
            
            file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER    
                
            
          }# end loop model replicate
        } # end of if model  
        
      ############################################################################################################### 
        if (BAG == "NO") {}else{
          
          for (i in 1:n.rep)    
          {
            cat("> BAG MODEL run ",i,"predict  ...", "\n",append = FALSE)
            
            mod.tech <- "BAG"
            
            creat.subDir (paste(output,"/",proj[d],"_projection",sep=""),mod.tech)
            
            outpath <- paste(output,"/",proj[d],"_projection/",mod.tech,sep="")
            
            r <- stack( paste(output,"/",proj[d],"_projection/stack.rast.grd",sep="") )
            
            load(paste(output,"/",mod.tech,"/",mod.tech,"_model_",i,sep=''))
            
            pred <- predict(r, bag.tr ) 
            
            # unscale the data
            pred.mean <- par.resc[1,1]
            pred.sd <- par.resc[2,1]
            
            pred <- (pred * pred.sd) + pred.mean 
            
            # back transformed
            f <- as.function(transf.y) # define function
            pred <- calc(pred, fun=f)
            
            writeRaster(pred, filename = paste(outpath,"/",mod.tech,"_model_run_",i,".tif",sep=""), 
                        datatype="FLT8S", overwrite=TRUE)
            
            file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER    
            
          }# end loop model replicate
        } # end of if model 
        
      ############################################################################################################### 
        if (EGB == "NO") {}else{
          
          for (i in 1:n.rep)    # actually don't work for large raster on perso computer
          {
            cat("> EGB MODEL run ",i,"predict  ...", "\n",append = FALSE)
            
            mod.tech <- "EGB"
            
            creat.subDir (paste(output,"/",proj[d],"_projection",sep=""),mod.tech)
            
            outpath <- paste(output,"/",proj[d],"_projection/",mod.tech,sep="")
            
            r <- stack( paste(output,"/",proj[d],"_projection/stack.rast.grd",sep="") )
            
            # relaod model 
            file.name <- paste(output,"/",mod.tech,"/",mod.tech,"_model_",i,sep='')
            bst <- xgb.load(file.name)
            
            pred <- predict.mat.egb(r,model=bst, cut=5000) 

		pred <- check.res(pred,res=res(r),ext=extent(r))
            
            # unscale the data
            pred.mean <- par.resc[1,1]
            pred.sd <- par.resc[2,1]
            
            pred <- (pred * pred.sd) + pred.mean 
            
            # back transformed
            f <- as.function(transf.y) # define function
            pred <- calc(pred, fun=f)
            
            writeRaster(pred, filename = paste(outpath,"/",mod.tech,"_model_run_",i,".tif",sep=""), 
                        datatype="FLT8S", overwrite=TRUE)
            
            file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER    
            
          }# end loop model replicate
        } # end of if model 
        
      ###############################################################################################################  
       
         if (GBM == "NO") {}else{
          
          for (i in 1:n.rep)    
          {
            cat("> GBM MODEL run ",i,"predict  ...", "\n",append = FALSE)
            
            mod.tech <- "GBM"
            
            creat.subDir (paste(output,"/",proj[d],"_projection",sep=""),mod.tech)
            
            outpath <- paste(output,"/",proj[d],"_projection/",mod.tech,sep="")
            
            r <- stack( paste(output,"/",proj[d],"_projection/stack.rast.grd",sep="") )
            
            load(paste(output,"/",mod.tech,"/",mod.tech,"_model_",i,sep=''))
            
            pred <- predict(r, gbm.tr , n.trees=nb.tree, type="response",na.rm = TRUE) 
            
            # unscale the data
            pred.mean <- par.resc[1,1]
            pred.sd <- par.resc[2,1]
            
            pred <- (pred * pred.sd) + pred.mean 
            
            # back transformed
            f <- as.function(transf.y) # define function
            pred <- calc(pred, fun=f)
            
            writeRaster(pred, filename = paste(outpath,"/",mod.tech,"_model_run_",i,".tif",sep=""), 
                        datatype="FLT8S", overwrite=TRUE)
            
            file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER    
            
          }# end loop model replicate
        } # end of if model 
        
        ############################################################################################################### 
        
        if (GLM == "NO") {}else{
          
          for (i in 1:n.rep)    
          {
            cat("> GLM MODEL run ",i,"predict  ...", "\n",append = FALSE)
            
            mod.tech <- "GLM"
            
            creat.subDir (paste(output,"/",proj[d],"_projection",sep=""),mod.tech)
            
            outpath <- paste(output,"/",proj[d],"_projection/",mod.tech,sep="")
            
            r <- stack( paste(output,"/",proj[d],"_projection/stack.rast.grd",sep="") )
            
            load(paste(output,"/",mod.tech,"/",mod.tech,"_model_",i,sep=''))
            
            pred <- predict(r, glm.tmp,type="response",na.rm = TRUE) 
            
            # unscale the data
            pred.mean <- par.resc[1,1]
            pred.sd <- par.resc[2,1]
            
            pred <- (pred * pred.sd) + pred.mean 
            
            # back transformed
            f <- as.function(transf.y) # define function
            pred <- calc(pred, fun=f)
            
            writeRaster(pred, filename = paste(outpath,"/",mod.tech,"_model_run_",i,".tif",sep=""), 
                        datatype="FLT8S", overwrite=TRUE)
            
            file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER    
            
          }# end loop model replicate
        } # end of if model 
        ############################################################################################################### 
        
        if (GLM_sw == "NO") {}else{
          
          for (i in 1:n.rep)    
          {
            cat("> GLM_sw MODEL run ",i,"predict  ...", "\n",append = FALSE)
            
            mod.tech <- "GLM_sw"
            
            creat.subDir (paste(output,"/",proj[d],"_projection",sep=""),mod.tech)
            
            outpath <- paste(output,"/",proj[d],"_projection/",mod.tech,sep="")
            
            r <- stack( paste(output,"/",proj[d],"_projection/stack.rast.grd",sep="") )
            
            load(paste(output,"/",mod.tech,"/",mod.tech,"_model_",i,sep=''))
            
            pred <- predict(r, glm.tmp.step, type="response",na.rm = TRUE) 
            
            # unscale the data
            pred.mean <- par.resc[1,1 ]
            pred.sd <- par.resc[2,1 ]
            
            pred <- (pred * pred.sd) + pred.mean 
            
            # back transformed
            f <- as.function(transf.y) # define function
            pred <- calc(pred, fun=f)
            
            writeRaster(pred, filename = paste(outpath,"/",mod.tech,"_model_run_",i,".tif",sep=""), 
                        datatype="FLT8S", overwrite=TRUE)
            
            file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER    
            
          }# end loop model replicate
        } # end of if model 
      ############################################################################################################### 
        
        if (GAM == "NO") {}else{
          
          for (i in 1:n.rep)    
          {
            cat("> GAM MODEL run ",i,"predict  ...", "\n",append = FALSE)
            
            mod.tech <- "GAM"
            
            creat.subDir (paste(output,"/",proj[d],"_projection",sep=""),mod.tech)
            
            outpath <- paste(output,"/",proj[d],"_projection/",mod.tech,sep="")
            
            r <- stack( paste(output,"/",proj[d],"_projection/stack.rast.grd",sep="") )
            
            load(paste(output,"/",mod.tech,"/",mod.tech,"_model_",i,sep=''))
            
            pred <- predict(r, gam.tmp, type="response",na.rm = TRUE ) 
            
            # unscale the data
            pred.mean <- par.resc[1,1 ]
            pred.sd <- par.resc[2,1 ]
            
            pred <- (pred * pred.sd) + pred.mean 
            
            # back transformed
            f <- as.function(transf.y) # define function
            pred <- calc(pred, fun=f)
            
            writeRaster(pred, filename = paste(outpath,"/",mod.tech,"_model_run_",i,".tif",sep=""), 
                        datatype="FLT8S", overwrite=TRUE)
            
            file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER  
            
          }# end loop model replicate
        } # end of if model 
        
      } # end of loop projection
      
} # end loop of element 
      #################################################################################################################################
      #################################################################################################################################

endCluster() # END OF MULTICORE CALCULATION

# proc.time() - ptm # check time
       
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#### COFFE TIMES
#################################################################################################################################
#################################################################################################################################




