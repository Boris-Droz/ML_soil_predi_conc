####################################################################################################################################
####################################################################################################################################
###					
### RELOAD PRODUCTED RASTER 
###  CREAT AVERAGE AND SD FOR EACH TECHNIQUE
### ==============================================
####################################################################################################################################
####################################################################################################################################
## v1.0 Mai 2018 - Boris DROZ - Project Copper
####################################################################################################################################
## DESCRIPTION
###############     #####################################
#                   --- CHECK PARAMETER UNTIL LINE 50 ---
##                   #####################################

##    Tested under R 3.1 and 3.4
################################
## Average and creat the sensmble for all run of predicted model
## Need to run the script Model_pred_v9.6 before
##
# Input: - output from the Model_pred_v9.6 script
#       
# Output: 
#=======
# - produced a average for each techique, the sd for each technique and the ENSEMBLE
#
####################################################################################################################################

## Set library path
#.libPaths("C:/Users/bdroz/Documents/R/win-library/3.4")

## library list
###############

# spatial analysis
library(raster)

####################################################################################################################################
## SCRIPT PARAMETER  --> NEED TO GO TROUH AND MAKE THE APPROPRIATE MODIFICATION !!!
#####################

## =========
inpath <-"D:/copper/input/pred"
output <- "D:/copper/output/20180330_model_v18_1000"

outpath <- output ## DON'T MODIFIED ##

## set your work space ## DON' MODIFIE ##
setwd(inpath) ## DON'T MODIFIED ##

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

###################################################################################################################################
####################################################################################################################################
## SCRIPT START HERE
####################################################################################################################################

# folder list of diff projection
proj <- list.dirs(path = inpath, full.names = FALSE, recursive = FALSE)
fn <- list.dirs(path = inpath, full.names = TRUE, recursive = FALSE)

# folder list of element
elem <- list.dirs(path = output, full.names = FALSE, recursive = FALSE)
fn.out <- list.dirs(path = output, full.names = TRUE, recursive = FALSE)

# 1) START LOOP FOR DIFFERENT ELEMENT
########################################

for (j in 1:length(fn.out))
{ 
  cat("> ",elem[j]," projection is running now...", "\n",append = FALSE)
  cat("##############################################", "\n",append = FALSE)
  
  output <- fn.out[j]
  
  ### 2) START LOOP FOR DIFFERENT PROJECTION (SCENARIO)
  #####################################################
  for (d in 1:length(proj))
  { 

    cat("> COMPUTE AVERAGE for projection ",proj[d],"...", "\n",append = FALSE)
    
    creat.subDir (paste(output,"/",proj[d],"_projection",sep=""),"average")
    
    # folder list --> list technic
    folder <- list.dirs(path = paste(output,"/",proj[d],"_projection",sep=""), full.names = TRUE, recursive = FALSE)
    
    mod.tech <- list.dirs(path = paste(output,"/",proj[d],"_projection",sep=""), full.names = FALSE, recursive = FALSE)
    
    # keep all path of mean and sd
    f.mean <- NULL
    f.sd <- NULL
    
    for (q in 2:length(folder)) # don't account average folder (it is in first)
    {
      cat( "--> ",mod.tech[q], "averaging ....", "\n",append = FALSE)
      
      fns <- list.files(folder[q],pattern=".tif$",full.names = TRUE) # list file for each technic
      
      ### A) COMPUTE AVERAGING FOR ONE TECHNIQUE 
      ##########################################
      
      f.mean <- c(f.mean, paste(folder[1],"/",mod.tech[q],"_average.tif",sep=""))
      f.sd <- c(f.sd, paste(folder[1],"/",mod.tech[q],"_sd.tif",sep=""))
      
      for (n in 1:length(fns)) # first summ all prediction
      {
        cat( "--> SUMM LOOP ",n, "....", "\n",append = FALSE)
      
          if (n == 1) { output.raster <- stack(fns[n])
          
          writeRaster(output.raster, filename = paste(folder[1],"/",mod.tech[q],"_average.tif",sep=""),
                      datatype="FLT8S", overwrite=TRUE) 
          } else {
            
            #reload raster, merge with the new on, rewrite and remove temp...
            output.raster <- stack( paste(folder[1],"/",mod.tech[q],"_average.tif",sep="") )
            output.raster <- calc( stack(fns[n],output.raster) , fun=sum )
            
            writeRaster(output.raster, filename = paste(folder[1],"/",mod.tech[q],"_average.tif",sep=""),datatype="FLT8S", overwrite=TRUE)
            
            file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER
          } 
      } # end sum raster
      
      # reload summ and average it
      output.raster <- raster( paste(folder[1],"/",mod.tech[q],"_average.tif",sep="") )/length(fns)
      
      writeRaster(output.raster, filename = paste(folder[1],"/",mod.tech[q],"_average.tif",sep=""),datatype="FLT8S", overwrite=TRUE)
      
      ### B) COMPUTE SD FOR ONE TECHNIQUE 
      ##########################################
      
      for (n in 1:length(fns)) # first creat sum of diff with mean
      {
        cat( "--> SD LOOP ",n, "....", "\n",append = FALSE)
        
        if (n == 1) { 
          mean <- raster( paste(folder[1],"/",mod.tech[q],"_average.tif",sep="") )
          
          output.raster <- (raster(fns[n])- mean)^2
        
          writeRaster(output.raster, filename = paste(folder[1],"/",mod.tech[q],"_sd.tif",sep=""),
                      datatype="FLT8S", overwrite=TRUE) 
          
          file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER
        } else {
          mean <- raster( paste(folder[1],"/",mod.tech[q],"_average.tif",sep="") )
          
          #reload raster, merge with the new on, rewrite and remove temp...
          output.raster <- stack( paste(folder[1],"/",mod.tech[q],"_sd.tif",sep="") )
          output.raster <- calc(stack( output.raster,(raster(fns[n])- mean)^2 ), fun= sum)
          
          writeRaster(output.raster, filename = paste(folder[1],"/",mod.tech[q],"_sd.tif",sep=""),datatype="FLT8S", overwrite=TRUE)
          
          file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER
          
        } 
      } # end sum raster
      
      output.raster <- sqrt( raster( paste(folder[1],"/",mod.tech[q],"_sd.tif",sep="") )/(length(fns)-1) )
      
      writeRaster(output.raster, filename = paste(folder[1],"/",mod.tech[q],"_sd.tif",sep=""),datatype="FLT8S", overwrite=TRUE)
      
      file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER    
      
    }# end loop model technic
    
    cat( "--> ENSEMBLE averaging ....", "\n",append = FALSE)
    
    # compute the average for the ensemble
    #########################################
    ens.mean <- calc(stack(f.mean), fun=sum)/(length(folder))
    
    writeRaster(ens.mean, filename = paste(folder[1],"/ENSEMBLE_average.tif",sep=""), 
                datatype="FLT8S", overwrite=TRUE)
    
    file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER    
    
    # compute the sd of the ensemble model 
    ######################################
    # creat a list of file for into all folder
    fns <- NULL
    for (n  in 2:length (folder))   
    {fns <- c(fns,list.files(folder[n],pattern=".tif$",full.names = TRUE))}
    # print(fns) # all file 
    
    for (n in 1:length(fns)) # first creat sum of diff with mean
    {
      cat( "--> SD ENSEMBLE LOOP ",n, "....", "\n",append = FALSE)
      
      if (n == 1) { 
        # reload ensemble raster
        ens.mean <- raster( paste(folder[1],"/ENSEMBLE_average.tif",sep="") )
        
        output.raster <- (raster(fns[n])- ens.mean)^2
      
        writeRaster(output.raster, filename = paste(folder[1],"/ENSEMBLE_sd.tif",sep=""),
                    datatype="FLT8S", overwrite=TRUE) 
      } else {
        
        # reload ensemble raster
        ens.mean <- raster(paste(folder[1],"/ENSEMBLE_average.tif",sep="")) 
        
        #reload raster, merge with the new on, rewrite and remove temp...
        output.raster <- stack( paste(folder[1],"/ENSEMBLE_sd.tif",sep="") )
        output.raster <- calc(stack( output.raster,(raster(fns[n])- ens.mean)^2 ), fun= sum)
        
        writeRaster(output.raster, filename = paste(folder[1],"/ENSEMBLE_sd.tif",sep=""),datatype="FLT8S", overwrite=TRUE)
        
        file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER
        
      } 
    } # end sum raster
    
    output.raster <- sqrt( raster( paste(folder[1],"/ENSEMBLE_sd.tif",sep="") )/(length(fns)-1) )
    
    writeRaster(output.raster, filename = paste(folder[1],"/ENSEMBLE_sd.tif",sep=""),datatype="FLT8S", overwrite=TRUE)
    
    file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER
    
    
  } # end loop of the projection  
} # end loop of element    