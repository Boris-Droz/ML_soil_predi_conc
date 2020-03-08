####################################################################################################################################
####################################################################################################################################
###					
###  VARIABLE SELECTION
### ===========================
###					
####################################################################################################################################
####################################################################################################################################
## HISTORIC:
## --------
##          v1.0, Decembre 2015 - Boris DROZ & Gerrad Jones, ETHZ & EAWAG
##                Projet Se around the world
##          v2-4 modified December 2016 - B.Droz copper project - unistra
##		      v5.7 July 2017 - - choose analysis and divers para.
## 
####################################################################################################################################
## DESCRIPTION
###############      ######################################
#                   --- CHECK PARAMETER UNTIL LINE 190 ---
##                   #####################################

##    Tested under R 3.1
#########################
# 
# Transform all continuous var. --> Normalized !!!
#
# Performed divers analysis in order to select the best var. for further model
##  you have possibilities to made or not all methods
#           i)    comparaison between data on raster and extracted data
#		        ii) Moran analysis --> spatial interpolation
#           iii)   correlation plot between variables
#           iv)  cluster analysis
#           v)   principal components analysis (PCA) 
#           vi)    Neural network & cancelation procedure
#           vii)   GLM - univarite
#           viii)  GLM - Lasso procedure
#           ix) Random Forest - mean decrese accurancy

# Input: 
#=======  i) data : var. to predict and cooordinate for each data point 
#                  *txt tab file with 3 header"X_WGS84","Y_WGS84" and "Var_to_pred"
#         ii) path.pred : # folder with all *.geotiff "predictive variables"
##                        should have same extent, coordinate and resolution.

# Output: 
#=======
# Boxplot of all data between resample and remain
# Correlation predict between all variable
# PCA analyis for the data set
# NNET:
#   - Model performance: i) plot OBS vs PRED --> slope, Intercept, R2 
#                            and Root of mean square error between the outputs and targets (RMSE)
#                         ii) MOdel --> AIC, AICc and ratio between nb iteration and nb iteration max (ratio_CONV)
# creat a average folder
####################################################################################################################################
## Load library
##################

## Set library path if necessary
#.libPaths("C:/Users/Public/Documents/R/win-library/3.1")

# stat analysis
library(ade4)
library(boot)
library(MASS)
library(pracma)
library(qpcR)
library(pvclust)
library(nFactors)
library(ape)

# data read store or plot
library(foreign)
library(ggplot2)
library(gplots)
library(plyr)
library(Hmisc)
library(corrplot)
library(stringr)

# spatial analysis
library(raster)
library(gdalUtils)

## Regression Modeling library
library(rms)  
library(gam)
library(nnet)
library(neuralnet) 
library(NeuralNetTools)# importance and sensitivity analysis of nnet
library(glmnet)
library(randomForest)
library(AICcmodavg)

## library parrallel core
library(snowfall)
library(parallel)

####################################################################################################################################
## SCRIPT PARAMETER
#####################
####################################################################################################################################
## SCRIPT PARAMETER  
#####################

# Set folder 
## =========
inpath <- "D:/copper/input/pred" # input
outpath <-"D:/copper/output/20180307_var_selection" # output

path.pred <- paste(inpath,"/varinput",sep="") # folder with all *.geotiff "predictive variables"

## set your work space ## DON'T MODIFIE ##
setwd(inpath)

# *.txt File name of data points
data <- "20170613_Data_Cu_Vine_EU"

# MASK OF THE PROJECTION --> only in case of i) COMPARE RASTER DATA AND THE EXTRACT VALUES
##########################
mask.proj <- "YES" ## YES or NO

mask.r <- "maskvine_wgs84_250m.tif"

## Choose which analysis to performed
#######################################
COMP <- "YES" # comprison between rasster and extract data
MI <- "NO" # Moran index spatial analysis
CORR <- "NO" # correlation plot
CLUST <- "NO" # cluster analysis
PCA <- "NO" # principal compound analysis
NNET <- "YES" # Neural network & cancelation procedure
GLM <- "NO" # glm univariate
LASSO <- "NO" # lasso glm procedure 
RF <- "NO" # random forest test

## TRANSFORMATION
## ==============  
## list of transformation to satisfied normality of the data
# Transformation list same number as predictive variable 
## # alist(x=,log10(abs(x))) ,alist(x=,sqrt(abs(x)))) # EXEMPLE OF TRANSF.
#==============================================================

transf.y <- alist(x=,log10(abs(x))) # dep. variable

## selected predictive variable - exemple
#transf <- list(  alist(x=,log10(abs(x))), alist(x=,(x)), alist(x=,(x)), 
 #               alist(x=,(1/x)),alist(x=,(x)),  alist(x=,(1/x)), 
  #              alist(x=,log10(abs(x))), alist(x=,(x)), alist(x=,(x)) )

transf <- list( alist(x=,log10(abs(x))), alist(x=,(x)), alist(x=,(x)),alist(x=,(x)),alist(x=,(x)),
		alist(x=,(x)),alist(x=,(x)),alist(x=,(x)),alist(x=,sqrt(abs(x))),alist(x=,(1/x)),
		 alist(x=,log10(abs(x))), alist(x=,sqrt(abs(x))),alist(x=,sqrt(abs(x))),
		alist(x=,log10(abs(x))), alist(x=,log10(abs(x))), alist(x=,(x)), alist(x=,log10(abs(x))), alist(x=,log10(abs(x))),
		alist(x=,(x)), alist(x=,(x)), alist(x=,log10(abs(x))), alist(x=,(1/x)), alist(x=,log10(abs(x))), 
		alist(x=,log10(abs(x))), alist(x=,(x)),alist(x=,(x)), alist(x=,(x)), alist(x=,(x)), alist(x=,log10(abs(x))), alist(x=,log10(abs(x))) )

length(transf)

## NUMBER OF MODEL BUILDING
# ========================
n.rep <- 1
## RESAMPLING WITHIN A EQUAL SCALE INTO THE GRID CELL
#====================================================
resampling <- "NO" ## YES or NO

## reSampling paramater if resampling
# Define the degree grid cell resolution (Proj. WGS_1984)
siz.cell <- 2.5

## FUNCTION PARAMETER
## ==================
## glm-gam family type of distribution
fam <- "gaussian" # not actif!

## NNET PARAMETER
##################
## Weight decay Folow recomendation of B. D. Ripley: "Pattern Recognition and Neural Networks", Cambridge, 1996.
UNIT.Max <- 10 # UMBER OF HIDDEN UNITS --> between the number of input nodes and number of output nodes 
            	
it.max <- 1000 # Need to be 10'000 to be MONTE-CARLO Permutation --> on final model

## RANDOM FOREST
################
nb.tree <- 100

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
##########################
##########################
##########################
##                        ##
##    ADD-IN FUNCTIONS    ##
##                        ##
##########################
########################################################################################################################################
#########################################################
########################################################
## Neuronal Network input cancellation FUNCTION
## ============================================
##
## v2.0 - June 2015 -  Boris Droz &G errad Jones ETHZ & EAWAG
## 
## Procedure recommaneded by Héctor F. Satizábal M. (2007)
###########################################################
## Parameter
## =========
## n : number of hidden nods of the nnet
## w : number of weight decay

nnet.cancel <- function (x.train,y.train,n,w,it.max= 10000)
{
  
  require(nnet);
  require(MASS);
  require(qpcR);
  
  # INITIATE A NULL TABLE
  rmse.table <- NULL;
  
  x.tr <- x.train # temporary x data set for training nnet
  
  for (n in 1:(ncol(x.train)-1))
  {
    
    # TRAIN NEURAL NETS
    data.in <- cbind(y.train,x.tr)
    
    net <- nnet(x.tr,y.train, data=data.in, size = n, 
                linout = TRUE, maxit = it.max, decay = w,
                 trace = FALSE)
    
    # PREDICT THE OBS
    pred.y <- predict(net, newdata=x.tr,type="raw",na.rm = TRUE)
    
    # PERFORMANCE OF THE MODEL
    lin.corr <- lm(y.train[,1] ~ pred.y) #Matrix ."slope","Int","R2","RMSE" y.in[,1]
    
    if (nrow(summary(lin.corr)$coefficients)==1) # Fuck up check
    { 
      # TRAIN NEURAL NETS
      net <- nnet(x.tr,y.train, data=data.in, size = n, 
                  linout = TRUE, maxit = it.max, decay = w,
                  trace = FALSE)
      
      # PREDICT THE OBS
      pred.y <- predict(net, newdata=x.tr,type="raw",na.rm = TRUE)
      
      # PERFORMANCE OF THE MODEL
      lin.corr <- lm(y.train ~ pred.y)  
      
    }else{ }
    
    slope <- summary(lin.corr)$coefficients[2,1] #slope
    std.slope <- summary(lin.corr)$coefficients[2,2]
    int. <- summary(lin.corr)$coefficients[1,1] #Intercept
    std.int <- summary(lin.corr)$coefficients[1,2]
    rsquare <- summary(lin.corr)$r.squared #r squared
    RMSE <- sqrt(mean((y.train - pred.y)^2,na.rm = TRUE)) #RMSE
    
    perf.model <- c(slope=slope, std.slope=std.slope, int.=int., std.int=std.int, rsquare=rsquare, RMSE=RMSE)
    
    train.rmse <- perf.model[6]
    
    # NAME OF VARIABLE USED
    var.sel <- paste((dimnames(x.tr)[[2]]),"/",sep="",collapse="")
    
    # KEEP RESULT TO SEE EVOLUTION OF THE PERFOMANCE
    if (n == 1) { rmse.table <- c(Var.sel = var.sel, perf.model) }else{  
      rmse.table <- rbind( rmse.table, c(Var.sel = var.sel, perf.model))}
    
    diff.table <- NULL
    
    for (m in 1:ncol(x.tr))
    {
      # REPLACE ONE VARIABLE BY THAN AVERAGE
      x.temp <- x.tr
      x.temp [,m] <- mean(x.tr[,m])
      
      # CALCULATE RMSE FOR NEW SET 
      new.rmse <- sqrt(mean( (y.train - predict(net, newdata=x.temp,type="raw" ))^2 ))
      
      diff <- train.rmse - new.rmse
      
      # APPEND EACH diff TO A VECTOR
      if (m == 1)  diff.table <- diff else  diff.table <- rbind( diff.table, diff);
      
    }
    rownames(diff.table) <- colnames(x.tr) # RE-NAME NEW TABLE
    
    ## FIND THE LESS RELEVANT INPUT
    l.r <- order(diff.table[,1],decreasing = TRUE)
    
    ## DISCARD THE LESS RELEVANT VARIABLE
    x.tr <- x.tr[,-l.r[1]]
  }
  
  rownames(rmse.table) <- rev(seq (from=1, to= (ncol(x.train)-1), by=1))   
  rmse.table <- as.data.frame(rmse.table)
  
  rmse.table <- rmse.table[order(rmse.table$RMSE,decreasing = FALSE),]
  
  return(rmse.table)
}

########################################################################################################################################
#################
## Adjusted D2 ##
#################

adj.D2.glm <- function(glmobj)
        {
          go <- glmobj
          D2 <- (go$null.deviance - go$deviance)/go$null.deviance
          p <- length(go$coefficients)
          n <- length(go$fitted)
          adj.D2 <- 1 - ((n - 1)/(n - p)) * (1 - D2)
          if (adj.D2 < 0) 
              {
                adj.D2 <- 0
                return(adj.D2)
              } else 
              {
                return(adj.D2)
              }
        }

#################################################################################################################################
#################
## Rescaled D2 ##
#################

R2.glm.rsadjnag <- function(glm.obj) 
        {
          n.s <- length(glm.obj$y)
          p <- length(glm.obj$coefficients)
          R2.rsnag <- (1-(exp(-glm.obj$null.deviance/2)/exp(-glm.obj$deviance/2))^(2/n.s))/
            (1-exp(-glm.obj$null.deviance/2)^(2/n.s))        	
          1 - ((n.s-1)/(n.s-p))*(1-R2.rsnag)      	
        }

########################################################################################################################################################
###################################################################################
## VARIABALE SLECTION BASED ON RANDOM TREE ANALYSIS 
### =================================================
##
## Procedure from Carpita et al. 2014 "Football Mining with R" 
## in Data Mining Applications with R
##
## Mai 2015 - Boris DROZ & Gerrad Jones, ETHZ & EAWAG
###################################################################################

RF.var.sel <- function (x.train,y.train,nb.tree=1000)
                {
                  
                  require(randomForest);
                  # Generate a matrix Zf of pseudo-covariates by randomly permuting
                  Zf <- x.train[sample(nrow(x.train)),] 
                  # Covariate and pseudo-covariate matrices
                  dtset.pseudo <- data.frame(cbind(x.train, Zf, y.train)) 
                  # rename var. dependant!!
                  dimnames(dtset.pseudo)[[2]][length(dimnames(dtset.pseudo)[[2]])] <- "Y"
                  # Generated tree analysis 
                  rf <- randomForest(Y ~ ., data=dtset.pseudo, ntree = nb.tree, importance=TRUE) 
                  # Reporte the total decrease in node impurities (TDNI) and increase error                  
                  VIMs <- rf$importance
                  
                  return(VIMs)           
                }     
#########################################################
#########################################################
#################################################
## TUNE NNET FUNCTION
## ===================
## 
## Find out the optimal parameter (hidden nods and weight decay)
##  for the nnet function
##
## Folow recomendation of B. D. Ripley: "Pattern Recognition and Neural Networks", Cambridge, 1996.
##
## max.nods need to be higher than the number of variable, however no rule of tumbs are develloped for the real values.
## Input are data. table containing Y and X.

## Mai 2015 - Boris DROZ & Gerrad Jones, ETHZ & EAWAG
## mod. 17.6.2015

###############################################################

nnet.para.opt <- function (x.train,y.train,max.nods=30,it.max= 10000)
  {

  require(nnet);
  require(MASS);
  
  # INITIATE A NULL TABLE
  sse.table <- NULL;
  
  # SEARCH FOR OPTIMAL WEIGHT DECAY
  # RANGE OF WEIGHT DECAYS SUGGESTED BY B. RIPLEY
  for (w in c(0.0001, 0.001, 0.01))
  {
    # SEARCH FOR OPTIMAL NUMBER OF HIDDEN UNITS
    for (n in 1:max.nods)
    {
      # UNITIATE A NULL VECTOR
      sse <- NULL;
      # FOR EACH SETTING, RUN NEURAL NET MULTIPLE TIMES
      for (i in 1:10)
      {
        # INITIATE THE RANDOM STATE FOR EACH NET
        set.seed(i);
        # TRAIN NEURAL NETS
        net <- nnet(x.train,y.train, data=cbind(y.train,x.train), size = n, 
                    linout = TRUE, maxit = it.max, decay = w, trace = FALSE); #print (net)
        
        # CALCULATE SSE FOR TESTING SET
        #test.sse <- sum((test.set$Y - predict(net, test.set))^2);
        
        # CALCULATE RMSE FOR TESTING SET 
        test.sse <- sqrt(mean((y.train - predict(net, newdata=x.train,type="raw" ))^2))
        
        # APPEND EACH SSE TO A VECTOR
        if (i == 1) sse <- test.sse else sse <- rbind(sse, test.sse);
      }
      # APPEND AVERAGED SSE WITH RELATED PARAMETERS TO A TABLE
      sse.table <- rbind(sse.table, c(WT = w, UNIT = n, RMSE = mean(sse), SD_RMSE =sd(sse)));
    }
  }     
  sse.table <- as.data.frame(sse.table)
  sse.table <- sse.table[order(sse.table$RMSE),]
  
  return(sse.table)
}

###############################################################################################################
##############################
## CV-NNET ##
##############################
#
## Droz. B - 12.6.2015
#
## Adapted function from Ecospat package for continuous data
## all data with same weights are conbsidered
#############################################################

CV.NNET. <- function(NNET.obj, data.cv, K=10, cv.lim = 10, name.sp)
{
  
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 )
  
  # CONTROL IF THE NUMBER OF OBSERVATIONS IS SUFFICIENT
  
  if ((K > n) || (K <= 1))
    stop("K outside allowable range")
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  K <- K
  K.lst <- round(K)
  kvals <- unique(round(n/(1:floor(n/2))))
  temp <- abs(kvals - K.lst)
  if (!any(temp == 0))
  {
    K.lst <- kvals[temp == min(temp)][1]
  }
  if (K.lst != K)
  {
    cat("K has been set to", K.lst, "\n",append = F)
  }
  
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  K.lst <- K.lst
  
  f <- ceiling(n/K.lst)
  s <- sample(rep(1:K.lst,f), n)
  
  # RESPONSE PREDICTION
  
  for (i in 1:K.lst)
  {
    j.out <- id[(s == i)]
    j.out <- sort(j.out)
    j.in <- id[(s != i)]
    j.in <- sort(j.in)
    
    data.cal.cv <- data.cv[j.in,]
    data.test.cv<- data.cv[j.out,]
    
    nnet.cal <- update(NNET.obj, data = data.cal.cv)
    
    nnet.val <- predict(NNET.obj, data.test.cv , type = "raw")
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(nnet.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(nnet.val), after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0) ,predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}
#################################################################################################################################################
# FUNCTION check and produced subDir folder
###########################################
#
creat.subDir <- function (mainDir,subDir)
{
  if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
  } else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
  }
}
#################################################################################################################################################

####################################################################################################################################
####################################################################################################################################
## SCRIPT STAR HERE
####################################################################################################################################
##############################################
# -- LOADING DATA SET -----
##############################################
setwd(inpath)

## pred. var file (batch all file in input) 
fns <- list.files(path.pred,pattern=".tif$",full.names = TRUE); print(fns)

# predicitve variable name list
var.list <- strsplit(gsub(".tif", "", list.files(path.pred,pattern=".tif$",full.names = FALSE)),split="_"); print(var.list)

# Open the data set and the unknow point around the world
data.in <- read.table(paste(data,".txt",sep=""),header = TRUE ,na.strings = "NaN"); head(data.in)

##################################
## --- FANCY RESAMPLING OR NOT---  ##
##################################

## equal.scale - average degree cell
##--------------------------------

# pos dep. var after x and y coord
dep.var <- max( match("Y_WGS84",names(data.in)), 
                match("X_WGS84",names(data.in)) ) + 1

# position of the last original column
ncol.data0 <- ncol(data.in)

if (resampling=="YES") 
{
  cat("> AVERAGE BY CELL GRID ...", "\n",append = FALSE)
  
  data.in$X_WGS84 <- round( data.in$X_WGS84/siz.cell )*siz.cell
  
  data.in$Y_WGS84 <-  round( data.in$Y_WGS84/siz.cell )*siz.cell
  
  cell.long <- unique(data.in$X_WGS84)
  cell.lat <- unique(data.in$Y_WGS84)
  
  d.samp <- matrix(data = NA, nrow = 0, ncol= ncol(data.in), dimnames = list(c(),dimnames(data.in)[[2]]))
  
  nb <- 0 # counter
  
  ## mean of all data in the same cell
  for (s in 1:length(cell.long))  
  {
    for (k in 1:length(cell.lat))
    {
      if (table( data.in$X_WGS84==cell.long[s] & 
                   data.in$Y_WGS84==cell.lat[k])["FALSE"]<=nrow(data.in)-1){
        
        nb <- nb + 1
        
        #if (min(x_coor,y_coor,dep.var)==1) { nb <- NULL} 
        nb <- rep(nb, (min( match("Y_WGS84",names(data.in)), 
                            match("X_WGS84",names(data.in)))-1) ) 
      
        d.samp <- rbind(d.samp,
                      c(nb,cell.long[s],cell.lat[k],
                        mean(data.in[data.in$X_WGS84==cell.long[s] & data.in$Y_WGS84==cell.lat[k],
                                     (dep.var+1):ncol(data.in)])))
      
    }else{ d.samp <- d.samp} 
  }
}
d.samp <- as.data.frame (d.samp)

# Control same projection, spatial extent and resolution.
# coord.ref similar ....
ext <- cbind(d.samp$X_WGS84,d.samp$Y_WGS84)

#colnames(r.proj) <- c("Longitude","Latitude",dep.var)

for (z in 1:length(fns))
{
  cat("> Extract Pred.", var.list[[z]][1]," ...", "\n",append = FALSE)
  
  pred.var <- raster(fns[z])
  
  pos <- ncol(d.samp)
  
  # AGREGATE THE DATA OR NOT IN FUNCTION OF THE RESOLUTION                    
  if (siz.cell/res(pred.var)[1] > 1) {
    
    pred.var <- aggregate(pred.var,fact=siz.cell/res(pred.var)[1],fun=mean)
    
    ## APPLIED DATA MODIFICATION
    f <- as.function(transf[[z]]) # define function
    pred.var <- calc(pred.var, fun=f) # applied modification
    
    ## RESCALE THE PREDICTOR with normal rescaling
    pred.mean <- cellStats(pred.var, stat='mean')
    pred.sd <- cellStats(pred.var, stat='sd')
    
    pred.var <- (pred.var-pred.mean)/pred.sd # function rescaling
    
  }else{ }  
  
  d.samp <- cbind(d.samp,extract(pred.var,ext))  
  
  #UNIFIED NAME OF THE RASTER
  colnames(d.samp)[pos+1] <- var.list[[z]][1]
  
}

} else {
  d.samp <- data.in 
  
  # coord.ref similar ....
  ext <- cbind(d.samp["X_WGS84"],d.samp["Y_WGS84"])
  
  for (z in 1:length(fns))
  {
    cat("> Extract Pred.", var.list[[z]][1]," ...", "\n",append = FALSE)
    
    pred.var <- raster(fns[z])
    
    ## APPLIED DATA MODIFICATION
    f <- as.function(transf[[z]]) # define function
    pred.var <- calc(pred.var, fun=f) # applied modification
    
    ## RESCALE THE PREDICTOR with normal rescaling
    pred.mean <- cellStats(pred.var, stat='mean')
    pred.sd <- cellStats(pred.var, stat='sd')
    
    pred.var <- (pred.var-pred.mean)/pred.sd # function rescaling
    
    pos <- ncol(d.samp)
    
    d.samp <- cbind(d.samp,extract(pred.var,ext))  
    
    #UNIFIED NAME OF THE RASTER
    colnames(d.samp)[pos+1] <- var.list[[z]][1]
  }  
  
}

file.remove(dir(tmpDir(),full.name=TRUE)) #REMOVE TEMPORARY RASTER

write.table(d.samp, file=paste(outpath,"/data_sampled.txt",sep=""),sep="\t",
            append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)

######################################################################################################################################
## extract variable check --
#############################

for (z in (ncol.data0+1):ncol(d.samp))
  {
	print(names(d.samp)[z])
	print(length(d.samp[,z]))
	print(length(na.omit(d.samp[,z])))
}

#####################################################################################################################################
## ORGINIZED DATA SET FOR VAR: SELECTION
##======================================

#############################################################
data.in <- cbind(d.samp$X_WGS84,d.samp$Y_WGS84, d.samp[,dep.var],d.samp[,(ncol.data0+1):ncol(d.samp)])

names(data.in)[3] <- names(d.samp)[dep.var] # reload the name of dep.var

# modified the dep. var ##
##########################
### APPLIED DATA MODIFICATION for the dep.var
f <- as.function(transf.y) # define function 
data.in[,3] <- f(data.in[,3]) # applied modification

data.in <- na.omit(data.in) # delet row with NA value

# separe data and coordonnate
coor.xy <- data.in [,1:2]
data.in <- data.in [,3:ncol(data.in)]

## RESCALE THE PREDICTOR with normal rescaling
pred.mean <- mean(data.in[,1])
pred.sd <- sd(data.in[,1])

data.in[,1] <- (data.in[,1]-pred.mean)/pred.sd # function rescaling

# rewrite names of header
names(data.in) <- c(names(d.samp)[dep.var],names(d.samp)[(ncol.data0+1):ncol(d.samp)])
names(coor.xy) <- c("X_WGS84","Y_WGS84")

# creat x and y table
y.in <- data.in[,1,drop = FALSE]
x.in <- data.in[,2:ncol(data.in)]
####################################################################
################################################################################################################################
## VARIABLE SELECTION
##=====================
############################################################  
## IGNITITA LIST AND VECTOR --> used if loop greater than 1
########################################################### 

## NNET model performance
modl.perf <- matrix(0,n.rep,9,dimnames=list(c(1:n.rep),c("CV_slope","CV_b","CV_R2","CV_ME",
                                                         "CV_MAE","CV_RMSE",
                                                         "AIC","AICc","CONV")))

modl.perf.l <- list(modl.perf,modl.perf,modl.perf,modl.perf)

## weight of nnet by garson method
garson.plot.l <- list(list(),list(),list(),list())

# Prepare output matrix - variable
mat.var.pca <- matrix(0,n.rep,(length(var.list)-2),dimnames=list(c(1:n.rep),
              paste("Var.Axe_",seq(from=1,to=(length(var.list)-2),by=1),sep="")))

mat.var.pca.l <- list (mat.var.pca, mat.var.pca ,mat.var.pca, mat.var.pca)

# PCA LIST
# Prepare output matrix - all 
mat.all.pca <- matrix(0,n.rep,(length(var.list)-2),dimnames=list(c(1:n.rep),
              paste("All.Axe_",seq(from=1,to=(length(var.list)-2),by=1),sep="")))

mat.all.pca.l <- list (mat.all.pca,mat.all.pca,mat.all.pca,mat.all.pca)

# CORPLOT LIST
df.cor.pred.var.l <- list(list(),list(),list(),list())

df.cor.pred.all.l <- list(list(),list(),list(),list())

## initialise the list data --- GLM performance / Optimum

mod.perf.glm <- list(list(NULL,NULL,NULL,NULL),list(NULL,NULL,NULL,NULL),
                     list(NULL,NULL,NULL,NULL),list(NULL,NULL,NULL,NULL))

mod.opt.glm <- list(NULL,NULL,NULL,NULL)

# RANDOM FOREST - mean decrease accurancy (MDA)
RF.list <- list(NULL,NULL,NULL,NULL)

## NNET optimal parameter
WT.opt <- list(NULL,NULL,NULL,NULL)
UNIT.opt <- list(NULL,NULL,NULL,NULL)

## --> end ignitia temporary list
#################################################################################################################################

## START THE VARIABLE SELECTION
## ============================
  
for (i in 1:n.rep)
    {
      ## SAMPLING
      cat("> Replication ",i," is starting now...", "\n",append = FALSE)
  	  cat("#############################################", "\n",append = FALSE)
      
      t <- 1 # used before with several sampling technique   
      
#################################################################################################################################     

## i) COMPARE RASTER DATA AND THE EXTRACT VALUES
## =============================================

if (COMP=="NO"){ }else{

Raster.data <- NULL # keep raster data
Raster.var <- NULL # keep var name
      
for (q in 1:length(fns))
      {
  
      cat("Comp. raster VS extract for ", var.list[[q]][1], "\n",append = FALSE)
      
      # read the raster 
      pred.var <- raster(fns[q])
      
      # re-extract data
      ext.data <- extract(pred.var,ext)
      
      if (i == 1) { # plot only the graph for rep. run 1
      jpeg(paste(outpath,"/Pred.var_plot_",var.list[[q]][1],".jpg",sep=""),width = 35, height = 14,units="cm",res=150)
      
        plot(pred.var, main = var.list[[q]][1])
      
      dev.off() } else {}
      
        if (mask.proj == "YES")
            {
            pred.var <- mask(pred.var,raster(paste(inpath,"/",mask.r,sep="") ) )
          
            }else {}
        
        ## Calculate summary of the raster data set 
        pred.mean <- cellStats(pred.var, stat='mean')
        pred.sd <- cellStats(pred.var, stat='sd')
        pred.min <- cellStats(pred.var, stat='min')
        pred.max <- cellStats(pred.var, stat='max')
        
        ## Calculate summary of the extract data set # remove the NA
        ext.mean <- mean(ext.data, na.rm = TRUE)
        ext.sd <- sd(ext.data, na.rm = TRUE)
        ext.min <- min(ext.data, na.rm = TRUE)
        ext.max <- max(ext.data, na.rm = TRUE)
        
        # keep summary data
        Raster.data <- rbind( Raster.data,c(pred.min, pred.mean, pred.sd, pred.max, 
                                            ext.min, ext.mean, ext.sd, ext.max) ) 
        ## APPLIED DATA MODIFICATION
        f <- as.function(transf[[q]]) # define function
        pred.var <- calc(pred.var, fun=f) #
      
        #RESCALE THE PREDICTOR with normal rescaling
        pred.mean <- cellStats(pred.var, stat='mean')
        pred.sd <- cellStats(pred.var, stat='sd')
        
        pred.var <- (pred.var-pred.mean)/pred.sd # function rescaling
        
        names(pred.var) <- var.list[[q]][1]
      
      ## Plot comparison
      ##################
      jpeg(paste(outpath,"/rastVSextr_plot_",var.list[[q]][1],"rep_",i,".jpg",sep=""),width = 35, height = 14,units="cm",res=150)
      
            # compare data and the extracted data
          # layout boxplot is at the bottom 
          nf <- layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(3,1))
          par(mar=c(3.1, 3.1, 1.1, 2.1))
          hist(pred.var,xlim=c(-4,4), col = "grey" )
          boxplot(x.in[,q], horizontal=TRUE,  outline=TRUE,
                  ylim=c(-4,4), frame=F, col = "green1", width = 10)
          
      dev.off()
      
      # iginite a sampling for the raster
      s <- sample( 1:ncell(pred.var), 50000 , replace=TRUE )
      
      p <- sample( 1:nrow(x.in), 50000 , replace=TRUE)
      
      ## compute a pair t-test
      t_comp <- t.test(x.in[p,q],pred.var[s] ,paired = TRUE, var.equal = TRUE)
      
      if (q==1){
          out.table <- cbind(t=t_comp$statistic,p=t_comp$p.value)
          row.names(out.table) <- var.list[[q]][1]
        }else{
          out.table <- rbind(out.table,c(t=t_comp$statistic,p=t_comp$p.value))
          row.names(out.table)[q] <- var.list[[q]][1]
        }
      file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

      } # end of loop pred. var 

write.table(out.table, file=paste(outpath,"/pred.var_t_test_",i,".txt",sep=""),sep="\t",
            append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)

row.names(Raster.data) <- row.names(out.table)

colnames(Raster.data) <- c("pred_min","pred_mean","pred_sd","pred_max", 
                            "ext_min","ext_mean","ext_sd","ext_max")

write.table(Raster.data, file=paste(outpath,"/pred.var_summary_",i,".txt",sep=""),sep="\t",
            append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)

} # end technique
#################################################################################################################################     

## i) MORAN ANALYSIS --> SPATIAL INTERPOLATION
## =============================================

if (MI=="NO"){ }else{

Raster.var <- NULL # keep var name
Moran.index <- NULL

# calculate the distance matrix for the coordonate
mat.dist <- as.matrix(dist(coor.xy))

mat.dist.inv <- 1/mat.dist
diag(mat.dist.inv) <- 0

# compute Moran index for all variable
# follow Gittleman and Kot (1990)
# p> 0.05 to reject the null hypothesis of spatila autocorrelation if p > 0.05 no autocorrelation!!

for (im in 1:ncol(data.in))
{
 test <- Moran.I(data.in[,im], mat.dist.inv, scaled = TRUE, na.rm = TRUE, alternative = "greater")

 Moran.index <-rbind(Moran.index,c(test$observed,test$expected,test$sd,test$p.value))

if (im==1) { Raster.var <- names(d.samp)[dep.var] 
	}else{Raster.var <- c(Raster.var,var.list[[(im-1)]][1]) }
	
}

colnames(Moran.index) <- c("observed","expected","sd","p.value")
rownames(Moran.index) <- Raster.var

write.table(Moran.index, file=paste(outpath,"/Moran.index_",i,".txt",sep=""),sep="\t",
            append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)

} # end technique
#################################################################################################################################

## ii) CORR PLOT BETWEEN PREDICTOR
## ===============================

if (CORR=="NO"){ }else{

df.cor.pred.var.l[[t]][[i]] <-data.frame(round(cor(x.in),3))

# Save correlations in your workspace
write.table(df.cor.pred.var.l[[t]][[i]],file=paste(outpath,"/Cor_Pred_var_rep_",i,".txt",sep=""),
            sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)

df.cor.pred.all.l[[t]][[i]] <-data.frame(round(cor(data.in),3))

# Save correlations in your workspace
write.table(df.cor.pred.all.l[[t]][[i]],file=paste(outpath,"/Cor_Pred_all_rep_",i,".txt",sep=""),
            sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)

} # end technique
#################################################################################################################################

## iii) CLUSTER ANALYSIS of the pred.variable
## ======================
# ref: Suzuki, R. and Shimodaira, H. (2006)

if (CLUST=="NO"){ }else{

cat("-> Cluster analysis...", "\n",append = FALSE)

        # Ward Hierarchical Clustering with Bootstrapped p values
        fit <- pvclust(x.in, method.hclust="ward",
                       method.dist="euclidean",nboot=10000)

        jpeg(paste(outpath,"/Pred.var_cluster_rep_",i,".jpg",sep=""),
             width = 30, height = 15,units="cm",res=150)
        
            plot(fit) # dendogram with p values
            # add rectangles around groups highly supported by the data
            pvrect(fit, alpha=.95) 
        
        dev.off()

} # end technique
#################################################################################################################################

## iV) PCA
## ========

if (PCA=="NO"){ }else{

      ## PCA - variable ##
      ######################
      pca.pred.var <- dudi.pca(x.in,scannf = FALSE, nf= ncol(x.in) )
      
      # Variance explained by all axe 
      mat.var.pca.l[[t]][i,] <- pca.pred.var$eig[1:(length(var.list)-2)]/sum(pca.pred.var$eig)
      
      ## PCA - all ##
      ######################
      pca.pred.all <- dudi.pca(data.in,scannf = FALSE, nf= ncol(data.in) )
      
      # Variance explained by all axe in ratio
      mat.all.pca.l[[t]][i,] <- pca.pred.all$eig[1:(length(var.list)-2)]/sum(pca.pred.all$eig)
      
      # Plot PCA correlation circle #
      ###############################
      #if (i == 1) { # write a plot only if first run
        jpeg(paste(outpath,"/PCA_CorCircle_rep_",i,".jpeg",sep=""),width = 30, height = 15,units="cm",res=150)
          par(mfrow=c(1,2))
        
            s.corcircle(pca.pred.var$co,clabel=.7, cgrid = 2, full = FALSE, sub = "var.indep", 
                        csub = 2.5, possub = "bottomleft", box = TRUE)
            
            s.corcircle(pca.pred.all$co,clabel=.7, cgrid = 2, full = FALSE, sub = "all.var", csub = 2.5, 
                        possub = "bottomleft", box = TRUE)
        
        dev.off()
      
      
       pdf(paste(outpath,"/PCA_CorCircle_rep_",i,".pdf",sep=""),width = 11, height = 5.5)
           par(mfrow=c(1,2))
            
            s.corcircle(pca.pred.var$co,clabel=.7, cgrid = 2, full = FALSE, sub = "var.indep", 
                        csub = 2.5, possub = "bottomleft", box = TRUE)
            
            s.corcircle(pca.pred.all$co,clabel=.7, cgrid = 2, full = FALSE, sub = "all.var", csub = 2.5, 
                        possub = "bottomleft", box = TRUE)
        
        dev.off()
      #}else{}
      
      ## Compute root mean square loading factor --> pca var. importance
      ##########################################
      mydata <- x.in
      
      # Determine Number of Factors to Extract based on Kaiser (1960) equation
      #library(nFactors)
      ev <- eigen(cor(mydata)) # get eigenvalues
      ap <- parallel(subject=nrow(mydata),var=ncol(mydata),
                     rep=100,cent=.05)
      nf <- nScree(x=ev$values, aparallel=ap$eigen$qevpea) #plotnScree(nf) 
      nf <- nf$Components$nkaiser # extract Kaiser method of selecting PCA var number
      
      load.f <- pca.pred.var$co[,1:nf] # loading factor
      
      f <- function(x){(sum(x^2))^0.5}
      
      rms.load <- apply(load.f,1,f) # root mean square of important var.
      
      if (i==1) { mat.rms.load <- rms.load  
            }else{ mat.rms.load <- rbind(mat.rms.load,rms.load) 
                  row.names(mat.rms.load) <- seq(from=1,to=i, by=1) }
      
      # Save correlations in your workspace
      write.table(mat.rms.load,file=paste(outpath,"/pca_var_imp_",i,".txt",sep=""),
                  sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)

## compute euclidien distance (on PCA1 and PCA2)
##########################################

fun <- function(x){sqrt(sum(x[1]^2,x[2]^2))} # distance of vector two

euc.dist <- apply(pca.pred.var$co[,1:2],1, fun)

if (i==1) { mat.euc.dist<- euc.dist  
            }else{ mat.euc.dist <- rbind(mat.euc.dist,euc.dist) 
			row.names(mat.euc.dist) <- seq(from=1,to=i, by=1) }

      ## Plot PCA point classes 
      ###########################

      ## classified data point per cluster analysis
      ##############################################
      ## Zonal analysis 
      ############################################
      # Ward Hierarchical Clustering without Bootstrapped 
      # simple analysis by 5 classes
  
      d <- dist(mydata, method = "euclidean") # distance matrix
      fit <- hclust(d, method="ward")

      groups <- cutree(fit, k=5) # cut tree into 5 clusters
      
      jpeg(paste(outpath,"/PCA_clustering_rep_",i,".jpg",sep=""),width = 30, height = 15,units="cm",res=150)
          par(mfrow=c(1,2))  
          
          plot(fit)
          # draw dendogram with red borders around the 5 clusters
          rect.hclust(fit, k=5, border="red") 

          # data point
          s.class(pca.pred.var$li, cgrid = 2, fac= as.factor(groups),cstar = 0, 
              cellipse = 0, col=c("#a5cf84","#f052b1","#3d82eb","#e86519","#545454"))

      dev.off()

} # end technique
#################################################################################################################################

## V) NETWORK ANALYSIS SELECTION
## ===============================

if (NNET=="NO"){ }else{

      ## OPTIMIZED PARAMETER PERFORMANCE OF NNET
      ##########################################
      
      if (i ==1) {  cat("> NNET Parameter optimisation sampling ...", "\n",append = FALSE)
                    
                    opt.para <- nnet.para.opt(x.in, y.in, max.nods = UNIT.Max, it.max = it.max)

			  write.table(opt.para, file=paste(outpath,"/NNET_all_para_test.txt",sep=""),sep="\t",
                  	append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)

                    WT.opt[[t]] <- opt.para$WT[1]; UNIT.opt[[t]] <- opt.para$UNIT[1]
                      
                    WT <- WT.opt[[t]]; UNIT <- UNIT.opt[[t]]
        
			  write.table(c(WT_opt=WT.opt[[t]],Unit_opt = UNIT.opt[[t]]),
			              file=paste(outpath,"/NNET_Opt_para.txt",sep=""),
			              sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
        
			              n_train <- nrow(data.in) # number data used in var. selection

			  write.table(n_train,
			              file=paste(outpath,"/Number_data_var.sel.txt",sep=""),
			              sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
      
      } else {  WT <- WT.opt[[t]]; UNIT <- UNIT.opt[[t]]}

      cat("> NNET MODEL run ",i,"sampling  ...", "\n",append = FALSE)
      
      net.tr <- nnet(x.in,y.in,data=data.in, decay = WT, size = UNIT, linout = TRUE, 
                     maxit = it.max, trace= FALSE, Hess = TRUE)

      # PREDICT THE OBS
      pred.y <- predict(net.tr, newdata=x.in,type="raw",na.rm = TRUE)

      # PERFORMANCE OF THE MODEL -- follow Pineiro 2008
      lin.corr <- lm(y.in[,1] ~ as.vector(pred.y)) #Matrix ."slope","Int","R2","RMSE"

      if (nrow(summary(lin.corr)$coefficients)==1) # Fuck up check
      { 
        # TRAIN NEURAL NETS
        net.tr <- nnet(x.in,y.in, data=data.in, size = n, 
                    linout = TRUE, maxit = it.max, decay = w,
                    trace = FALSE)
        
      }else{ }
	
	    net.tr0 <- net.tr # keep old model for sensitivity analysis    
  
      # Cancelation variable procedure
      ################################
	    var.can <- nnet.cancel(x.train = x.in, y.train = y.in, n= UNIT,w = WT, it.max = it.max)
  
	    write.table(var.can,
	            file=paste(outpath,"/NNET_var_cancel_",i,".txt",sep=""),
	            sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
        
      ## Plot result
      ##############
    	jpeg(paste(outpath,"/NNET_plot_rep_",i,".jpg",sep=""),width = 30, height = 30,units="cm",res=150)
      
    	    plotnet(net.tr0)
      
    	dev.off()
      
      ## VAR. IMPORTANCE
      ######################
      ## weight algorithm (Garson 1991)
      ##-------------------------------
      
      g  <- garson(net.tr, "Log_Se", bar_plot = TRUE, x_lab = NULL,
                   y_lab = "Rel.importance", wts_only = FALSE)
      
      garson.plot.l [[t]][[i]] <- g
      
      par(mar=c(3,4,1,1),family='serif'); g
      
      ggsave(filename = paste(outpath,"/NNET_var.imp_rep_",i,".jpeg",sep=""), plot = last_plot(), width = 30, height = 15, units= "cm")
            
      ## MODEL PERFORMANCE
      ####################
      ## "train_slope","train_b","train_R2","train_RMSE", "test_slope","test_b","test_R2","test_RMSE", "AIC","AICc","CONV")))
	    ## CROSS - Validate model
    	CV.data <- CV.NNET.(net.tr, data.in, K=10, cv.lim = 10, name.sp=dep.var)
	 
	    obs <- CV.data$obs
      pred <- CV.data$predictions
      
      ## OBS VS PRED - CV set --> follow Pineiro 2008
      ##------------------------
      
      lin.corr <- lm(obs ~ pred) # Matrix ."slope","Int","R2","RMSE"
      modl.perf.l[[t]][i,1] <- summary(lin.corr)$coefficients[2,1] # slope
      modl.perf.l[[t]][i,2] <- summary(lin.corr)$coefficients[1,1] # Intercept
      modl.perf.l[[t]][i,3] <- summary(lin.corr)$r.squared # r squared
      modl.perf.l[[t]][i,4] <- mean((obs - pred),na.rm = TRUE) #ME
      modl.perf.l[[t]][i,5] <- mean(abs((obs - pred)),na.rm = TRUE) # MAE
      modl.perf.l[[t]][i,6] <- sqrt(mean((obs - pred)^2,na.rm = TRUE)) # RMSE

	    jpeg(paste(outpath,"/NNET_Train_OBSvsPRED_rep_",i,".jpg",sep=""),width = 30, height = 15,units="cm",res=150)
      
        ## Plot observed VS predicted for train
        plot (obs,pred, xlab="Observed_log_Se", ylab="Predict_log_Se")
        
    	dev.off()
      
    	## AIC and AICC of the Best model
    	################################
    	RSS <- sum((y.in - predict(net.tr, newdata=x.in, type="raw"))^2) 
    	
    	modl.perf.l[[t]][i,7] <- 2*sum(net.tr$wts!=0) - length(y.in)*log(RSS/length(y.in)) # AIC
    	
    	modl.perf.l[[t]][i,8] <- modl.perf.l[[t]][i,5] + (2*sum(net.tr$wts!=0)+(sum(net.tr$wts!=0)+1))/
    	  (length(y.in) -  sum(net.tr$wts!=0)-1) #AICc
  
   
      ## CONVERGENCE
      ##------------
      out <- capture.output(nnet(x.in,y.in,data=data.in, decay = WT, size = UNIT,linout = TRUE, maxit = it.max, Hess = TRUE))
     
      if (sum(str_count(out,"converged")) == 1){
        
          nb_conv <- as.numeric(str_extract_all(out[sum(str_count(out,"iter"))+2],"[[:print:]]{4}") [[1]][2])
          
          }else{
          nb_conv <- it.max
        }
      
      modl.perf.l[[t]][i,9] <- nb_conv/it.max
 } # end technique   
#################################################################################################################################

## VI) UNIVARIATE GLM
## ==================

if (GLM=="NO"){ }else{

      data.cal <- data.in    
      
      v.pred <- var.list
      
      mat.model.perf <- matrix(0,length(v.pred),4,dimnames=list(v.pred,c("glm.quad.d2","glm.quad.AIC","glm.quad.AICc","glm.quad_CV_error")))
      
      mat.model.optim <- matrix(0,length(v.pred),1,dimnames=list(v.pred,c("glm.quad.optim")))
      
      for (j in c(2:(length(v.pred)+1)))
        {                
        cat("-> Calibrating GLMs for variable ",names(data.cal[j]),"...", "\n",append = F)
        
        ## MODELS CALIBRATION ##
        ########################
        
        df.tmp.input <- na.omit(data.frame(data.cal[,c(1,j)]))
       
        row.names(df.tmp.input) <- 1:dim(df.tmp.input)[1]
        	
        # GLM with Quadratic term #
        ###########################
        
        glm.tmp.univ.quad <<- glm(eval(parse(text = paste(paste(names(df.tmp.input[1])), "~ poly(", names(df.tmp.input[2]),",2)", collapse = ""))),
                                  data=df.tmp.input,family=gaussian,maxit = it.max) 
        
        # MODEL FIT - DEVIANCE 
        #--------------------                   
        
        # Model Fit = Adj. resc.Nagelkerke R2        
        mat.model.perf[j-1,1] <- round(adj.D2.glm(glm.tmp.univ.quad),3)
       
        # MODEL PREDICTIVE POWER 
        #-----------------------                         
        
        df.out.cv.glm <- cv.glm(df.tmp.input, glm.tmp.univ.quad, K=10)
        
        glm.val <- predict.glm(glm.tmp.univ.quad , newdata= df.tmp.input[2], type = "response")
        
        #CV estimated predicto error and adj. cross vali estimation
        mat.model.perf[j-1,4] <- df.out.cv.glm$delta[1]
        
        # MAX-AIC AICc
        mat.model.perf[j-1,2] <- AIC(glm.tmp.univ.quad) 
        
        mat.model.perf[j-1,3] <- AICc(glm.tmp.univ.quad)
                                     
        ###########################################
        ## PLOTS OF RESPONSE CURVES AND OPTIMUMS ##	
        ###########################################
        
        jpeg(paste(outpath,"/GLM_Resp.curves_",names(df.tmp.input[2]),"_rep_",i,".jpg",sep=""),
             width = 30, height = 15,units="cm",res=300)
              
              xmin <- min(df.tmp.input[,2])	
              xmax <- max(df.tmp.input[,2])
              ymax <- max(glm.tmp.univ.quad$fitted)	
              ymin <- min(glm.tmp.univ.quad$fitted)
              
              ###################
              ## GLM QUADRATIC ##
              ###################					   		   
              
              plot(df.tmp.input[,2],glm.tmp.univ.quad$fitted,xlab = names(df.tmp.input[2]),ylab = paste("log(",dep.var,")",sep=""),
                   main = paste("GLM"," ~ pol(",names(df.tmp.input[2]),",2)",sep=""),xlim = c(xmin,xmax), ylim = c(ymin,ymax),type="n")
              
              df.glm.optim <- data.frame(VAR=df.tmp.input[,2],PRED=glm.tmp.univ.quad$fitted)
              
              # lines(df.tmp.pred.glm.tmp.univ.quad[order(df.tmp.pred.glm.tmp.univ.quad[,2]),],col="red",lwd=4)
              
              abline(h=0,col="black")
              abline(h=1,col="black")  
              
              points(df.tmp.input[,2],glm.tmp.univ.quad$fitted,cex=.1,col="darkgrey")
              points(df.tmp.input[df.tmp.input[,1]==1,2],rep(1,length(df.tmp.input[df.tmp.input[,1]==1,1])),col="green",pch=4,cex=.7)
              points(df.tmp.input[df.tmp.input[,1]==0,2],rep(0,length(df.tmp.input[df.tmp.input[,1]==0,1])),col="red",pch=1,cex=.7)
              
              if (length(unique(df.glm.optim[df.glm.optim$PRED==max(df.glm.optim$PRED),"VAR"]))<2)
              {
                mat.model.optim[j-1,1] <- unique(df.glm.optim[df.glm.optim$PRED==max(df.glm.optim$PRED),"VAR"])
                abline(v=mat.model.optim[j-1,1],col="red")
              } else {
                
                mat.model.optim[j-1,1] <- mean(unique(df.glm.optim[df.glm.optim$PRED==max(df.glm.optim$PRED),"VAR"]))
                abline(v=mat.model.optim[j-1,1],col="red")
              }
              
              text(xmin,0.95,labels=paste("Adj.D2 = ",mat.model.perf[j-1,1],sep=""),pos=4)	  
                      
        dev.off()	
        
        
        #############################################################################
        # End of loop predictors
      }		
 
      # Export output #
      #################
      
      write.table(mat.model.perf, file=paste(outpath,"/GLM_Model_perf_rep_",i,".txt",sep=""),
                  sep="\t",append=F,row.names=T,col.names=T,quote=F)
      write.table(mat.model.optim, file=paste(outpath,"/GLM_Model_optim_rep_",i,".txt",sep=""),
                  sep="\t",append=F,row.names=T,col.names=T,quote=F)
      #write.table(mat.model.permut.glm, file=paste(outpath,"/",s.tech[t],"/Model_permut_rep_",i,".txt",sep=""),sep="\t",append=F,row.names=T,col.names=T,quote=F)

      if (i==1) mod.perf.glm[[t]][[1]] <- mat.model.perf[,1] else mod.perf.glm[[t]][[1]] <- cbind(mod.perf.glm[[t]][[1]],mat.model.perf[,1]);
      if (i==1) mod.perf.glm[[t]][[2]] <- mat.model.perf[,2] else mod.perf.glm[[t]][[2]] <- cbind(mod.perf.glm[[t]][[2]],mat.model.perf[,2]);
      if (i==1) mod.perf.glm[[t]][[3]] <- mat.model.perf[,3] else mod.perf.glm[[t]][[3]] <- cbind(mod.perf.glm[[t]][[3]],mat.model.perf[,3]);
      if (i==1) mod.perf.glm[[t]][[4]] <- mat.model.perf[,4] else mod.perf.glm[[t]][[4]] <- cbind(mod.perf.glm[[t]][[4]],mat.model.perf[,4]);

      if (i==1) mod.opt.glm[[t]] <- mat.model.optim else mod.opt.glm[[t]] <- cbind(mod.opt.glm[[t]],mat.model.optim);

 } # end technique     
#################################################################################################################################

## VII) LASSO GLM PROCEDURE
## ============================

if (LASSO=="NO"){ }else{

  cat("-> LASSO PROCEDURE...", "\n",append = FALSE)
      
      x.in <- data.in[,2:(length(var.list)+1)]
      y.in <- data.in[,1]
      
      pdf(paste(outpath,"/Lasso_procedure_rep_",i,".pdf",sep=""),width = 11, height = 5.5)
          
      par(mfrow=c(1,2))
      
          fit.cv <- cv.glmnet(as.matrix(x.in), y.in , family="gaussian",nfolds=10)
          
          plot(fit.cv,cex.lab=1.5,cex.axis=1.5, las=1)
          
          headings <- dimnames(x.in)[[2]]
      
          color <- rainbow(length(headings), s = 1, v = 1, start = 0, 
                          end = max(1, length(headings) - 1)/length(headings), alpha = 1)
          
          plot(fit.cv$glmnet.fit,"norm",label=TRUE,
               col=color, lw=2,cex.lab=1.5,cex.axis=1.5, las=1)
        
          legend("bottomleft", headings,col=color, lty=1:length(headings),lw=2)
          
      dev.off()

} # end technique    
#################################################################################################################################

## VIII) RANDOM FOREST TREE ANALYSIS
## ===============================

if (RF=="NO"){ }else{

      cat("--> RANDOM FOREST TREE...", "\n",append = FALSE)

      RF.table <- RF.var.sel (x.in,y.in,nb.tree = nb.tree)
      
      if (i==1) RF.list[[t]] <- RF.table else RF.list[[t]] <- cbind(RF.list[[t]],RF.table);
      
    	#################################################################################################################################
    	# End of loop t sampling method
    	#}				
    	
  } # end technique

}# End of loop n rep

################################################################################################################################
#######################   
# Export output means #
######################

creat.subDir(mainDir= outpath, subDir = "average")

  #for ( t in 1:length(s.tech))
       # {
    t <- 1
        ####################
        ## II) NNET VAR. SEL
        ####################
        
        ## NNET - Model PerformANCE table
        #--------------------------------
        Model_perf <- cbind(colMeans(modl.perf.l[[t]], na.rm = TRUE), apply(modl.perf.l[[t]],2,sd))
        
        dimnames(Model_perf) <- list(dimnames(modl.perf.l[[t]])[[2]],c("Average","SD"))
        
        write.table(modl.perf.l[[t]], file=paste(outpath,"/NNET_Model_perf_x-run.txt",sep=""),
                    sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
        
        write.table(Model_perf, file=paste(outpath,"/average/NNET_Model_perf.txt",sep=""),
                    sep="\t", append=FALSE, quote=FALSE)
        
        #####################
        ### Variable imp NNET
        #####################
        garson.plot <- garson.plot.l[[t]]
        
        # Re-order data by row names
        order <- as.numeric(row.names(garson.plot[[1]]$data))
        sort <- sort(order, index.return=TRUE)
        garson.plot[[1]]$data[,1] <- garson.plot[[1]]$data[sort$ix,1] 
        garson.plot[[1]]$data[,2] <- garson.plot[[1]]$data[sort$ix,2] 
        
        mean.plot <- garson.plot[[1]]$data
        row.names(mean.plot) <- sort$x
        
        sd.plot <- list(mean.plot)
        
        for (n in 2:n.rep) # loop to calculate mean
            {
              # Re-order data by row names for each plot
              sort <- sort(as.numeric(row.names(garson.plot[[n]]$data)), index.return=TRUE)
              garson.plot[[n]]$data[,1] <- garson.plot[[n]]$data[sort$ix,1] 
              garson.plot[[n]]$data[,2] <- garson.plot[[n]]$data[sort$ix,2]
              
              mean.plot[,1] <- garson.plot[[n]]$data[,1] + mean.plot[,1]
              sd.plot[[n]] <- garson.plot[[n]]$data
            }
        
        mean.plot[,1] <- mean.plot[,1]/n.rep
        
        std <- c()
        
        for (q in 1:nrow(sd.plot[[1]])) # loop for sd calculation
            { 
              temp.pred.var <- c()
              for (d in 1:n.rep)
                  {
                    temp.pred.var <- c(temp.pred.var, sd.plot[[d]][q,1])
                  }
              std <- c(std,sd(temp.pred.var))    
            }
        
        # Re-order data by rel_imp --> IDEA Toi be similar than input to have everything on the right order --> don't work properly
        sort <- sort(mean.plot[,1], index.return=TRUE)
        mean.plot[,1] <- mean.plot[sort$ix,1] 
        mean.plot[,2] <- mean.plot[sort$ix,2] 
        std <- std[sort$ix]
  
        ## Write mean data into a table
        row.names(mean.plot) <- sort$ix
        
        g.table <- cbind(mean.plot[,1],std)
        
        dimnames(g.table) <- list(mean.plot[,2],c("Average","SD"))
        
        write.table(g.table,file=paste(outpath,"/average/NNET_var.imp_mean.txt",sep=""))
        
        ## Plot nice bar plot
        ###########################
        
        pdf(paste(outpath,"/average/NNET_var.imp_mean.pdf",sep=""),width = 5.5, height = 5.5)
        
            bp <- barplot(mean.plot[,1],col=c("white"),beside=TRUE, axisnames=TRUE, space=c(0,0.5),
                          ylab= "Rel. importance",
                          ylim=c(0,(max(mean.plot[,1])+max(std))),
                          names.arg = mean.plot[,2],
                          cex.axis=1.5, cex.lab=1.5,las=2)            
            
            # the error bars
            ##################
            # Plot the vertical lines of the error bars
            # The vertical bars are plotted at the midpoints
            dg1 <- mean.plot[,1]
            stDevs <- std
            segments(bp, dg1 - 0, bp, dg1 + stDevs, lwd=1)
            # Now plot the horizontal bounds for the error bars
            # 2. The upper bar
            segments(bp - 0.2, dg1 + stDevs, bp + 0.2, dg1 + stDevs, lwd=1)
        
        dev.off()
        
        ##########################
        ## III) PCA and Cor. PLot
        ##########################
        
        var <- cbind(colMeans(mat.var.pca.l[[t]], na.rm = TRUE), apply(mat.var.pca.l[[t]],2,sd))
        all <- cbind(colMeans(mat.all.pca.l[[t]], na.rm = TRUE), apply(mat.all.pca.l[[t]],2,sd))
        
        dimnames(var) <- list(paste("Var._",seq(from=1,to=(length(var.list)-2),by=1)),c("Average","SD"))
  
        dimnames(all) <- list(paste("Var._",seq(from=1,to=(length(var.list)-2),by=1)),c("Average","SD"))
  
        write.table(var, file=paste(outpath,"/average/PCA_var.txt",sep=""),
                    sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
  
      	write.table(all, file=paste(outpath,"/average/PCA_all.txt",sep=""),
                    sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
              
        ## COR.Predict.table
        ##-------------------
        df.cor.pred.var <- df.cor.pred.var.l[[t]]
        df.cor.pred.all <- df.cor.pred.all.l[[t]]
        
        ## extract data and mean for cor. pred
        p <- nrow(df.cor.pred.all[[1]])
        names <- dimnames(df.cor.pred.all[[1]])
        
        std.cor.var <- matrix(0, p, p, dimnames = names)
        std.cor.all <- matrix(0, p, p, dimnames = names)
        
        for (q in 1:p)
            {
              for (r in 1:p)
                  {
                    temp.pred.var <- c()
                    temp.all <- c()
                    for (d in 1:n.rep)
                        {
                          temp.pred.var <- c(temp.pred.var, df.cor.pred.var[[d]][q,r])
                          temp.all <- c(temp.all,df.cor.pred.all[[d]][q,r])
                          
                        }
                        std.cor.var[q,r] <- sd(temp.pred.var)
                        std.cor.all[q,r] <- sd(temp.all)    
                  }  
            }
        
        mean.cor.var <- matrix(0, p, p, dimnames = names)
        mean.cor.all <- matrix(0, p, p, dimnames = names)
        
        for (n in 1:n.rep)
            {
              mean.cor.var <- df.cor.pred.var[[n]] + mean.cor.var
              mean.cor.all <- df.cor.pred.all[[n]] + mean.cor.all 
            }
        
        mean.cor.var <- mean.cor.var/n.rep
        mean.cor.all <- mean.cor.all/n.rep
        
        ## write mean and std cor predict
        write.table(std.cor.var,file=paste(outpath,"/average/cor.var.std.cor.txt",sep=""),
                    sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
        write.table(std.cor.all,file=paste(outpath,"/average/cor.all.std.txt",sep=""),
                    sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
        
        write.table(mean.cor.var,file=paste(outpath,"/average/cor.var.mean.txt",sep=""),
                    sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
        write.table(mean.cor.all,file=paste(outpath,"/average/cor.all.mean.txt",sep=""),
                    sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
        
        ## root mean square load of PCA
        write.table(mat.rms.load,file=paste(outpath,"/PCA_rms_loading_x-run.txt",sep=""),
                    sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
        
        # compute average
        av.rms.load <- apply(mat.rms.load,2,mean)
        
        write.table(av.rms.load,file=paste(outpath,"/average/PCA_mean_rms_loading.txt",sep=""),
                    sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)

	#compute euclidienne distance
	write.table(mat.euc.dist,file=paste(outpath,"/PCA_euc.dist_x-run.txt",sep=""),
                    sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
        
        # compute average
        av.euc.dist <- apply(mat.euc.dist,2,mean)
        
        write.table(av.euc.dist,file=paste(outpath,"/average/PCA_mean_euc.dist.txt",sep=""),
                    sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)

        
        # std of the avearge
        sd.rms.load <- apply(mat.rms.load,2,sd)
        
        write.table(sd.rms.load,file=paste(outpath,"/average/PCA_std_rms_loading.txt",sep=""),
                    sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
        
        ## COR.Plot
        #----------
        ## need to define as matrix
        var <- as.matrix(mean.cor.var)
        all <- as.matrix(mean.cor.all)
          
        ## Plot option  
        plot.method<-"color"     # options= "circle", "square", "ellipse", "number", "shade", "color", "pie"
        plot.type<-"lower"         # options= "full", "lower", "upper"
        plot.order<-"original"     # options= "original", "AOE", "FPC", "hclust", "alphabet"
              
        col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                                   "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))  
        
        pdf(paste(outpath,"/average/corr_plot.pdf",sep=""),width = 11, height = 5.5)
            par(mfrow=c(1,2))
            
            corrplot(var,method=plot.method,type=plot.type, order=plot.order,outline=T,tl.col="black",
                     tl.offset=0.5,cl.length=10,sig.level=0.05,insig="pch",pch=1,pch.cex=23,col=col2(20))
            
            corrplot(all,method=plot.method,type=plot.type, order=plot.order,outline=T,tl.col="black",
                     tl.offset=0.5,cl.length=10,sig.level=0.05,insig="pch",pch=1,pch.cex=23,col=col2(20))
        
        dev.off()
      
      ######################
      ## IV) GLM UNIVARIATE
      #####################
      
      Av.perf.glm <- cbind(rowMeans(mod.perf.glm[[t]][[1]], na.rm = TRUE),apply(mod.perf.glm[[t]][[1]],1,sd),
                           rowMeans(mod.perf.glm[[t]][[2]], na.rm = TRUE),apply(mod.perf.glm[[t]][[2]],1,sd),
                           rowMeans(mod.perf.glm[[t]][[3]], na.rm = TRUE),apply(mod.perf.glm[[t]][[3]],1,sd),
                           rowMeans(mod.perf.glm[[t]][[4]], na.rm = TRUE),apply(mod.perf.glm[[t]][[4]],1,sd))
      
      dimnames(Av.perf.glm) <- list(dimnames(mat.model.perf)[[1]],
                                    paste(c("","SD_"),rep(dimnames(mat.model.perf)[[2]], each =2)))
      
      write.table(Av.perf.glm, file=paste(outpath,"/average/GLM.mod.perf.txt",sep=""),sep="\t",
                  append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
      
      av.opt.glm <- cbind(rowMeans(mod.opt.glm[[t]], na.rm = TRUE),apply(mod.opt.glm[[t]],1,sd))
      
      dimnames(av.opt.glm) <- list(dimnames(av.opt.glm)[[1]], c("Average","SD_")) 
                  
      
      write.table(av.opt.glm, file=paste(outpath,"/average/GLM.mod.opt.txt",sep=""),sep="\t",
                        append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
      
      #####################
      ## VI) RANDOM FOREST
      #####################
      write.table(RF.list[[t]],file=paste(outpath,"/random_forest_xrun.txt",sep=""),sep="\t",
                  append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
     
      RF.average <- cbind(rowMeans(RF.list[[t]], na.rm = TRUE), apply(RF.list[[t]],1,sd))
     
   	  dimnames(RF.average) <- list(dimnames(RF.average)[[1]], c("Average","SD_"))   	
      
      write.table(RF.average,file=paste(outpath,"/average/random_forest.txt",sep=""),sep="\t",
                  append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)      
   # }

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#### COFFE TIMES
#################################################################################################################################
#################################################################################################################################




