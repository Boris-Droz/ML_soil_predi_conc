####################################################################################################################################
####################################################################################################################################
###					                ###
###  MODEL PROJECTION		    ###   
###                         ###					
####################################################################################################################################
####################################################################################################################################
## Historic
## --------
## v7.0 September 2016 - Boris DROZ & Gerrad Jones, ETHZ & EAWAG
## v8.0 Januar 2017 - Boris DROZ --> more technic implemented
## v9.0 include xboost and ranger
## v10.0 include extraction outlier (option)
## v12.0 January 2018 - inclu stop criteria --> decrease over-fitting
####################################################################################################################################
## DESCRIPTION
################    #####################################
#                   --- CHECK PARAMETER UNTIL LINE 225 ---
##                   #####################################

##    Tested under R 3.1
#########################

## Sampling by equal scale for raster and a point data set
## Rescaling all variable 
## Performed separately several modelling technics (choose by user)
##        Machine learning model
##        ----------------------
##                  - Neural Networks (nnet -- nnet package v.7.3-12)
##                  - Stuttgart Neural Network Simulator (SNNS -- RSNNS package v.0.4-7) 
##                  - Extrem learning machine (ELM -- elmNN package v. 1.0)
##                  - Random Forest (RF -- randomForest package v.4.6-12) 
##                  - Case-specific Random Forest (CSRF --- ranger package v 0.6.0)
##
##                  - Bagging trre method (BAG -- ipred package v.0.9-5)
##                  - Extreme Gradient Boosting (EGB --- xgboost v.0.6-4)
##                  - Generalized Boosted Regression Models (GBM -- gbm package v.2.1.1)
##            
##        Simple linear model
##        -------------------
##                  - Generalized Linear Models (GLM)
##                  - GLM step-wise (GLM_sw)
##                  - Generalized Additive Models (GAM --package gam v. 1.14)
##
# Input: 
#=======  i) data : var. to predict and cooordinate for each data point 
#                  *txt tab file with 3 header"X_WGS84","Y_WGS84" and "Var_to_pred"
#         ii) path.pred : # folder with all *.geotiff "predictive variables"
##                        should have same extent, coordinate and resolution.
#
# Output: - model projection
#=======  - result from the cross validation
#         - Importance are compute using garson for nnet and rsnns, mean decrease in node impurity for tree method
#             and permutation procedure for linear model
# 
####################################################################################################################################
## Load library
##################

## Set library path
#.libPaths("C:/Users/Public/Documents/R/win-library/3.1")

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

## library parrallel core
library(snowfall)
library(parallel)

####################################################################################################################################
## SCRIPT PARAMETER  --> NEED TO GO TROUH AND MAKE THE APPROPRIATE MODIFICATION !!!
#####################

# Set folder 
## =========
inpath <- "D:/copper/input/pred"
output <-"D:/copper/output/20180325_model_v18"

outpath <- output ## DON'T MODIFIED ##

path.pred <- paste(inpath,"/current",sep="") 

## set your work space ## DON' MODIFIE ##
setwd(inpath) ## DON'T MODIFIED ##

## DATA Point
##===========
## X Y coordonnate should be nominate : X_WGS84  and Y_WGS84
## col names should be identical to var.y and var.x
# File name of data points 
data <- "20170613_Data_Cu_Vine_EU"

Raster.data <- "YES" # specify the pred var if no include predicitve variable in data above

# in case of Raster.data <- "NO" should be specify: 
var.x <-  c("ZBRD", "ZlogAI", "ZET","ZlogPrecip", "ZlogIrrigation", "ZlogTOC_GEMAS", 
            "ZEmiss_PCA1", "ZEmiss_PCA2", "ZSoilChem_PCA1","ZSoilPhys_PCA1", "ZVeg_PCA1")

## RESAMPLING WITHIN A EQUAL SCALE INTO THE GRID CELL
#====================================================
resampling <- "NO" ## YES or NO

## reSampling paramater if resampling
# Define the degree grid cell resolution (Proj. WGS_1984)
siz.cell <- 0.002083333

# REMOVE OUTLIER ??
#####################
rem.out <- "NO" ## YES or NO

per.out <- 0.95 # which percentil is considered as outlier classic 0.95 based on z score
replace <- "YES" ## YES or NO if yes replacement by mean

## TRANSFORMATION
## ==============  
var.y <- c("CUt_mgkg") ; length(var.y)

transf.y <- list(alist(x=,log10(abs(x))) )

length (transf.y)# dep. variable transf must simi length as y.var

# Transformation list same number as predictive variable 
## # alist(x=,log10(abs(x))) ,alist(x=,sqrt(abs(x))) # EXEMPLE OF TRANSF.
#==============================================================

transf <- list( alist(x=,log10(abs(x))) , alist(x=,(x)) , alist(x=,(x)),  alist(x=,(1/x)), 
			 alist(x=,(x)) , alist(x=,(x)) , 
			 alist(x=,log10(abs(x))), alist(x=,log10(abs(x))),  alist(x=,(x)), alist(x=,(x)) )

length(transf)

# RESCALE VARIABLE --> rescale to Z-score to normalise data !!!
res.q <- "YES"

## NUMBER OF MODEL BUILDING
# ========================
n.rep <- 10

## -- CHOOSE MODELLING TECHN:
## YES or NO
#############################
# machine learning method
NNET <- "YES"
RSNNS <-"NO"
ELM <- "NO"

RF <- "YES"
CSRF <- "YES"

BAG <- "YES"
EGB <- "NO"
GBM <- "NO"

# linear additive
GLM <- "NO"
GLM_sw <- "NO"
GAM <- "NO"

## FUNCTION PARAMETER
## ==================
## PERMUTATION FOR var. opt. Calc. (glm, gam, gbm model)
n.perm <- 100

# xx fold cross validation calibration
fold.cv <- 10

## NNET PARAMETER
##################
## Weight decay Folow recomendation of B. D. Ripley: "Pattern Recognition and Neural Networks", Cambridge, 1996.
## between 0.1 and 0.01
deca <- 0.1

UNIT.Max <- 15 # NUMBER OF HIDDEN UNITS --> between the number of input nodes and number of output nodes 

it.max <- 10000 # Need to be 10'000 to be MONTE-CARLO Permutation

# RANDOM FOREST PARAMETER
## number of tree --> RF and BAG, GBM model (100 good compromise)
################
nb.tree <- 100

overfit="YES" # additionnal restriction used for RF and CSRF to avoid overfitting
alpha.lim = 0.1
##### .....

####################################################################################################################################
####################################################################################################################################
## ADD IN FUNCTION
#########################################################
#######################################################################################################
## lm paramter table
####################

asses.lm. <- function (x,y)
          {
            lin.corr <- lm(y ~ x) #Matrix ."slope","Int","R2","RMSE"
            slope <- summary(lin.corr)$coefficients[2,1] #slope
            std.slope <- summary(lin.corr)$coefficients[2,2]
            int. <- summary(lin.corr)$coefficients[1,1] #Intercept
            std.int <- summary(lin.corr)$coefficients[1,2]
            rsquare <- summary(lin.corr)$r.squared #r squared
            MAE <- mean( abs((y - x)),na.rm = TRUE )# mean abs error
            MSE <- mean((y - x)^2,na.rm = TRUE) # mean square error
            RMSE <- sqrt(mean((y - x)^2,na.rm = TRUE)) #root mean square error
            
            lm.vector <- c(slope=slope, std.slope=std.slope, int.=int., std.int=std.int, 
                           rsquare=rsquare, MAE=MAE, MSE=MSE, RMSE=RMSE)
            
            return(lm.vector) 
          }  
  
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

#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
##########   ---- CROSS VALIDATION FUNCTIONS   ##################################################################################################
#################################################################################################################################################
## All fucntion are adapted from Ecospat package for continuous data
## All data with same weights are considered

## Historic
## --------
## CV-NNET Droz. B - 12.6.2015 -modified 7.9.2015 / 10.2.2017
## CV-RSNNS Droz. B - February 2017
## CV-ELM Droz. B - February 2017
## CV- RF Droz. B - 12.6.2015
## CV - CSRF Droz B. 17.3.2017
## CV- BAG Droz. B - February 2017
## CV- GBM Droz. B - February 2017

## CV- GLM Droz. B - 12.6.2015 --> work for GAM too
## CV- GLM step Droz. B - 12.6.2015
## 

#################################################################################################################################################
## CV-NNET ##
CV.NNET. <- function(data.cv, nnet.unit, nnet.WT=0.01, it.max =1000, K=10, cv.lim = 10, name.sp)
{
  
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 ) 
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  
  K.lst <- K
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  
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
    
    x.cal <- data.cal.cv[,colnames(data.cal.cv)!=name.sp]
    y.cal <- data.cal.cv[,colnames(data.cal.cv)==name.sp]
    
    x.test <- data.test.cv[,colnames(data.test.cv)!=name.sp]
    
    nnet.cal <-  nnet(x.cal,y.cal, data=data.cal.cv, decay = nnet.WT, size = nnet.unit, 
                      linout = TRUE, maxit = it.max, Hess = TRUE, trace = FALSE)  
    
    nnet.val <- predict(nnet.cal, newdata = x.test , type = "raw",na.rm = TRUE)
    
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
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0),
                                    predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
   return(df.res)    
  }
################################################################################################################################################# 
## CV-RSNNS ##
CV.RSNNS. <- function(data.cv, nral.hidden= round(ncol(data.cv)*2/3) , it.max =1000, K=10, cv.lim = 10, name.sp)
{
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 ) 
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  
  K.lst <- K
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  
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
    
    x.cal <- data.cal.cv[,colnames(data.cal.cv)!=name.sp]
    y.cal <- data.cal.cv[,colnames(data.cal.cv)==name.sp]
    
    x.test <- data.test.cv[,colnames(data.test.cv)!=name.sp]
    
    nnet.cal <- mlp(x.cal,y.cal, size= nral.hidden, maxit= it.max, learnFunc = "Std_Backpropagation")
    
    nnet.val <- predict(nnet.cal,x.test)
    
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
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0),
                                    predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}
#################################################################################################################################################
CV.ELM. <- function(data.cv, elm.obj, K=10, cv.lim = 10, name.sp)
{
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 ) 
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  
  K.lst <- K
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  
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
    
    x.cal <- data.cal.cv[,colnames(data.cal.cv)!=name.sp]
    y.cal <- data.cal.cv[,colnames(data.cal.cv)==name.sp]
    
    x.test <- data.test.cv[,colnames(data.test.cv)!=name.sp]
    
    elm.cal <- update(elm.obj, data = data.cal.cv, weights=rep(1,nrow(data.cal.cv))) 
    
    elm.val <- predict(elm.cal,x.test)
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(elm.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(elm.val), after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0),
                                    predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}

#################################################################################################################################################
## CV-RF ##
CV.RF. <- function(RF.obj, data.cv, K=10, cv.lim = 10, name.sp)
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
    
    rf.cal <- update(RF.obj, data = data.cal.cv, weights=rep(1,nrow(data.cal.cv)))
    
    rf.val <- predict(rf.cal, data.test.cv , type = "response",weights=rep(1,nrow(data.test.cv)))
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(rf.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(rf.val), after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0) ,predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}

################################################################################################################################################
## CV-CSRF ##
CV.CSRF. <- function(CSRF.obj, data.cv, K=10, cv.lim = 10, name.sp)
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
    
    csrf.cal <- update(CSRF.obj, data = data.cal.cv) 
    csrf.val <- predict(csrf.cal, data.test.cv , type = "response") 
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(csrf.val)$predictions
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(csrf.val)$predictions, after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0) , predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}

#################################################################################################################################################
## CV-BAG ##
CV.BAG. <- function(data.cv, K=10, cv.lim = 10, nb.tree= 100, name.sp)
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
    
    x.cal <- data.cal.cv[,colnames(data.cal.cv)!=name.sp]
    y.cal <- data.cal.cv[,colnames(data.cal.cv)==name.sp]
    
    bag.cal <- ipredbagg (y.cal , x.cal, nbag=nb.tree)
    
    bag.val <- predict(bag.cal, newdata = data.test.cv , type = "raw" , na.rm = TRUE)
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(bag.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(bag.val), after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0) ,predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}

#################################################################################################################################################
## CV-EGB ##
CV.EGB. <- function(data.cv, K=10, cv.lim = 10, it.max= 1000, name.sp)
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
    
    # CREAT DATA SET FOR x-FOLD CV
    data.cal.cv <- data.cv[j.in,]
    data.test.cv<- data.cv[j.out,]
    
    x.cal <- data.cal.cv[,colnames(data.cal.cv)!=name.sp]
    y.cal <- data.cal.cv[,colnames(data.cal.cv)==name.sp]
    
    x.test <- data.test.cv[,colnames(data.test.cv)!=name.sp]
    
    ## -- MODEL --
    # set the parameter of the booster  (default)
    para <- list(booster= "gbtree", max_depth = 6, eta = 0.3, silent = 1, nthread = 2,
                 objective = "reg:linear", eval_metric = "rmse")
    
    # run model
    egb.cal <- xgboost(data = as.matrix(x.cal), label = t(as.vector(y.cal)), missing = NA, weight = NULL,
                         params = para , nrounds = it.max, verbose = 0 )
    
    egb.val <-  predict(egb.cal, as.matrix(x.test)) 
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- egb.val
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, egb.val, after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0) , predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}
#################################################################################################################################################
## CV-GBM ##
CV.GBM. <- function(gbm.obj, data.cv, nb.tree= 100, K=10, cv.lim = 10, name.sp)
{
  data.cv <- as.data.frame(data.cv)
  
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 ) 
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  
  K.lst <- K
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  
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
    
    x.cal <- data.cal.cv[,colnames(data.cal.cv)!=name.sp]
    y.cal <- data.cal.cv[,colnames(data.cal.cv)==name.sp]
    
    x.test <- data.test.cv[,colnames(data.test.cv)!=name.sp]
    
    gbm.cal <- update(gbm.obj, data = data.cal.cv, 
                      distribution = "gaussian") 
    
    gbm.val <- predict.gbm(gbm.cal, newdata=as.data.frame(x.test), 
                           n.trees=nb.tree, type="response",na.rm = TRUE)
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(gbm.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(gbm.val), after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0),
                                    predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}

###############################################################################################################
## CV-GLM   ##
CV.glm. <- function(glm.obj, K=10, cv.lim = 10, name.sp)
{
  
  data.cv <<- glm.obj$data
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
    
    glm.cal <- update(glm.obj, data = data.cal.cv, weights=rep(1,nrow(data.cal.cv)))
    
    glm.val <- predict(glm.cal, data.test.cv , type = "response",weights=rep(1,nrow(data.test.cv)))
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(glm.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(glm.val), after=length(vect.predicted))
    }	
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0) ,predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}

###############################################################################################################
## CV-GLM stepwise
CV.glm.step <- function(glm.obj, K=10, cv.lim = 10, name.sp)
{
  
  data.cv <<- glm.obj$data
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
    
    #w.cal <- rep(1,nrow(data.cal.cv))
    #w.test <- rep(1,nrow(data.test.cv))
    
    glm.cal <- update(glm.obj, data = data.cal.cv ) #, weights= w.cal)
    
    glm.val <- predict(glm.cal, data.test.cv , type = "response") #,weights= w.test)
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(glm.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(glm.val), after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0) ,predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}

#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
##########   ---- PERMUTATION VARIABLE FUNCTIONS   ##################################################################################################
##             --> calculate the relative importance of each variable !!!!!
#################################################################################################################################################

## Historic
## --------
## Original function write by C.Randin modified by B.Droz!
## var.imp.rf Droz. B - 12.6.2015
## var.imp.glm.step Droz. B - 12.6.2015 --> work for glm too

## Input :
## =======
### model :  model object
### cal : calibration dataset used to build the model
### names.pred: name of the predicting variables
### nperm : number of time each variable is permutated

####################################################################################################################################3333333
#########################
## VarImp RandomForest ##
#########################
var.imp.rf <- function(model,cal,names.pred,nperm=100)
{
	ref<-as.numeric(as.character(predict(model,cal)))
	VarImp<-vector()
	
	for (i in 1:length(names.pred))
	{
		print(names.pred[i])
		refi<-vector()
		
		for (j in 1:nperm)
		{
			if (j%%100==0)
			{
				cat("> Permutation ",j, "\n",append = FALSE)
			}
			
			cali<-cal
			cali[,names.pred[i]]<-cali[sample(1:nrow(cali),nrow(cali)),names.pred[i]]

			refi<-c(refi,1-cor(ref,as.numeric(as.character(predict(model,cali)))))
		}
		
		VarImp<-c(VarImp,round(mean(refi),3))
	}
	
	names(VarImp)<-names.pred
	return<-VarImp
}

## VarImp GLM, GLM Step, GAM
#############################
var.imp.glm.step <- function(model,cal,names.pred,nperm=100)
{
	
	ref<-predict(model,cal)
	VarImp<-vector()
	
	for (i in 1:length(names.pred))
	{
		print(names.pred[i])
		refi<-vector()
		
		for (j in 1:nperm)
		{
			if (j%%100==0)
			{
				cat("> Permutation ",j, "\n",append = F)
			}
			
			cali<-cal
			cali[,names.pred[i]]<-cali[sample(1:nrow(cali),nrow(cali)),names.pred[i]]
			refi<-c(refi,1-cor(ref,predict(model,cali)))
		}
		
		VarImp<-c(VarImp,round(mean(refi),3))
	}
	
	names(VarImp)<-names.pred
	return<-VarImp
}

####################################################################################################################################
####################################################################################################################################
## SCRIPT START HERE
####################################################################################################################################
## PARRALLEL CORE SETUP
##--------------------

ptm <- proc.time()# ignite timer

beginCluster(detectCores()) # ACTIVATE THIS MULTI CORE CALCULATION 

###################################################################################################################################
# -- LOADING DATA SET -----
##############################################

# Open the data set and the unknow point around the world
data.in <- read.table(paste(data,".txt",sep=""),header = TRUE ,na.strings = "NaN"); head(data.in)

if (Raster.data == "NO")
    {
      # if  you used data.in directly:  
      to.match <- c("X_WGS84","Y_WGS84",var.y, var.x)
      
      as.matrix(data.in[ ,match(to.match,names(data.in)) ] ) -> d.samp
      var.list <- var.x # keep a list of variable
      
      # position of the before the first x data set
      ncol.data0 <- min(match(var.x,colnames(d.samp) ) ) -1
      
      RESCAL <- NULL
      
    }else{
      
      ## pred. var file (batch all file in input) 
      fns <- list.files(path.pred,pattern=".tif$",full.names = TRUE); print(fns)
      
      # predicitve variable name list
      var.list <- gsub(".tif", "", list.files(path.pred,pattern=".tif$",full.names = FALSE)); print(var.list)
      
      var.x <- var.list
  
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
    
    RESCAL <- NULL
    
    ## RESAMPLING or NOT
    ####################
    if (resampling=="YES") 
    {
      ## RESAMPLE FIELD DATA
      ######################
      
      cat("> AVERAGE BY CELL GRID ...", "\n",append = FALSE)
      
    # define the avearge coordonnate
      data.in$X_WGS84 <- round( data.in$X_WGS84/siz.cell )*siz.cell
      data.in$Y_WGS84 <-  round( data.in$Y_WGS84/siz.cell )*siz.cell
      
      #cell.long <- unique(data.in$X_WGS84)
      #cell.lat <- unique(data.in$Y_WGS84)
    
    	# get unique row of new coordonate
    	cell <- unique(cbind(data.in$X_WGS84,data.in$Y_WGS84))
      
      d.samp <- matrix(data = NA, nrow = 0, ncol= ncol(data.in), dimnames = list(c(),dimnames(data.in)[[2]]))
      
      nb <- 0 # counter
      
      ## mean of all data in the same cell
    for ( s in 1: nrow(cell))
    {
    nb <- nb + 1
    
    d.temp <- data.in[data.in$X_WGS84==cell[s,1] & data.in$Y_WGS84==cell[s,2],]
    
    if (nrow(d.temp)==1) { d.samp <- rbind(d.samp,d.temp)
     }else{ d.samp <- rbind(d.samp, apply(d.temp,2,mean, na.rm=TRUE) )
    }
     
    }
    
      d.samp <- as.data.frame (d.samp)
      
      # Control same projection, spatial extent and resolution.
      # coord.ref similar ....
      ext <- cbind(d.samp$X_WGS84,d.samp$Y_WGS84)
      
      #colnames(r.proj) <- c("Longitude","Latitude",dep.var)
      
      ## RESAMPLE RASTER PRED VAR IF NECESSARY
      ########################################
      
      for (z in 1:length(fns))
      {
        cat("> Extract Pred.", var.list[[z]][1]," ...", "\n",append = FALSE)
        
        pred.var <- raster(fns[z])
        
        pos <- ncol(d.samp)
        
        # AGREGATE THE DATA OR NOT IN FUNCTION OF THE RESOLUTION                    
        if (1 < siz.cell/res(pred.var)[1] & 1==siz.cell/res(pred.var)[1] ) {
          
          pred.var <- aggregate(pred.var,fact=siz.cell/res(pred.var)[1],fun=mean)
          
        }else{ }  
        
        ## APPLIED DATA MODIFICATION
        f <- as.function(transf[[z]]) # define function
        pred.var <- calc(pred.var, fun=f) # applied modification
        
        ## RESCALE THE PREDICTOR with normal rescaling
        pred.mean <- cellStats(pred.var, stat='mean')
        pred.sd <- cellStats(pred.var, stat='sd')
        
        if (pred.mean<Inf & pred.mean>-Inf){}else{
          pred.mean <- mean(extract(pred.var,ext), na.rm=TRUE)
          pred.sd <- sd(extract(pred.var,ext), na.rm=TRUE)
        }
        
        RESCAL <- cbind( RESCAL,c(pred.mean,pred.sd) ) # keep rescalling data
        
        pred.var <- (pred.var-pred.mean)/pred.sd # function rescaling
        
        d.samp <- cbind(d.samp,extract(pred.var,ext))  
        
        #UNIFIED NAME OF THE RASTER
        colnames(d.samp)[pos+1] <- var.list[[z]][1]
        
      }
      
    } else {
      d.samp <- data.in 
      
      # coord.ref similar ....
      ext <- cbind(d.samp["X_WGS84"],d.samp["Y_WGS84"])
      
      # EXTRACT PREDICTING VARIABLE
      #############################
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
    
        if (pred.mean<Inf & pred.mean>-Inf){}else{
          	pred.mean <- mean(extract(pred.var,ext), na.rm=TRUE)
          	pred.sd <- sd(extract(pred.var,ext), na.rm=TRUE)
          	}
    
        RESCAL <- cbind( RESCAL,c(pred.mean,pred.sd) ) # keep rescalling data
        
        pred.var <- (pred.var-pred.mean)/pred.sd # function rescaling
        
        pos <- ncol(d.samp)
        
        d.samp <- cbind(d.samp,extract(pred.var,ext))  
        
        #UNIFIED NAME OF THE RASTER
        colnames(d.samp)[pos+1] <- var.list[[z]][1]
      }  
    }
    
    # extract outlier or not?
    #########################
    if (rem.out == "NO") {d.samp <- d.samp 
        }else{
          dat.cont <- d.samp[, (ncol(d.samp)-length(transf)+1) :ncol(d.samp)] # considered data
          
          ## replace the dat or not??
          if (replace == "NO"){
            fun.out <- function(x) {!scores(x, type="z", prob=per.out)} # function
            pos.out <-  apply(apply(dat.cont, 2, fun.out),1,sum)==9 # pos of outlier row
            d.samp <- d.samp[pos.out,]; nrow(d.samp) # delet the outlier data
          }else{
            fun.out <- function(x) {scores(x, type="z", prob=per.out)} # function
            dat.cont[apply(dat.cont, 2, fun.out)] <-0 # zero etannt mean because already z score normed
            d.samp[, (ncol(d.samp)-length(transf)+1) :ncol(d.samp)] <- dat.cont # replace data
          }
        }
    
    # data table of all values
    write.table(d.samp, file=paste(outpath,"/data_calib.txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
}# end loop selecting predicting variables

##########################################################################################################
## ---  BUILD MODEL FOR DIVERS DEPENDENT VARIABLE ---
#######################################################

# keep original output
output.0 <- output

for (tu in 1:length(var.y))
  {
  
  # creat a new folder for the considered dep. var
  mod.tech <- var.y[tu]
  outpath <- paste(output.0,"/",mod.tech,sep="")
  creat.subDir (output.0,mod.tech)
  
  output <- outpath # define the new output folder for the current dep.var.
  
  # redifine the de.var pos for the run 
  dep.var <- match( var.y[tu], colnames(d.samp) )
  x.var.pos <-match( var.x, colnames(d.samp) )

#####################################################################################################################################
## ORGINIZED DATA SET 
##======================================

data.in <- cbind(d.samp[,dep.var],d.samp[, x.var.pos])

colnames(data.in)[1] <- colnames(d.samp)[dep.var] # reload the name of dep.var

# modified the dep. var 
########################
### APPLIED DATA MODIFICATION for the dep.var
f <- as.function(transf.y[[tu]]) # define function 
data.in[,1] <- f(data.in[,1]) # applied modification

if (res.q == "NO") 
    {
    data.in <- na.omit(data.in) # delet row with NA value
  
    } else {## RESCALE THE PREDICTOR with normal rescaling

	  data.in <- na.omit(data.in) # delet row with NA value
      
    pred.mean <- mean(data.in[,1])
    pred.sd <- sd(data.in[,1]) 
    
    RESCAL <- cbind( c(pred.mean,pred.sd), RESCAL ) # keep rescalling data
    
    data.in[,1] <- (data.in[,1]-pred.mean)/pred.sd # function rescaling
    
    data.in <- na.omit(data.in) # delet row with NA value
    
    # rewrite names of header
    names(data.in) <- c(names(d.samp)[dep.var],names(d.samp)[(ncol.data0+1):ncol(d.samp)])
    colnames(RESCAL) <- c(names(d.samp)[dep.var],names(d.samp)[(ncol.data0+1):ncol(d.samp)])
    row.names(RESCAL) <- c("mean","sd")
    
    # write rescall data
    write.table(RESCAL, file=paste(outpath,"/rescal_para.txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    }

# creat x and y table
y.in <- data.in[,1,drop = FALSE]
x.in <- data.in[,2:ncol(data.in)]

#################################################################################################################################
#################################################################################################################################
## initialise the list data

# --- model performance ---
modl.perf.nnet <- NULL
modl.perf.rsnns <- NULL
modl.perf.elm <- NULL
modl.perf.rf <- NULL
modl.perf.csrf <- NULL

modl.perf.bag <- NULL
modl.perf.egb <- NULL
modl.perf.gbm <- NULL

modl.perf.glm <-NULL
modl.perf.glmstep <- NULL
modl.perf.gam <- NULL

# --- model importance ---
varimp.nnet <- NULL
varimp.rsnns <- NULL
varimp.elm <- NULL
varimp.rf <- NULL
varimp.csrf <- NULL

varimp.bag <- NULL
varimp.egb <- NULL
varimp.gbm <- NULL

varimp.glm <- NULL
varimp.glmstep <- NULL
varimp.gam <- NULL

# --- AVERAGE TABLE ---
mod.name <- NULL
modl.perf.av <- NULL
varimp.av <- NULL

count.while <- 0
pred.st.list <-list(NULL)

#####################################################################################################################################
#####################################################################################################################################      
## -------------------------   MODEL CALIBRATION ------------------------------------------
#####################################################################################################################################   
if (NNET == "NO") {}else{
    
  for (i in 1:n.rep)    
    {
    mod.tech <- "NNET"
    
    outpath <- paste(output,"/",mod.tech,sep="")
    
    creat.subDir (output,mod.tech)
    
    ######################
    ### NETWORK ANALYSIS
    #####################                       
    ## OPTIMIZED PARAMETER PERFORMANCE OF NNET
    ##########################################
 
      if (i ==1) {  cat("> NNET Parameter one layer ...", "\n",append = FALSE)
        
        ##(from the FAQ for a commercial neural network software company)  
        ## Number of inputs + outputs) * (2/3) -- Heaton 2008
        nb.hid1 <- round (ncol(data.in) *2/3)
        
        nb.hid2 <- round ( nb.hid1*2/3 )
        
        nb.hidden =c(nb.hid1,nb.hid2)
        
        WT.opt <- deca
        
      } else {}
    
    cat("> NNET MODEL run ",i,"sampling  ...", "\n",append = FALSE)
    
    net.tr <- nnet(x.in,y.in,data=data.in, decay = WT.opt, size = nb.hid1, 
                   linout = TRUE, maxit = it.max, Hess = TRUE) 
    
    # Save Model - R file
    file.name<-paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('net.tr'),file=file.name)
    
    ######################
    ## Variable importance
    ######################
    ## weight algorithm (Garson 1991)
    ##-------------------------------
    g  <- garson(net.tr, colnames(y.in), bar_plot = TRUE, x_lab = NULL,
                 y_lab = "Rel.importance", wts_only = FALSE)
  
  	    # Re-order data by row names
        order <- as.numeric(row.names(g$data))
        sort <- sort(order, index.return=TRUE)
        g$data[,1] <- g$data[sort$ix,1] 
        g$data[,2] <- g$data[sort$ix,2] 
  
  	varimp.nnet <- rbind(varimp.nnet,g$data[,1])
  
  	colnames(varimp.nnet) <- g$data[,2]
  
    ## MODEL PERFORMANCE
    ####################
    ## PREDICT THE OBS
  	pred.y <- predict(net.tr, newdata=x.in,type="raw",na.rm = TRUE)
  	
  	# PERFORMANCE OF THE MODEL
  #	lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
  	M.perf <- asses.lm.(pred.y, y.in[,1] )
  	
  	names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    
    ## CROSS - Validate model
    CV.data  <- CV.NNET.(data.in, nnet.unit= nb.hid1, nnet.WT= WT.opt, it.max= it.max, 
                         K=fold.cv, cv.lim = 10, name.sp = colnames(y.in))
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, height = 15,units="cm",res=150)
      
        ## Plot observed VS predicted for CV
        plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
      
    } else {}
    
    ## AIC and AICC of the Best model
    ################################
    RSS <- sum((y.in - predict(net.tr, newdata=x.in, type="raw"))^2) 
    
    aic.temp <- 2*sum(net.tr$wts!=0) - length(y.in)*log(RSS/length(y.in)) # AIC
    
    aicc.temp <- aic.temp + (2*sum(net.tr$wts!=0)+(sum(net.tr$wts!=0)+1))/
                               (length(y.in) -  sum(net.tr$wts!=0)-1) #AICc
    
    modl.perf.nnet <- rbind(modl.perf.nnet, c(M.perf,CV.perf, RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    
  } # end of loop repetition
  
  # write all model perf run
  write.table(modl.perf.nnet, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  #write table for all run
  write.table(varimp.nnet, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.nnet,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.nnet,2,mean) ) 
  
  # name technique update
  mod.name  <- c(mod.name, mod.tech)
  
  row.names(modl.perf.av) <- mod.name
  row.names(varimp.av) <- mod.name
  
} # end loop modeling technique
  
##########################################################################################################################################################
if (RSNNS == "NO") {}else{
  
  for (i in 1:n.rep)    
  { 
    mod.tech <- "RSNNS"
    
    outpath <- paste(output,"/",mod.tech,sep="")
    
    creat.subDir (output,mod.tech)
  
    ## Neural Networks (RSNNS) -- RSNNS package
    ##############################################
    ## OPTIMIZED PARAMETER PERFORMANCE OF RSNNS
    ##########################################
    
    if (i ==1) {  cat("> RSNNS Parameter two layer ...", "\n",append = FALSE)
      
      ##(from the FAQ for a commercial neural network software company)  
      ## Number of inputs + outputs) * (2/3) -- Heaton 2008
      nb.hid1 <- round (ncol(data.in) *2/3)
      
      nb.hid2 <- round (nb.hid1*2/3)
      
      nb.hidden =c(nb.hid1,nb.hid2)
      
    } else {}
    
    cat("> RSNNS MODEL run ",i,"sampling  ...", "\n",append = FALSE)
      
      rsnns.tr <- mlp(x.in,y.in, size= nb.hidden, maxit= it.max, learnFunc = "Std_Backpropagation")
      
    # Save Model - R file
    file.name <- paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('rsnns.tr'),file=file.name)
    
    ######################
    ## Variable importance
    ######################
    ## weight algorithm (Garson 1991)
    ##-------------------------------
    # compute a simple one layer neural network (var.importance not possible on two layer)
    rsnns.g <- mlp(x.in,y.in, size= nb.hidden[1], maxit= it.max, learnFunc = "Std_Backpropagation")
    
    g  <- garson(rsnns.g)
    
    # Re-order data by row names
    order <- as.numeric(row.names(g$data))
    sort <- sort(order, index.return=TRUE)
    g$data[,1] <- g$data[sort$ix,1] 
    g$data[,2] <- g$data[sort$ix,2] 
    
    varimp.rsnns <- rbind(varimp.rsnns,g$data[,1])
    
    colnames(varimp.rsnns) <- g$data[,2]
    
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict(rsnns.tr, x.in) 
    
    # PERFORMANCE OF THE MODEL
    #lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in[,1])
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    
    ## CROSS - Validate model
    CV.data  <- CV.RSNNS.(data.in, nral.hidden= nb.hidden , it.max = it.max, 
                           K=fold.cv, cv.lim = 10, name.sp = colnames(y.in))
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i ==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, height = 15,units="cm",res=150)
      
        ## Plot observed VS predicted for CV
        plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
    }else{}
    
    ## AIC and AICC of the Best model
    ################################
    RSS <- sum((y.in - predict(rsnns.tr, newdata=x.in, type="raw"))^2)
    
    aic.temp <- 2*sum( weightMatrix(rsnns.tr)!=0 ) - length(y.in)*log(RSS/length(y.in)) # AIC
    
    aicc.temp <- aic.temp + (2*sum( weightMatrix(rsnns.tr)!=0 )+(sum( weightMatrix(rsnns.tr)!=0 )+1))/
      (length(y.in) -  sum( weightMatrix(rsnns.tr)!=0 )-1) #AICc
    
    modl.perf.rsnns <- rbind(modl.perf.rsnns, c(M.perf, CV.perf, RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    
  } # end of loop repetition
  
  # write all model perf run
  write.table(modl.perf.rsnns, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  #write table for all run
  write.table(varimp.rsnns, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.rsnns,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.rsnns,2,mean) )
  
  # name technique update
  mod.name  <- c(mod.name, mod.tech)
  
  row.names(modl.perf.av) <- mod.name
  row.names(varimp.av) <- mod.name
  
} # end loop modeling technique 
  
##########################################################################################################################################################  
if (ELM == "NO") {}else{
    
    for (i in 1:n.rep)    
    {
    
    mod.tech <- "ELM"
    
    outpath <- paste(output,"/",mod.tech,sep="")
    
    creat.subDir (output,mod.tech)
    
    ######################
    ### NETWORK ANALYSIS
    #####################                       
    ## OPTIMIZED PARAMETER PERFORMANCE OF ELM
    ##########################################
    
    if (i ==1) {  cat("> ELM Parameter one layer ...", "\n",append = FALSE)
      
      ##(from the FAQ for a commercial neural network software company)  
      ## Number of inputs + outputs) * (2/3) -- Heaton 2008
      nb.hid1 <- round (ncol(data.in) *2/3)
      
      nb.hid2 <- round (( nb.hid1+ncol(y.in))*2/3 )
      
      nb.hidden =c(nb.hid1,nb.hid2)
      
      WT.opt <- 0.01
      
    } else {}
    
    cat("> ELM MODEL run ",i,"sampling  ...", "\n",append = FALSE)
    
    if (i ==1) {  cat("> ELM choose algorithm ...", "\n",append = FALSE)
    
              algo <- c("sig","sin","radbas","hardlim","hardlims","satlins",
                      "tansig","tansig","tribas","poslin","purelin")
    
              for (s in 1:length(algo)) # test the best algorythm
                  {
                    # train the model
                    elm.tr <- elmtrain(x=x.in, y=y.in,nhid=nb.hid1, actfun=algo[s] )
                    
                    # PREDICT THE OBS
                    pred.y <- predict(elm.tr, newdata=as.data.frame(data.in)[,2:ncol(data.in)],type="raw")
                    
                    # PERFORMANCE OF THE MODEL
                    lin.corr <- asses.lm.(pred.y,y.in[,1]) 
                    
                    if (s==1){rse <-lin.corr[6] } else { rse <- c(rse,lin.corr[6]) }
                    
                  }
            # min rse = best algorithm
            b.algo <-  which(rse == min(rse), arr.ind = TRUE) 
            
            cat("> ", algo[b.algo], " is selected ...", "\n",append = FALSE)
            } else{}
   
    # train elm
    elm.tr <- elmtrain(x=x.in, y=y.in,nhid=nb.hid1, actfun=algo[b.algo] )
    
    # Save Model - R file
    file.name<-paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('elm.tr'),file=file.name)
    
    ######################
    ## Variable importance
    ######################
    varimp <- var.imp.rf(elm.tr,x.in,var.list,nperm=n.perm)
    
    varimp.elm <- rbind(varimp.elm,(varimp/sum(varimp)) )
    
    ## MODEL PERFORMANCE
    ####################
    ####################
    # PREDICT THE OBS
    pred.y <- predict( elm.tr, newdata=as.data.frame(data.in)[,2:ncol(data.in)],type="raw" ) 
    
    # PERFORMANCE OF THE MODEL
    #lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in[,1])
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    
    ## CROSS - Validate model
    CV.data <- CV.ELM. (data.in, elm.tr, K=fold.cv, cv.lim = 10, name.sp = colnames(y.in) )
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, 
          height = 15,units="cm",res=150)
      
        ## Plot observed VS predicted for CV
        plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
      
    } else {}
    
    ## AIC and AICC of the Best model
    ################################
    RSS <- sum((y.in - predict(elm.tr, newdata=x.in, type="raw"))^2) 
    
    aic.temp <- 2*sum(elm.tr$outweight!=0) - length(y.in)*log(RSS/length(y.in)) # AIC
    
    aicc.temp <- aic.temp + (2*sum(elm.tr$outweight!=0)+(sum(elm.tr$outweight!=0)+1))/
      (length(y.in) -  sum(elm.tr$outweight!=0)-1) #AICc
    
    modl.perf.elm <- rbind(modl.perf.elm, c(M.perf, CV.perf, RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    
  } # end of loop repetition

  # write all model perf run
  write.table(modl.perf.elm, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

  #write table for all run
  write.table(varimp.elm, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.elm,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.elm,2,mean) )
  
  # name technique update
  mod.name  <- c(mod.name, mod.tech)
  
  row.names(modl.perf.av) <- mod.name
  row.names(varimp.av) <- mod.name
  
} # end loop modeling technique
  
#####################################################################################################################################
if (RF == "NO") {}else{
  
  for (i in 1:n.rep)    
  {
    mod.tech <- "RF"
    
    outpath <- paste(output,"/",mod.tech,sep="")
    
    creat.subDir (output,mod.tech)  
    
    #####################
    ## Random Forest
    #####################
    cat("> RANDOM FOREST MODEL run ",i,"sampling  ...", "\n",append = FALSE)
    
    # First we tune the parameter
    if (i==1) { rfpa <- tuneRF (x = x.in,y = t(y.in),data = data.in, ntree = nb.tree, na.action=na.omit)
          }else{}
    
    if (overfit=="YES") {nb.hid1 <- round (ncol(data.in) *2/3)
                  rf <- randomForest(eval(parse(text= paste(colnames(y.in),"~.",sep=""))), data =data.in, 
                       ntree = nb.tree, mtry = rfpa,importance=TRUE, na.action=na.omit, maxnodes=nb.hid1)
        } else { 
                  rf <- randomForest(eval(parse(text= paste(colnames(y.in),"~.",sep=""))), data =data.in, 
                      ntree = nb.tree, mtry = rfpa,importance=TRUE, na.action=na.omit )      
        }
    # Save Model - R file
    file.name <- paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('rf'),file=file.name) 
    
    ######################
    ## Variable importance type 2 mean decrease in node impurity
    ######################
    varimp <- rf$importance[,2]
    
    varimp.rf <- rbind(varimp.rf,(varimp/sum(varimp)) ) 
    
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict( rf, newdata=x.in ) 
    
    # PERFORMANCE OF THE MODEL
    #lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in[,1]) ; # print (M.perf)
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    
    ###CROSS - Validate model
    CV.data  <- CV.RF.(rf, data.in,  K=fold.cv, cv.lim = 10, name.sp = colnames(y.in))  
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i ==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, 
          height = 15,units="cm",res=150)
      
        ## Plot observed VS predicted for CV
        plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
    }else{} 
    
    ## AIC and AICC of the Best model
    ################################
    RSS <- sum((y.in - predict(rf, newdata=x.in))^2) 
    
    aic.temp <- 2*mean(treesize(rf)) - length(y.in)*log(RSS/length(y.in)) # AIC
    
    aicc.temp <- aic.temp + (2*mean(treesize(rf))+(mean(treesize(rf))+1))/
      (length(y.in) -  mean(treesize(rf))-1) #AICc
    
    modl.perf.rf <- rbind(modl.perf.rf, c(M.perf, CV.perf, RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    
  } # end of loop repetition
  
  # write all model perf run
  write.table(modl.perf.rf, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  #write table for all run
  write.table(varimp.rf, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
      
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.rf,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.rf,2,mean) )
  
  # name technique update
  mod.name  <- c(mod.name, mod.tech)
  
  row.names(modl.perf.av) <- mod.name
  row.names(varimp.av) <- mod.name
  
} # end loop modeling technique

#####################################################################################################################################
if (CSRF == "NO") {}else{
  
  for (i in 1:n.rep)    
  {
    mod.tech <- "CSRF"
    
    outpath <- paste(output,"/",mod.tech,sep="")
    
    creat.subDir (output,mod.tech)  
    
    ##############################
    ## Case specific Random forest
    ##############################
    cat("> Case specific RANDOM FOREST MODEL run ",i,"sampling  ...", "\n",append = FALSE)
    
    # First we tune the parameter
    if (i==1) { rfpa <- tuneRF (x = x.in,y = t(y.in),data = data.in, ntree = nb.tree, na.action=na.omit)
                rfpa <- rfpa[match( max(rfpa[,2]), rfpa[,2] ),1]
          }else{}
    
    if (overfit=="YES") {
        nb.hid1 <- round (ncol(data.in) *2/3)
    
        csrf <- ranger(eval(parse(text= paste(colnames(y.in),"~.",sep=""))), data =as.data.frame(data.in), 
                   num.trees = nb.tree,mtry= rfpa, importance="impurity", 
                   min.node.size = nb.hid1, splitrule = "maxstat", alpha = alpha.lim)
      }else{
        csrf <- ranger(eval(parse(text= paste(colnames(y.in),"~.",sep=""))), data =as.data.frame(data.in), 
                     num.trees = nb.tree,mtry= rfpa, importance="impurity")
      }
    
    # Save Model - R file
    file.name <- paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('csrf'),file=file.name) 
    
    ######################
    ## Variable importance
    ######################
    varimp <- csrf$variable.importance
    
    varimp.csrf <- rbind(varimp.csrf,(varimp/sum(varimp)) ) 
    
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict(csrf, dat=x.in) 
    
    # PERFORMANCE OF THE MODEL
    #lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y$predictions, y.in[,1])
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    
    ###CROSS - Validate model
    CV.data  <- CV.CSRF.(csrf, as.data.frame(data.in),  K=fold.cv, cv.lim = 10, name.sp = colnames(y.in))  
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i ==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, 
          height = 15,units="cm",res=150)
      
      ## Plot observed VS predicted for CV
      plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
    }else{} 
    
    ## AIC and AICC of the Best model
    ################################
    RSS <- sum((y.in - pred.y$predictions )^2) 
    
    aic.temp <- NA # AIC
    
    aicc.temp <- NA #AICc
    
    modl.perf.csrf <- rbind(modl.perf.csrf, c(M.perf, CV.perf, RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    
  } # end of loop repetition
  
  # write all model perf run
  write.table(modl.perf.csrf, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  #write table for all run
  write.table(varimp.csrf, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.csrf,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.csrf,2,mean) )
  
  # name technique update
  mod.name  <- c(mod.name, mod.tech)
  
  row.names(modl.perf.av) <- mod.name
  row.names(varimp.av) <- mod.name
  
} # end loop modeling technique
  
#####################################################################################################################################
if (BAG == "NO") {}else{
    
  for (i in 1:n.rep)    
  {
    mod.tech <- "BAG"
    
    outpath <- paste(output,"/",mod.tech,sep="")
    
    creat.subDir (output,mod.tech)
    
    ######################
    ### BAGGING ANALYSIS -- adabag
    #####################                       
    ## OPTIMIZED PARAMETER PERFORMANCE OF BAG
    ##########################################
    cat("> BAG MODEL run ",i,"sampling  ...", "\n",append = FALSE)
    
    bag.tr <- ipredbagg (y.in[,1], x.in, nbag=nb.tree)
    
    # Save Model - R file
    file.name<-paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('bag.tr'),file=file.name)
    
    ######################
    ## Variable importance --> node impurity (similar as permutation)
    ######################
    for (n in 1: nb.tree)
      { 
        varimp.t <- bag.tr$mtrees[[n]]$btree$variable.importance # extract the var. imp
        sort <- sort(names(varimp.t), index.return=TRUE) # sort by alphabetic order
        varimp.t <- varimp.t[sort$ix]
        varimp.t <- varimp.t[names(varimp.t)!=colnames(y.in)] # remove the var. dependent
       
         if (n==1){ varimp <- varimp.t} else { varimp <- rbind( varimp, varimp.t )}
      }
    
    varimp.bag <- rbind( varimp.bag, apply(varimp, 2, sum)/ sum(varimp) )
    
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict(bag.tr, newdata=data.frame(data.in), type="raw",na.rm = TRUE)
    
    # PERFORMANCE OF THE MODEL
    #lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in[,1])
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    
    ## CROSS - Validate model
    CV.data  <- CV.BAG. (data.in, K= fold.cv, cv.lim = 10, nb.tree = nb.tree, name.sp = colnames(y.in))
      
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, height = 15,units="cm",res=150)
      
        ## Plot observed VS predicted for CV
        plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
      
    } else {}
    
    ## AIC and AICC of the Best model
    ################################
    RSS <- sum((y.in - predict(bag.tr, newdata= data.in, type = "raw"))^2) 
    
    for (n in 1: nb.tree)
        { 
          if (n==1){ nb.weight <- sum(bag.tr$mtrees[[n]]$btree$frame$wt!=0) }
            else{ nb.weight <- c( nb.weight, sum(bag.tr$mtrees[[n]]$btree$frame$wt!=0) )}
        }
    nb.weight <- mean(nb.weight)
    
    aic.temp <- 2*nb.weight - length(y.in)*log(RSS/length(y.in)) # AIC
    
    aicc.temp <- aic.temp + (2*nb.weight +(nb.weight+1) )/
      (length(y.in) -  nb.weight-1) #AICc
    
    modl.perf.bag <- rbind(modl.perf.bag, c(M.perf, CV.perf, RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    
  } # end of loop repetition
  
  # write all model perf run
  write.table(modl.perf.bag, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  #write table for all run
  write.table(varimp.bag, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.bag,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.bag,2,mean) )
  
  # name technique update
  mod.name  <- c(mod.name, mod.tech)
  
  row.names(modl.perf.av) <- mod.name
  row.names(varimp.av) <- mod.name
     
} # end loop modeling technique

#####################################################################################################################################
if (EGB == "NO") {}else{
  
  for (i in 1:n.rep)    
  {
    mod.tech <- "EGB"
    
    outpath <- paste(output,"/",mod.tech,sep="")
    
    creat.subDir (output,mod.tech)  
    
    #####################
    ## EGB -- START --> get info on http://xgboost.readthedocs.io/en/latest/get_started
    #####################
    cat("> Extreme Gradient Boosting run ",i,"sampling  ...", "\n",append = FALSE)
    
    # set the parameter of the booster  (default)
    para <- list(booster= "gbtree", max_depth = 6, eta = 0.3, silent = 1, nthread = 2,
                 objective = "reg:linear", eval_metric = "rmse")
    
    # run model
    egb.model <- xgboost(data = as.matrix(x.in), label = t(as.vector(y.in)), missing = NA, weight = NULL,
                params = para , nrounds = it.max, verbose = 0 )
            
    # Save Model - R file
    file.name <- paste(outpath,"/",mod.tech,"_model_",i,sep='')
    xgb.save(egb.model, file.name) 
    
    ######################
    ## Variable importance
    ######################
    varimp <- xgb.importance(colnames(x.in), model = egb.model)
    
    # weird but need the plot to calculate the importance overwise don't do
    xgb.plot.importance (varimp, rel_to_first = TRUE, xlab = "Relative importance")
    
    #re-organized to have alphabetic order
    order <- varimp$Feature
    sort <- sort(order, index.return=TRUE)
    varimp <- varimp$Importance[sort$ix]
    
    varimp.egb <- rbind(varimp.egb, varimp ) 
    
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict(egb.model, as.matrix(x.in)) 
    
    # PERFORMANCE OF THE MODEL
    #lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in[,1])
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    
    ###CROSS - Validate model
    CV.data  <- CV.EGB.(data.in,  K=fold.cv, cv.lim = 10, it.max= it.max, name.sp = colnames(y.in))  
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i ==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, 
          height = 15,units="cm",res=150)
      
      ## Plot observed VS predicted for CV
      plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
    }else{} 
    
    ## AIC and AICC of the Best model
    ################################
    RSS <- sum((y.in - pred.y)^2) 
    
    aic.temp <- NA # AIC
    
    aicc.temp <- NA #AICc
    
    modl.perf.egb <- rbind(modl.perf.egb, c(M.perf, CV.perf, RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    
  } # end of loop repetition
  
  # write all model perf run
  write.table(modl.perf.egb, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  #write table for all run
  write.table(varimp.egb, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.egb,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.egb,2,mean) )
  
  # name technique update
  mod.name  <- c(mod.name, mod.tech)
  
  row.names(modl.perf.av) <- mod.name
  row.names(varimp.av) <- mod.name
  
} # end loop modeling technique

#####################################################################################################################################
if (GBM == "NO") {}else{
  
  for (i in 1:n.rep)    
  {
    mod.tech <- "GBM"
    
    outpath <- paste(output,"/",mod.tech,sep="")
    
    creat.subDir (output,mod.tech)
    
    ## OPTIMIZED PARAMETER PERFORMANCE OF GBM
    ##########################################
    
    form <- as.formula(paste(colnames(y.in), "~", paste(var.list[!var.list %in% colnames(y.in)],
                                                        collapse = " + ")))
    
    cat("> GBM MODEL run ",i,"sampling  ...", "\n",append = FALSE)
    
    gbm.tr <- gbm(form, data=as.data.frame(data.in), n.trees=nb.tree, distribution = "gaussian" )
    
    # Save Model - R file
    file.name<-paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('gbm.tr'),file=file.name)
    
    ######################
    ## Variable importance -- permutation method
    ######################
    gbm.imp <- relative.influence(gbm.tr, n.trees= nb.tree, scale.=FALSE, sort.=FALSE)

    varimp.gbm <- rbind(varimp.gbm,(gbm.imp /sum(gbm.imp)) ) # relative importance to 1
    
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict.gbm(gbm.tr, newdata=as.data.frame(x.in), 
                        n.trees=nb.tree, type="response",na.rm = TRUE)
    
    # PERFORMANCE OF THE MODEL
    #lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in[,1])
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    
    ## CROSS - Validate model
    CV.data  <- CV.GBM.(gbm.tr, data.in, nb.tree= nb.tree, K=fold.cv, cv.lim = 10, name.sp = colnames(y.in))
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, height = 15,units="cm",res=150)
      
      ## Plot observed VS predicted for CV
      plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
      
    } else {}
    
    ## AIC and AICC of the Best model
    ################################
    RSS <- sum((y.in - predict.gbm(gbm.tr, newdata=as.data.frame(x.in), 
                                   n.trees=nb.tree, type="response",na.rm = TRUE))^2) 
    
    for (n in 1: nb.tree)
      { 
        if (n==1){ nb.weight <- sum(unlist(gbm.tr$trees[[n]])!=0) }
        else{ nb.weight <- c( nb.weight, sum(unlist(gbm.tr$trees[[n]])!=0) )}
      }
    nb.weight <- mean(nb.weight)
   
    aic.temp <- 2*nb.weight - length(y.in)*log(RSS/length(y.in)) # AIC
    
    aicc.temp <- aic.temp + (2*nb.weight +(nb.weight+1) )/
      (length(y.in) -  nb.weight-1) #AICc
    
    modl.perf.gbm <- rbind(modl.perf.gbm, c(M.perf, CV.perf, RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    
  } # end of loop repetition

  # write all model perf run
  write.table(modl.perf.gbm, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  #write table for all run
  write.table(varimp.gbm, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.gbm,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.gbm,2,mean) )
  
  # name technique update
  mod.name  <- c(mod.name, mod.tech)
  
  row.names(modl.perf.av) <- mod.name
  row.names(varimp.av) <- mod.name
  
} # end loop modeling technique
  
##################################################################################################################################### 
if (GLM == "NO") {}else{
    
    for (i in 1:n.rep)    
    {
    mod.tech <- "GLM"
    
    outpath <- paste(output,"/",mod.tech,sep="")
    
    creat.subDir (output,mod.tech) 
  
	  cat("> GLM MODEL run ",i,"sampling  ...", "\n",append = FALSE)

    #####################
    ## GLM - Multivariate
    #####################
    data.cal <- data.frame(data.in)    
    
    names.pred <- var.list
    
    name.sp <- colnames(y.in)
    
    for (p in 1:length(names.pred)) 
        { 
          if (p == 1) 
          {
            fmula.pred.glm <- as.vector(paste("poly(",names.pred[p],",2)",sep=""),mode="character")  
          }
          else  
          {
            fmula.pred.glm <- paste(fmula.pred.glm," + poly(",names.pred[p],",2)",sep="") 
          }
        }		 	
    
    ###########
    # GLM FIT #
    ###########
    glm.tmp <- glm(eval(parse(text = paste(paste(name.sp), "~",fmula.pred.glm, collapse = ""))),
                             data=data.cal,family=gaussian,maxit = 100)    
  
    # Save Model - R file
    file.name <- paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('glm.tmp'),file=file.name) 
    
    ######################
    ## Variable importance
    ######################
    varimp <- var.imp.glm.step(glm.tmp,data.cal,names.pred,nperm = n.perm)
    
    varimp.glm <- rbind(varimp.glm, (varimp /sum(varimp)) )
  
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict(glm.tmp, newdata = data.cal, type="response",na.rm = TRUE)
    
    # PERFORMANCE OF THE MODEL
    #lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in[,1])
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    
    ## CROSS - Validate model
    CV.data  <- CV.glm.(glm.tmp, K=fold.cv, cv.lim = 10, name.sp = colnames(y.in))  
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i ==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, 
          height = 15,units="cm",res=150)
      
        ## Plot observed VS predicted for CV
        plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
    }else{}
    
    ## AIC and AICC of the Best model
    ################################
    ## AIC and AICC of the Best model
    ################################
    RSS <- sum((y.in - predict(glm.tmp, newdata = data.cal, type="response",na.rm = TRUE))^2) 
    
    aic.temp <- AIC(glm.tmp) # AIC
    aicc.temp <- AICc(glm.tmp) #AICc
    
    modl.perf.glm <- rbind(modl.perf.glm, c(M.perf, CV.perf, RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    
    } # end of loop repetition
  
  # write all model perf run
  write.table(modl.perf.glm, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  #write table for all run
  write.table(varimp.glm, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.glm,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.glm,2,mean) )
  
  # name technique update
  mod.name  <- c(mod.name, mod.tech)
  
  row.names(modl.perf.av) <- mod.name
  row.names(varimp.av) <- mod.name
  
} # end loop modeling technique  
  
##################################################################################################################################### 
if (GLM_sw == "NO") {}else{
  
    for (i in 1:n.rep)    
    {
    mod.tech <- "GLM_sw"
    
    outpath <- paste(output,"/",mod.tech,sep="")
    
    creat.subDir (output,mod.tech)     
    
    ###############
    ## GLM -STEP ##
    ###############
    data.cal <- data.frame(data.in)    
    
    names.pred <- var.list
    
    name.sp <- colnames(y.in)
    
    glm.tmp.step <- step(glm(eval(parse(text = paste(paste(name.sp), "~",fmula.pred.glm, collapse = ""))),
                             data=data.cal,family=gaussian,maxit = 100),trace = F,direction="both")
    # Save Model - R file
    file.name <- paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('glm.tmp.step'),file=file.name) 
    
    ######################
    ## Variable importance
    ######################
    varimp <- var.imp.glm.step(glm.tmp.step,data.cal,names.pred,nperm = n.perm)
    
    varimp.glmstep <- rbind(varimp.glmstep, (varimp/sum(varimp)) )
    
    ####################
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict(glm.tmp.step, newdata = data.cal, type="response",na.rm = TRUE)
    
    # PERFORMANCE OF THE MODEL
    #lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in[,1])
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    
    ## CROSS - Validate model
    CV.data  <- CV.glm.step(glm.tmp.step, K=fold.cv, cv.lim = 10, name.sp = colnames(y.in)) 
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i ==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, height = 15,units="cm",res=150)
      
          ## Plot observed VS predicted for CV
          plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
    }else{}  
    
    ## AIC and AICC of the Best model
    ################################
    RSS <- sum((y.in - predict(glm.tmp.step, newdata = data.cal, type="response",na.rm = TRUE))^2) 
    
    aic.temp <- AIC(glm.tmp.step) # AIC
    aicc.temp <- AICc(glm.tmp.step) #AICc
    
    modl.perf.glmstep <- rbind(modl.perf.glmstep, c(M.perf, CV.perf, RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    
  } # end of loop repetition

# write all model perf run
write.table(modl.perf.glmstep, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

#write table for all run
write.table(varimp.glmstep, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

# calculate avearge values for model perf and variable importance
modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.glmstep,2,mean) )
varimp.av <- rbind( varimp.av, apply(varimp.glmstep,2,mean) )

# name technique update
mod.name  <- c(mod.name, mod.tech)

row.names(modl.perf.av) <- mod.name
row.names(varimp.av) <- mod.name

} # end loop modeling technique  
  
##################################################################################################################################### 
if (GAM == "NO") {}else{
    
  for (i in 1:n.rep)    
  {
    mod.tech <- "GAM"
    
    outpath <- paste(output,"/",mod.tech,sep="")
    
    creat.subDir (output,mod.tech) 
    
    cat("> GAM MODEL run ",i,"sampling  ...", "\n",append = FALSE)
    
    #####################
    ## GAM - Multivariate
    #####################
    data.cal <- data.frame(data.in)    
    
    names.pred <- var.list
    
    name.sp <- colnames(y.in)
    
    for (p in 1:length(names.pred)) 
      { 
        if (p == 1) 
        {
          fmula.pred.glm <- as.vector(paste("poly(",names.pred[p],",2)",sep=""),mode="character")  
        }
        else  
        {
          fmula.pred.glm <- paste(fmula.pred.glm," + poly(",names.pred[p],",2)",sep="") 
        }
      }		 	
    
    ###########
    # GAM FIT #
    ###########
    gam.tmp <- gam(eval(parse(text = paste(paste(name.sp), "~",fmula.pred.glm, collapse = ""))),
                   data=data.cal,family=gaussian,maxit = 100)    
    
    # Save Model - R file
    file.name <- paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('gam.tmp'),file=file.name) 
    
    ######################
    ## Variable importance
    ######################
    varimp <- var.imp.glm.step(gam.tmp,data.cal,names.pred,nperm = n.perm)
      
    varimp.gam <- rbind(varimp.gam, (varimp/sum(varimp)) )
    
    ## MODEL PERFORMANCE
    ####################
    ####################
    # PREDICT THE OBS
    pred.y <- predict(gam.tmp, newdata = data.cal, type="response",na.rm = TRUE)
    
    # PERFORMANCE OF THE MODEL
    #lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in[,1])
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    
    ## CROSS - Validate model
    CV.data  <- CV.glm.(gam.tmp, K=fold.cv, cv.lim = 10, name.sp = colnames(y.in))  
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i ==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, 
          height = 15,units="cm",res=150)
      
        ## Plot observed VS predicted for CV
        plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
    }else{}
    
    ## AIC and AICC of the Best model
    ################################
    RSS <- sum((y.in - predict(gam.tmp, newdata = data.cal, type="response",na.rm = TRUE))^2) 
    
    aic.temp <- AIC(gam.tmp) # AIC
    aicc.temp <- AICc(gam.tmp) #AICc
    
    modl.perf.gam <- rbind(modl.perf.gam, c(M.perf, CV.perf, RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    
  } # end of loop repetition

# write all model perf run
write.table(modl.perf.gam, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

#write table for all run
write.table(varimp.gam, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

# calculate avearge values for model perf and variable importance
modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.gam,2,mean) )
varimp.av <- rbind( varimp.av, apply(varimp.gam,2,mean) )

# name technique update
mod.name  <- c(mod.name, mod.tech)

row.names(modl.perf.av) <- mod.name
row.names(varimp.av) <- mod.name

} # end loop modeling technique  
  
#####################################################################################################################################
##################################################################################################################################### 
 ## Write average values
# define folder
outpath <- paste(output,"/average",sep="")

# check if folder exist
creat.subDir (output,"average")

#write table for average perf and var imp
write.table(modl.perf.av, file=paste(outpath,"/model_perf_mean.txt",sep=""),sep="\t", 
            append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

write.table(varimp.av, file=paste(outpath,"/var_imp_mean.txt",sep=""),sep="\t", 
            append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

} # end loop diff dep variable
 
endCluster() # END OF MULTICORE CALCULATION

proc.time() - ptm # check time

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#### END SCRIPT ---- COFFE TIMES ----
#################################################################################################################################
#################################################################################################################################




