library(R6)
library(data.tree)
library(survival)
library(dplyr)


UnivariateSurvivalTree = R6Class(
  
  
  classname = "UnivariateSurvivalTree",
  
  public = list(
    
    # MAIN DATA
    traindata = NULL,
    originalrownames = NULL,
    covariates = NULL,
    time = NULL,
    censor = NULL,
    
    # TREE DATA
    nsize = NULL, # terminal node size threshold
    
    # INITIALIZER
    initialize = function(
    traindata, time, censor, covariates = c(), nsize = 10
    ){
      if (!is.data.frame(traindata) || nrow(traindata) == 0) {
        stop("\"traindata\" variable must be a nonempty data frame")
      }
      # Rename rows of "traindata" to a sequence of consecutive integers
      # Store the original row names in case they are required later
      self$originalrownames = rownames(traindata)
      rownames(traindata) = seq(1:nrow(traindata))
      # Make sure the survival and covariate variable names are actually valid
      for (name in c(time, censor, covariates)) {
        if (!(name %in% names(traindata))) {
          stop(
            paste(name, "is not a variable in the data", sep = " ")
          )
        }
      }
      # Bulk initiation of various variables
      self$traindata = traindata
      self$covariates = covariates
      self$time = time
      self$censor = censor
      self$nsize = nsize
    },
    #Find the best split
    bestsplit = function(data){
      ncovas=length(self$covariates)
      opt.split.bycov = data.frame(covariate=NA,bestsplitval=NA,lrstat=NA)
      for (i in 1:ncovas) {
        if (length(unique(data[self$covariates[i]][[1]]))==1)
        {opt.split.bycov[i,1]=self$covariates[i]
        opt.split.bycov[i,2]=0
        opt.split.bycov[i,3]=0}
        else{
          data.sort<- data %>% arrange(data[self$covariates[i]])
          n = nrow(data)
          lrstat.old=0
          k=1
          while (k<=n-1){
            while ((data.sort[self$covariates[i]][k,]==data.sort[self$covariates[i]][k+1,]) & k<n){
              k=k+1
            }
            if(k<n){
              ind=c(rep(1,k),rep(0,n-k))
              Y<-Surv(data.sort[self$time][[1]],data.sort[self$censor][[1]]==1)
              lrstat=survdiff(Y~ind)[[5]]
              if (lrstat>lrstat.old){
                lrstat.max=lrstat
                split.pos=k 
              }
              lrstat.old=lrstat
              k=k+1
            }
          }
          opt.split.bycov[i,1]=self$covariates[i]
          opt.split.bycov[i,2]=data.sort[self$covariates[i]][split.pos,1]
          opt.split.bycov[i,3]=lrstat.max
        }
      }
      logr=max(opt.split.bycov[,3])
      max=which.max(opt.split.bycov[,3])
      split=opt.split.bycov[max,]
      v=c(-split[[2]],rep(0,length(self$covariates)))
      v[max+1]=1
      mylist=list(v,logr)
      return(mylist)
    },
    #Add children function
    createtreer = function(nodeb){
      n = nrow(nodeb$data)
      if (n > nsize){
        X=nodeb$data
        subsetX=strtoi(rownames(X))
        newdata=self$traindata[subsetX,]
        bsplit <- self$bestsplit(newdata)
        v = bsplit[[1]]
        nodeb$optv = v
        lrnodes <-self$lrnode.make(nodeb$data,v)
        nl = nrow(lrnodes[[1]])
        nr = nrow(lrnodes[[2]])
        if (nl>0 && nr>0){
          splits <- lrnodes[[3]]
          lrstat <- self$lr.compute(nodeb$data,splits)
          nodeb$lrstat <- lrstat
          lcnodename = paste0(nodeb$name,"l"); rcnodename = paste0(nodeb$name,"r")
          nodeb$AddChild(lcnodename,data=lrnodes[[1]]) 
          nodeb$AddChild(rcnodename,data=lrnodes[[2]])
          self$createtreer(nodeb$children[[1]])
          self$createtreer(nodeb$children[[2]])
        } else {
          nodeb$KMestimate <- self$KM.compute(nodeb$data)
        }
      } else {
        nodeb$KMestimate <- self$KM.compute(nodeb$data)
      }
    },
    
    #Make left and right children function 
    lrnode.make = function(X,v){
      Z=as.matrix(cbind(1,X))
      splits=Z%*%v
      XL=X[splits<=0,]; XR=X[splits>0,]
      return(list(XL,XR,splits))
    },
    
    #Compute log rank statistic function 
    lr.compute = function(X,splits,data1=self$traindata,time1=self$time,censor1=self$censor){
      subsetX = strtoi(rownames(X))
      dataflr=data1[subsetX,]
      Y<-Surv(dataflr[time1][[1]],dataflr[censor1][[1]]==1)
      lrstat=survdiff(Y~splits<=0)[[5]]
      return(lrstat)
    },
    #Compute the KM censoring estimate
    KM.compute = function(X){
      subsetX = strtoi(rownames(X))
      datasubset = self$traindata[subsetX,]
      Y<-Surv(datasubset[self$time][[1]],datasubset[self$censor][[1]]==1)
      kmfit=survfit(Y~1)
      YC<-Surv(datasubset[self$time][[1]],datasubset[self$censor][[1]]==0)
      kmfitc=survfit(YC~1)
      return(list(kmfit,kmfitc)) 
    },
    createtree = function(subset) {
      subX = self$traindata[subset, self$covariates]
      subtree = Node$new("Node0", data = subX)
      self$createtreer(subtree)
      #subtree$Do(function(node) node$KMest<-self$KM.compute(node), filterFun = isLeaf)
      return(subtree)
    }
  ) 
) 
      
      
      
    