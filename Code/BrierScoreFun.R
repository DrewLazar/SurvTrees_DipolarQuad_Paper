#install.packages("SurvMetrics")
library(data.tree)
library(survival)
library(SurvMetrics)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Load the Remission data
load("../../../Data/Remission.rda")

alldata=Remission
#dipolartreep <- kfoldPruning(dipolartree,Remission,3)


#Test of Brier Score with no censoring
#survt <- c(1,2,3,4)
#status <- c(1,1,1,1)
#df <- data.frame(survt, status)
#Survobj = Surv(df$survt,df$status==1)
#pre_sp = c(.1,.25,.80,.95)
#Brier(Survobj,pre_sp,2.5)
#should have?
#(1/4)*(.1^2 + .25^2 + .2^2 + .05^2)

#Test of Integrated Brier Score with no censoring
#survt <- c(1,2,3,4)
#status <- c(1,1,1,1)
#object = Surv(df$survt,df$status==1)
#IBSrange = c(2.5, 3.5)
#sp_matrix = cbind(c(.1,.25,.80,.95),c(.07,.25,.76,.90))
#IBS(object, sp_matrix, IBSrange)

Brierscore = function(treetoassess,testdata,IBSrange,time,censor){

#2
sendtesttree <- function(treetoaddtest,testX){
  splittestatnode <- function(nodein){
    v = nodein$optv
    Z = as.matrix(cbind(1, nodein$testX))
    splits = Z %*% v
    nodein$children[[1]]$testX = nodein$testX[splits < 0, ]
    nodein$children[[2]]$testX = nodein$testX[splits >= 0, ]
  }
  treetoaddtest$testX<-testX
  treetoaddtest$Do(function(node) splittestatnode(node), filterFun = isNotLeaf)
}

#3
timesandspreds <- function(node,IBSR){
  times = strtoi(rownames(node$data))
  surpreds  = summary(node$KMest[[1]],times=IBSR,extend=TRUE)
  return(list(times,surpreds[[6]]))
}
treetoassess$Do(function(node) node$timespreds<-timesandspreds(node,IBSrange), filterFun = isLeaf)
a = treetoassess$Get('timespreds', filterFun = isLeaf)
cola = ncol(a)
times = a[[1,1]]
for (j in 1:(cola-1)){
  times = c(times,a[[1,j+1]])
}
times = matrix(times,nrow=nrow(testdata),byrow=TRUE)
predsc = rep(a[[2,1]],length(a[[1,1]]))
preds = matrix(predsc, nrow = length(a[[1,1]]), byrow = TRUE)
predsm= preds
for (i in 2:cola){
  predsc = rep(a[[2,i]],length(a[[1,i]]))
preds = matrix(predsc, nrow = length(a[[1,i]]), byrow = TRUE)
predsm = rbind(predsm,preds)
}
sp_matrix = cbind(times,predsm)
sp_matrix = sp_matrix[order(sp_matrix[,1],decreasing=FALSE),]
sp_matrix = sp_matrix[,-1]

object<-Surv(testdata[time][[1]],testdata[censor][[1]]==1)
IBS(object, sp_matrix, IBSrange)
}

#dipolartreep=kfoldPruning(dipolartree,Remission,3)
IBSrange = c(2, 5,8, 12, 16, 26)
BS=Brierscore(dipolartree,Remission,IBSrange,'survt','status')
