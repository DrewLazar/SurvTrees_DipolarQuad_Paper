
bootstrapPruning <- function(intree,dipolarmodel,bts = c(),splitrule) {

#bts[1] = number of bootstrap samples, bts[2] = percentage size of bootstrap relative to sample
#splitrule==0 means use Z %*% v < 0 for splitting rule, splitrule==1 means use 
# Z %*% v <= 0 for splitting rule. 

#I) Get list of pruned trees from all data (T1 > T2 > ... > Tm)  tree and list of pruning stats (alpha_1 < ... < alpha_m)

#Add pruning criteria numbers function 
pruningstat<- function(node) {
  GTh<-sum(node$Get("lrstat",filterFun=isNotLeaf))
  Sh <- node$totalCount - node$leafCount
  gh <- GTh/Sh
}
alldata=dipolarmodel$traindata 
#function for pruning all children from node 
prunenodesfun = function(opttree, prunenode){
  prunenodel = paste0(prunenode,'l')
  Prune(opttree, pruneFun = function(x) x$name!=prunenodel)
  prunenoder = paste0(prunenode,'r')
  Prune(opttree, pruneFun = function(x) x$name!=prunenoder)
}


#function for getting list of pruned trees and pruning stats while pruning off 
#sequence of weakest links 
lrpruningfun = function(treetoprune) {
  pstats = c()
  opttree = list()
  opttree[[1]] = Clone(treetoprune)
  opttree[[1]]$Do(function(node) node$prnstat<- pruningstat(node), filterFun = isNotLeaf)  
  pstats[1] = min(opttree[[1]]$Get("prnstat", filterFun = isNotLeaf)) 
  i=0
  while (opttree[[i+1]]$totalCount != 1){
    i = i+1 
    prunenode=names(which.min(opttree[[i]]$Get("prnstat", filterFun = isNotLeaf)))
    prunenodesfun(opttree[[i]],prunenode)
    opttree[[i+1]] = Clone(opttree[[i]])
    opttree[[i+1]]$Do(function(node) node$prnstat<- pruningstat(node), filterFun = isNotLeaf) 
    if (opttree[[i]]$totalCount != 1){
      pstats[i+1] = min(opttree[[i+1]]$Get("prnstat", filterFun = isNotLeaf))
    }
  }
  opttree=opttree[-(i+1)]
  a = list(pstats,opttree)
  return(a)
} 
listprune = lrpruningfun(intree)
listprunestat = listprune[[1]]
listprunetree = listprune[[2]]


#II) Get list of geometric means of pruning stat 

#function for creating geometric means of sequence of pairs
gmsfun = function(cforgms){
  n=length(cforgms); a=c()
  for (i in 1:(n-1)){
    a[i]=sqrt(cforgms[i]*cforgms[i+1])
  }
  return(a)
}

geomeans = gmsfun(listprunestat)

#III) Get bootstrap samples
#sample range lies between 1 to 5
btsamps = matrix(0,bts[1],floor(bts[2]*nrow(alldata)))
for (i in 1:bts[1]){
btsamps[i,]<- order(sample(1:nrow(alldata),floor(bts[2]*nrow(alldata)),replace=F)) 
}

#IV) Build trees on bootstrap samples using geometric means

#Function for getting list of indices of greatest number in list less than 
#geometric means 
vecfromgms <- function(prunestatlist) {
  lg = length(geomeans)
  retlist = c()
  for (i in 1:lg){
    if (sum(prunestatlist <= geomeans[[i]]) == 0){
      retlist[i]=1
    } else {
      retlist[i] = max(which(prunestatlist <= geomeans[[i]]))
    } 
  }
  return(retlist)
}

#Create tree for each bootstrap sample
treeforfolds = function(){
  f = bts[1]
  lg = length(geomeans)
  lfoldtrees = list()
  for (i in 1:bts[1]){
    dipolartree2<-dipolarmodel$createtree(btsamps[i,]) 
    list.fold.prune = lrpruningfun(dipolartree2)
    list.fold.prunestat = list.fold.prune[[1]]
    list.fold.prunetree = list.fold.prune[[2]]
    a = vecfromgms(list.fold.prunestat)
    trainfoldtree = list.fold.prunetree[a]
    lfoldtrees[[i]] = trainfoldtree
  }
  return(lfoldtrees)
}
lfoldtrees = treeforfolds() 



#V) Finding best pruned tree from list of pruned trees from whole training set 


#function for sending a test set down a tree, subsets get attached to nodes 
sendtesttree <- function(treetoaddtest,testX){
  splittestatnode <- function(nodein){
    v = nodein$optv
    Z = as.matrix(cbind(1, nodein$testX))
    splits = Z %*% v
    if (splitrule == 0){
    nodein$children[[1]]$testX = nodein$testX[splits < 0, ]
    nodein$children[[2]]$testX = nodein$testX[splits >= 0, ]
    } else {
      nodein$children[[1]]$testX = nodein$testX[splits <= 0, ]
      nodein$children[[2]]$testX = nodein$testX[splits > 0, ]
    }
  }
  treetoaddtest$testX<-testX
  treetoaddtest$Do(function(node) splittestatnode(node), filterFun = isNotLeaf)
}

#Function for computing log rank statistic with X sent down tree
logrank.compute <- function(traintree){
  logrank.node.compute = function(trainnode){
    subnames = strtoi(rownames(trainnode$testX))
    subtestdata = alldata[subnames,]
    Y = Surv(
      subtestdata[time][[1]], subtestdata[censor][[1]] == 1
    )
    v = trainnode$optv
    Z = as.matrix(cbind(1, trainnode$testX))
    splits = Z %*% v
    if (splitrule == 0){
    lrstat = survdiff(Y ~ splits < 0)[[5]]
    } else {
      lrstat = survdiff(Y ~ splits <= 0)[[5]]
    }
    return(lrstat)
  } 
  traintree$Do(function(node) node$lrtest<-logrank.node.compute(node), filterFun = isNotLeaf)
}

#Function for computing goodness of fit stat for data X_a sent down tree T
GFSfun <- function(testX, tree){
  sendtesttree(tree,testX)
  logrank.compute(tree)
  GFS<-sum(tree$Get("lrtest",filterFun=isNotLeaf)) 
} 
    
X = intree$data
f = bts[1]
lg = length(geomeans)
Gam = matrix(nrow = f, ncol = lg)
for (i in 1:f){
  for (j in 1:lg){
    bstree = lfoldtrees[[i]][[j]]
    Ga = GFSfun(X,bstree) - sum(bstree$Get("lrstat",filterFun=isNotLeaf)) 
    Gam[i,j] = Ga
  }
}

obks =  (1/bts[1])*colMeans(Gam)


statformax <- function(alpha){
a=c() 
for (i in 1:lg){
  treeformax = listprunetree[[i]]
  a[i] = sum(treeformax$Get("lrstat",filterFun=isNotLeaf)) + obks[i] - 
          alpha*(treeformax$totalCount - treeformax$leafCount) 
}
return(a)
}

bestgeomean = geomeans[which.max(statformax(1))] 

FinalPrunedtree = listprunetree[max(which(listprunestat <= bestgeomean))][[1]]
}
