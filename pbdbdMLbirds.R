library(PBD)

k = 1
load("pbdbd-birdsdata.rda")	
newlist = GenusTreeList[[k]]
resultsTab = expand.grid(tree=k,initpar=1:6,clade=1:length(newlist),cladesize=NA,b=NA,mu_1=NA,lambda_1=NA,mu_2=NA,loglik=NA,df=NA,conv=NA,la=NA,mu=NA,loglikbd=NA,dfbd=NA,convbd=NA)
resultsTabshort = NULL
initparsopttable = rbind(c(0.1,0.01,0.1),c(0.1,1,0.1),c(0.1,0.01,1),c(1,0.01,1),c(1,1,1))

for(i in 1:length(newlist))
{
	tree = newlist[[i]][[1]]
	missingspec = newlist[[i]][[2]]
	branch.times = branching.times(tree)

  if(length(branch.times) >= 9)
  {
   	for(it in 1:5)
    {
	    focrow = which(resultsTab$clade==i & resultsTab$initpar==it)
      res1 = pbd_ML(branch.times, initparsopt = initparsopttable[it,], exteq = 1, missnumspec = missingspec)
      res2 = bd_ML(branch.times, initparsopt = initparsopttable[it,1:2], missnumspec = missingspec)
	    resultsTab[focrow,4:16] = c(length(branch.times)+1,res1,res2[c(1:2,5:7)])
    }
    focrows = which(resultsTab$clade==i & resultsTab$initpar < 6)
    focrow = which(resultsTab[focrows,8] == max(resultsTab[focrows,8]))
    resultsTab[which(resultsTab$clade==i & resultsTab$initpar==6),4:11] = resultsTab[focrows[focrow[1]],4:11]
    focrow = which(resultsTab[focrows,9] == max(resultsTab[focrows,9]))
    resultsTab[which(resultsTab$clade==i & resultsTab$initpar==6),12:16] = resultsTab[focrows[focrow[1]],12:16]
    resultsTabshort = rbind(resultsTabshort,resultsTab[which(resultsTab$clade==i & resultsTab$initpar==6),])
  	print(i)
  }
}
save(resultsTab,resultsTabshort,file = 'pbdbd-birdsout.RData')
