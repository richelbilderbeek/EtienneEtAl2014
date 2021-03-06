# Do maximum likelihood parameter estimation on 'pbd_sim_dens2_' files,
# which contains the branching times of simulated PBD phylogenies.
# Write the parameter estimates to 'pbdbd.Rout'

library(laser)
library(TESS)
library(DDD)
library(PBD)
library(subplex)
muvec = c(0,0.1,0.2)
la2vec = c(0.1,0.3,1)
for(i in 1:3) { 
  for(j in 1:3) { 
    for(k in 1:1000) {
      datafile = paste('pbd_sim_dens2_',i,'-',j,'-3.out',sep = '')
      data = scan(file = datafile,skip = k - 1, nlines = 1)
      flush.console()
      pars1 = c(0.5,muvec[i],la2vec[j],muvec[i])
      out1 = bd_ML(brts = data,initparsopt = c(pars1[1]*pars1[3]/(pars1[3]+pars1[4]),pars1[2] + 0.01),cond = 1,soc = 2,btorph = 0)
      write(c(pars1,k,length(data)+1,as.numeric(unlist(out1))),file = 'pbdbd.Rout',ncolumns = 13, append = T, sep = ',')
    }
  }
}