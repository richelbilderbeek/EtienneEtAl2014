library(ape)
library(PBD)

# Fabre tree
datatree = read.tree('http://datadryad.org/bitstream/handle/10255/dryad.8494/PRIMATES-Chronogram_273_taxa.tre')
brts = rev(sort(as.numeric(branching.times(datatree))))
missnumspec = 103
pars2 = c(1,1,0)
endmc = 1000
result = pbd_bootstrap(brts, endmc = endmc,missnumspec = missnumspec)
save(result,file = 'Primates_Fabre_bootstrap.Rdata')

# Springer trees
datatree = read.nexus('Primates_Springer.tre')
outgroup = c("Galeopterus_variegatus","Cynocephalus_volans","Lagomorpha","Tupaia_glis","Tupaia_minor")
result = list()
for(i in 1:4)
{
   datatree[[i]] = drop.tip(datatree[[i]],outgroup)
   brts = 100*rev(sort(as.numeric(branching.times(datatree[[i]]))))
   missnumspec = 9
   pars2 = c(1,1,0)
   endmc = 1000
   MLout = pbd_ML(brts,initparsopt = c(0.66,0.63,0.99),idparsopt = 1:3,exteq = 1,missnumspec = missnumspec)
   MLpars = as.numeric(unlist(MLout[1:4]))
   exp_durspec = pbd_durspec_mean(c(MLpars[1],MLpars[4],MLpars[3]))
   median_durspec = pbd_durspec_quantile(c(MLpars[1],MLpars[3],MLpars[4]),0.5)
   if(i == 1)
   {
       result[[1]] = unlist(cbind(ntips = length(brts) + 1,MLout,exp_durspec,median_durspec))
   } else {
       result[[1]] = rbind(result[[1]],unlist(cbind(ntips = length(brts) + 1,MLout,exp_durspec,median_durspec)))
   }
}
save(result,file = 'Primates_Springer.Rdata')
