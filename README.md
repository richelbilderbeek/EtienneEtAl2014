# EtienneEtAl2014

Non-official repository investigating the scripts from 'Estimating the duration of speciatition from phylogenies' by Etienne, Morlon &amp; Lambert, 2014.

These are just the scripts downloaded from the Dryad repository (DOI: `10.5061/dryad.js88n`) to learn from and experiment with.

As the PBD package is licensed under GPL-2, I assume this code uses the same license.

## File descriptions

From a `readme.txt` file:

```
PBD_0.8.tar.gz contains the R package PBD needed to estimate PBD model parameters and properties of the duration of speciation.

pbd_sim_dens2.m contains the Matlab code to simulate, using the Gillespie algorithm, phylogenies under the protracted speciation model. The function pbd_sim_cpp in the PBD package does this more efficiently, but this was not used for generating the simulated data of the paper.

pbd_sim_dens2*.out contain the phylogenetic branching times resulting from calling pbd_sim_dens2.m with various parameters (which are in the *-part).

pbd_sim_dens2_ML.out contains the resulting maximum likelihood parameters and likelihood.

pbd_sim_recstr_full.m is auxiliary to pbd_sim_dens2.m. This was written by Bart Haegeman.

pbd_sim_step2a.m is auxiliary to pbd_sim_dens2.m. This was written by Bart Haegeman.

pbdbd.R contains R code to estimate the parameters of the constant-rate birth-death model for the phylogenies in pbd_sim_dens2*.out

pbdbd.Rout contains the result of pbdbd.R

pbdbd-birdsdata.RData contains the phylogenetic branching times of the bird clades. This was provided by Alex Pigot.

pbdbd-birdsout.RData contains the output after running pbdbdMLbirds.R on pbdbd-birdsdata.RData

pbdbdMLbirds.R contains R code to estimate parameters of the protracted speciation model and the constant-rate birth-death model for the avian clades in pbdbd-birdsdata.RData.

pbdMLprimates.R contains R code to estimate parameters of the protracted speciation model for the primate clades.

Primates_Fabre.tre contains the Newick tree of the primates according to Fabre et al 2009. See the main text for the Dryad link.

Primates_Fabre_bootstrap.Rdata contains the results of running pbd_bootstrap on the Fabre et al 2009 primate tree.

Primates_Springer.Rdata contains the results of running pbd_ML on the Springer et al 2012 primate tree.

Primates_Springer.tre contains the Newick trees of the primates according to Springer et al 2012. See the main text for the Dryad link.
```
