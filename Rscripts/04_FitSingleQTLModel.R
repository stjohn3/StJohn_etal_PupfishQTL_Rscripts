# Script 04
# Fit a single-QTL model at putative QTL position

#Variables for cross and map
mycross <- convert2cross2(YOUR_CROSS_HERE)
map <- insert_pseudomarkers(mycross$gmap, step=1)

#genotype probability 
#Calculating genotype probs
pr<-calc_genoprob(
  mycross,
  map,
  error_prob = 0.0001,
  map_function ="kosambi",
  lowmem = FALSE,
  cores = 0)
pr<-clean(pr)

# use only if applicable
#Get covariate
Xcovar <- get_x_covar(mycross)

#scan1
#Scan1 HK
out <- scan1(pr, mycross$pheno, Xcovar=Xcovar)

#Find position with max LOD score
max_pos <- max(out, map,  lodcolumn = 1)

#Pull out genotype probabilities 
pr_max <- pull_genoprobpos(pr, map, max_pos$chr, max_pos$pos)  

#Fit a single QTL model
out_fit1 <- fit1(pr_max, mycross$pheno[,1], addcovar=Xcovar, blup=F) 

#calculate % variance explained
LOD<-(out_fit1$lod)
PVE<-1-(10^(-(2/length(out_fit1$fitted))*out_fit1$lod))