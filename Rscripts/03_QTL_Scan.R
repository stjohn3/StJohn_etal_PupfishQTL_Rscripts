# Script 03
# Scan for QTL peaks

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
  #kinship matrix
  #Make kinship matrix
  kinship <- calc_kinship(pr)
  kinship_loco <- calc_kinship(pr, "loco")
  
  # use only if applicable
  #Get covariate
  Xcovar <- get_x_covar(mycross)
  
  #scan1
  #Scan1 HK
  out <- scan1(pr, mycross$pheno, Xcovar=Xcovar)
  out <- clean(out)
  #Scan1 LM
  out_pg <- scan1(pr,mycross$pheno, kinship, Xcovar=Xcovar)
  out_pg <- clean(out_pg)
  #Scan1 LOCO
  out_pg_loco <- scan1(pr,mycross$pheno, kinship_loco, Xcovar=Xcovar)
  out_pg_loco <- clean(out_pg_loco)
  
  #permutations
  #HK permutations
  operm <- scan1perm(pr, mycross$pheno, Xcovar=Xcovar, n_perm=1000,cores=((detectCores(all.tests = FALSE, logical = TRUE))/2))
  
  thr = summary(operm,alpha=c(manual.alpha))
  QTLPeak_DF_HK<-find_peaks(scan1_output = out,
                                map = map,
                                threshold = thr, 
                                prob = 0.95, expand2markers = TRUE)%>%
    mutate(method=paste("HK"))
  
  #LM
  operm <- scan1perm(pr, mycross$pheno, kinship=kinship, Xcovar=Xcovar, n_perm=1000,cores=((detectCores(all.tests = FALSE, logical = TRUE))/2))
  
  thr = summary(operm,alpha=c(manual.alpha))
  QTLPeak_DF_LM<-find_peaks(scan1_output = out_pg,
                                map = map,
                                threshold = thr,
                                prob = 0.95, expand2markers = TRUE)%>%
    mutate(method=paste("LM"))
  #LOCO
  
  operm <- scan1perm(pr, mycross$pheno, kinship=kinship_loco, Xcovar=Xcovar, n_perm=1000,cores=((detectCores(all.tests = FALSE, logical = TRUE))/2))
  
  thr = summary(operm,alpha=c(manual.alpha))
  QTLPeak_DF_LOCO<-find_peaks(scan1_output = out_pg_loco,
                                  map = map,
                                  threshold = thr,
                                  prob = 0.95,
                                  expand2markers = TRUE)%>%
    mutate(method=paste("LOCO"))
  

  #Calculate Bayes credible intervals
  Interval_peak<-bayes_int(out, #changed based on scan type 
                           map, 
                           lodcolumn= 1, 
                           chr=3,
                           prob=0.95)