# Script 00
# Load libraries and raw genotype and phenotype files for use with Rqtl, Rqtl2, and ASMap

#### Load Libraries ####
library("librarian")
shelf(qtl,qtl2,magrittr,dplyr,tidyverse,parallel,LinkageMapView,
ASMap,qtl2convert,IRanges,GenomicRanges,data.table,plyr,
ggfortify,circlize,Hmisc,psych,qtlDesign,qtl2pleio,cowplot)

#### Set working directory ####
setwd(PLACE YOUR WORKING DIRECTORY HERE)

#### Load in genotype and phenotypes files for Rqtl and Rqtl2 ####

Crescent.Pond<-read.cross("csvs", dir=".",
                                          genfile="./Genotype_Files/CRP_MxP_QTL_pupfish_F2_geno.csv",
                                          phefile="./Phenotype_Files/CRPHybrid_Phenotypes.csv",
                                          genotypes=c("a","h","b"), 
                                          na.strings=c("-"))
Crescent.Pond.mapthis<-Crescent.Pond
Little.Lake<-read.cross("csvs", dir=".",
                                        genfile="LL_MxP_QTL_pupfish_F2_geno.csv",
                                        phefile="LLHybrid_Phenotypes.csv",
                                        genotypes=c("a","h","b"), 
                                        na.strings=c("-"))
Little.Lake.mapthis<-Little.Lake

#### Convert rqtl format for use with ASMap package ####
CRP.ASM<-convert2bcsft(Crescent.Pond, BC.gen=0, F.gen=2, estimate.map=FALSE)
LL.ASM<-convert2bcsft(Little.Lake, BC.gen=0, F.gen=2, estimate.map=FALSE)