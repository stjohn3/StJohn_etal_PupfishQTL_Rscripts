# Script 06
# Find out which genes, that are withing a QTL regions, also show signs of being under selection in the wild

#load in files that contain genes under selection 
P_Sweeps<-read.table("YOUR_DIRECTORY/MP_95th_20k_P_sweep_gene_annotations.txt",header = T)%>%
  mutate(Gene=tolower(Gene))%>%
  mutate(batch=paste("P_Sweeps"))


M_Sweeps<-read.table("YOUR_DIRECTORY/QTL_CRP_LL/MP_95th_20k_M_sweep_gene_annotations.txt",header = T)%>%
  mutate(Gene=tolower(Gene))%>%
  mutate(batch=paste("M_Sweeps"))

All.Swept.fixed.genes<-rbind(P_Sweeps, M_Sweeps)


# Make dataframe of genes that are under selection and overlap with QTL regions
left_join(genes, #data frame of genes within QTL ranges
          All.Swept.fixed.genes, #dataframe of genes under selection (change based on which one you want to use)
          by=c("gene"="Gene"))%>%
  na.omit()%>%
  distinct()->Swept.Genes.WithinQTL

#Add in information about whether these genes are SGV, Introgressed, or de novo
#read in files containing this information
M_local<-read.table("YOUR_DIRECTORY/M_sweep_snps_gene_localities.txt", header=T)%>%
  mutate(species="M")

P_local<-read.table("YOUR_DIRECTORY/P_sweep_snps_gene_localities.txt", header=T)%>%
  mutate(species="P")

All_Local<-rbind(M_local,P_local)%>%
  mutate(Gene=tolower(Gene))

genes.per.qtl.localities<-left_join(Swept.Genes.WithinQTL, #dataframe of genes with signatures of selection
                                    All_Local, #dataframe with SGV, introgressed, and De novo info
                                    by = c("gene"="Gene"))%>%
  select(gene,Intro_pop,Local,Dist.Cat,species)%>%
  distinct()