# Script 02.B
# Assembling Linkage Map for Crescent Pond

#### Estimate initial linkage map ####
CRP.ASM.Map<-mstmap.cross(CRP.ASM.Pull,
                          bychr = F,
                          trace = TRUE,
                          dist.fun ="kosambi",
                          p.value = 1e-20,
                          id="ID")

#### Push Back some of the pulled markers ####
CRP_push<-CRP.ASM.Map

CRP_push<-pushCross(CRP_push, type = "co.located")
CRP_push<-pushCross(CRP_push, type = "seg.distortion", pars =list(seg.ratio =".3:.4:.3"))

#### Estimate recombination fractions and form linkage groups ####
CRP.ASM.Map.df<-CRP_push
CRP.ASM.Map.df<-est.rf(CRP.ASM.Map.df)
CRP.ASM.Map.df <- formLinkageGroups(CRP.ASM.Map.df,
                                    max.rf=0.35,
                                    min.lod=5,
                                    reorgMarkers=TRUE)
#### Identify problem markers ####
dropone <- droponemarker(CRP.ASM.Map.df,error.prob=0.01)
badmar <- rownames(summary(dropone, lod.column=2))
CRP.ASM.Map.drop  <- drop.markers(CRP.ASM.Map.df, badmar)

# for loop to generate list of problematic markers
list.markers.to.drop<-NULL
for(i in 28:length(CRP.ASM.Map.drop$geno)){
  drop.df<-CRP.ASM.Map.drop$geno[[i]]
  drop.df$map
  list.markers.to.drop<-append(list.markers.to.drop, drop.df$map)
}

badmar<-rownames(as.data.frame(list.markers.to.drop))

CRP.ASM.Map.drop  <- drop.markers(CRP.ASM.Map.drop, badmar)

#### Merge an linkage groups that were split up (only use if needed) ####
CRP.ASM.Map.drop.merge<-mergeCross(CRP.ASM.Map.drop,
                                   merge = list(
                                     "14" = c("14","26", "27")))



CRP.ASM.Map.drop.merge<-mergeCross(CRP.ASM.Map.drop.merge,
                                   merge = list(
                                     "8" = c("8", "25")))

heatMap(CRP.ASM.Map.drop.merge, lmax = 70)

#### Write out Linkage Map ####
write.cross(CRP.ASM.Map.drop.merge,"csvs", "YOURDIRECTORYHERE")

