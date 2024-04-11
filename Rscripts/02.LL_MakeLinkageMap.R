# Script 02.A
# Assembling Linkage Map for Little Lake

#### Estimate initial linkage map ####
LL.ASM.Map<-mstmap.cross(LL.ASM.Pull,
                         bychr = F,
                         trace = TRUE,
                         dist.fun ="kosambi",
                         p.value = 1e-14,
                         id="ID")

#### Push back some of the pulled markers ####
LL_push<-LL.ASM.Map

LL_push<-pushCross(LL_push, type = "co.located")
LL_push<-pushCross(LL_push, type = "seg.distortion", pars =list(seg.ratio =".3:.4:.3"))

#### Estimate recombination fractions and form linkage groups ####
LL.ASM.Map.df<-LL_push
LL.ASM.Map.df<-est.rf(LL.ASM.Map.df)
LL.ASM.Map.df <- formLinkageGroups(LL.ASM.Map.df,
                                   max.rf=0.35,
                                   min.lod=5,
                                   reorgMarkers=TRUE)

#### Identify problem markers ####
dropone <- droponemarker(LL.ASM.Map.df,error.prob=0.01)
badmar <- rownames(summary(dropone, lod.column=2))
LL.ASM.Map.drop  <- drop.markers(LL.ASM.Map.df, badmar)

#for loop to generate list of problematic markers
list.markers.to.drop<-NULL
for(i in 25:length(LL.ASM.Map.drop$geno)){
  drop.df<-LL.ASM.Map.drop$geno[[i]]
  print(drop.df$map)
  list.markers.to.drop<-append(list.markers.to.drop, drop.df$map)
}

badmar<-rownames(as.data.frame(list.markers.to.drop))

LL.ASM.Map.drop  <- drop.markers(LL.ASM.Map.drop, badmar)


#### Merge any linkage groups that were split up (only use if needed) ####
LL.ASM.Map.drop.merge<-mergeCross(LL.ASM.Map.drop,
                                  merge = list(
                                    "5" = c("5.1", "5.2"),
                                    "7" =c("7.1", "7.2"),
                                    "18"=c("18.1","18.2"),
                                    "23"=c("23.1","23.2")
                                  ))
LL.ASM.Map.drop.merge_24<-LL.ASM.Map.drop.merge

#### Write out Linkage Map ####
write.cross(LL.ASM.Map.drop.merge_24,"csvs", "YOURDIRECTORYHERE")

