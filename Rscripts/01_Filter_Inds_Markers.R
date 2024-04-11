# Script 01
# Filtering raw data before analysis

#### Removing individuals with 100% missing data ####
LL.ASM <- subset(LL.ASM, ind=(ntyped(LL.ASM)>1))
CRP.ASM <- subset(CRP.ASM, ind=(ntyped(CRP.ASM)>1))

#### Ensuring correct individuals are in Little Lake genotype data ####

LL.ASM<-subset(LL.ASM, ind = c("LLPM001","LLPM002","LLPM003","LLPM004","LLPM006","LLPM007","LLPM010","LLPM014","LLPM036","LLPM048","LLPM049","LLPM051","LLPM052","LLPM053","LLPM054",
                               "LLPM055","LLPM057","LLPM058","LLPM060","LLPM061","LLPM062","LLPM063","LLPM064","LLPM065","LLPM066","LLPM067","LLPM068","LLPM069","LLPM070","LLPM071",
                               "LLPM072","LLPM073","LLPM074","LLPM075","LLPM076","LLPM077","LLPM079","LLPM080","LLPM081","LLPM082","LLPM083","LLPM084","LLPM085","LLPM086","LLPM087",
                               "LLPM088","LLPM089","LLPM090","LLPM091","LLPM092","LLPM093","LLPM094","LLPM095","LLPM102","LLPM104","LLPM108","LLPM110","LLPM111","LLPM112","LLPM113",
                               "LLPM114","LLPM117","LLPM119","LLPM120","LLPM122","LLPM123","LLPM125","LLPM127","LLPM131","LLPM135","LLPM138","LLPM146","LLPM147","LLPM152","LLPM153",
                               "LLPM155","LLPM158","LLPM159","LLPM172","LLPM178","LLPM181","LLPM183","LLPM184","LLPM186","LLPM187","LLPM189","LLPM190","LLPM193","LLPM199","LLPM200",
                               "LLPM203","LLPM216","LLPM217","LLPM220","LLPM221","LLPM228","LLPM230","LLPM231","LLPM234","LLPM235","LLPM241","LLPM244","LLPM246","LLPM249","LLPM250",
                               "LLPM259","LLPM265","LLPM266","LLPM267","LLPM268","LLPM271","LLPM272","LLPM273","LLPM275","LLPM280","LLPM283","LLPM285","LLPM287","LLPM289","LLPM292",
                               "LLPM293","LLPM296","LLPM298","LLPM303","LLPM304","LLPM309","LLPM310","LLPM315","LLPM319","LLPM320","LLPM321","LLPM325","LLPM326","LLPM329","LLPM331",
                               "LLPM332","LLPM335","LLPM336","LLPM338","LLPM339","LLPM342","LLPM343","LLPM344","LLPM345","LLPM346","LLPM347","LLPM350","LLPM351","LLPM352","LLPM354",
                               "LLPM355","LLPM356","LLPM357","LLPM358","LLPM359","LLPM360","LLPM362","LLPM363","LLPM364","LLPM365","LLPM366","LLPM367","LLPM370","LLPM371","LLPM372",
                               "LLPM373","LLPM374","LLPM380","LLPM382","LLPM389","LLPM391","LLPM392","LLPM394","LLPM400","LLPM406","LLPM409","LLPM410","LLPM420","LLPM422","LLPM424",
                               "LLPM425","LLPM426","LLPM431","LLPM433","LLPM434","LLPM436","LLPM437","LLPM438","LLPM439","LLPM440","LLPM441","LLPM442","LLPM443","LLPM444","LLPM445",
                               "LLPM447","LLPM448","LLPM451","LLPM455","LLPM459","LLPM460","LLPM461","LLPM462","LLPM463","LLPM464","LLPM465","LLPM466","LLPM468","LLPM469","LLPM475",
                               "LLPM476","LLPM477","LLPM478","LLPM479","LLPM481","LLPM482","LLPM483","LLPM484","LLPM486","LLPM487","LLPM488","LLPM490","LLPM491","LLPM492","LLPM494",
                               "LLPM496","LLPM498","LLPM500","LLPM501","LLPM502","LLPM503","LLPM504","LLPM505","LLPM506","LLPM507","LLPM508","LLPM509","LLPM511","LLPM512","LLPM513",
                               "LLPM514","LLPM517","LLPM518","LLPM520","LLPM522","LLPM524","LLPM525","LLPM526","LLPM527","LLPM528","LLPM529","LLPM530","LLPM531","LLPM532","LLPM533",
                               "LLPM535","LLPM537","LLPM538","LLPM540","LLPM541","LLPM542","LLPM543","LLPM544","LLPM545","LLPM546","LLPM547","LLPM548","LLPM549","LLPM550","LLPM551",
                               "LLPM552","LLPM553","LLPM554","LLPM555","LLPM556","LLPM557","LLPM558","LLPM559","LLPM560","LLPM561","LLPM563","LLPM564"))

#### Remove Markers with Extremely High or Low Heterozygosity ####
Heterozygote.markers.CRP<-statMark(CRP.ASM, stat.type = "marker")$marker$AB
CRP.ASM<-drop.markers(CRP.ASM, c(markernames(CRP.ASM)[Heterozygote.markers.CRP > 0.98], markernames(CRP.ASM)[Heterozygote.markers.CRP < 0.1])) 

Heterozygote.markers.LL<-statMark(LL.ASM, stat.type = "marker")$marker$AB
LL.ASM<-drop.markers(LL.ASM, c(markernames(LL.ASM)[Heterozygote.markers.LL > 0.98], markernames(LL.ASM)[Heterozygote.markers.LL < 0.1])) 

#### Repeat removing individuals with 100% missing data since we removed markers ####
LL.ASM <- subset(LL.ASM, ind=(ntyped(LL.ASM)>1))
CRP.ASM <- subset(CRP.ASM, ind=(ntyped(CRP.ASM)>1))

#### Pulling out markers that suffer from segregation distortion ####
# will see if some of these markers can be retained in future step
LL.ASM%>%
  pullCross(., type = "missing", pars = list(miss.thresh =0.75))%>%
  pullCross(., type = "seg.distortion", pars =list(seg.ratio =".25:.5:.25"))%>%
  pullCross(., type = "co.located") ->LL.ASM.Pull

LL.ASM.Pull


CRP.ASM%>%
  pullCross(., type = "missing", pars = list(miss.thresh =0.72))%>%
  pullCross(., type = "seg.distortion", pars =list(seg.ratio =".25:.5:.25"))%>%
  pullCross(., type = "co.located") ->CRP.ASM.Pull

CRP.ASM.Pull