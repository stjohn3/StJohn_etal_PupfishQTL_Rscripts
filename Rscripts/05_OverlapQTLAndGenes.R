# Script 05
# Overlap QTL ranges with genes (of known putative function)

#Make data frame of QTL ranges to use when overlapping intervals 
#set important variables
working.map<-map
Trait<-(names(max_pos)[3])
Chr<-as.numeric(((max_pos)[1]))

#get scaffolds that are within QTL peak range
scaffolds<-working.map%>%
  dplyr::filter(chr==Chr)%>%
  dplyr::mutate(pos=as.character(pos),
                pos=as.numeric(pos))%>%
  dplyr::filter(pos<Interval_peak[3] +1, # hi
                pos>(Interval_peak[1]-1)) #lo

#Final input for interval overlap
QTL_Range<-data.frame(Trait=as.character(Trait),
           Chr=(scaffolds))# select column with scaffold names
  
# load in gff file with gene intervals
Gene_list<-read.table("YOUR_DIRECTORY_HERE/c_brontotheroides.known_putative_function.genes_only.gff",header=T)

#update column names
colnames(Gene_list)<-c("scaffold", "start","end","gene")

# remove NAs
Gene_list%<>%
  na.omit()

#Overlap QTL ranges with gene ranges 
gr1<-with(QTL_Range,GRanges(reference.chr, #column name for scaffold ID
                                  IRanges(start=low,
                                          end=high,
                                          names=Trait)))
gr2<-with(Gene_list,GRanges(scaffold,
                            IRanges(start=start,
                                    end=end,
                                    names=gene)))

type1<-findOverlaps(query = gr2, subject = gr1, type = 'within')
type2<-findOverlaps(query = gr2, subject = gr1, type = 'any')

type1.df<-as.data.frame(type1)
type2.df<-as.data.frame(type2)

genes<-rbind(type1.df,type2.df)%>%
  na.omit()
