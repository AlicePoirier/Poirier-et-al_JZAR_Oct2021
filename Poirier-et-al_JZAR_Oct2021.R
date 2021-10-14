#
#### R code for Poirier et al., Journal of Zoo and Aquarium Research, submitted October 2021
#### Olfactory communication in tamarins: a comparative analysis of scents from wild and captive populations
#



### Datasets ----
# Poirier-et-al_JZAR_Dataset_Oct2021.csv
library(dplyr)

Dataset<-read.csv(file.choose(), header=T, na.strings="NA", sep=",")
Dataset[,c(1:12)]<- lapply(Dataset[,c(1:12)], factor)
Dataset<- Dataset %>% arrange(Sample, Area)
head(Dataset)


# Standardisation, then arcsin-transformation to reduce bottom/ceiling effect, and then log to meet distribution requirements
Dataset.std<- data.frame(Dataset %>%
                     group_by(Sample) %>%
                     summarise(Area=Area, SArea=sum(Area),
                               RPA=round((Area/sum(Area))*100, digits=3),
                               logRPA=round(log(((Area/sum(Area))*100)+1), digits=3),
                               asinRPA=round((log(asin(sqrt((Area/sum(Area)*100)/100))+0.01)), digits=5)) %>%
                     add_tally(name="Npeak") )
Dataset.std<- Dataset.std %>% arrange(Sample, Area)

Dataset.all<- cbind(Dataset, Dataset.std[,-(1:2)])
head(Dataset.all)


# Sample.list with number of compounds in each sample
library(reshape2)

Sample.list<- melt(xtabs(~Site+Group+Type+Sample, data=Dataset.all, addNA=T))
Sample.list<- Sample.list[which(Sample.list$value>0),]
head(Sample.list)




### Descriptive stats ----
min(Sample.list$value)
max(Sample.list$value)
mean(Sample.list$value)
sd(Sample.list$value)


Sample.list %>% group_by(Site) %>%  
  summarise(min=min(value), max=max(value), mean=mean(value), sd=sd(value), med=median(value),
               se=round(sd(value)/sqrt(sum(n())), digits=3))

Sample.list %>% group_by(Type) %>%  
  summarise(min=min(value), max=max(value), mean=mean(value), sd=sd(value), med=median(value),
            se=round(sd(value)/sqrt(sum(n())), digits=3))




### Tests on number of compounds ----
# Material & Methods - test difference between body and gland
wilcox.test(value~Type, data=Sample.list)
wilcox.test(value~Type, data=Sample.list[which(Sample.list$Site=="Captive"),])
wilcox.test(value~Type, data=Sample.list[which(Sample.list$Site=="Wild"),])

# Material & Methods - test difference between wild groups
wilcox.test(value~Group, data=Sample.list[which(Sample.list$Site=="Wild"),])

# Results - comparison captive/wild
wilcox.test(value~Site, data=Sample.list)




### Non-metric Multidimentional Scaling ----
library(vegan)

head(Dataset.all)
A.cla<-arrange(unique(Dataset.all[,c(1:11)]), Sample)
mamaA<- xtabs(logRPA~Sample+Peak, Dataset.all)
mamaA[1:6,1:6]
vdimaA<-vegdist(mamaA, method="bray")

NMDSA2d<- metaMDS(vdimaA, k=2, autotransform=F, noshare=F, wascores=F, expand=T, trymax=200, trace=F)
nmdsA2d<- scores(NMDSA2d)
colnames(nmdsA2d)<- c("NMDS1.2d","NMDS2.2d")
NMDSA2d$stress

Dataset.NMDS<- cbind(A.cla, nmdsA2d)
head(Dataset.NMDS)


# ANOSIM
# Material & Methods - gland vs body
anosim(vdimaA, strata=A.cla$Site, grouping=A.cla$Type, permutations=999, distance="bray")

# Results - comparison captive/wild
anosim(vdimaA, A.cla$Site, permutations=999, distance="bray")




### Figure 2 ----
library(ggplot2)
library(ggsignif)
library(cowplot)

# Number of compounds
Site.Npeak<-
  ggplot(Sample.list, aes(x=Site, y=value, fill=Site))+
  geom_boxplot(color="black")+
  geom_point(shape=1, size=2)+
  scale_fill_manual(values=c("#E69F00","#009E73"))+   # colorblind dark yellow + green
  xlab("Study condition")+
  ylab("Number of compounds detected")+  
  labs(title=expression(paste(bold("a."))))+
  guides(fill="none")+
  geom_signif(comparisons = list(c("Captive", "Wild")), test="wilcox.test", map_signif_level=TRUE)+
  theme_classic(base_size=11)+ 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))


# NMDS chemical composition
Site.nmds<-
  ggplot(Dataset.NMDS, mapping=aes(x=NMDS1.2d, y=NMDS2.2d, shape=Site, color=Site, linetype=Site))+
  geom_point(size=2)+
  scale_color_manual(values=c("#E69F00","#009E73"))+    # colorblind dark yellow + green
  scale_shape_manual(values=c(15,16))+
  stat_ellipse(aes(x=NMDS1.2d, y=NMDS2.2d), level=0.95)+
  scale_linetype_manual(values=c("longdash","solid"))+
  xlab("NMDS 1") + ylab("NMDS 2")+ 
  labs(title=expression(paste(bold("b."))))+
  guides(colour="none", shape="none", linetype="none", size="none")+
  annotate("text", x=0.2, y=-0.7, label="2D stress = 0.055")+
  theme_classic(base_size=11)+ 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))


# Plot
#dev.new()
Fig2<- plot_grid(Site.Npeak, Site.nmds, rel_heights=c(1,1), ncol=2, align="h")
ggsave(plot=Fig2, filename="enter destination file here", 
       device="tiff", scale=1.4, dpi=300, units="cm", width=15.5, height=7)




### Tests on asinRPA, only peaks found in both sites ----
CommonPeaks<- c("P07","P11","P13","P17","P21","P34","P49","P53")
Dataset.commonPeaks<- filter(Dataset.all, Peak %in% CommonPeaks)

#
for (i in 1:length(CommonPeaks)){
  Peaki=CommonPeaks[i]
  print(Peaki)
  print(wilcox.test(asinRPA~Site, data=Dataset.commonPeaks[which(Dataset.commonPeaks$Peak==Peaki),]))
}
# Signif: P07 and P21 only


