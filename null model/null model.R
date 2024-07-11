raup_crick_abundance_one_comparison = function(null.one.use,null.two.use,spXsite, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=9999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE){  
  ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended.  The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).  Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  If ties are split (as we recommend) the dissimilarity (default) and similarity (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work. The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.    
  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model    
  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }  
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))  
  ##make the spXsite matrix into a new, pres/abs. matrix:
  ceiling(spXsite/max(spXsite))->spXsite.inc  
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite.inc, MARGIN=2, FUN=sum)  
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(spXsite, MARGIN=2, FUN=sum)
  ##make_null:
  ##looping over each pairwise community combination:
  for(null.one in null.one.use){
    for(null.two in null.two.use){     
      null_bray_curtis<-NULL
      for(i in 1:reps){
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(spXsite[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');      
        ##same for com2:
        com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(spXsite[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        null.spXsite = rbind(com1,com2); # null.spXsite;        
        ##calculate null bray curtis
        null_bray_curtis[i] = distance(null.spXsite,method='bray-curtis');        
      }; # end reps loop
      ## empirically observed bray curtis
      obs.bray = distance(spXsite[c(null.one,null.two),],method='bray-curtis');
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);    
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);     
      rc = (num_less_than_in_null )/reps; # rc;     
      if(split_ties){       
        rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      };     
      if(!classic_metric){          
        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1          
        rc = (rc-.5)*2
      };
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix      
      results.out = as.data.frame(results[null.two,null.one])
      colnames(results.out) = colnames(as.matrix(results))[null.one]
      rownames(results.out) = rownames(as.matrix(results))[null.two]
      print(c(null.one,null.two,date()));     
    }; ## end null.two loop    
  }; ## end null.one loop
  if(as.distance.matrix){ ## return as distance matrix if so desired
    results<-as.dist(results)
  }    
  return(results.out)
}; ## end function
setwd("D:/703/王文银/0C-李飞数据/分长期人工草地/零模型")
otu<-read.csv("bacteria.csv",row.names = 1,header = T)
otu<-t(otu)
dim(otu)
phylo<-read.tree("rep_set_tree_bacteria.tre")
match.phylo.otu = match.phylo.data(phylo, t(otu))
str(match.phylo.otu)
beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T))
dim(beta.mntd.weighted)
beta.mntd.weighted[1:5,1:5]
write.csv(beta.mntd.weighted,'spe.MNTD_weighted_bacteria.csv',quote=F)
identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted))
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted))
beta.reps = 999
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps))
dim(rand.weighted.bMNTD.comp)
for (rep in 1:beta.reps) {  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F)) 
  print(c(date(),rep))  
}
weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data))
dim(weighted.bNTI)
for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,]
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals)
    rm("rand.vals")    
  }
}
rownames(weighted.bNTI) = colnames(match.phylo.otu$data)
colnames(weighted.bNTI) = colnames(match.phylo.otu$data)
weighted.bNTI
write.csv(weighted.bNTI,"spe.weighted_bNTI_bacteria.csv",quote=F)
pdf("spe.weighted_bNTI_Histogram_bacteria.pdf")
hist(weighted.bNTI)
dev.off()

otu = read.csv("bacteria.csv",row.names = 1,header = T)
otu<-t(otu)
library(permute)
library(gee)
library(vegan)
library(ape)
library(picante)
library(ecodist)
#########
######## calculate Bray-Curtis########
#########
library(ecodist)
bray.out = as.matrix(distance(otu,method = 'bray-curtis'))
write.csv(bray.out,"bray_weighted.csv");
####################
########## for Raup-Crick ##########
####################
metric = 'RC'
abund.for.names = 'weighted'
colnames(otu) = gsub(pattern = "X",replacement = "",x = colnames(otu))
print(dim(otu))
no.of.samples = nrow(otu)
metric = 'RC'
abund = T
abund.for.names = 'weighted'
source("Raup_Crick_Abundance_One_Comparison.r")
for (i in 1:(no.of.samples - 1)){
  for (j in (i + 1):(no.of.samples)){
    raup.crick.out = raup_crick_abundance_one_comparison(null.one.use = i,null.two.use = j,otu, plot_names_in_col1=FALSE, classic_metric=FALSE, split_ties=TRUE, reps=999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE);
    print(raup.crick.out)
    write.csv(raup.crick.out,paste('RC_',abund.for.names,'_Diss_comm1_',i,'_comm2_',j,'.csv',sep=""),row.names=T,quote=F);
  }
}
## find RC reps that still need to be done
all.files = list.files(pattern = 'RC_')
all.RC = all.files[grep("RC_weighted_Diss_comm1_",all.files)]
length(all.RC)
bray = read.csv('bray_weighted.csv',row.names=1);
raup.crick.out = bray
raup.crick.out[,] = -999
for (i in 1:length(all.RC)) {
  RC.temp = read.csv(paste(all.RC[i],sep=""),row.names=1)
  if (nrow(RC.temp) != 1 ) {print(c(i,"ERROR"))} 
  if (ncol(RC.temp) != 1 ) {print(c(i,"ERROR"))} 
  raup.crick.out[which(rownames(raup.crick.out) == rownames(RC.temp)),which(colnames(raup.crick.out) == colnames(RC.temp))] = RC.temp[1,1]
  print(RC.temp)
}
raup.crick.out = as.data.frame(as.matrix(as.dist(raup.crick.out)))
head(raup.crick.out)
raup.crick.out[1:5,1:5]
write.csv(raup.crick.out,'RC_weighted_bacteria.csv',quote=F);
pdf(paste('RC.',abund.for.names,'_Histogram_bacteria.pdf',sep=""))
raup.crick.out.for.hist = as.dist(raup.crick.out)
raup.crick.out.for.hist = raup.crick.out.for.hist[raup.crick.out.for.hist != -999]
hist(raup.crick.out.for.hist,breaks=40)
abline(v=c(-0.95,0.95),col=2,lwd=2,lty=2)
dev.off()


bNTI<-read.table("spe.weighted_bNTI_wetland.csv",header=T,sep=",",row.names=1)
rc<-read.table("RC_weighted_wetland.csv",header=T,sep=",",row.names=1)
# 4 # influences of different processes
## match names in bNTI and RC
rcc=rc[match(rownames(bNTI),rownames(rc)),match(colnames(bNTI),colnames(rc))]
bNTI.v<-as.vector(as.dist(bNTI))
rc.v<-as.vector(as.dist(rcc))
id.selectna<-(bNTI.v<=2&bNTI.v>=(-2))
num.pair<-length(bNTI.v)
select.h<-sum(bNTI.v>2)/num.pair#??��ѡ??
select.l<-sum(bNTI.v<(-2))/num.pair#???ʻ?ѡ??
disper.h<-sum(rc.v[id.selectna]>0.95)/num.pair
disper.l<-sum(rc.v[id.selectna]<(-0.95))/num.pair
drift<-sum(rc.v[id.selectna]<=0.95&rc.v[id.selectna]>=(-0.95))/num.pair
res=data.frame(select.h,select.l,disper.h,disper.l,drift,num.pair)
write.csv(res,"gen_Processes_wetland.csv")

bNTI<-read.table("spe.weighted_bNTI_alpine.csv",header=T,sep=",",row.names=1)
rc<-read.table("RC_weighted_alpine.csv",header=T,sep=",",row.names=1)
# 4 # influences of different processes
## match names in bNTI and RC
rcc=rc[match(rownames(bNTI),rownames(rc)),match(colnames(bNTI),colnames(rc))]
bNTI.v<-as.vector(as.dist(bNTI))
rc.v<-as.vector(as.dist(rcc))
id.selectna<-(bNTI.v<=2&bNTI.v>=(-2))
num.pair<-length(bNTI.v)
select.h<-sum(bNTI.v>2)/num.pair#??��ѡ??
select.l<-sum(bNTI.v<(-2))/num.pair#???ʻ?ѡ??
disper.h<-sum(rc.v[id.selectna]>0.95)/num.pair
disper.l<-sum(rc.v[id.selectna]<(-0.95))/num.pair
drift<-sum(rc.v[id.selectna]<=0.95&rc.v[id.selectna]>=(-0.95))/num.pair
res=data.frame(select.h,select.l,disper.h,disper.l,drift,num.pair)
write.csv(res,"gen_Processes_alpine.csv")

bNTI<-read.table("spe.weighted_bNTI_HD.csv",header=T,sep=",",row.names=1)
rc<-read.table("RC_weighted_HD.csv",header=T,sep=",",row.names=1)
# 4 # influences of different processes
## match names in bNTI and RC
rcc=rc[match(rownames(bNTI),rownames(rc)),match(colnames(bNTI),colnames(rc))]
bNTI.v<-as.vector(as.dist(bNTI))
rc.v<-as.vector(as.dist(rcc))
id.selectna<-(bNTI.v<=2&bNTI.v>=(-2))
num.pair<-length(bNTI.v)
select.h<-sum(bNTI.v>2)/num.pair#??��ѡ??
select.l<-sum(bNTI.v<(-2))/num.pair#???ʻ?ѡ??
disper.h<-sum(rc.v[id.selectna]>0.95)/num.pair
disper.l<-sum(rc.v[id.selectna]<(-0.95))/num.pair
drift<-sum(rc.v[id.selectna]<=0.95&rc.v[id.selectna]>=(-0.95))/num.pair
res=data.frame(select.h,select.l,disper.h,disper.l,drift,num.pair)
write.csv(res,"gen_Processes_HD.csv")

bNTI<-read.table("spe.weighted_bNTI_sown.csv",header=T,sep=",",row.names=1)
rc<-read.table("RC_weighted_sown.csv",header=T,sep=",",row.names=1)
# 4 # influences of different processes
## match names in bNTI and RC
rcc=rc[match(rownames(bNTI),rownames(rc)),match(colnames(bNTI),colnames(rc))]
bNTI.v<-as.vector(as.dist(bNTI))
rc.v<-as.vector(as.dist(rcc))
id.selectna<-(bNTI.v<=2&bNTI.v>=(-2))
num.pair<-length(bNTI.v)
select.h<-sum(bNTI.v>2)/num.pair#??��ѡ??
select.l<-sum(bNTI.v<(-2))/num.pair#???ʻ?ѡ??
disper.h<-sum(rc.v[id.selectna]>0.95)/num.pair
disper.l<-sum(rc.v[id.selectna]<(-0.95))/num.pair
drift<-sum(rc.v[id.selectna]<=0.95&rc.v[id.selectna]>=(-0.95))/num.pair
res=data.frame(select.h,select.l,disper.h,disper.l,drift,num.pair)
write.csv(res,"gen_Processes_sown.csv")

bNTI<-read.table("spe.weighted_bNTI_LS.csv",header=T,sep=",",row.names=1)
rc<-read.table("RC_weighted_LS.csv",header=T,sep=",",row.names=1)
# 4 # influences of different processes
## match names in bNTI and RC
rcc=rc[match(rownames(bNTI),rownames(rc)),match(colnames(bNTI),colnames(rc))]
bNTI.v<-as.vector(as.dist(bNTI))
rc.v<-as.vector(as.dist(rcc))
id.selectna<-(bNTI.v<=2&bNTI.v>=(-2))
num.pair<-length(bNTI.v)
select.h<-sum(bNTI.v>2)/num.pair#??��ѡ??
select.l<-sum(bNTI.v<(-2))/num.pair#???ʻ?ѡ??
disper.h<-sum(rc.v[id.selectna]>0.95)/num.pair
disper.l<-sum(rc.v[id.selectna]<(-0.95))/num.pair
drift<-sum(rc.v[id.selectna]<=0.95&rc.v[id.selectna]>=(-0.95))/num.pair
res=data.frame(select.h,select.l,disper.h,disper.l,drift,num.pair)
write.csv(res,"gen_Processes_LS.csv")

setwd("D:/703/wwy/0C-李飞数据/分长期人工草地/零模型")
####( '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', 'gray')
dat<-read.csv("不同类型草地中零模型分类.csv", header=TRUE, row.names = 1)
dat$Taxonomy <- factor(rownames(dat), levels = rev(rownames(dat)))
dat<- melt(dat, id = 'Taxonomy')

#group <- read.delim('group.txt', sep = '\t', stringsAsFactors = FALSE)
#names(group)[1] <- 'variable'
#phylum_top10 <- merge(phylum_top10, group, by = 'variable')

library(reshape2)
library(ggplot2)

###c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', 'gray')
dat$variable<-factor(dat$variable,levels= c("ASM", "AKM", "SDM", "ARS", "ARL" ))
p1 <- ggplot(dat, aes(variable, 100 * value, fill = Taxonomy)) +
  geom_col(position = 'stack', width = 0.6) +
  #facet_wrap(~group, scales = 'free_x', ncol = 2) +
  scale_fill_manual(values =  rev(c('#BEBADA', '#FDB462','#BC80BD', '#CCEBC5', 'gray'))) +
  labs(x = '', y = 'Relative importance (%)') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.title = element_blank(), legend.text = element_text(size = 12),legend.position = "top" )
p1
ggsave('不同类型草地中零模型分类_bacterial.pdf', p1, width = 11, height = 10)
ggsave('不同类型草地中零模型分类_bacterial.png', p1, width = 11, height = 10)

dat<-read.csv("不同类型草地中零模型分类_fungi.csv",header=TRUE,row.names = 1)
dat$Taxonomy <- factor(rownames(dat), levels = rev(rownames(dat)))
dat<- melt(dat, id = 'Taxonomy')

#group <- read.delim('group.txt', sep = '\t', stringsAsFactors = FALSE)
#names(group)[1] <- 'variable'
#phylum_top10 <- merge(phylum_top10, group, by = 'variable')

library(reshape2)
library(ggplot2)

###c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', 'gray')
dat$variable<-factor(dat$variable,levels= c("ASM", "AKM", "SDM", "ARS", "ARL" ))
p2 <- ggplot(dat, aes(variable, 100 * value, fill = Taxonomy)) +
  geom_col(position = 'stack', width = 0.6) +
  #facet_wrap(~group, scales = 'free_x', ncol = 2) +
  scale_fill_manual(values =  rev(c('#BEBADA', '#BC80BD', '#CCEBC5', 'gray'))) +
  labs(x = '', y = 'Relative importance (%)') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.title = element_blank(), legend.text = element_text(size = 12),legend.position = "top")
  
p2
ggsave('不同类型草地中零模型分类_fungi.pdf', p2, width = 11, height = 10)
ggsave('不同类型草地中零模型分类_fungi.png', p2, width = 11, height = 10)



setwd("D:/703/wwy/0C-李飞数据/分长期人工草地/零模型")
data<-read.csv("bNTI.csv",row.names = 1,header = T)
dt = factor(data$group, levels=c('ASM','AKM','SDM','ARS','ARL'))
library(ggplot2)
#windowsFonts(TNM = windowsFont("Times New Roman"))
p3<-ggplot(data, aes(x=dt,y=bNTI))+geom_boxplot(aes(fill=group))+
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ 
  geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(fill=group),width =0.2,shape = 21,size=2.5)+
  scale_color_manual(values=c("black","black","black","black","black","black"))+
  geom_hline(aes(yintercept=2), color="#AA0000",linetype="dashed")+
  geom_hline(aes(yintercept=-2),color="#AA0000",linetype="dashed")+
  theme(legend.position="none", #不需要图例
        axis.text.x=element_text(colour="black",size=14), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(size=14,face="bold"), #设置x轴刻度标签的字体属性
        axis.title.y=element_text(size = 14,face="bold"), #设置y轴的标题的字体属性
        axis.title.x=element_text(size = 14,face="bold"), #设置x轴的标题的字体属性
        plot.title = element_text(size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank(),
  )+
  theme_bw(base_size = 18)+
  theme(panel.grid=element_blank())+
  ylab("βNTI")+xlab("")
#p+ facet_wrap(~ group, scales="free")
p3
ggsave('βNTI_bacteria.pdf', p3, width = 8, height = 6)
ggsave('βNTI-bacteria.png', p3, width = 8, height = 6)



setwd("D:/703/wwy/0C-李飞数据/分长期人工草地/零模型")
data<-read.csv("bNTI_fungi.csv",row.names = 1,header = T)
dt = factor(data$group, levels=c('ASM','AKM','SDM','ARS','ARL'))
library(ggplot2)
#windowsFonts(TNM = windowsFont("Times New Roman"))
p4<-ggplot(data, aes(x=dt,y=bNTI))+geom_boxplot(aes(fill=group))+
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ 
  geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(fill=group),width =0.2,shape = 21,size=2.5)+
  scale_color_manual(values=c("black","black","black","black","black","black"))+
  geom_hline(aes(yintercept=2), color="#AA0000",linetype="dashed")+
  geom_hline(aes(yintercept=-2),color="#AA0000",linetype="dashed")+
  theme(legend.position="none", #不需要图例
        axis.text.x=element_text(colour="black",size=14), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(size=14,face="bold"), #设置x轴刻度标签的字体属性
        axis.title.y=element_text(size = 14,face="bold"), #设置y轴的标题的字体属性
        axis.title.x=element_text(size = 14,face="bold"), #设置x轴的标题的字体属性
        plot.title = element_text(size=15,face="bold",hjust = 0.5)#, #设置总标题的字体属性
        # panel.grid.major = element_blank(), #不显示网格线
        #panel.grid.minor = element_blank()
  )+
  ylab("βNTI")+xlab("")
#p+ facet_wrap(~ group, scales="free")
p4
ggsave('βNTI_fungi.pdf', p4, width = 10, height = 6)
ggsave('βNTI_fungi.png', p4, width = 8, height = 6)

p5<-ggarrange(p3, p4)
p5
ggsave('βNTI汇总.pdf', p5, width = 10, height = 6)
ggsave('βNTI汇总.png', p5, width = 10, height = 6)
p6<-ggarrange(p1, p2)
p6
ggsave('随机过程和确定过程.pdf', p6, width = 10, height = 6)
ggsave('随机过程和确定过程.png', p6, width = 10, height = 6)



setwd("D:/703/wwy/0C-李飞数据/分长期人工草地/零模型")
data<-read.csv("bNTI细菌+真菌.csv",row.names = 1,header = T)
dt = factor(data$group, levels=c('ASM','AKM','SDM','ARS','ARL'))
library(ggplot2)
#windowsFonts(TNM = windowsFont("Times New Roman"))
p<-ggplot(data, aes(x=dt,y=bNTI))+
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ 
  geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(shape=16, position = position_jitter(0.2))+
  scale_color_manual(values=c("black","black","black","black","black","black"))+
  geom_hline(aes(yintercept=2), color="#AA0000",linetype="dashed")+
  geom_hline(aes(yintercept=-2),color="#AA0000",linetype="dashed")+
  theme_bw(base_size = 18)+
  theme(panel.grid=element_blank())+
  theme(axis.text = element_text(colour = 'black'),
        legend.position = "")+
  ylab("βNTI")+xlab("")+
  facet_wrap(~ type, scales="free")
p
ggsave('βNTI汇总.pdf', p, width = 8, height = 5)
ggsave('βNTI汇总.png', p, width = 8, height = 5)
