
###Functions
fisher.fun = function(datt){
  if (is.null(colnames(datt))){
    colnames(datt) = paste0("var",1:ncol(datt))
  }
  x2all = NULL
  xx = as.factor(datt[,1])
  for (jj in 2:ncol(datt)){        
    yy = as.factor(datt[,jj])
    cnames = colnames(datt)
    namey = cnames[jj]
    namex = cnames[1]
    
    out0 = fisher.test(xx,yy)
    x21 = matrix(out0$p.value);
    colnames(x21) = paste0(namex,"__with__",namey)
    nlevy = length(levels(yy))
    nlevx = length(levels(xx))
    
    if ((nlevx==2) & (nlevy==2)){
      x2all = cbind(x2all,x21)
    } else {
      mtable<-table(xx, yy)
      out1 = fisher.multcomp(mtable)
      x1 = out1$p.value
      x2 = matrix(x1,nrow=1)
      if (nlevy ==2 ){
        colnames(x2) = paste0(paste0("multcomp-",namex,":",rep(rownames(x1),ncol(x1)),":",rep(colnames(x1),each=nrow(x1))),"__with__",paste0(namey,":",paste0(noquote(levels(yy)), collapse = ":")))
      } else {
        colnames(x2) = paste0(paste0("multcomp-",namex,":",rep(rownames(x1),ncol(x1))),"__with__",paste0(namey,":",rep(colnames(x1),each=nrow(x1)))) }
      x2all = cbind(x2all,x21,x2)
    }
  }
  x2all = x2all[,is.na(x2all[1,])==FALSE]
  return(x2all)
}

kd.fun = function(xx,yymat,mtd){
  xx = as.factor(xx)
  ngns = ncol(yymat)
  if ( is.null(colnames(yymat)) ==TRUE ) {
    ngnsn = paste0("Column_",1:ngns)
  } else {
    ngnsn = colnames(yymat)
  }
  
  fun.tmp = function(xx,yy,mtd){
    tst1 = kruskal.test(yy~factor(xx))
    pval1 = tst1$p.value
    tst2 = posthoc.kruskal.dunn.test(x=yy, 
                                     g=xx,
                                     p.adjust.method=mtd)
    
    pval2= tst2$p.value
    pval.v2.tmp = c(pval2)
    pval.v.tmp = c(pval1,pval.v2.tmp[is.na(pval.v2.tmp)==FALSE])
    pval.v = matrix(pval.v.tmp,nrow=1)
    ngrps = ncol(pval2)
    gnamesc = colnames(pval2)
    gnamesr = rownames(pval2)
    ind.use = which(is.na(pval2)==FALSE,arr.ind=TRUE, useNames =FALSE)
    tmp1 = gnamesr[ind.use[,1]]
    tmp2 = gnamesc[ind.use[,2]]
    colnames(pval.v) = c("kruskal.test_pval",paste0("Dunn: ",tmp2," v.s. ",tmp1))
    return(pval.v) 
  }
  pval.mat=fun.tmp(xx,yymat[,1],mtd)
  if (ngns>1){
    for (jj in 2:ngns){
      yy=yymat[,jj]
      pval.mat1=fun.tmp(xx,yy,mtd)
      pval.mat=rbind(pval.mat,pval.mat1)
    }
  }
  rownames(pval.mat)=ngnsn
  return(pval.mat)
}



####Figure 1####

Species<-data[c(1,12,948:1465,1539:1595)]
vec.spec=apply(Species[,-c(1:2)]>0,2,mean)
hist(vec.spec)
vec.spec2=vec.spec
vec.spec2[which(vec.spec<=0.2)]="Shared with equal or less than 20% of the samples"
vec.spec2[which((vec.spec<=0.5)&(vec.spec>0.2))]="Shared with 20-50% of the samples"
vec.spec2[which(vec.spec>=0.5)]="Shared with at least 50% of the samples"
table(vec.spec2)
shared<-as.data.frame(vec.spec2)
mat.spec=NULL
for (jj in 1:nrow(Species)){
  out1=vec.spec2[which(Species[jj,-c(1:2)]>0)] 
  out2=c(sum(out1=="Shared with equal or less than 20% of the samples"),sum(out1=="Shared with 20-50% of the samples"),sum(out1=="Shared with at least 50% of the samples"))
  mat.spec=rbind(mat.spec,out2)
}
Tot=apply(mat.spec,1,sum)
hist(Tot)
mat.spec2=mat.spec[order(Tot),]
Species2=Species[order(Tot),]

dat=data.frame(ID=rep(1:nrow(Species2),3),ID2=rep(tmp7,3),Type=rep(Species2[,2],3),
               Value=c(mat.spec2[,1],mat.spec2[,2],mat.spec2[,3]),
               Type2=rep(c("Shared with equal or less than 20% of the samples","Shared with 20-50% of the samples","Shared with at least 50% of the samples"),each=nrow(Species2)) )

dat$Type3=factor(dat$Type2,levels=c("Shared with equal or less than 20% of the samples","Shared with 20-50% of the samples","Shared with at least 50% of the samples"))

ggplot(dat, aes(fill=Type3, y=Value, x=-ID)) + 
  geom_bar(position="stack", stat="identity", color="black") +
  theme_test()+
  scale_fill_manual(values=c( "indianred2","paleturquoise3","pink4"))+
  scale_y_continuous(name="Number of species")+
  theme(text = element_text(colour = "black", size=34),
        axis.text=element_text(colour="black", size=34),
        axis.title.y = element_text(vjust=5),
        panel.border = element_rect(size=2.4),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 20),
        strip.background = element_rect(size=2),
        axis.title.x = element_blank(),
        legend.justification = "top",
        legend.text = element_text(  size=26),
        axis.ticks = element_line(color = "black", size=1.2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1.5,"cm"),
        plot.title = element_blank())+
  theme(plot.margin=unit(c(1,0.5,1.5,1),"cm"))

dat <- dat[order(dat$Value) , ]

ggplot(dat, aes(fill=Type3, y=Value, x=-ID)) + 
  geom_bar(position="fill", stat="identity", color="black") +
  theme_test()+
  scale_fill_manual(values=c( "indianred2","paleturquoise3","pink4"))+
  scale_y_continuous(name="Proportion")+
  theme(text = element_text(colour = "black", size=34),
        axis.text=element_text(colour="black", size=34),
        axis.title.y = element_text(vjust=5),
        panel.border = element_rect(size=2.4),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 20),
        strip.background = element_rect(size=2),
        axis.title.x = element_blank(),
        legend.justification = "top",
        legend.text = element_text(  size=26),
        axis.ticks = element_line(color = "black", size=1.2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1.5,"cm"),
        plot.title = element_blank())+
  theme(plot.margin=unit(c(1,0.5,1.5,1),"cm"))


####Figure 2####
################

####Figure 2A###


b2<-ggplot(data)+
  aes(data$Group,  data$InvSimpsonVirome)+
  geom_boxplot(width=0.6,lwd=1,outlier.colour = NA,aes(fill=Group))+
  scale_x_discrete(limits=c("Control", "PBC", "NAFLD1", "NAFLD2"))+
  scale_y_continuous(name="Inverse Simpson index")+
  scale_fill_manual(limits=c("Control", "PBC", "NAFLD1", "NAFLD2"),
                    values=c("gray93" ,"lightblue3", "rosybrown3","indianred3"))+
  geom_jitter(width = 0.2, alpha=0.5,size=3, shape=21,colour="black", 
              aes(fill=Group))+
  theme_test()+
  ggtitle(label="Viral diversity")+
  theme(text = element_text(colour = "black", size=30),
        axis.text=element_text(colour="black", size=30),
        axis.title.y = element_text(size=30, vjust=5),
        axis.title.x = element_blank(),
        panel.border = element_rect(size=2),
        axis.text.x = element_blank(),
        plot.title = element_text(colour="black", size=34,
                                  hjust=0.5,vjust=5),
        legend.justification = "top",
        legend.text = element_text(  size=28),
        axis.ticks = element_line(color = "black", size=1.2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.8, "cm"),
        legend.key.width = unit(1,"cm"))+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
b2
hist(data$InvSimpsonVirome)
plant.lm <- lm(data$InvSimpsonVirome~ data$Main_Diagnosis, data = data)
plant.av <- aov(plant.lm)
summary(plant.av)
tukey.test <- TukeyHSD(plant.av)
tukey.test


####Figure 2B###

Taxa<-data[c(1,12,878:880)]
Taxa_mean <- aggregate(Taxa[,3:5], 
                       by=list(Group=Taxa$Group), FUN=mean)
Taxa_mean1 <- melt(Taxa_mean, id.vars=c("Group"), 
                   variable.name = "Taxa", value.name="Relative_Abundance")
names(Taxa_mean1)[2]<-"Taxa"
names(Taxa_mean1)[3]<-"Relative_Abundance"
cp <- coord_polar(theta = "y")
cp$is_free <- function() TRUE
levels(Taxa_mean1$Group) <- c("Control",  "NAS 0-4","NAS 5-8/LCI" , "Mild PBC")
Taxa_mean1$Group<-ordered(Taxa_mean1$Group, levels=c("Control","Mild PBC",  "NAS 0-4","NAS 5-8/LCI" ))


ggplot(Taxa_mean1,aes(Group,Relative_Abundance,fill=Taxa,decreasing=T))+
  geom_bar(position="fill", stat="identity", alpha=1)+
  scale_y_continuous(name="Relative Abundance")+
  scale_x_discrete(name="")+
  theme_test()+
  coord_polar("y")+
  cp +
  facet_wrap(~Group, scales = "free",nrow = 3) +
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = c("steelblue4",
                               "lightgoldenrod1",
                               "seagreen4"),
                    breaks=c("Bacteriophages", "Mammalian",
                             "Plant"),
                    labels=c("Bacteriophages", "Mammalian viruses", "Other viruses, including plant/food derived viruses"),
                    guide=guide_legend(ncol =1))+
  theme(text = element_text(colour = "black", size=20),
        axis.title.y = element_blank(),
        panel.border = element_rect(size=2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        axis.title.x = element_blank(),
        legend.justification = "top",
        strip.text.x = element_text(size = 24),
        strip.background = element_rect(size=2),
        legend.title = element_blank(),
        legend.text=element_text(size = 22),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1.2,"cm"),
        plot.title = element_blank())+
  theme(plot.margin=unit(c(0,0,0,0),"cm"))

hist(data$Bacteriophages)
hist(data$Mammalian)
hist(data$Plant)
xx=factor(data$Group)
ks.table<-kd.fun(xx,data[,878:880],"fdr")
ks.table<-as.data.frame(ks.table)
ks.table$AdjPvalue<-p.adjust(ks.table$kruskal.test_pval ,method="fdr")


####Figure 2C####

names(data1)
data2<-data1[c(12,23,24,30,31,33,36,76,78,79,80:83,85,91,96,116,118,858,859,861,1694:2211)]
set.seed(1234)
rf <- randomForest(data2$Group ~ ., data=data2,importance=TRUE,proximity=TRUE,na.action = na.roughfix)
varImpPlot(rf)
importance<-as.data.frame(rf$importance)
importance<-importance[order(-importance$MeanDecreaseGini),]
importance1<-importance[c(1:40),]
row.names(importance1)
print(importance1$Var1)
importance1<-as.data.frame(importance1)
importance1<- importance1 %>% rownames_to_column("Var1")
write.csv(importance1,file= "importanceRF.csv")

importance1_Group1 <- read_excel("~/importanceRF_edit.xlsx")

ggplot(data = importance1_Group1,
       aes(x =  reorder(Var1,MeanDecreaseGini),
           y =MeanDecreaseGini))+
  geom_bar(stat = "identity", color="black", fill="gray80")+
  coord_flip(clip = "off")+
  
  scale_y_continuous(name="Mean decrease in Gini",breaks=seq(0,5, by = 1))+
  theme_test()+
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1.5),
        panel.background = element_rect(),
        axis.text.x = element_text(size=34, colour="black", vjust=-6),
        axis.title.x = element_text(size=34, colour="black", vjust=-18),
        axis.text.y = element_text(size=26, colour="black"),
        axis.title.y = element_blank(),
        plot.title = element_text(size=30, vjust = 10, hjust = 0.5),
        axis.line.y = element_line(color = "black", size=0.9),
        axis.line.x = element_line(color = "black", size=0.9),
        axis.ticks = element_line(color = "black", size=0.9),
        axis.ticks.length = unit(0.4, "cm"),
        strip.text.y = element_text(hjust=0,vjust = 0.5,angle=180,face="bold", size=30, 
                                    colour="black" ),
        plot.margin=unit(c(2,2,2,2),"cm"))


####Figure 2D####

row.names(importance1)
RFtaxa<-data1[c("Group1",                                 
                "Species_phages_Lactococcus_phage_949_p" ,      
                "Species_phages_Lactococcus_phage_50101_p",     
                "Species_phages_Lactococcus_phage_BK5_T_p" ,    
                "Species_phages_Lactococcus_phage_Tuc2009_p"  , 
                "Species_phages_Lactococcus_phage_phiLC3_p"    ,
                "Species_mammalian_Parvovirus_NIH_CQV_p"       ,
                "Species_phages_Lactococcus_phage_63301_p"    , 
                "Species_phages_Lactococcus_phage_D4412_p"     ,
                "Species_phages_Lactococcus_phage_TP901_1_p"   ,
                "Species_phages_unclassified_Punalikevirus_p"  ,
                "Species_phages_Streptococcus_phage_TP_778L_p", 
                "Species_phages_Lactococcus_phage_BM13_p"      ,
                "Species_phages_Lactococcus_phage_M5938_p"     ,
                "Species_phages_Streptococcus_phage_20617_p"   ,
                "Species_mammalian_Avian_gyrovirus_2_p"        ,
                "Species_phages_Streptococcus_phage_TP_J34_p"  ,
                "Species_phages_Salmonella_phage_SJ46_p"       ,
                "Species_phages_Streptococcus_phage_YMC_2011_p",
                "Species_phages_Lactococcus_phage_phiL47_p"    ,
                "Species_phages_Lactococcus_phage_bIL67_p"    )]

Species<-RFtaxa[2:ncol(RFtaxa)]
clinical<-data1[c(1,12)]
comb<-cbind(clinical,Species)
comb2 <- gather(data =comb, key = Class, value = Abundance,-c(1:2))
ggplot(data =comb2, mapping = aes(x = ID,
                                  y = Class,
                                  fill = factor(Abundance), decreasing=T)) +
  geom_tile( colour = "white", size = 0.1) +
  facet_grid(.~Group1, scales = "free_x",space = "free_x")+
  scale_fill_manual(values=c("dodgerblue4", "firebrick3"),labels=c("absent", "present"), name="")+
  theme_test()+
  theme(text = element_text(colour = "black", size=25),
        axis.title = element_blank(),
        axis.text.y = element_text(colour="black", size=16),
        axis.text.x = element_blank(),
        panel.border = element_rect(size=2),
        axis.title.x = element_text(vjust=-5, face="bold"),
        legend.justification = "top",
        legend.text = element_text(size=22),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1,"cm"),
        strip.text.x = element_text(size = 20),
        strip.background = element_rect(size=2),
        plot.title = element_blank())+
  theme(plot.margin=unit(c(1,0.5,1.5,1),"cm"))

tmp1=fisher.fun(RFtaxa)
tmp1<-as.data.frame(tmp1)
summary(glm(data1$Group1~PPI+data1$Species_mammalian_Avian_gyrovirus_2_p, family="binomial", data=data1))
###Same for other taxa


####Figure 3####
################

data<-data1[c(1,12,1466:1498)]
data<-data[, c(1:2, 2+which(colSums(data[3:ncol(data)]) >0))]
names(data)
Spec<-data[3:27]
Spec<-Spec[,order(colSums(Spec),decreasing=TRUE)]
N<-10
taxa_list<-colnames(Spec)[1:N]
N<-length(taxa_list)
new_x<-data.frame(Spec[,colnames(Spec) %in% taxa_list],Others=rowSums(Spec[,!colnames(Spec) %in% taxa_list]))
mean1=apply(new_x,2,mean)
order1=order(mean1,decreasing=T)
new_x=new_x[,order1]
colnames(new_x)
new_x<-data.frame(Spec[,colnames(Spec) %in% taxa_list],Others=rowSums(Spec[,!colnames(Spec) %in% taxa_list]))
ls1=c( "Escherichia_phages",
       "Enterobacteria_phages",
       "Leuconostoc_phages", 
       "Lactobacillus_phages",
       "Streptococcus_phages",
       "Pseudomonas_phages", 
       "Shigella_phages")
tmp=colnames(new_x)
ind1=tmp%in%ls1
order1=c(which(tmp==ls1[1]),which(tmp==ls1[2]),which(tmp==ls1[3]),which(tmp==ls1[4]),which(tmp==ls1[5]),which(tmp==ls1[6]),which(tmp==ls1[7]),which(ind1==FALSE) )
new_x=new_x[,order1]
colnames(new_x)
tmp2=colnames(new_x)
tmp3=factor(tmp2,levels=tmp2)
order2=order(new_x[,1],new_x[,2],new_x[,3],new_x[,4],new_x[,5],new_x[,6],new_x[,7],new_x[,8],new_x[,9],new_x[,10],decreasing=T)
new_x=new_x[order2,]
tmp4=data1$Group1[order2] 

df2 = data.frame(Taxa=rep(tmp3,each=nrow(new_x)),
                 Value=unlist(new_x),Type=rep(1:nrow(new_x), ncol(new_x)), 
                 Type2=rep(tmp4,ncol(new_x))
)

tmp5=c()
tmp5[tmp4=="NAS04"]=1:sum(tmp4=="NAS04")
tmp5[tmp4=="NAS58LCI"]=1:sum(tmp4=="NAS58LCI")
df4 = data.frame(Taxa=rep(tmp3,each=nrow(new_x)),
                 Value=unlist(new_x),Type=rep(1:nrow(new_x), ncol(new_x)),
                 Type2=rep(tmp4,ncol(new_x)),
                 Type3=rep(tmp5,ncol(new_x))
)

levels(df4$Type2) <- c("NAS 0-4", "NAS 5-8/LCI")
df4$Type3Num<-as.numeric(df4$Type3)
ggplot(df4,aes(Type3Num,Value,fill=Taxa))+
  geom_area() +
  theme_test()+
  facet_grid(.~Type2, scales = "free_x",space = "free_x")+
  scale_y_continuous(name="Proportion")+
  theme(text = element_text(colour = "black", size=34),
        axis.text=element_text(colour="black", size=34),
        axis.title.y = element_text(vjust=5),
        panel.border = element_rect(size=2.4),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.justification = "top",
        legend.text = element_text(  size=26),
        axis.ticks = element_line(color = "black", size=1.2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1.5,"cm"),
        strip.background = element_rect(size=2),
        plot.title = element_blank())+
  theme(plot.margin=unit(c(1,0.5,1.5,1),"cm"))+
  scale_fill_manual(values = c( "orange",
                                "seagreen4",
                                "springgreen1",
                                "tomato3",
                                "yellow3",
                                "palevioletred3",
                                "lightgoldenrod1",
                                "steelblue4",
                                "steelblue3",
                                "mediumpurple4",
                                "cyan4"),
                    guide=guide_legend(ncol =1))


####Figure 4####
################

###Figure 4A


Taxa<-data[c(1,7,878:880)]
Taxa_mean <- aggregate(Taxa[,3:5], 
                       by=list(Group=Taxa$Group2), FUN=mean)
Taxa_mean1 <- melt(Taxa_mean, id.vars=c("Group2"), 
                   variable.name = "Taxa", value.name="Relative_Abundance")
names(Taxa_mean1)[2]<-"Taxa"
names(Taxa_mean1)[3]<-"Relative_Abundance"
cp <- coord_polar(theta = "y")
cp$is_free <- function() TRUE

ggplot(Taxa_mean1,aes(Group2,Relative_Abundance,fill=Taxa,decreasing=T))+
  geom_bar(position="fill", stat="identity", alpha=1)+
  scale_y_continuous(name="Relative Abundance")+
  scale_x_discrete(name="")+
  theme_test()+
  coord_polar("y")+
  cp +
  facet_wrap(~Group2, scales = "free",nrow = 2) +
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = c("steelblue4",
                               "lightgoldenrod1",
                               "seagreen4"),
                    breaks=c("Bacteriophages", "Mammalian",
                             "Plant"),
                    labels=c("Bacteriophages", "Mammalian viruses", "Other viruses, including plant/food derived viruses"),
                    guide=guide_legend(ncol =1))+
  theme(text = element_text(colour = "black", size=20),
        axis.title.y = element_blank(),
        panel.border = element_rect(size=2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        axis.title.x = element_blank(),
        legend.justification = "top",
        strip.text.x = element_text(size = 24),
        strip.background = element_rect(size=2),
        legend.title = element_blank(),
        legend.text=element_text(size = 22),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1.2,"cm"),
        plot.title = element_blank())+
  theme(plot.margin=unit(c(0,0,0,0),"cm"))

xx=factor(data1$Group2)
ks.table<-kd.fun(xx,data1[,878:880],"fdr")
ks.table<-as.data.frame(ks.table)
ks.table$AdjPvalue<-p.adjust(ks.table$kruskal.test_pval ,method="fdr")
View(ks.table)


####Figure 4B####

a1<-ggplot(data1)+
  aes(factor(data1$Group2), data1$ShannonVirome)+
  geom_boxplot(width=0.6,lwd=1,outlier.colour = NA,aes(fill=Group2))+
  scale_fill_manual(values=c( "lemonchiffon2","lightgoldenrod3"))+
  geom_jitter(width = 0.2, alpha=0.7, size=3, shape=21,colour="black", 
              aes(fill=Group2))+
  theme_test()+
  stat_compare_means()+
  ggtitle(label="Shannon index")+
  theme(text = element_text(colour = "black", size=30),
        axis.text=element_text(colour="black", size=30),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(size=3),
        axis.text.x = element_text( size=30),
        plot.title = element_text(colour="black", size=34,
                                  hjust=0.5,vjust=5),
        legend.justification = "top",
        legend.text = element_text(  size=24),
        axis.ticks = element_line(color = "black", size=1.2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.8, "cm"),
        legend.key.width = unit(1,"cm"))+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

a1

a2<-ggplot(data1)+
  aes(data1$Group2, data1$InvSimpsonVirome)+
  geom_boxplot(width=0.6,lwd=1,outlier.colour = NA,aes(fill=Group2))+
  scale_fill_manual(values=c( "lemonchiffon2","lightgoldenrod3"))+
  geom_jitter(width = 0.2, alpha=0.7, size=3, shape=21,colour="black", 
              aes(fill=Group2))+
  theme_test()+
  stat_compare_means()+
  ggtitle(label="Inverse Simpson index")+
  theme(text = element_text(colour = "black", size=30),
        axis.text=element_text(colour="black", size=30),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(size=3),
        axis.text.x = element_text( size=30),
        plot.title = element_text(colour="black", size=34,
                                  hjust=0.5,vjust=5),
        legend.justification = "top",
        legend.text = element_text(  size=24),
        axis.ticks = element_line(color = "black", size=1.2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.8, "cm"),
        legend.key.width = unit(1,"cm"))+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

a2
ggarrange(a1,a2, ncol=2)



####Figure 4C#### 

data<-data1[c(7,23,24,30,31,33,36,76,78,79,80:83,85,91,96,116,118,858,859,861,1694:2211)]
data$Group2<-data$Fibrosis_F2_F4_LCI
droplevels(data)
df1 <- droplevels(data[!grepl("^\\s*$", data$Fibrosis_F2_F4_LCI),,drop=FALSE] )
data<-df1

set.seed(12212)
rff2 <- randomForest(data$Group2 ~ .,data=data,importance=TRUE,proximity=TRUE,na.action = na.roughfix)
varImpPlot(rff2)
importancefib2<-as.data.frame(rff2$importance)
importancefib2<-importancefib2[order(-importancefib2$MeanDecreaseGini),]
importancefib12<-importancefib2[c(1:40),]
write.csv(importancefib12, file="importancefib12.csv")
importancefib12_edit <- read_excel("~/importancefib12_edit.xlsx")

ggplot(data = importancefib12New,
       aes(x =  reorder(Var1,MeanDecreaseGini),
           y =MeanDecreaseGini))+
  geom_bar(stat = "identity", color="black", fill="gray80")+
  coord_flip(clip = "off")+
  scale_y_continuous(name="Mean decrease in Gini",breaks=seq(0,5, by = 1))+
  theme_test()+
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1.5),
        panel.background = element_rect(),
        axis.text.x = element_text(size=34, colour="black", vjust=-6),
        axis.title.x = element_text(size=34, colour="black", vjust=-18),
        axis.text.y = element_text(size=26, colour="black"),
        axis.title.y = element_blank(),
        plot.title = element_text(size=30, vjust = 10, hjust = 0.5),
        axis.line.y = element_line(color = "black", size=0.9),
        axis.line.x = element_line(color = "black", size=0.9),
        axis.ticks = element_line(color = "black", size=0.9),
        axis.ticks.length = unit(0.4, "cm"),
        strip.text.y = element_text(hjust=0,vjust = 0.5,angle=180,face="bold", size=30, 
                                    colour="black" ),
        plot.margin=unit(c(2,2,2,2),"cm"))


#####Figure 4D####


row.names(importancefib12)
RFtaxafib2<-data1[c("Group2",   
                    "Species_phages_Lactococcus_phage_Tuc2009_p" ,  
                    "Species_phages_Lactococcus_phage_phiLC3_p",    
                    "Species_phages_Lactococcus_phage_KSY1_p"    ,   
                    "Species_phages_Lactococcus_phage_949_p"    ,   
                    "Species_phages_Leuconostoc_phage_phiLN12_p",  
                    "Species_phages_Lactococcus_phage_TP901_1_p"  , 
                    "Species_phages_Lactococcus_phage_63301_p" ,                    
                    "Species_phages_Lactococcus_phage_ul36_p"  ,    
                    "Species_phages_Lactococcus_phage_50101_p"  ,   
                    "Species_phages_Leuconostoc_phage_phiLN04_p" ,  
                    "Species_phages_Lactococcus_phage_r1t_p"     ,  
                    "Species_phages_Bacteroides_phage_B40_8_p"  ,   
                    "Species_phages_Lactobacillus_phage_phiAT3_p"  ,
                    "Species_phages_Streptococcus_phage_YMC_2011_p",
                    "Species_phages_Leuconostoc_phage_phiLN03_p" ,  
                    "Species_phages_Lactococcus_phage_BK5_T_p"   ,  
                    "Species_phages_Streptococcus_phage_5093_p"  ,  
                    "Species_phages_Lactobacillus_phage_J_1_p"  ,   
                    "Species_phages_Lactococcus_phage_M6165_p" )]


spec<-RFtaxafib2[2:ncol(RFtaxafib2)]
for(j in colnames(spec)[1:ncol(spec)]){
  set(spec, i=NULL, j=j, value= as.numeric(spec[[j]]))
}

clinical<-data1[c(1,7)]
comb<-cbind(clinical,spec)
spec.tab<- gather(data = comb, key = Class, value = Abundance,-c(1:2))

ggplot(data = spec.tab, mapping = aes(x = ID,
                                      y = Class,
                                      fill = factor(Abundance), decreasing=T)) +
  geom_tile( colour = "white", size = 0.1) +
  facet_grid(.~spec.tab$Group2, scales = "free_x",space = "free_x")+
  scale_fill_manual(values=c("dodgerblue4", "firebrick3"),labels=c("absent", "present"), name="")+
  theme_test()+
  theme(text = element_text(colour = "black", size=25),
        axis.title = element_blank(),
        axis.text.y = element_text(colour="black", size=16),
        axis.text.x = element_blank(),
        panel.border = element_rect(size=2),
        axis.title.x = element_text(vjust=-5, face="bold"),
        legend.justification = "top",
        legend.text = element_text(size=22),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1,"cm"),
        
        strip.text.x = element_text(size = 20),
        strip.background = element_rect(size=2),
        plot.title = element_blank())+
  theme(plot.margin=unit(c(1,0.5,1.5,1),"cm"))

tmp1=fisher.fun(RFtaxa)
tmp1<-as.data.frame(tmp1)
summary(glm(data1$Group2~PPI+data1$Species_phages_Bacteroides_phage_B40_8_p, family="binomial", data=data1))
#same for other taxa


#####Figure 5#####
##################

###Figure 5A######


datnew1<-data1[c("Group1", "Age", "AST_U..l")]
test1<-glm(Group1~.,data=datnew1,family="binomial") 
summary(test1)
pred1<-predict(test1, type="response")
roc1<-roc(data1$Group1~pred1,plot=TRUE,grid=TRUE, print.auc=TRUE)
cv.glm(datnew1,test1)$delta 

datnew1<-data1[c("Group1",
                 "Age", 
                 "AST_U..l",
                 "InvSimpson_16S")]
test2<-glm(Group1~.,data=datnew1, family="binomial")
summary(test2)
pred2<-predict(test2, type="response")
roc2<-roc(data1$Group1~pred2,plot=TRUE,  grid=TRUE, print.auc=TRUE)
cv.glm(datnew1,test2)$delta

datnew1<-data1[c("Group1",
                 "Age",
                 "AST_U..l",
                 "InvSimpsonVirome")]

test3<-glm(Group1~.,data=  datnew1, family="binomial")
summary(test3)
pred3<-predict(test3, type="response")
roc3<-roc(data1$Group1~pred3,plot=TRUE,  grid=TRUE, print.auc=TRUE)
cv.glm(datnew1,test3)$delta 

d <- lrm(Group1 ~ Age+AST_U..l, data=data1)     
e <- lrm(Group1 ~ Age+AST_U..l+InvSimpsonVirome, data=data1)
lrtest(d,e)

data1$pred1<-predict(test1, type="response")
data1$pred2<-predict(test2, type="response")
data1$pred3<-predict(test3, type="response")

ggroct<-ggplot(data1, aes(d =data1$Group1 )) + 
  geom_roc(data=data1, aes(d = data1$Group1, m = pred1),
           n.cuts = 0, size=1.7, colour="gray54", linetype=1)+
  geom_roc(data=data1, aes(d = data1$Group1, m = pred2),
           n.cuts = 0, size=1.7, colour="gray29", linetype=1)+
  geom_roc(data=data1, aes(d = data1$Group1, m = pred3),
           n.cuts = 0, size=1.7, colour="dodgerblue3", linetype=1)+
  theme_test()+   
  scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .2))+
  scale_y_continuous("Sensitivity", breaks = seq(0, 1, by = .2))+
  theme(text = element_text(colour = "black", size=32),
        axis.text=element_text(colour="black", size=32),
        axis.title.y = element_text(vjust=5),
        axis.title.x = element_text(vjust=-5),
        panel.border = element_rect(size=2.5),
        axis.text.x = element_text( size=32),
        plot.title = element_text(colour="black", size=32,
                                  hjust=0.5,vjust=5),
        legend.justification = "top",
        legend.text = element_text(  size=24),
        axis.ticks = element_line(color = "black", size=1.2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.8, "cm"),
        legend.key.width = unit(1,"cm"))+
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50")+
  theme(plot.margin=unit(c(1,1,1,1),"cm")) 

ggroct


###Figure 5B####


datnew1<-data1[c("Group2","Age", 
                 "AST_U..l",
                 "Platelets_x1E9..l")]

test1fib<-glm(Group2~.,data=  datnew1, family="binomial")
summary(test1fib)
pred1fib<-predict(test1fib, type="response")
rocfib1<-roc(data1$Group2~pred1fib,plot=TRUE,  grid=TRUE, print.auc=TRUE)
cv.glm(datnew1,test1fib)$delta 

datnew1<-data1[c("Group2",
                 "Age", 
                 "AST_U..l",
                 "Platelets_x1E9..l",
                 "InvSimpson_16S"  )]

test2fib<-glm(Group2~.,data=datnew1,family="binomial")
summary(test2fib)
pred2fib<-predict(test2fib, type="response")
rocfib2<-roc(data1$Group2~pred2fib,plot=TRUE,grid=TRUE, print.auc=TRUE)
cv.glm(datnew1,test2fib)$delta 


datnew1<-data1[c("Group2",
                 "Age", 
                 "AST_U..l",
                 "Platelets_x1E9..l",
                 "InvSimpsonVirome")]

test3fib<-glm(Group2~.,data=datnew1,family="binomial")  
summary(test3fib)
pred3fib<-predict(test3fib, type="response")
rocfib3<-roc(data1$Fibrosis_F2_F4_LCI~pred3fib,plot=TRUE,  grid=TRUE, print.auc=TRUE)
cv.glm(datnew1,test3fib)$delta

a <- lrm(Group2 ~ Age+AST_U..l+Platelets_x1E9..l, data=data1)     
b <- lrm(Group2 ~ Age+AST_U..l+Platelets_x1E9..l+InvSimpsonVirome, data=data1)
lrtest(a,b)

data1$pred1fib<-predict(test1fib, type="response")
data1$pred2fib<-predict(test2fib, type="response")
data1$pred3fib<-predict(test3fib, type="response")

ggroct2<-ggplot(data1) + 
  geom_roc(data=data1, aes(d =data1$Group2, m = pred1fib),
           n.cuts = 0, size=1.7, colour="gray54", linetype=1)+
  geom_roc(data=data1, aes(d =data1$Group2, m = pred2fib),
           n.cuts = 0, size=1.7, colour="gray29", linetype=1)+
  geom_roc(data=data1, aes(d = data1$Group2, m = data1$pred3fib),
           n.cuts = 0, size=1.7, colour="dodgerblue3", linetype=1)+
  theme_test()+  
  scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .2))+
  scale_y_continuous("Sensitivity", breaks = seq(0, 1, by = .2))+
  theme(text = element_text(colour = "black", size=32),
        axis.text=element_text(colour="black", size=32),
        axis.title.y = element_text(vjust=5),
        axis.title.x = element_text(vjust=-5),
        panel.border = element_rect(size=2.5),
        axis.text.x = element_text( size=32),
        plot.title = element_text(colour="black", size=32,
                                  hjust=0.5,vjust=5),
        legend.justification = "top",
        legend.text = element_text(  size=24),
        axis.ticks = element_line(color = "black", size=1.2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.8, "cm"),
        legend.key.width = unit(1,"cm"))+
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50")+
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

ggroct2


#####Figure 5C####

t1<-ggbarplot( data1, x = "LBP_Steatosis_Grade", y = "pred3",
               add = "mean_se", color="black", fill="LBP_Steatosis_Grade", size=1,
               add.params = list(size = 1)) +
  coord_cartesian(clip = "off")+
  scale_fill_brewer(palette = "Blues")+
  scale_y_continuous(name="Model (viral diversity + AST + Age)")+
  geom_jitter(width = 0.2, alpha=0.5, size=3, shape=21,colour="black",
              aes(fill=factor(LBP_Steatosis_Grade)))+
  theme_test()+
  scale_x_discrete(name="Grade of steatosis ")+
  stat_compare_means()+
  theme(text = element_text(colour = "black", size=30),
        axis.text=element_text(colour="black", size=30),
        axis.title.y = element_text( size=30, vjust=5),
        axis.title.x = element_text( size=30, vjust=-5),
        panel.border = element_rect(size=1.2),
        axis.text.x = element_text( size=30),
        axis.text.y = element_text( size=30),
        plot.title = element_text(colour="black", size=30,
                                  hjust=0.5,vjust=5),
        legend.justification = "none",
        legend.text = element_text(  size=24),
        axis.ticks = element_line(color = "black", size=1.2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.8, "cm"),
        legend.key.width = unit(1,"cm"))+
  theme(plot.margin=unit(c(0.5,0.5,1,0.5),"cm"))
t1

t1f<-ggbarplot(data1, x = "LBP_Fibrosis_Stage", y = "pred3fib",
               add = "mean_se", color="black", fill="LBP_Fibrosis_Stage", size=1,
               add.params = list(size = 1)) +
  coord_cartesian(clip = "off")+
  scale_fill_brewer(palette = "Blues")+
  scale_y_continuous(name="Model (viral diversity + AST + Age + Plt)")+
  geom_jitter(width = 0.2, alpha=0.5, size=3, shape=21,colour="black",
              aes(fill=factor(LBP_Fibrosis_Stage)))+
  theme_test()+
  stat_compare_means()+
  scale_x_discrete(name="Stage of fibrosis ")+
  theme(text = element_text(colour = "black", size=30),
        axis.text=element_text(colour="black", size=30),
        axis.title.y = element_text( size=30, vjust=5),
        axis.title.x = element_text( size=30, vjust=-5),
        panel.border = element_rect(size=1.2),
        axis.text.x = element_text( size=30),
        axis.text.y = element_text( size=30),
        plot.title = element_text(colour="black", size=30,
                                  hjust=0.5,vjust=5),
        legend.justification = "none",
        legend.text = element_text(  size=24),
        axis.ticks = element_line(color = "black", size=1.2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.8, "cm"),
        legend.key.width = unit(1,"cm"))+
  theme(plot.margin=unit(c(0.5,0.5,1,0.5),"cm"))
t1f

t2<-ggbarplot(data1, x = "LBP_Inflammation_Grade", y = "pred3",
              add = "mean_se", color="black", fill="LBP_Inflammation_Grade", size=1,
              add.params = list(size = 1)) +
  coord_cartesian(clip = "off")+
  scale_fill_brewer(palette = "Blues")+
  scale_y_continuous(name="Model (viral diversity + AST + Age)")+
  geom_jitter(width = 0.2, alpha=0.5, size=3, shape=21,colour="black",
              aes(fill=factor(LBP_Inflammation_Grade)))+
  theme_test()+
  stat_compare_means()+
  scale_x_discrete(name="Grade of inflammation ")+
  theme(text = element_text(colour = "black", size=30),
        axis.text=element_text(colour="black", size=30),
        axis.title.y = element_text( size=30, vjust=5),
        axis.title.x = element_text( size=30, vjust=-5),
        panel.border = element_rect(size=1.2),
        axis.text.x = element_text( size=30),
        axis.text.y = element_text( size=30),
        plot.title = element_text(colour="black", size=30,
                                  hjust=0.5,vjust=5),
        legend.justification = "none",
        legend.text = element_text(  size=24),
        axis.ticks = element_line(color = "black", size=1.2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.8, "cm"),
        legend.key.width = unit(1,"cm"))+
  theme(plot.margin=unit(c(0.5,0.5,1,0.5),"cm"))
t2

t3<-ggbarplot(data1, x = "LBP_Ballooning", y = "pred3",
              add = "mean_se", color="black", fill="LBP_Ballooning", size=1,
              add.params = list(size = 1)) +
  coord_cartesian(clip = "off")+
  scale_fill_brewer(palette = "Blues")+
  scale_y_continuous(name="Model (viral diversity + AST + Age)")+
  geom_jitter(width = 0.2, alpha=0.5, size=3, shape=21,colour="black",
              aes(fill=factor(LBP_Ballooning)))+
  theme_test()+
  stat_compare_means()+
  scale_x_discrete(name="Grade of ballooning ")+
  theme(text = element_text(colour = "black", size=30),
        axis.text=element_text(colour="black", size=30),
        axis.title.y = element_text( size=30, vjust=5),
        axis.title.x = element_text( size=30, vjust=-5),
        panel.border = element_rect(size=1.2),
        axis.text.x = element_text( size=30),
        axis.text.y = element_text( size=30),
        plot.title = element_text(colour="black", size=30,
                                  hjust=0.5,vjust=5),
        legend.justification = "none",
        legend.text = element_text(  size=24),
        axis.ticks = element_line(color = "black", size=1.2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.8, "cm"),
        legend.key.width = unit(1,"cm"))+
  theme(plot.margin=unit(c(0.5,0.5,1,0.5),"cm"))
t3

ggarrange(t1,t2,t3,t1f,widths = c(1,1.1,1,1.2), ncol=4)

####Figure 6####
################


datab<-data1[c(1,12,1466:1498)]
datab<-datab[, c(1:2, 2+which(colSums(datab[3:ncol(datab)]) >0))]
Generab<-datab[3:27]
Generab<-Generab[,order(colSums(Generab),decreasing=TRUE)]
N<-20
taxa_listb<-colnames(Generab)[1:N]
N<-length(taxa_listb)
new_xb<-data.frame(Generab[,colnames(Generab) %in% taxa_listb],Others=rowSums(Generab[,!colnames(Generab) %in% taxa_listb]))
dfb<-NULL
for (i in 1:dim(new_xb)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_xb),Taxa=rep(colnames(new_xb)[i],dim(new_xb)[1]),Value=new_xb[,i],Type=c(data1$ID))
  if(i==1){dfb<-tmp} else {dfb<-rbind(dfb,tmp)}
}


data<-data1[c(1,12,650:856)]
spec<-data[3:209]
spec<-spec[,order(colSums(spec),decreasing=TRUE)]
N<-200
taxa_list<-colnames(spec)[1:N]
N<-length(taxa_list)
new_x<-data.frame(spec[,colnames(spec) %in% taxa_list],Others=rowSums(spec[,!colnames(spec) %in% taxa_list]))


df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Type=c(data1$ID))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}

dftogether<-rbind(dfb,df)
ggplot(data = dftogether, mapping = aes(x = Type,
                                        y = Taxa,
                                        fill = Value)) +
  geom_tile( colour = "white", size = 0.2) +
  xlab(label = "Sample")+
  theme_test()+
  scale_fill_gradient2(low="darkseagreen3", mid="deepskyblue2", 
                       high="red",trans="log",  
                       na.value = "darkseagreen2",
                       breaks = c(0.002,0.05,1,30), name="Log10 rel. abundance")+
  theme(text = element_text(colour = "black", size=25),
        axis.text=element_text(colour="black", size=15),
        axis.text.y=element_text(size=24, color=c("red3", "black")),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(size=2),
        axis.title.x = element_text(vjust=-5, face="bold"),
        legend.justification = "top",
        legend.text = element_text(size=22),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1,"cm"),
        strip.text.x = element_text(size = 24),
        strip.background = element_rect(size=2),
        plot.title = element_blank())+
  theme(plot.margin=unit(c(1,0.5,1.5,1),"cm"))+
  scale_y_discrete(limits=c(
    "Genus_Cronobacter",
    "Cronobacter_phages",
    "Genus_Brochothrix",
    "Brochothrix_phages",
    "Genus_Varibaculum",
    "Rhodoferax_phages",
    "Genus_Clostridium",
    "Clostridium_phages",
    "Genus_Yersinia",
    "Yersinia_phages",
    "Genus_Klebsiella",
    "Klebsiella_phages",
    "Genus_Enterococcus",
    "Enterococcus_phages",
    "Genus_Bacteroides",
    "Bacteroides_phages",
    "Genus_Staphylococcus",
    "Staphylococcus_phages",
    "Genus_Actinobaculum",
    "Edwardsiella_phages",
    "Genus_Propionibacterium",
    "Propionibacterium_phages",
    "Genus_Pseudomonas",
    "Pseudomonas_phages",
    "Genus_Shigella",
    "Shigella_phages",
    "Genus_Salmonella",
    "Salmonella_phages",
    "Genus_Streptococcus",
    "Streptococcus_phages",
    "Genus_Lactobacillus",
    "Lactobacillus_phages",
    "Genus_Scardovia",
    "Enterobacteria_phages",
    "Genus_Leuconostoc",
    "Leuconostoc_phages",
    "Genus_Escherichia",
    "Escherichia_phages",
    "Genus_Lactococcus",
    "Lactococcus_phages"),
    labels=c(
      bquote(italic("Cronobacter")),
      bquote(italic("Cronobacter")~"phages"),
      bquote(italic("Brochothrix")),
      bquote(italic("Brochothrix")~"phages"),
      bquote(italic("Rhodoferax")),
      bquote(italic("Rhodoferax")~"phages"),          
      bquote(italic("Clostridium")),
      bquote(italic("Clostridium")~"phages"),  
      bquote(italic("Yersinia")),
      bquote(italic("Yersinia")~"phages"),  
      bquote(italic("Klebsiella")),
      bquote(italic("Klebsiella")~"phages"),  
      bquote(italic("Enterococcus")),
      bquote(italic("Enterococcus")~"phages"),  
      bquote(italic("Bacteroides")),
      bquote(italic("Bacteroides")~"phages"),  
      bquote(italic("Staphylococcus")),
      bquote(italic("Staphylococcus")~"phages"),  
      bquote(italic("Edwardsiella")),
      bquote(italic("Edwardsiella")~"phages"),  
      bquote(italic("Propionibacterium")),
      bquote(italic("Propionibacterium")~"phages"), 
      bquote(italic("Pseudomonas")),
      bquote(italic("Pseudomonas")~"phages"),  
      bquote(italic("Shigella")),
      bquote(italic("Shigella")~"phages"),  
      bquote(italic("Salmonella")),
      bquote(italic("Salmonella")~"phages"),  
      bquote(italic( "Streptococcus")),
      bquote(italic("Streptococcus")~"phages"),  
      bquote(italic( "Lactobacillus")),
      bquote(italic("Lactobacillus")~"phages"),  
      bquote(italic("Enterobacteria")),
      bquote(italic("Enterobacteria")~"phages"),
      bquote(italic( "Leuconostoc")),
      bquote(italic("Leuconostoc")~"phages"),  
      bquote(italic( "Escherichia")),
      bquote(italic("Escherichia")~"phages"),  
      bquote(italic( "Lactococcus")),
      bquote(italic("Lactococcus")~"phages")  
    ))


####Supplementary Figure 1####
##############################


heatmapdata<-data1[c(12,7,23,24,30,31,32,33,36,42,47,54,71,62,63,64,65,66,76,78:85,
                     89:91,96:98,116:119,111,101,126:128,132,857:861,863:1595)]
heatmapdata<-heatmapdata[,c(1:49,49+which(colMeans(heatmapdata[50:ncol(heatmapdata)]) >0))]
dat.tmp<-heatmapdata
ext.mat = dat.tmp[,-c(1:49)]>0
nzcts = apply(ext.mat,2,sum)  
sum(nzcts>0)
prop1 = 0.1
data2 = dat.tmp[,c(1:49,which(nzcts>(nrow(ext.mat)*prop1))+49)]
heatmapdata<-data2
names(heatmapdata)
Phages<-heatmapdata[c(65:177)]
clinical<-heatmapdata[c(9,31:33,13,43,38,39,40:42)]
names(heatmapdata)

A <- as.matrix(Phages)
B <- as.matrix(clinical)

corrrr<- rcorr(A, B, type="spearman")

corP<-corrrr$P
corR<-corrrr$r
corN<-corrrr$n

row.names(corR)

C=corP[c(114:124),c(1:113)]

D=corR[c(114:124),c(1:113)]

N=corN[c(114:124),c(1:113)]

cormeltP <- melt(C)
cormeltN <- melt(N)
cormeltR <- melt(D)

cormeltP$value<-round(cormeltP$value, digits=3)
cormeltN$value<-round(cormeltN$value, digits=3)
cormeltR$value<-round(cormeltR$value, digits=3)

combicor1<-cbind(cormeltN, cormeltR, cormeltP)

names(combicor1)[3]<-"number"
names(combicor1)[6]<-"rvalue"
names(combicor1)[9]<-"pvalue"

combicorclear1<-combicor1[c(1:3,6,9)]

newdata1 <- combicorclear1[which(combicorclear1$pvalue<0.05),]

DR1<-newdata1[c(1,2,4)]
CP1<-newdata1[c(1,2,5)]

remeltP1<-dcast(CP1, Var1~Var2)
remeltR1<-dcast(DR1, Var1~Var2)

remeltPP1<-column_to_rownames(remeltP1, var = "Var1")
remeltRR1<-column_to_rownames(remeltR1, var = "Var1")

remeltPP1<-as.matrix(remeltPP1)
remeltRR1<-as.matrix(remeltRR1)

stars1 <-cut(remeltPP1, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

test<-melt(remeltRR1)

ggplot(test, aes(Var1, Var2, fill=value)) +
  ggtitle("Title")+
  geom_tile(colour = "white", size = 0.2) +
  scale_fill_gradient2(low="dodgerblue4", mid="white", high="firebrick3") +
  theme_minimal()+  
  labs(x="",y="",fill="R") + 
  geom_text(aes(label=stars1), color="black", size=8, hjust=0.5, vjust=0.8)+
  theme_test()+
  theme(text = element_text(colour = "black", size=10),
        axis.text=element_text(colour="black", size=10, face="italic"),
        axis.title.y = element_text(vjust=5, face="bold"),
        panel.border = element_rect(size=2),
        axis.text.x = element_text(colour="black", size=10, 
                                   angle=45, vjust=1, hjust=1),
        axis.title.x = element_text(vjust=-5, face="bold"),
        legend.justification = "top",
        legend.text = element_text( size=10),
        axis.ticks = element_line(color = "black", size=1.2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(),
        legend.key.size = unit(0.8, "cm"),
        legend.key.width = unit(1,"cm"),
        plot.title = element_blank())+
  theme(plot.margin=unit(c(1,0.5,1.5,1),"cm"))

summary(test$Var2)

Phages2<-data1[,c( 'Genus_mammalian_Protoparvovirus',
                   'Species_phages_Escherichia_virus_HK97',
                   'Species_phages_Lactococcus_phage_ASCC191',
                   'Species_phages_Bacteroides_phage_B40_8',
                   'Lactobacillus_phages',
                   'Species_phages_Lactococcus_phage_phi7',
                   'Species_phages_Escherichia_phage_TL_2011b',
                   'Species_phages_unclassified_Uetakevirus',
                   'Species_phages_Enterobacteria_phage_BP_4795',
                   'Species_phages_Bacteroides_phage_B124_14',
                   'Species_phages_Enterobacteria_phage_Sf101',
                   'Species_phages_Enterobacteria_phage_mEp460',
                   'Bacteroides_phages',
                   'Species_phages_Lactococcus_phage_M6165',
                   'Species_phages_Leuconostoc_phage_phiLNTR3',
                   'Species_phages_Lactococcus_phage_28201',
                   'Species_phages_Lactococcus_phage_BM13',
                   'Species_phages_Enterobacteria_phage_HK633',
                   'Species_phages_Lactococcus_phage_phi7',
                   'Species_phages_Stx2_converting_phage_1717',
                   'Species_phages_Leuconostoc_phage_phiLN6B',
                   'Species_phages_Leuconostoc_phage_phiLNTR3',
                   'Species_phages_Lactococcus_phage_BK5_T',
                   'Species_phages_Leuconostoc_phage_phiLN25',
                   'Species_phages_Streptococcus_phage_DT1',
                   'Species_phages_Lactococcus_phage_KSY1',
                   'Species_phages_Streptococcus_phage_Abc2',
                   'Species_phages_Lactococcus_phage_50101',
                   'Species_phages_Lactococcus_phage_D4412',
                   'Species_phages_Lactococcus_phage_phiLC3',
                   'Species_phages_Leuconostoc_phage_phiLN12',
                   'Species_phages_Leuconostoc_phage_phiLN03',
                   'Species_phages_Lactococcus_phage_M6165',
                   'Species_phages_Leuconostoc_phage_phiLN25',
                   'Species_phages_Leuconostoc_phage_phiLN6B',
                   'Species_phages_Lactococcus_phage_63301',
                   'Species_phages_Lactococcus_phage_c2',
                   'Leuconostoc_phages',
                   'Species_phages_Leuconostoc_phage_phiLN04',
                   'Species_phages_Lactococcus_phage_bIL67',
                   'Leuconostoc_phages',
                   'Species_phages_Lactococcus_phage_WRP3',
                   'Species_phages_Lactococcus_phage_D4410',
                   'Species_phages_Leuconostoc_phage_Ln_8',
                   'Species_phages_Leuconostoc_phage_phiLN03',
                   'Species_phages_Leuconostoc_phage_phiLN12',
                   'Species_phages_Lactococcus_phage_M5938',
                   'Species_phages_Lactococcus_phage_M6162',
                   'Species_phages_Leuconostoc_phage_Lmd1',
                   'Species_phages_Leuconostoc_phage_phiLN04',
                   'Lactobacillus_phages',
                   'Pseudomonas_phages' )]

clinical2<-heatmapdata[c(9,31:33,13,43,38,39,40:42)]
heatmap_clinical<-cbind(clinical2, Phages2)
write.csv(heatmap_clinical, file="heatmap_clinical2.csv")
heatmap_clinical_import <- read_excel("~/heatmap_clinical2_edit.xlsx")
names(heatmap_clinical_import)
clinical2<-heatmap_clinical_import[c(1:11)]
Phages2<-heatmap_clinical_import[c(12:53)]

A <- as.matrix(Phages2)
B <- as.matrix(clinical2)
corrrr<- rcorr(A, B, type="spearman")

corP<-corrrr$P
corR<-corrrr$r
corN<-corrrr$n

row.names(corR)

C=corP[c(43:53),c(1:42)]

D=corR[c(43:53),c(1:42)]

N=corN[c(43:53),c(1:42)]

test<-melt(D)

stars <- cut(C, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

ggplot(melt(D), aes(Var1, Var2, fill=value)) +
  geom_tile(colour = "white", size = 0.2) +
  scale_fill_gradient2(low="dodgerblue4", mid="white", high="firebrick3") +
  theme_minimal()+  
  labs(x="",y="",fill="R") + 
  geom_text(aes(label=stars), color="black", size=12, hjust=0.5, vjust=0.8)+
  theme_test()+ 
  theme(text = element_text(colour = "black", size=20),
        axis.text=element_text(colour="black", size=20),
        axis.title.y = element_text(vjust=5, face="bold"),
        panel.border = element_rect(size=2),
        axis.text.x = element_text(colour="black", size=20,
                                   angle=45, vjust=1, hjust=1),
        axis.title.x = element_text(vjust=-5, face="bold"),
        legend.justification = "top",
        legend.text = element_text( size=24),
        axis.ticks = element_line(color = "black", size=1.2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(),
        legend.key.size = unit(0.8, "cm"),
        legend.key.width = unit(1,"cm"),
        plot.title = element_blank())+
  theme(plot.margin=unit(c(1,0.5,1.5,1),"cm"))


####Supplementary Figure 2######
################################


levels(NAS58LCI$PPI) <- c("NAS 5-8/LCI, no PPI use","NAS 5-8/LCI, daily PPI use" )

CC<-NAS58LCI[c(1,54,878:880)]
CC_mean <- aggregate(CC[,3:5], 
                     by=list(PPI=CC$PPI), FUN=mean)
CC_mean1 <- melt(CC_mean, id.vars=c("PPI"), 
                 variable.name = "Taxa", value.name="Relative_Abundance")
names(CC_mean1)[2]<-"Taxa"
names(CC_mean1)[3]<-"Relative_Abundance"
cp <- coord_polar(theta = "y")
cp$is_free <- function() TRUE
names(CC_mean1$Taxa)

ggplot(CC_mean1,aes(PPI,Relative_Abundance,fill=Taxa,decreasing=T))+
  geom_bar(position="fill", stat="identity", alpha=1)+
  scale_y_continuous(name="Relative Abundance")+
  scale_x_discrete(name="")+
  theme_test()+
  coord_polar("y")+
  cp +
  facet_wrap(~PPI, scales = "free",nrow = 3) +
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = c("steelblue4",
                               "lightgoldenrod1",
                               "seagreen4"),
                    breaks=c("Bacteriophages", "Mammalian",
                             "Plant"),
                    labels=c("Bacteriophages", "Mammalian viruses", "Other viruses, including plant/food derived viruses"),
                    guide=guide_legend(ncol =1))+
  
  theme(text = element_text(colour = "black", size=20),
        axis.title.y = element_blank(),
        panel.border = element_rect(size=2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        axis.title.x = element_blank(),
        legend.justification = "top",
        strip.text.x = element_text(size = 16),
        strip.background = element_rect(size=2),
        legend.title = element_blank(),
        legend.text=element_text(size = 22),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1.2,"cm"),
        plot.title = element_blank())+
  theme(plot.margin=unit(c(0,0,0,0),"cm"))

xx=factor(NAS58LCI$PPI)
ks.table<-kd.fun(xx,NAS58LCI[,878:880],"fdr")
ks.table<-as.data.frame(ks.table)
View(ks.table)




levels(F2_F4$PPI) <- c("F2-F4, no PPI use","F2-F4, daily PPI use" )

CC<-F2_F4[c(1,54,878:880)]
CC_mean <- aggregate(CC[,3:5], 
                     by=list(PPI=CC$PPI), FUN=mean)
CC_mean1 <- melt(CC_mean, id.vars=c("PPI"), 
                 variable.name = "Taxa", value.name="Relative_Abundance")
names(CC_mean1)[2]<-"Taxa"
names(CC_mean1)[3]<-"Relative_Abundance"
cp <- coord_polar(theta = "y")
cp$is_free <- function() TRUE
names(CC_mean1$Taxa)

ggplot(CC_mean1,aes(PPI,Relative_Abundance,fill=Taxa,decreasing=T))+
  geom_bar(position="fill", stat="identity", alpha=1)+
  scale_y_continuous(name="Relative Abundance")+
  scale_x_discrete(name="")+
  theme_test()+
  coord_polar("y")+
  cp +
  facet_wrap(~PPI, scales = "free",nrow = 3) +
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = c("steelblue4",
                               "lightgoldenrod1",
                               "seagreen4"),
                    breaks=c("Bacteriophages", "Mammalian",
                             "Plant"),
                    labels=c("Bacteriophages", "Mammalian viruses", "Other viruses, including plant/food derived viruses"),
                    guide=guide_legend(ncol =1))+
  
  theme(text = element_text(colour = "black", size=20),
        axis.title.y = element_blank(),
        panel.border = element_rect(size=2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        axis.title.x = element_blank(),
        legend.justification = "top",
        strip.text.x = element_text(size = 16),
        strip.background = element_rect(size=2),
        legend.title = element_blank(),
        legend.text=element_text(size = 22),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1.2,"cm"),
        plot.title = element_blank())+
  theme(plot.margin=unit(c(0,0,0,0),"cm"))

xx=factor(F2_F4$PPI)
ks.table<-kd.fun(xx,F2_F4[,878:880],"fdr")
ks.table<-as.data.frame(ks.table)
View(ks.table)

summary(glm(data1$Group1~PPI+data1$InvSimpsonVirome, family="binomial", data=data1))
summary(glm(data1$Group2~PPI+data1$InvSimpsonVirome, family="binomial", data=data1))


