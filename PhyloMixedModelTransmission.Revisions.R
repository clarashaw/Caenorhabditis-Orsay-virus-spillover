setwd("C:/Users/cls6630/OneDrive - The Pennsylvania State University/Spillover")

library(dplyr)
library(tidyr)
library(ape)
library(MCMCglmm)
library(coda)

#Tree is from https://zenodo.org/record/3571457#.YJvu76hKg2x
#Stevens et. al. 2020 Current Biology
#Use CTree.ultra2 from PhylogMixedModelRevisions.R

#Load my infection data
data<-read.csv("TransmissionData.csv") #Get this from running "FBRTransmission

#Make a column called "Tree.match" that has the same species names as in the tree
data$Tree.match[data$Species=="castelli"]<-"CCAST"
data$Tree.match[data$Species=="tropicalis"]<-"CTROP"
data$Tree.match[data$Species=="wallacei"]<-"CWALL"
data$Tree.match[data$Species=="doughertyi"]<-"CDOUG"
data$Tree.match[data$Species=="latens"]<-"CLATE"
data$Tree.match[data$Species=="elegans"]<-"CELEG"
data$Tree.match[data$Species=="waitukubuli"]<-"CWAIT"
data$Tree.match[data$Species=="macrosperma"]<-"CMACR"
data$Tree.match[data$Species=="sulstoni"]<-"CSULS"

#Assign transmission abilities
#Make an transmission column that gets 0 if undetectable in exposure plate,
#1 if passage 1 is undetectable,
#2 if passage 1 has something but it eventually dies out before passage 5
#3 if the passage makes it through passage 5

#Also make a column for average time on the plate
data.time<-data%>%group_by(Group, Strain, Species, Tree.match)%>%summarise(avg.time=mean(Days.on.plate))
data1<-select(data, Group, Strain, Species, Tree.match,Passage, Ct)
datawide = data1 %>% spread(Passage, Ct)

datawide<-full_join(datawide,data.time, by=c("Group","Strain","Species","Tree.match"))
colnames(datawide)<-c("Group","Strain","Species","Tree.match","P0","P1","P2","P3","P4","P5","average.time")

datawide$transmission[datawide$P0=="Inf"]<-0
datawide$transmission[datawide$P1=="Inf"]<-1
datawide$transmission[datawide$P5>0 & datawide$P5!="Inf"]<-3
datawide$transmission[is.na(datawide$transmission)]<-2

#filter just for species in the tree 
datawide$P0[datawide$P0==Inf]<-40 #assign this so that it can be used in the virus amplification on primary exposure plates model
tree.data<-filter(datawide, !is.na(Tree.match))

#Use CTree.small from Susceptibility analysis
plot(CTree.small)

#Calculate distances from C. elegans
Distances<-cophenetic.phylo(CTree.small)
D.eleg<-as.data.frame(Distances[,"CELEG"])
colnames(D.eleg)<-"dist.from.elegans"
D.eleg$Tree.match<-rownames(D.eleg)

#Calculate the inverse relatedness matrix (phylogenetic covariance matrix)
class(CTree.ultra2)<-"phylo"
INphylo <- inverseA(CTree.small, nodes="ALL",scale=TRUE)

tree.data<-as.data.frame(tree.data) #for some reason, this is necessary

#Put transmission data together with distance data
tree.data<-left_join(tree.data, D.eleg, by=c("Tree.match"))

#Get correlation between P0 and distance from C.elegans
cor(tree.data$dist.from.elegans, tree.data$P0) #0.4770969

ggplot(tree.data, aes(dist.from.elegans, transmission))+geom_point(alpha=0.2)
ggplot(tree.data, aes(P0, transmission))+geom_point(alpha=0.5)+xlab('Ct on exposure plate')+
  ylab("Transmission ability")+theme_bw()

##########################################################################
#MCMC Generalized linear mixed effects models

#Each random effect needs a prior G
#Residuals get a prior R
#Use default prior for fixed (normal distribution with mean 0 and variance of 10^8)

#with phylogeny and species
priors1 <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                         G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                         G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)),
        R = list(V = 1, nu = 0.002)) 

model1.T <- MCMCglmm(transmission ~ dist.from.elegans+P0, 
                  random = ~Tree.match+Species+Strain,
                  ginverse = list(Tree.match=INphylo$Ainv), 
                  prior=priors1,nitt = 100000000,thin=5000,burnin = 30000,
                  data = tree.data, verbose = FALSE, family="gaussian")
summary(model1.T) #DIC=66.23027  
heidel.diag(model1.T$Sol) #Passed
heidel.diag(model1.T$VCV) #Passed

#determine median predicition for each fixed effect
#distance
dist<-seq(0,0.686613, length.out=30)
m<-model1.T$Sol[,1]
b<-model1.T$Sol[,2]
b1<-model1.T$Sol[,3]
mean.P0<-mean(tree.data$P0)
QT<- NULL
for (j in dist){
  POINTS<-NULL
  distance<-j
  for (i in 1:19994){
    point<-m[i]+(b[i])*distance+(b1[i]*mean.P0)
    POINTS<-c(POINTS, point)
  }
  Q<-quantile(POINTS, c(0.025, 0.5, 0.975))
  qt<-c(distance[1], Q[1], Q[2],Q[3])
  QT<-rbind(QT,qt)
}

QT<-as.data.frame(QT)
colnames(QT)<-c("dist.from.elegans","inf.025","inf.5","inf.975")

xlabel<-expression(paste("Phylogenetic distance from ",italic("C. elegans")))
Fig4a<-ggplot(tree.data, aes(dist.from.elegans, transmission))+
  geom_jitter(alpha=0.5, width=0.01, height=0.1)+
  geom_line(data=QT, aes(dist.from.elegans, inf.5), color="red")+
  geom_line(data=QT, aes(dist.from.elegans, inf.025), color="red", lty="dashed")+
  geom_line(data=QT, aes(dist.from.elegans, inf.975), color="red", lty="dashed")+
  theme_bw()+
  labs(x=xlabel,
       y=("Transmission ability")
  )+
  theme(axis.text = element_text(color="black"))+
  scale_y_continuous(limits=c(-0.5,4), breaks=c(0,1,2,3,4))

#amplification in primary exposure populations
P0.seq<-seq(11,40, length.out=30)
m<-model1.T$Sol[,1]
b<-model1.T$Sol[,2]
b1<-model1.T$Sol[,3]
mean.dist<-mean(tree.data$dist.from.elegans)
QT2<- NULL
for (j in P0.seq){
  POINTS<-NULL
  P0<-j
  for (i in 1:19994){
    point<-m[i]+(b[i])*mean.dist+(b1[i]*P0)
    POINTS<-c(POINTS, point)
  }
  Q<-quantile(POINTS, c(0.025, 0.5, 0.975))
  qt<-c(P0[1], Q[1], Q[2],Q[3])
  QT2<-rbind(QT2,qt)
}

QT2<-as.data.frame(QT2)
colnames(QT2)<-c("P0","inf.025","inf.5","inf.975")

Fig4b<-ggplot(tree.data, aes(P0, transmission))+
  geom_point(alpha=0.5)+
  geom_line(data=QT2, aes(P0, inf.5), color="red")+
  geom_line(data=QT2, aes(P0, inf.025), color="red", lty="dashed")+
  geom_line(data=QT2, aes(P0, inf.975), color="red", lty="dashed")+
  theme_bw()+
  labs(x="Viral amplification in \nprimary exposure populations (Ct)",
       y=("Transmission ability")
  )+
  theme(axis.text = element_text(color="black"))+
  scale_y_continuous(limits=c(-0.5,4), breaks=c(0,1,2,3,4))
##############################################

model2.T <- MCMCglmm(transmission ~ dist.from.elegans, 
                   random = ~Tree.match+Species+Strain,
                   ginverse = list(Tree.match=INphylo$Ainv), 
                   prior=priors1,nitt = 100000000,thin=5000,burnin = 30000,
                   data = tree.data, verbose = FALSE, family="gaussian")
summary(model2.T) #DIC=70.24553 
heidel.diag(model2.T$Sol) #passes
heidel.diag(model2.T$VCV) #passes

model3.T <- MCMCglmm(transmission ~ P0, 
                   random = ~Tree.match+Species+Strain,
                   ginverse = list(Tree.match=INphylo$Ainv), 
                   prior=priors1,nitt = 100000000,thin=5000,burnin = 30000,
                   data = tree.data, verbose = FALSE, family="gaussian")
summary(model3.T) #DIC=66.86288
heidel.diag(model3.T$Sol) #passed
heidel.diag(model3.T$VCV) #passed


model4.T <- MCMCglmm(transmission ~ 1, 
                   random = ~Tree.match+Species+Strain,
                   ginverse = list(Tree.match=INphylo$Ainv), 
                   prior=priors1,nitt = 100000000,thin=5000,burnin = 30000,
                   data = tree.data, verbose = FALSE, family="gaussian")
summary(model4.T) #DIC=70.43542 
heidel.diag(model4.T$Sol) #passed
heidel.diag(model4.T$VCV) #passed


priors2 <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                         G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)),
                R = list(V = 1, nu = 0.002)) 

model6.T <- MCMCglmm(transmission ~ dist.from.elegans+P0, 
                   random = ~Species+Strain,
                   ginverse = list(Tree.match=INphylo$Ainv), 
                   prior=priors2,nitt = 100000000,thin=5000,burnin = 30000,
                   data = tree.data, verbose = FALSE, family="gaussian")
summary(model6.T) #66.74811
heidel.diag(model6.T$Sol) #passed
heidel.diag(model6.T$VCV) #passed


model13.T <- MCMCglmm(transmission ~ dist.from.elegans, 
                    random = ~Species+Strain,
                    ginverse = list(Tree.match=INphylo$Ainv), 
                    prior=priors2,nitt = 100000000,thin=5000,burnin = 30000,
                    data = tree.data, verbose = FALSE, family="gaussian")
summary(model13.T) #70.39674
heidel.diag(model13.T$Sol) #passed
heidel.diag(model13.T$VCV) #passed

model14.T <- MCMCglmm(transmission ~ P0, 
                    random = ~Species+Strain,
                    ginverse = list(Tree.match=INphylo$Ainv), 
                    prior=priors2,nitt = 100000000,thin=5000,burnin = 30000,
                    data = tree.data, verbose = FALSE, family="gaussian")
summary(model14.T) #67.13791 
heidel.diag(model14.T$Sol)
heidel.diag(model14.T$VCV)


model16.T <- MCMCglmm(transmission ~ 1, 
                    random = ~Species+Strain,
                    ginverse = list(Tree.match=INphylo$Ainv), 
                    prior=priors2,nitt = 100000000,thin=5000,burnin = 30000,
                    data = tree.data, verbose = FALSE, family="gaussian")
summary(model16.T) #70.43529  
heidel.diag(model16.T$Sol)
heidel.diag(model16.T$VCV)


M1DIC<-model1.T$DIC 
M2DIC<-model2.T$DIC 
M3DIC<-model3.T$DIC 
M4DIC<-model4.T$DIC 
M6DIC<-model6.T$DIC 
M13DIC<-model13.T$DIC 
M14DIC<-model14.T$DIC 
M16DIC<-model16.T$DIC 

M1DIC-M1DIC
M3DIC-M1DIC
M6DIC-M1DIC
M14DIC-M1DIC
M2DIC-M1DIC
M16DIC-M1DIC
M4DIC-M1DIC
M13DIC-M1DIC

#DIC weights
#caluclate deltaDIC
del1<-(M1DIC-M1DIC)/2
del2<-(M2DIC-M1DIC)/2
del3<-(M3DIC-M1DIC)/2
del4<-(M4DIC-M1DIC)/2
del6<-(M6DIC-M1DIC)/2
del13<-(M13DIC-M1DIC)/2
del14<-(M14DIC-M1DIC)/2
del16<-(M16DIC-M1DIC)/2

sumall<-exp(-1*del1)+exp(-1*del2)+exp(-1*del3)+exp(-1*del4)+exp(-1*del6)+exp(-1*del13)+exp(-1*del14)+exp(-1*del16)

#model weights
model1.weight<-exp(-1*del1)/sumall
model3.weight<-exp(-1*del3)/sumall
model6.weight<-exp(-1*del6)/sumall
model14.weight<-exp(-1*del14)/sumall
model2.weight<-exp(-1*del2)/sumall
model16.weight<-exp(-1*del16)/sumall
model4.weight<-exp(-1*del4)/sumall
model13.weight<-exp(-1*del13)/sumall

#models that include distance are 1,2,6,13
dist.weight<-(exp(-1*del1)+exp(-1*del2)+exp(-1*del6)+exp(-1*del13))/sumall
#0.5580364
#models that include P0 are 1,3,6,14
P0.weight<-(exp(-1*del1)+exp(-1*del3)+exp(-1*del6)+exp(-1*del14))/sumall
# 0.8617425
#models that include phylogeny are 1,2,3,4,
phyl.weight<-(exp(-1*del1)+exp(-1*del2)+exp(-1*del3)+exp(-1*del4))/sumall
#0.5455507
#models that include distance or phylogeny are 1,2,3,4,6,13
phyldist.weight<-(exp(-1*del1)+exp(-1*del2)+exp(-1*del3)+exp(-1*del4)+exp(-1*del6)+exp(-1*del13))/sumall
# 0.791883

###########################################################
#Get correlation with susceptibility experiment
dat.small<-filter(data, Passage==0)
View(dat.small)
dat.small$Ct[dat.small$Ct==Inf]<-40
dat.avg<-dat.small%>%group_by(Strain)%>%summarise(mean.ct=mean(Ct))

data.susc<-data.exposed #From PhyloMixedModel.R
data.susc$Ct[data.susc$Ct==Inf]<-40
dat.susc.2<-data.susc%>%group_by(Strain)%>%summarise(mean.ct.susc=mean(Ct))

both<-inner_join(dat.avg, dat.susc.2, by="Strain")
ggplot(both, aes(mean.ct.susc, mean.ct))+geom_point()
cor(both$mean.ct, both$mean.ct.susc)

######################################################
#Calculate variance explained by fixed effects in best model (model1.T) for
#Figure 4C
#distance from C. elegans
meandist<-mean(tree.data$dist.from.elegans)
chain.t<-as.numeric(model1.T$Sol[,'dist.from.elegans'])
OUT<-NULL
for (j in 1:length(chain.t)){
  for (i in 1:length(tree.data$dist.from.elegans)){
    dif<-chain.t[j]*meandist-chain.t[j]*tree.data$dist.from.elegans[i]
    dif2<-dif^2
    output<-c(i, j, dif2)
    OUT<-rbind(OUT, output)
  }
}
OUT<-as.data.frame(OUT)
colnames(OUT)<-c("i","j","dif2")
sumi<-OUT%>%group_by(j)%>%summarise(sumdif<-sum(dif2))
marg<-sum(sumi)/length(sumi$j)
sumi<-as.data.frame(sumi)
colnames(sumi)<-c("j", "sumi")

#amplification in primary exposure populations
meanP0<-mean(tree.data$P0)
chain.p<-as.numeric(model1.T$Sol[,'P0'])
OUT1<-NULL
for (j in 1:length(chain.p)){
  for (i in 1:length(tree.data$P0)){
    dif<-chain.p[j]*meanP0-chain.p[j]*tree.data$P0[i]
    dif2<-dif^2
    output<-c(i, j, dif2)
    OUT1<-rbind(OUT1, output)
  }
}
OUT1<-as.data.frame(OUT1)
colnames(OUT1)<-c("i","j","dif2")
sumi1<-OUT1%>%group_by(j)%>%summarise(sumdif<-sum(dif2))
marg1<-sum(sumi1)/length(sumi1$j)
sumi1<-as.data.frame(sumi1)
colnames(sumi1)<-c("j", "sumi")

lambda.dist<- sumi$sumi/(model1.T$VCV[,'Tree.match']+model1.T$VCV[,'Species']+model1.T$VCV[,'Strain']+model1.T$VCV[,'units']+sumi$sumi+sumi1$sumi)
mean(lambda.dist) #0.4659902
HPDinterval(lambda.dist) #6.049574e-09 0.8904562

lambda.P0<- sumi1$sumi/(model1.T$VCV[,'Tree.match']+model1.T$VCV[,'Species']+model1.T$VCV[,'Strain']+model1.T$VCV[,'units']+sumi$sumi+sumi1$sumi)
mean(lambda.P0) #0.4477158
HPDinterval(lambda.P0) #6.591477e-08 0.883488

#posterior probability of the phylogenetic signal
lambda1 <- (model1.T$VCV[,'Tree.match'])/
  (model1.T$VCV[,'Tree.match']+model1.T$VCV[,'Species']+model1.T$VCV[,'Strain']+model1.T$VCV[,'units']+sumi$sumi+sumi1$sumi)
mean(lambda1) #0.04341552
HPDinterval(lambda1) #2.204891e-12 0.1711036

lambda2 <- (model1.T$VCV[,'Species'])/
  (model1.T$VCV[,'Tree.match']+model1.T$VCV[,'Species']+model1.T$VCV[,'Strain']+model1.T$VCV[,'units']+sumi$sumi+sumi1$sumi)
mean(lambda2) #0.01382678
HPDinterval(lambda2) #6.181967e-11 0.05616818

lambda3 <- (model1.T$VCV[,'Strain'])/
  (model1.T$VCV[,'Tree.match']+model1.T$VCV[,'Species']+model1.T$VCV[,'Strain']+model1.T$VCV[,'units']+sumi$sumi+sumi1$sumi)
mean(lambda3) #0.005120081
HPDinterval(lambda3) #1.324389e-11 0.02070458

lambda4 <- (model1.T$VCV[,'units'])/
  (model1.T$VCV[,'Tree.match']+model1.T$VCV[,'Species']+model1.T$VCV[,'Strain']+model1.T$VCV[,'units']+sumi$sumi+sumi1$sumi)
mean(lambda4) #0.02393171
HPDinterval(lambda4)# 0.002267975 0.05516267

lambda.both <- (model1.T$VCV[,'Tree.match']+sumi$sumi)/
  (model1.T$VCV[,'Tree.match']+model1.T$VCV[,'Species']+model1.T$VCV[,'Strain']+model1.T$VCV[,'units']+sumi$sumi+sumi1$sumi)
mean(lambda.both) #0.5094057
HPDinterval(lambda.both)#0.01188442 0.9304754

#pull lambdas together for the figure
Figure.T<-data.frame(lambda=c(mean(lambda.dist),mean(lambda.P0),mean(lambda1), mean(lambda2), mean(lambda3), mean(lambda4), mean(lambda.both),NA),
                   low=c(HPDinterval(lambda.dist)[1],HPDinterval(lambda.P0)[1], HPDinterval(lambda1)[1], HPDinterval(lambda2)[1], HPDinterval(lambda3)[1], HPDinterval(lambda4)[1], HPDinterval(lambda.both)[1],NA),
                   high=c(HPDinterval(lambda.dist)[2],HPDinterval(lambda.P0)[2], HPDinterval(lambda1)[2], HPDinterval(lambda2)[2], HPDinterval(lambda3)[2], HPDinterval(lambda4)[2], HPDinterval(lambda.both)[2],NA),
                   label=c("phylogenetic distance","virus on P0","Phylogeny","Species","Strain","Residual","both","bar"))
Figure.T$label <- factor(Figure.T$label, levels = c("virus on P0","phylogenetic distance","Phylogeny", "Species", "Strain", "Residual","bar","both"))

write.csv(Figure.T, "var.explained.transmission.csv")

Fig4c<-ggplot(Figure.T, aes(x=label, y=lambda))+geom_point()+
  geom_errorbar(aes(ymin=low, ymax=high), width=0.1)+
  theme_bw()+ylab("Variance explained")+
  theme(axis.text = element_text(color="black"))+
  scale_x_discrete(labels=c("Ct of primary\nexposure \npopulations\n(fixed)","phyl. dist.\nfrom \nC. elegans\n(fixed)","pairwise\nphyl. dist.\n(random)", "species\n(random)","strain\n(random)","residual","","both phylogenetic\neffects"))+
  geom_vline(aes(xintercept='bar'), lty="dashed")+
  theme(axis.text = element_text(color="black"))+
  xlab("")

Figure4<-plot_grid(Fig4b, Fig4a,Fig4c, labels=c("A","B","C"), align="h",axis="b",ncol=3, rel_widths = c(0.5,0.5,0.8))
save_plot("SpilloverFigure4.jpg",Figure4, base_width = 14)

save.image(file = "Transmission.RData")
#######
mean(tree.data$dist.from.elegans)
max(tree.data$dist.from.elegans)
min(tree.data$dist.from.elegans)
