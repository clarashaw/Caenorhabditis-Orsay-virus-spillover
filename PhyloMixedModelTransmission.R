setwd("C:/Users/cls6630/OneDrive - The Pennsylvania State University/Spillover")

library(dplyr)
library(tidyr)
library(ape)
library(MCMCglmm)
library(coda)

#Tree is from https://zenodo.org/record/3571457#.YJvu76hKg2x
#Stevens et. al. 2020 Current Biology
CTree<-read.tree("caenorhabditis_35species_1167orthos.treefile")
CTree$tip.label

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

#Work on the tree
plot(CTree) #The tree needs to be rooted

#root the tree. DPACH is the outgroup
CTree<-root(CTree, "DPACH", node = NULL, resolve.root = FALSE,
            interactive = FALSE, edgelabel = FALSE)

#The tree needs to be ultrameric
is.ultrametric(CTree) #FALSE

#make the tree ultrameric
CTree.ultra2<-chronopl(CTree, lambda=1) #lambda=0 makes it so each branch has different evolutionary rates

is.ultrametric(ultramericTree)

#drop the outgroup DPACH because there are two (DCORO is still an outgroup)
#drop tips not in the smaller set ()
CTree.ultra2<-drop.tip(CTree.ultra2, c("DPACH"), trim.internal = TRUE, subtree = FALSE, 
                root.edge = 0, rooted = is.rooted(CTree), collapse.singles = TRUE,
                interactive = FALSE)

#need unique node labels
CTree.ultra2$node.label<-as.character(c(1:33))

plot(CTree.ultra2)

#Calculate distances from C. elegans
Distances<-cophenetic.phylo(CTree.ultra2)
D.eleg<-as.data.frame(Distances[,"CELEG"])
colnames(D.eleg)<-"dist.from.elegans"
D.eleg$Tree.match<-rownames(D.eleg)

#Calculate the inverse relatedness matrix (phylogenetic covariance matrix)
INphylo <- inverseA(CTree.ultra2, nodes="TIPS",scale=TRUE)

tree.data<-as.data.frame(tree.data) #for some reason, this is necessary

#Put transmission data together with distance data
tree.data<-left_join(tree.data, D.eleg, by=c("Tree.match"))

#Get correlation between P0 and distance from C.elegans
cor(tree.data$dist.from.elegans, tree.data$P0) #0.4614948

##########################################################################
#Generalized linear mixed effects models

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
                  prior=priors1,nitt = 100000000,thin=10000,burnin = 300000,
                  data = tree.data, verbose = FALSE, family="gaussian")
summary(model1.T) #DIC=66.34947
heidel.diag(model1.T$Sol) #Passed
heidel.diag(model1.T$VCV)

model2.T <- MCMCglmm(transmission ~ dist.from.elegans, 
                   random = ~Tree.match+Species+Strain,
                   ginverse = list(Tree.match=INphylo$Ainv), 
                   prior=priors1,nitt = 100000000,thin=10000,burnin = 300000,
                   data = tree.data, verbose = FALSE, family="gaussian")
summary(model2.T) #DIC=69.88948
heidel.diag(model2.T$Sol) 
heidel.diag(model2.T$VCV)

model3.T <- MCMCglmm(transmission ~ P0, 
                   random = ~Tree.match+Species+Strain,
                   ginverse = list(Tree.match=INphylo$Ainv), 
                   prior=priors1,nitt = 100000000,thin=10000,burnin = 300000,
                   data = tree.data, verbose = FALSE, family="gaussian")
summary(model3.T) #DIC=66.83158
heidel.diag(model3.T$Sol) 
heidel.diag(model3.T$VCV)


model4.T <- MCMCglmm(transmission ~ 1, 
                   random = ~Tree.match+Species+Strain,
                   ginverse = list(Tree.match=INphylo$Ainv), 
                   prior=priors1,nitt = 100000000,thin=10000,burnin = 300000,
                   data = tree.data, verbose = FALSE, family="gaussian")
summary(model4.T) #DIC=70.44053
heidel.diag(model4.T$Sol) 
heidel.diag(model4.T$VCV)


priors2 <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                         G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)),
                R = list(V = 1, nu = 0.002)) 

model6.T <- MCMCglmm(transmission ~ dist.from.elegans+P0, 
                   random = ~Species+Strain,
                   ginverse = list(Tree.match=INphylo$Ainv), 
                   prior=priors2,nitt = 100000000,thin=10000,burnin = 300000,
                   data = tree.data, verbose = FALSE, family="gaussian")
summary(model6.T) #66.93487
heidel.diag(model6.T$Sol)
heidel.diag(model6.T$VCV)


model13.T <- MCMCglmm(transmission ~ dist.from.elegans, 
                    random = ~Species+Strain,
                    ginverse = list(Tree.match=INphylo$Ainv), 
                    prior=priors2,nitt = 100000000,thin=10000,burnin = 300000,
                    data = tree.data, verbose = FALSE, family="gaussian")
summary(model13.T) #69.83412
heidel.diag(model13.T$Sol)
heidel.diag(model13.T$VCV)


model14.T <- MCMCglmm(transmission ~ P0, 
                    random = ~Species+Strain,
                    ginverse = list(Tree.match=INphylo$Ainv), 
                    prior=priors2,nitt = 100000000,thin=10000,burnin = 300000,
                    data = tree.data, verbose = FALSE, family="gaussian")
summary(model14.T) #67.13872 
heidel.diag(model14.T$Sol)
heidel.diag(model14.T$VCV)


model16.T <- MCMCglmm(transmission ~ 1, 
                    random = ~Species+Strain,
                    ginverse = list(Tree.match=INphylo$Ainv), 
                    prior=priors2,nitt = 100000000,thin=10000,burnin = 300000,
                    data = tree.data, verbose = FALSE, family="gaussian")
summary(model16.T) #70.4342
heidel.diag(model16.T$Sol)
heidel.diag(model16.T$VCV)


M1DIC<-model1.T$DIC #66.34947
M2DIC<-model2.T$DIC #70.29106
M3DIC<-model3.T$DIC #66.88259
M4DIC<-model4.T$DIC #70.44053
M6DIC<-model6.T$DIC #66.93487
M13DIC<-model13.T$DIC #70.46175
M14DIC<-model14.T$DIC #67.13939
M16DIC<-model16.T$DIC #70.43571

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
#0.5424126
#models that include P0 are 1,3,6,14
P0.weight<-(exp(-1*del1)+exp(-1*del3)+exp(-1*del6)+exp(-1*del14))/sumall
#0.8582431
#models that include phylogeny are 1,2,3,4,
phyl.weight<-(exp(-1*del1)+exp(-1*del2)+exp(-1*del3)+exp(-1*del4))/sumall
#0.5481034
#models that include distance or phylogeny are 1,2,3,4,6,13
phyldist.weight<-(exp(-1*del1)+exp(-1*del2)+exp(-1*del3)+exp(-1*del4)+exp(-1*del6)+exp(-1*del13))/sumall
#0.7835963

###########################################################
#Get correlation coef
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
