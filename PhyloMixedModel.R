#Use this code to analyze spillover data for effects of phylogeny

setwd("C:/Users/cls6630/OneDrive - The Pennsylvania State University/Spillover")

library(dplyr)
library(tidyr)
library(ape)
library(MCMCglmm)
library(coda)

#Load the tree (from https://zenodo.org/record/3571457#.YJvu76hKg2x)
#this is the tree in Stevens et. al. 2020 Current Biology
CTree<-read.tree("caenorhabditis_35species_1167orthos.treefile")
CTree$tip.label

#Load my spillover data
data<-read.csv("spillover.blocks.csv")
data$Species[data$Species=="sp.43"]<-"vivipara"

#Make a column called "Tree.match" that has the same species names as in the tree
data$Tree.match[data$Species=="afra"]<-"CAFRA"
data$Tree.match[data$Species=="castelli"]<-"CCAST"
data$Tree.match[data$Species=="quiockensis"]<-"CQUIO"
data$Tree.match[data$Species=="virilis"]<-"CVIRI"
data$Tree.match[data$Species=="vivipara"]<-"CVIVI"
data$Tree.match[data$Species=="monodelphis"]<-"CMONO"
data$Tree.match[data$Species=="parvicauda"]<-"CPARV"
data$Tree.match[data$Species=="guadaloupensis"]<-"CGUAD"
data$Tree.match[data$Species=="uteleia"]<-"CUTEL"
data$Tree.match[data$Species=="brenneri"]<-"CBREN"
data$Tree.match[data$Species=="tropicalis"]<-"CTROP"
data$Tree.match[data$Species=="wallacei"]<-"CWALL"
data$Tree.match[data$Species=="doughertyi"]<-"CDOUG"
data$Tree.match[data$Species=="briggsae"]<-"CBRIG"
data$Tree.match[data$Species=="nigoni"]<-"CNIGO"
data$Tree.match[data$Species=="sinica"]<-"CSINI"
data$Tree.match[data$Species=="tribulationis"]<-"CTRIB"
data$Tree.match[data$Species=="zanzibari"]<-"CZANZ"
data$Tree.match[data$Species=="latens"]<-"CLATE"
data$Tree.match[data$Species=="remanei"]<-"CREMA"
data$Tree.match[data$Species=="elegans"]<-"CELEG"
data$Tree.match[data$Species=="inopinata"]<-"CINOP"
data$Tree.match[data$Species=="kamaaina"]<-"CKAMA"
data$Tree.match[data$Species=="becei"]<-"CBECE"
data$Tree.match[data$Species=="nouraguensis"]<-"CNOUR"
data$Tree.match[data$Species=="panamensis"]<-"CPANA"
data$Tree.match[data$Species=="waitukubuli"]<-"CWAIT"
data$Tree.match[data$Species=="macrosperma"]<-"CMACR"
data$Tree.match[data$Species=="sulstoni"]<-"CSULS"
data$Tree.match[data$Species=="monodelphis"]<-"CMONO"

#filter out N2s and JU1580s (+ and - control plates)
data.exposed<-filter(data, Strain != "N2")
data.exposed<-filter(data.exposed, Strain!="JU1580")

#Make an infection column that gets 1 if I deemed plate to be infected (Ct<29.54051); otherwise 0.
#infection cutoff = SD.line from SpilloverDataAnalysis
data.exposed$infection <- ifelse(data.exposed$Ct <= 29.54051, 1,0) 

#get rid of plates that had no Cts (were removed from experiment usually due to contamination)
data.exposed<-filter(data.exposed, !is.na(Ct))

#filter just for species in the tree 
tree.data<-filter(data.exposed, !is.na(Tree.match))

#Get fraction of plates infected and uninfected as well as average time on plate
tree.summary<-tree.data%>%group_by(Strain, Species, Tree.match) %>% 
  summarise(exposed=n(),inf=sum(infection), avg.days=mean(Days.on.plate))
tree.summary$n.inf<-tree.summary$exposed-tree.summary$inf

#Now look at the tree.
plot(CTree) #The tree needs to be rooted

#root the tree. DPACH is the outgroup
CTree<-root(CTree, "DPACH", node = NULL, resolve.root = FALSE,
            interactive = FALSE, edgelabel = FALSE)

#The tree needs to be ultrameric
is.ultrametric(CTree) #FALSE

#make the tree ultrameric
#chronopl() estimates the node ages of a tree using a semi-parametric method based on penalized
#likelihood (Sanderson 2002). The branch lengths of the input tree are interpreted as mean numbers
#of substitutions (i.e., per site).
CTree.ultra2<-chronopl(CTree, lambda=1) 

is.ultrametric(CTree.ultra2) #TRUE

#drop the outgroup DPACH because there are two (DCORO is still an outgroup)
CTree.ultra2<-drop.tip(CTree.ultra2, c("DPACH"), trim.internal = TRUE, subtree = FALSE, 
                root.edge = 0, rooted = is.rooted(CTree), collapse.singles = TRUE,
                interactive = FALSE)

#need unique node labels- add them in
CTree.ultra2$node.label<-as.character(c(1:33))

#Get distances from C. elegans
Distances<-cophenetic.phylo(CTree.ultra2)
D.eleg<-as.data.frame(Distances[,"CELEG"])
colnames(D.eleg)<-"dist.from.elegans"
D.eleg$Tree.match<-rownames(D.eleg)

#Calculate the inverse relatedness matrix (phylogenetic covariance matrix)
INphylo <- inverseA(CTree.ultra2, nodes="TIPS",scale=TRUE)

#Finish preparing the data - turn to a dataframe and join with distances dataframe caluclated above.
tree.summary<-as.data.frame(tree.summary) #for some reason, this is necessary
tree.summary<-full_join(tree.summary, D.eleg, by=c("Tree.match"))
tree.summary<-filter(tree.summary, !is.na(Species))#these are the ones in the tree that I don't have samples for.
tree.summary$prop.inf<-tree.summary$inf/tree.summary$exposed

#Find the min and max distance from C. elegans.
min(tree.summary$dist.from.elegans)
max(tree.summary$dist.from.elegans)
mean(tree.summary$dist.from.elegans)

#send this file to the super computer in PhyloMixedModel.super 1-5 files
write.csv(tree.summary, "tree.summary.csv", row.names = FALSE)

#Each random effect needs a prior G
#Residuals get a prior R
#Use default prior for fixed (normal distribution with mean 0 and variance of 10^8)

#model 1 has distance as a fixed effect and species as a random effect
priors1 <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)),                
                  R = list(V = 1, nu = 0.002))

model1.S <- MCMCglmm(cbind(inf, n.inf) ~ dist.from.elegans, 
                   random = ~Species,
                   prior=priors1,nitt = 100000000,thin=10000,burnin = 300000,
                   data = tree.summary, verbose = FALSE, family="multinomial2")
summary(model1.S) #DIC=51.84086  
heidel.diag(model1.S$Sol) 
heidel.diag(model1.S$VCV) 

#model2 includes distance with phylogeny and species in random effects
priors2 <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                         G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)), 
                R = list(V = 1, nu = 0.002))

model2.S <- MCMCglmm(cbind(inf, n.inf) ~ dist.from.elegans, 
                  random = ~Tree.match+Species,
                  start=list(QUASI=FALSE),
                  ginverse = list(Tree.match=INphylo$Ainv), 
                  prior=priors2,nitt = 100000000,thin=10000,burnin = 300000,
                  data = tree.summary, verbose = FALSE, family="multinomial2")
summary(model2.S) #DIC=50.72315
heidel.diag(model2.S$Sol) 
heidel.diag(model2.S$VCV) 

#model 3 does not include distance, but does include both random effects (use same priors as model 2)
model3.S <- MCMCglmm(cbind(inf, n.inf) ~ 1, 
                   random = ~Tree.match+Species,
                   start=list(QUASI=FALSE),
                   ginverse = list(Tree.match=INphylo$Ainv), 
                   prior=priors2,nitt = 100000000,thin=10000,burnin = 300000,
                   data = tree.summary, verbose = FALSE, family="multinomial2")
summary(model3.S) #DIC: 52.91207 
heidel.diag(model3.S$Sol)
heidel.diag(model3.S$VCV)

#without phylogeny, with species
priors4 <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)),
                R = list(V = 1, nu = 0.002))


model4.S <- MCMCglmm(cbind(inf, n.inf) ~ 1, 
                   random = ~Species,
                   prior=priors4,nitt = 100000000,thin=10000,burnin = 300000,
                   data = tree.summary, verbose = FALSE, family="multinomial2")
summary(model4.s) #DIC=54.4782 
heidel.diag(model4.S$Sol) 
heidel.diag(model4.S$VCV)

#DIC weights
#caluclate deltaDIC

#From super computer output 8/27/2021
model1DIC<-model1.S$DIC #51.84407
model2DIC<-model2.S$DIC #50.72315
model3DIC<-model3.S$DIC #52.91207 
model4DIC<-model4.S$DIC #54.48397

model1DIC-model2DIC
model3DIC-model2DIC
model4DIC-model2DIC

del1<-(model1DIC-model2DIC)/2
del2<-(model2DIC-model2DIC)/2
del3<-(model3DIC-model2DIC)/2
del4<-(model4DIC-model2DIC)/2


sumall<-exp(-1*del1)+exp(-1*del2)+exp(-1*del3)+exp(-1*del4)
  
dist.weight<-(exp(-1*del1)+exp(-1*del2))/sumall
phyl.weight<-(exp(-1*del2)+exp(-1*del3))/sumall
both.weight<-(exp(-1*del2)+exp(-1*del1)+exp(-1*del3))/sumall

model2.weight<-exp(-1*del2)/sumall
model1.weight<-exp(-1*del1)/sumall
model3.weight<-exp(-1*del3)/sumall
model4.weight<-exp(-1*del4)/sumall


############################################
#No longer using this analysis
#posterior probabilities
meandist<-mean(tree.summary$dist.from.elegans)
chain.susc<-as.numeric(model2$Sol[,'dist.from.elegans'])
OUT.susc<-NULL
for (j in 1:length(chain.susc)){
  for (i in 1:length(tree.summary$dist.from.elegans)){
  dif<-chain.susc[j]*meandist-chain.susc[j]*tree.summary$dist.from.elegans[i]
  dif2<-dif^2
  output<-c(i, j, dif2)
  OUT.susc<-rbind(OUT.susc, output)
  }
}
OUT.susc<-as.data.frame(OUT.susc)
colnames(OUT.susc)<-c("i","j","dif2")
sumi.susc<-OUT.susc%>%group_by(j)%>%summarise(sumdif<-sum(dif2))
sumi.susc<-as.data.frame(sumi.susc)
colnames(sumi.susc)<-c("j", "sumi")

dist.var<- sumi.susc$sumi/(model2$VCV[,'Tree.match']+model2$VCV[,'Species']+model2$VCV[,'units']+sumi.susc$sumi+(pi^2)/3)
mean(dist.var)
HPDinterval(dist.var)

#posterior probability of the phylogenetic signal
lambda1 <- (model2$VCV[,'Tree.match'])/
  (model2$VCV[,'Tree.match']+model2$VCV[,'Species']+model2$VCV[,'units']+sumi.susc$sumi+(pi^2)/3)
mean(lambda1) #0.1068
HPDinterval(lambda1) 

lambda2 <- (model2$VCV[,'Species'])/
  (model2$VCV[,'Tree.match']+model2$VCV[,'Species']+model2$VCV[,'units']+sumi.susc$sumi+(pi^2)/3)
mean(lambda2) 
HPDinterval(lambda2) 

lambda3 <- (model2$VCV[,'units'])/
  (model2$VCV[,'Tree.match']+model2$VCV[,'Species']+model2$VCV[,'units']+sumi.susc$sumi+(pi^2)/3)
mean(lambda3) 
HPDinterval(lambda3) 

lambda4 <- (model2$VCV[,'Tree.match']+model2$VCV[,'Species']+sumi.susc$sumi)/
  (model2$VCV[,'Tree.match']+model2$VCV[,'Species']+model2$VCV[,'units']+sumi.susc$sumi+(pi^2)/3)
mean(lambda4) 
HPDinterval(lambda4)

Figure<-data.frame(lambda=c(mean(lambda0), mean(lambda1), mean(lambda2), mean(lambda3)),
                   low=c(HPDinterval(lambda0)[1],HPDinterval(lambda1)[1], HPDinterval(lambda2)[1], HPDinterval(lambda3)[1]),
                   high=c(HPDinterval(lambda0)[2],HPDinterval(lambda1)[2], HPDinterval(lambda2)[2], HPDinterval(lambda3)[2]),
                   label=c("phylogenetic distance","phylogeny","species-level","residual"))
Figure$label <- factor(Figure$label, levels = c("phylogenetic distance", "phylogeny", "species-level", "residual"))

ggplot(Figure, aes(x=label, y=lambda))+geom_point()+
  geom_errorbar(aes(ymin=low, ymax=high))+theme_bw()+ylab("Variance explained")+
  xlab("")+theme(axis.text = element_text(color="black"))








































