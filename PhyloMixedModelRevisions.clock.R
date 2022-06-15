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

tree.summary<-tree.data%>%group_by(Strain, Species, Tree.match) %>% 
  summarise(exposed=n(),inf=sum(infection), avg.days=mean(Days.on.plate))
tree.summary$n.inf<-tree.summary$exposed-tree.summary$inf

#Now look at the tree.
plot(CTree) #The tree needs to be rooted

#Find species in the tree that I don't have
other<-setdiff(CTree$tip.label, tree.summary$Tree.match)

#root the tree. DPACH is the outgroup
CTree<-root(CTree, "DPACH", node = NULL, resolve.root = FALSE,
            interactive = FALSE, edgelabel = FALSE)

#The tree needs to be ultrametric
is.ultrametric(CTree) #FALSE

#make the tree ultrametric
#CTree.ultra2<-chronos(CTree, model="relaxed", lambda=0) #PHIIC=216.96

#CTree.ultra2<-chronos(CTree, model="discrete", lambda=0) #PHIIC = 123.25

#CTree.ultra2<-chronos(CTree, model="correlated", lambda=0) #PHIIC= 216.96 

CTree.ultra2<-chronos(CTree, model="clock") #PHIIC=87.25

is.ultrametric(CTree.ultra2) #TRUE

#Prune for our species
CTree.small<-drop.tip(CTree.ultra2, other, trim.internal = TRUE, subtree = FALSE, 
                      root.edge = 0, rooted = is.rooted(CTree), collapse.singles = TRUE,
                      interactive = FALSE)

#need unique node labels- add them in
CTree.small$node.label<-as.character(c(1:28))

#Get distances from C. elegans
Distances<-cophenetic.phylo(CTree.small)
D.eleg<-as.data.frame(Distances[,"CELEG"])
colnames(D.eleg)<-"dist.from.elegans"
D.eleg$Tree.match<-rownames(D.eleg)

#Calculate the inverse relatedness matrix (phylogenetic covariance matrix)
class(CTree.small)<-"phylo"
INphylo <- inverseA(CTree.small, nodes="TIPS",scale=TRUE)

#Finish preparing the data - turn to a dataframe and join with distances dataframe caluclated above.
tree.summary<-as.data.frame(tree.summary) 
tree.summary<-full_join(tree.summary, D.eleg, by=c("Tree.match"))
tree.summary$prop.inf<-tree.summary$inf/tree.summary$exposed

#Find the min and max distance from C. elegans.
min(tree.summary$dist.from.elegans) #0
max(tree.summary$dist.from.elegans) #1.062222
mean(tree.summary$dist.from.elegans) #0.3666099

#correct the average number of days for this one
tree.summary$avg.days[tree.summary$avg.days==-44379]<-30

write.csv(tree.summary, "tree.summary.revised.csv", row.names = FALSE)

#Each random effect needs a prior G
#Residuals get a prior R
#Use default prior for fixed (normal distribution with mean 0 and variance of 10^8)

#model 1 has distance as a fixed effect and species as a random effect
priors1 <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)),                
                  R = list(V = 1, nu = 0.002))

model1.S <- MCMCglmm(cbind(inf, n.inf) ~dist.from.elegans, 
                   random = ~Species,
                   prior=priors1,nitt = 2000000000,thin=100000,burnin = 300000,
                   data = tree.summary, verbose = FALSE, family="multinomial2")
summary(model1.S) #DIC=51.84581  
heidel.diag(model1.S$Sol) #passes
heidel.diag(model1.S$VCV) #passes

#model2 includes distance with phylogeny and species in random effects
priors2 <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                         G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)), 
                R = list(V = 1, nu = 0.002))

model2.S <- MCMCglmm(cbind(inf, n.inf) ~dist.from.elegans, 
                  random = ~Tree.match+Species,
                  start=list(QUASI=FALSE),
                  ginverse = list(Tree.match=INphylo$Ainv), 
                  prior=priors2,nitt = 500000000,thin=50000,burnin = 300000,
                  data = tree.summary, verbose = FALSE, family="multinomial2")
summary(model2.S) #DIC=50.1151
heidel.diag(model2.S$Sol) #passes!
heidel.diag(model2.S$VCV) #passes!

#Create figure 2 - what is the median prediction for effect of distance?
dist<-seq(0,1.062222, length.out=30)
m<-model2.S$Sol[,1]
b<-model2.S$Sol[,2]
QS<- NULL
for (j in dist){
  POINTS<-NULL
  distance<-j
  for (i in 1:9970){
  point<-m[i]+(b[i])*distance
  POINTS<-c(POINTS, point)
  }
  Q<-quantile(POINTS, c(0.025, 0.5, 0.975))
  qs<-c(distance[1], Q[1], Q[2],Q[3])
  QS<-rbind(QS,qs)
}

QS<-as.data.frame(QS)
colnames(QS)<-c("dist.from.elegans","inf.025","inf.5","inf.975")

require(boot)
QS$inv.inf.5<-inv.logit(QS$inf.5)
QS$inv.inf.025<-inv.logit(QS$inf.025)
QS$inv.inf.975<-inv.logit(QS$inf.975)
xlabel<-expression(paste("Phylogenetic distance from ",italic("C. elegans")))
Fig2a<-ggplot(tree.summary, aes(dist.from.elegans, prop.inf))+
  geom_jitter(alpha=0.5, width = 0.01, height=0.01)+
  geom_line(data=QS, aes(dist.from.elegans, inv.inf.5), color="red")+
  geom_line(data=QS, aes(dist.from.elegans, inv.inf.025), color="red", lty="dashed")+
  geom_line(data=QS, aes(dist.from.elegans, inv.inf.975), color="red", lty="dashed")+
  theme_bw()+
  labs(x=xlabel,
       y=("Proportion populations infected")
       )+
  theme(axis.text = element_text(color="black"))

#calculate variance explained by distance (fixed effect in model2.S)
meandist<-mean(tree.summary$dist.from.elegans)
chain.susc<-as.numeric(model2.S$Sol[,'dist.from.elegans'])
#This loop squares the deviations from the mean for each slope prediction in the chain
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
#caluclate the sum of squares for each value in the chain
sumi.susc<-OUT.susc%>%group_by(j)%>%summarise(sumdif<-sum(dif2))
sumi.susc<-as.data.frame(sumi.susc)
colnames(sumi.susc)<-c("j", "sumi")

write.csv(sumi.susc, "model2.S.variance.csv")

#R^2 of distance fixed effect
dist.var<- sumi.susc$sumi/(model2.S$VCV[,'Tree.match']+model2.S$VCV[,'Species']+model2.S$VCV[,'units']+sumi.susc$sumi+(pi^2)/3)
mean(dist.var) #0.890
HPDinterval(dist.var) #0.4866203 0.9963056

#posterior probability of pairwise phylogenetic distances
lambda1 <- (model2.S$VCV[,'Tree.match'])/
  (model2.S$VCV[,'Tree.match']+model2.S$VCV[,'Species']+model2.S$VCV[,'units']+sumi.susc$sumi+(pi^2)/3)
mean(lambda1) #0.05149834
HPDinterval(lambda1) #3.426104e-10 0.2199584

lambda.Species <- (model2.S$VCV[,'Species'])/
  (model2.S$VCV[,'Tree.match']+model2.S$VCV[,'Species']+model2.S$VCV[,'units']+sumi.susc$sumi+(pi^2)/3)
mean(lambda.Species) #0.04233608
HPDinterval(lambda.Species) #1.175506e-09 0.2048655

#phylogenetic signal (either part of phylogeny)
lambda <- (model2.S$VCV[,'Tree.match']+sumi.susc$sumi)/
  (model2.S$VCV[,'Tree.match']+model2.S$VCV[,'Species']+model2.S$VCV[,'units']+sumi.susc$sumi+(pi^2)/3)
mean(lambda) #0.9414986
HPDinterval(lambda)
# 0.7275703 0.9999469

#residual
lambda.res <- (model2.S$VCV[,'units'])/
  (model2.S$VCV[,'Tree.match']+model2.S$VCV[,'Species']+model2.S$VCV[,'units']+sumi.susc$sumi+(pi^2)/3)
mean(lambda.res) # 0.02153134
HPDinterval(lambda.res)
#4.520835e-07 0.09206744

#calculate heritability
heritability <- (model2.S$VCV[,'Tree.match'])/
  (model2.S$VCV[,'Tree.match']+model2.S$VCV[,'Species'])
mean(heritability) #0.5994237
HPDinterval(heritability)
#0.01335925     1

#pull together data for figure 2B
Figure<-data.frame(lambda=c(mean(dist.var),mean(lambda1), mean(lambda.Species), mean(lambda), mean(lambda.res),NA),
                   low=c(HPDinterval(dist.var)[1],HPDinterval(lambda1)[1], HPDinterval(lambda.Species)[1], HPDinterval(lambda)[1], HPDinterval(lambda.res)[1],NA),
                   high=c(HPDinterval(dist.var)[2],HPDinterval(lambda1)[2], HPDinterval(lambda.Species)[2], HPDinterval(lambda)[2], HPDinterval(lambda.res)[2],NA),
                   label=c("phylogenetic distance","pairwise","species-level","phylogeny (both)","residual","bar"))
Figure$label <- factor(Figure$label, levels = c("phylogenetic distance", "pairwise", "species-level", "residual", "bar","phylogeny (both)"))

write.csv(Figure, "df.for.Fig2.csv")

xlabels<-c(expression(paste("phyl. dist. from",italic("\nC. elegans\n"),"(fixed)")),
           "pairwise phyl. dist.\n(random)", "species\n(random)","residual","","both phylogenetic \neffects")

Fig2B<-ggplot(Figure, aes(x=label, y=lambda))+geom_point()+
  geom_errorbar(aes(ymin=low, ymax=high), width=0.1)+
  theme_bw()+ylab("Variance explained")+
  xlab("")+
  theme(axis.text = element_text(color="black"))+
  scale_x_discrete(labels=c("phyl. dist. \nfrom C. elegans\n(fixed)","pairwise\nphyl. dist.\n(random)", "species\n(random)","residual","","both phylogenetic \neffects"))+
  geom_vline(aes(xintercept='bar'), lty="dashed")

Figure2<-plot_grid(Fig2a, Fig2B, labels=c("A","B"), align="h",axis="b",ncol=2, rel_widths = c(0.7,1))
save_plot("SpilloverFigure2.jpg",Figure2, base_width = 9)

#model 3 does not include distance, but does include both random effects (use same priors as model 2)
model3.S <- MCMCglmm(cbind(inf, n.inf) ~ 1, 
                   random = ~Tree.match+Species,
                   start=list(QUASI=FALSE),
                   ginverse = list(Tree.match=INphylo$Ainv), 
                   prior=priors2,nitt = 100000000,thin=5000,burnin = 30000,
                   data = tree.summary, verbose = FALSE, family="multinomial2")
summary(model3.S) #DIC: 52.4848
heidel.diag(model3.S$Sol)
heidel.diag(model3.S$VCV)

#without phylogeny, with species
priors4 <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)),
                R = list(V = 1, nu = 0.002))


model4.S <- MCMCglmm(cbind(inf, n.inf) ~ 1, 
                   random = ~Species,
                   prior=priors4,nitt = 100000000,thin=5000,burnin = 30000,
                   data = tree.summary, verbose = FALSE, family="multinomial2")
summary(model4.S) #DIC=54.48263 
heidel.diag(model4.S$Sol) 
heidel.diag(model4.S$VCV)

#DIC weights
#caluclate deltaDIC
model1DIC<-model1.S$DIC 
model2DIC<-model2.S$DIC
model3DIC<-model3.S$DIC 
model4DIC<-model4.S$DIC 

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

save.image(file = "Susceptibility.RData")

mean(tree.summary$dist.from.elegans)
max(tree.summary$dist.from.elegans)
