setwd("C:/Users/cls6630/OneDrive - The Pennsylvania State University/Spillover")

library(ggplot2)
library(dplyr)
library(tidyr)
library(emdbook)

dat<-read.csv("2021_11_29_MJOrsayTM_FCK_FCH_GAA.manip.csv", header=TRUE, fill=TRUE, na.strings = "", stringsAsFactors = FALSE)
dat<-filter(dat, Experiment=="FCK")
dat$Ct<-as.numeric(dat$Ct)
dat$Dose<-as.numeric(dat$Dose)
dat$Ct[is.na(dat$Ct)]<-Inf

Fig<- ggplot(dat, aes(Dose, Ct, color=Strain))+geom_jitter(width = 0.02)+scale_y_reverse()+
  scale_x_log10()+xlab("Dilution of virus filtrate")+theme_bw()

save_plot("SpilloverSupplementFig1.jpg",Fig, base_height = 4.5)

#Figure out maximum likelihood estimate of dose of the stock (how many multiples of the TCID50 is it?)
dat$infected[dat$Ct==Inf]<-0
dat$infected[dat$Ct!=Inf]<-1
dat.group<-dat%>%group_by(Dose, Strain)%>%summarise(number.inf=sum(infected))
dat.group.JU1580<-filter(dat.group, Strain=="JU1580")

#The data of the number of infected at each diultion of the stock goes in here.
#The chance a plate will be infected 
SEQ<-seq(from=0.1, to=100000, by=2)
LIKE<-NULL
for (x in SEQ){
  like<-dbinom(dat.group.JU1580$number.inf[1],4,1-(.5)^(dat.group.JU1580$Dose[1]*x))*dbinom(dat.group.JU1580$number.inf[2],4,1-(.5)^(dat.group.JU1580$Dose[2]*x))*
    dbinom(dat.group.JU1580$number.inf[3],4,1-(.5)^(dat.group.JU1580$Dose[3]*x))*dbinom(dat.group.JU1580$number.inf[4],4,1-(.5)^(dat.group.JU1580$Dose[4]*x))*
    dbinom(dat.group.JU1580$number.inf[5],4,1-(.5)^(dat.group.JU1580$Dose[5]*x))*dbinom(dat.group.JU1580$number.inf[6],4,1-(.5)^(dat.group.JU1580$Dose[6]*x))*
    dbinom(dat.group.JU1580$number.inf[7],4,1-(.5)^(dat.group.JU1580$Dose[7]*x))*dbinom(dat.group.JU1580$number.inf[8],4,1-(.5)^(dat.group.JU1580$Dose[8]*x))
  output<-c(x, like)
  LIKE<-rbind(output, LIKE)
}

LIKE<-as.data.frame(LIKE)
colnames(LIKE)<-c("inf.dose.50","likelihood")
rownames(LIKE)<-NULL
LIKE$loglikelihood<-log(LIKE$likelihood)

Fig2<-ggplot(LIKE, aes(inf.dose.50, loglikelihood))+geom_line()+
  geom_hline(aes(yintercept=(max(loglikelihood)-1.92)))+
  xlab("TCID50 of per 20 microliters virul filtrate")+
  scale_x_continuous(limits=c(0,25000))+
  scale_y_continuous(limits=c(-20, 0))+theme_bw()

save_plot("SpilloverSupplementFig2.jpg",Fig2, base_height = 4.5)

which.max(LIKE$loglikelihood) #maximum likelihood estimate for the amount our stock is more concentrated than the ID50
#45719
LIKE$inf.dose.50[45719] #8562.1
max(LIKE$loglikelihood)-1.92

small<-filter(LIKE, loglikelihood>-4.99)
#19446
#3468