#Interactive effects of climate change, leaf age, and secondary metabolites 
#on plant growth, defense, and herbivory.

#LOADING LIBRARIES----------------------------------------------------
library(lme4)
library(ggplot2)
library(car)
library(multcomp)
library(plyr)
library(AICcmodavg)
library(betareg)
library(MuMIn)
library(car)
library(corrplot)

#LOADING DATASETS----------------------------------------------------

##GROWTH

#load growth data
grow <- read.csv(file="Piper_growth.csv",head=TRUE)
colnames(grow)[1] <- "Casa"

grow$total_gro <- grow$ht.2018.09.cm-grow$ht.2018.04.cm
grow$rel_gro<-grow$total_gro/grow$ht.2018.04.cm
grow$per_gro<-grow$rel_gro*100

grow <- grow[order(grow$Treatment),]
grow<-grow[-c(11:15),]#removing control (no chamber)

hist(grow$rel_gro)  #SW: not normal, this is a typical proportion
shapiro.test(grow$rel_gro)#LM: it is technically normal
hist(grow$total_gro)
shapiro.test(grow$total_gro)#LM: but so is this one. Both yield the same results, no treatment effect

#re-ordering factor levels so lm will compare everything to control
grow$Treatment <- factor(grow$Treatment, levels=c("control chamber", "CO2", "T°C", "T°C + CO2" ))


##CHEMISTRY

#Load standard curve data
ga <- read.csv(file = "GA_StandardCurve.csv", head=T)

ga.ab<-lm(ga$abs_avg~ga$ab_val_mg)
summary(ga.ab)
plot(ga$abs_avg~ga$ab_val_mg)
#y=mx+b, y=3.511354x+0.056570, R^2=0.9974

ga.conc<-lm(ga$abs_avg~ga$concen_mgml)
summary(ga.conc)
plot(ga$abs_avg~ga$concen_mgml)
#y=mx+b, y=0.702271 +0.056570 , R^2=0.9974

##--
##Load total phenolics data
phen <- read.csv(file = "Piper_phenolics.csv", head=T)

#delete blanks/negative controls
phen <- phen[order(phen$sample),]
phen<-phen[-c(1:4),]

#create new col for abs value in well
phen$ab_val_mg<-(phen$abs_avg-0.056570)/3.511354
phen$ab_val_mg<-as.numeric(phen$ab_val_mg)
phen$abs_avg<-as.numeric(phen$abs_avg)
phen$start_wt<-as.character(phen$start_wt)
phen$start_wt<-as.numeric(phen$start_wt)

#create column for absolute value/%dw
#ab_val/0.2mL (vol in well) = y/1.1mL (total volume in tube)
phen$pdw<-(((phen$ab_val_mg/0.2)*1.1)/(phen$start_wt))*100

#create new col for concentration
phen$concen<-(phen$abs_avg-0.056570)/0.702271

#change names for leaf stages
phen$stage<-as.character(phen$stage)
phen$stage[phen$stage=="y"]="Young"
phen$stage[phen$stage=="m"]="Mature"

phen<-phen[-c(1:30),]#removing control (no chamber)
phen <- phen[order(phen$treat),]#check removed correct rows
ggplot(phen, aes(x=treat, y=concen))+geom_boxplot()+geom_point()

#average triplicate readings for concentration
phen_ag<-aggregate(concen~treat+sample+stage+chamber, data=phen, FUN=mean)

#combine triplicate readings for %dw
phen_ag2<-aggregate(pdw~treat+sample+stage+chamber, data=phen, FUN=mean)
phen_ag$pdw<-(phen_ag$pdw)*100

#re-ordering factors
phen_ag2$treat <- factor(phen_ag2$treat, levels=c("control_chamber", "CO2", "TC", "TC+CO2" ))

##HERBIVORY

#load herbivory data
pg <- read.csv(file="Piper_herbivory.csv",head=TRUE)

#creating col for proportion herbivory
pg$percent_herbivory<-as.numeric(pg$percent_herbivory)
pg$prop_herb<-(pg$percent_herbivory/100)
pg$prop_herb<-as.numeric(pg$prop_herb)

pg <- pg[order(pg$treatment),]
pg1<-pg[-c(41:60),]#removing control (no chamber)
pg2<-aggregate(percent_herbivory~chamber + age + treatment,data=pg1,FUN=mean)#aggregate data

#get two datasets ready to combine
phen_ag2 <- phen_ag2[order(phen_ag2$sample),]
pg2 <- pg2[order(pg2$chamber),]

#combine aggregated herbivory and chemistry datasets
ph <- cbind(phen_ag2, percent_herbivory = pg2$percent_herbivory) 

#combining all data
ph <- ph[order(ph$stage),]
all.dat<-cbind(ph, growth = grow$total_gro) 
all.dat<-cbind(all.dat, per_gro = grow$per_gro)
#all.dat<-all.dat[order(all.dat$sample)]

#reordering
all.dat$treat <- factor(all.dat$treat, levels=c("control_chamber", "CO2", "TC", "TC+CO2" ))

#checking for colinearity
samp<-all.dat[,-c(1:4)]#cat vars
samp<-samp[,-4]#remove second growth var

L <- cor(samp)#correlation matrix
corrplot(L, method = "circle")
corrplot(L, method = "number")
#low colinearity

check <- lm(growth ~ percent_herbivory+pdw+stage+treat, data=all.dat)
summary(check)
vif(check)

#creating columns for proportions
all.dat$prop_herb<-all.dat$percent_herbivory/100
all.dat$prop_dw<-all.dat$pdw/100
all.dat$prop_gro<-all.dat$per_gro/100

all1<-lmer(prop_gro ~ treat + prop_herb + prop_dw + (1|chamber), 
		   data=all.dat, na.action = "na.fail")

b1<-betareg(prop_gro~treat+prop_herb+prop_dw, dat=all.dat)
summary(b1) #CO2 treatment (neg) and herbivory (pos) significant 
shapiro.test(resid(b1)) #normal!!

all.dat$prop_herb1<-all.dat$prop_herb+0.0001 #adding to avoid zero-inflation
b2<-betareg(prop_herb1~treat+prop_dw+prop_gro+stage, dat=all.dat)
summary(b2) #CO2+Temp treatment (pos), growth (pos), young lvs (neg) sig. phenolics marg sig
shapiro.test(resid(b2)) #normal!

b3<-betareg(prop_dw~treat+prop_herb+prop_gro+stage, dat=all.dat)
summary(b3)#no significance
shapiro.test(resid(b3)) #normal

#GROWTH + HERBIVORY
ggplot(all.dat, aes(prop_herb1, prop_gro))+
	geom_smooth(color="black",method = "lm")+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.30)+
	theme_classic()+
	theme(legend.position = "top",
		  text = element_text(size=15))+
	labs(x = "Proportion leaf herbivorized", y = "Relative growth in height")
