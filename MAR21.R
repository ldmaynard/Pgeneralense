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
summary(b1) #herbivory (neg) significant 
shapiro.test(resid(b1)) #normal!!

all.dat$prop_herb1<-all.dat$prop_herb+0.0001 #adding to avoid zero-inflation
b2<-betareg(prop_herb1~treat+prop_dw+prop_gro+stage, dat=all.dat)
summary(b2) #growth (neg), young lvs (neg) sig
shapiro.test(resid(b2)) #normal!

b3<-betareg(prop_dw~treat+prop_herb+prop_gro+stage, dat=all.dat)
summary(b3)#stage sig
shapiro.test(resid(b3)) #normal

##Plots

#GROWTH + HERBIVORY
ggplot(all.dat, aes(prop_herb1, prop_gro))+
	geom_smooth(color="black",method = "lm")+
	geom_jitter(position=position_jitter(width = 0.0), alpha=0.30)+
	theme_classic()+
	theme(legend.position = "top",
		  text = element_text(size=15))+
	labs(x = "Proportion leaf herbivorized", y = "Proportion change in height")
plot(all.dat$prop_gro~all.dat$prop_herb1)
summary(lm(all.dat$growth~all.dat$prop_herb1))
#R2=0.12 

#HERBIVORY + GROWTH
ggplot(all.dat, aes(prop_gro, prop_herb))+
	geom_smooth(color="black",method = "lm")+
	geom_jitter(position=position_jitter(width = 0.0), alpha=0.30)+
	theme_classic()+
	theme(legend.position = "top",
		  text = element_text(size=15))+
	labs(x = "Proportion change in height", y = "Proportion leaf herbivorized")
plot(all.dat$prop_herb~all.dat$prop_gro)
summary(lm(all.dat$prop_herb~all.dat$prop_gro))
#R2=0.10

#HERBIVORY + LEAF AGE
ggplot(all.dat, aes(stage, prop_herb))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.40)+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=15))+
	labs(x = "Leaf age", y = "Proportion leaf herbivorized")

#PHENOLICS + LEAF AGE
ggplot(all.dat, aes(stage, prop_dw))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.40)+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=15))+
	labs(x = "Leaf age", y = "Total phenolics (%dw in gallic acid equivalents)")

##spreading the data
library(data.table)
chem_herb_spread<-dcast(setDT(all.dat), chamber ~ stage, 
						value.var = c('treat', 'prop_dw', 'prop_herb'))
chem_herb_spread<-chem_herb_spread[,-2]#remove one of treatment cols??
colnames(chem_herb_spread)[2] <- "treat"

#combine with growth data
all.dat2 <- cbind(chem_herb_spread, growth = grow$rel_gro) 

#young leaves
b4<-betareg(growth~treat+prop_herb_Young+prop_dw_Young, dat=all.dat2)
summary(b4) #herb marg sig (neg)

all.dat2$prop_herb_Young1<-all.dat2$prop_herb_Young+0.00001
b5<-betareg(prop_herb_Young1~treat+growth+prop_dw_Young, dat=all.dat2)
summary(b5) #no sig

b6<-betareg(prop_dw_Young~growth+treat+prop_herb_Young, dat=all.dat2)
summary(b6) #no sig

#mature leaves
b7<-betareg(growth~treat+prop_herb_Mature+prop_dw_Mature, dat=all.dat2)
summary(b7) #herb sig (neg)

all.dat2$prop_herb_Mature1<-all.dat2$prop_herb_Mature+0.00001
b8<-betareg(prop_herb_Mature1~treat+growth+prop_dw_Mature, dat=all.dat2)
summary(b8) #T+CO2 treat (pos), growth (neg), and phenolics (neg) all significant

b9<-betareg(prop_dw_Mature~growth+treat+prop_herb_Mature, dat=all.dat2)
summary(b9) #no sig


##If separated, young and mature leaves still show decrease in growth with increased
##herbivory.
#However,  mature leaf group had sig effects with herbivory

#HERBIVORY + GROWTH
ggplot(all.dat2, aes(growth, prop_herb_Mature1))+
	geom_smooth(color="black",method = "logit")+
	geom_jitter(position=position_jitter(width = 0.0), alpha=0.30)+
	theme_classic()+
	theme(legend.position = "top",
		  text = element_text(size=15))+
	labs(x = "Proportion change in height", y = "Proportion leaf herbivorized")
plot(all.dat2$prop_herb_Mature1~all.dat2$growth)
summary(betareg(all.dat2$prop_herb_Mature1~all.dat2$growth))
#R2=0.13

#HERBIVORY + PHENOLICS
ggplot(all.dat2, aes(prop_dw_Mature, prop_herb_Mature1))+
	geom_smooth(color="black",method = "lm")+
	geom_jitter(position=position_jitter(width = 0.0), alpha=0.30)+
	theme_classic()+
	theme(legend.position = "top",
		  text = element_text(size=15))+
	labs(x = "Total phenolics (prop. dw in GAE)", y = "Proportion leaf herbivorized")
plot(all.dat2$prop_herb_Mature1~all.dat2$prop_dw_Mature)
summary(betareg(all.dat2$prop_herb_Mature1~all.dat2$prop_dw_Mature))


bmod<-betareg(prop_herb_Mature1~prop_dw_Mature, dat=all.dat2)
library(rcompanion)
plotPredy(data  = all.dat2,
		  y     = prop_herb_Mature1,
		  x     = prop_dw_Mature,
		  model = bmod,
		  xlab  = "Phenolics",
		  ylab  = "Proportion herb")

#HERBIVORY + TREATMENT
ggplot(all.dat2, aes(treat, prop_herb_Mature1))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.40)+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=15))+
	labs(x = "Treatment", y = "Proportion herbivory")+
	stat_summary(geom = 'text', label = c("A","AB","AB","B"),
				 fun = max, vjust = -1.5, size = 5.5)+
	scale_y_continuous(limits = c(0, 0.25))+
	theme(text = element_text(size=18))

ggplot(data=all.dat2, aes(x=treat, y=prop_herb_Mature1))+ 
	geom_point()+
	stat_summary(fun.data = "mean_se", colour="blue", size=1)+
	theme_classic()

#summary(glht(b8, linfct=mcp(treat="Tukey")))

library(emmeans)
d1<-emmeans(b8,pairwise~treat, type="response")
cld(d1$emmeans,  Letters ='ABCDEFGHIJKLMNOPQRS')

#HERBIVORY + GROWTH
ggplot(all.dat2, aes(growth, prop_herb_Mature1))+
	geom_smooth(color="black",method = "lm")+
	geom_jitter(position=position_jitter(width = 0.0), alpha=0.30)+
	theme_classic()+
	theme(legend.position = "top",
		  text = element_text(size=15))+
	labs(x = "Proportion change in height", y = "Proportion leaf herbivorized")

bmod2<-betareg(prop_herb_Mature1~growth, dat=all.dat2)
plotPredy(data  = all.dat2,
		  y     = prop_herb_Mature1,
		  x     = growth,
		  model = bmod2,
		  xlab  = "Phenolics",
		  ylab  = "Proportion herb")

#combining all

#aggregate herbivory by plant/chamber
herb_20<-aggregate(prop_herb~chamber+treatment,data=pg1,FUN=mean)

#aggregate phenolics by plant/chamber
phen_ag20<-aggregate(pdw~chamber+treat,data=phen_ag2,FUN=mean)

phen_ag20 <- phen_ag20[order(phen_ag20$chamber),]
herb_20 <- herb_20[order(herb_20$chamber),]
grow <- grow[order(grow$Casa),]
herb_gro1 <- cbind(herb_20, prop_gro = grow$rel_gro) 
all.dat3 <-cbind(herb_gro1, pdw = phen_ag20$pdw)
all.dat3$prop_dw<-all.dat3$pdw/100 #proportion phenolics

#all leaves
b10<-betareg(prop_gro~treatment+prop_herb+prop_dw, dat=all.dat3)
summary(b10) #herb sig (neg)

b11<-betareg(prop_herb~treatment+prop_gro+prop_dw, dat=all.dat3)
summary(b11) #CO2 treat (neg), growth (neg), and phenolics (neg) all significant

b12<-betareg(prop_dw~prop_gro+treatment+prop_herb, dat=all.dat3)
summary(b12) #no sig

#so same as mature leaves, kinda

#HERBIVORY + GROWTH
ggplot(all.dat3, aes(prop_gro, prop_herb))+
	geom_smooth(color="black",method = "lm")+
	geom_jitter(position=position_jitter(width = 0.0), alpha=0.30)+
	theme_classic()+
	theme(legend.position = "top",
		  text = element_text(size=15))+
	labs(x = "Proportion change in height", y = "Proportion leaf herbivorized")


#HERBIVORY + PHENOLICS
ggplot(all.dat3, aes(prop_dw, prop_herb))+
	geom_smooth(color="black",method = "lm")+
	geom_jitter(position=position_jitter(width = 0.0), alpha=0.30)+
	theme_classic()+
	theme(legend.position = "top",
		  text = element_text(size=15))+
	labs(x = "Total phenolics (prop. dw in GAE)", y = "Proportion leaf herbivorized")


bmod3<-betareg(prop_herb~prop_dw, dat=all.dat3)
plotPredy(data  = all.dat3,
		  y     = prop_herb,
		  x     = prop_dw,
		  model = bmod3,
		  xlab  = "Phenolics",
		  ylab  = "Proportion herb")

#HERBIVORY + TREATMENT
ggplot(all.dat3, aes(treatment, prop_herb))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.40)+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=15))+
	labs(x = "Treatment", y = "Proportion herbivory")+
	stat_summary(geom = 'text', label = c("AB","AB","A","B"),
				 fun = max, vjust = -1.5, size = 5.5)+
	scale_y_continuous(limits = c(0, 0.16))+
	theme(text = element_text(size=18))

ggplot(data=all.dat3, aes(x=treatment, y=prop_herb))+ 
	geom_point()+
	stat_summary(fun.data = "mean_se", colour="blue", size=1)+
	theme_classic()

#summary(glht(b11, linfct=mcp(treatment="Tukey")))

d2<-emmeans(b11,pairwise~treatment, type="response")
cld(d2$emmeans,  Letters ='ABCDEFGHIJKLMNOPQRS')
#doesn't make sense...

#HERBIVORY + GROWTH
ggplot(all.dat3, aes(prop_gro, prop_herb))+
	geom_smooth(color="black",method = "lm")+
	geom_jitter(position=position_jitter(width = 0.0), alpha=0.30)+
	theme_classic()+
	theme(legend.position = "top",
		  text = element_text(size=15))+
	labs(x = "Proportion change in height", y = "Proportion leaf herbivorized")

bmod4<-betareg(prop_herb~prop_gro, dat=all.dat3)
plotPredy(data  = all.dat3,
		  y     = prop_herb,
		  x     = prop_gro,
		  model = bmod4,
		  xlab  = "Phenolics",
		  ylab  = "Proportion herb")
