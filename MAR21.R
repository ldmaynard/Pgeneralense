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
library(viridis)

#LOADING & WRANGLING DATA----------------------------------------------------

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
grow$Casa<-as.factor(grow$Casa)
levels(grow$Casa)
grow <- grow[order(grow$Casa),]
all.dat<-cbind(ph, growth = grow$total_gro) 
all.dat<-cbind(all.dat, per_gro = grow$per_gro)


#creating columns for proportions
all.dat$prop_herb<-all.dat$percent_herbivory/100
all.dat$prop_dw<-all.dat$pdw/100
all.dat$prop_gro<-all.dat$per_gro/100

#write.csv(all.dat, "all.dat.csv")

#reordering
all.dat$treat <- factor(all.dat$treat, levels=c("control_chamber", "CO2", "TC", "TC+CO2" ))

##Spreading data — young and mature leaves = all.dat2
library(data.table)
chem_herb_spread<-dcast(setDT(all.dat), chamber ~ stage, 
						value.var = c('treat', 'prop_dw', 'prop_herb'))
chem_herb_spread<-chem_herb_spread[,-2]#remove one of treatment cols??
colnames(chem_herb_spread)[2] <- "treat"

#combine with growth data
all.dat2 <- cbind(chem_herb_spread, growth = grow$rel_gro) 

##Averaging all data = all.dat3
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

all.dat3$treatment<-as.factor(all.dat3$treatment)

#---

##CHECKING FOR COLINEARITY----
samp1<-all.dat[,-c(1:8)]#cat vars
samp1<-samp1[,-4]#remove second herb var

L1 <- cor(samp1)#correlation matrix
corrplot(L1, method = "circle")
corrplot(L1, method = "number")
#low colinearity

check1 <- lm(prop_gro ~ prop_herb+prop_dw+stage+treat, data=all.dat)
summary(check1)
vif(check1)


##ANALYSIS----


###GROWTH (using data averaged across all four leaves)----
b1<-betareg(prop_gro~treat+prop_herb+prop_dw, dat=all.dat)
summary(b1) #herbivory (neg) significant 
shapiro.test(resid(b1)) #normal!!

#GROWTH + HERBIVORY PLOT
ggplot(all.dat, aes(prop_herb, prop_gro))+
	geom_smooth(color="black",method = "lm")+
	geom_jitter(position=position_jitter(width = 0.0), alpha=0.30, size=2.5, aes(color=chamber))+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=19))+
	labs(x = "Proportion leaf herbivory", y = "Proportion change in height")+
	scale_color_viridis(discrete = T, option = "D")
#pseudo R^2=0.178

bmod<-betareg(prop_gro~prop_herb, dat=all.dat)
library(rcompanion)
plotPredy(data  = all.dat,
		  y     = prop_gro,
		  x     = prop_herb,
		  model = bmod,
		  xlab  = "Phenolics",
		  ylab  = "Proportion herb")


###PHENOLICS----
b3<-betareg(prop_dw~treat+stage, dat=all.dat)
summary(b3)#stage sig
shapiro.test(resid(b3)) #normal

#PHENOLICS + LEAF AGE PLOT
ggplot(data=all.dat, aes(x=stage, y=prop_dw))+ 
	geom_point(position=position_jitter(width = 0.0), alpha=0.30, aes(color=chamber), size=2.5)+
	stat_summary(fun.data = "mean_se", colour="black", size=1)+
	theme_classic()+
	labs(x="", y="Total phenolics (prop. dw in GAE)")+
	theme(text = element_text(size=18), legend.position = "none")+
	scale_color_viridis(discrete = T, option = "D")+
	geom_signif(comparisons = list(c("Mature", "Young")),
				map_signif_level = T, textsize=6)+
	scale_y_continuous(limits = c(0,.105))

#create cld
library(emmeans)
library(multcomp)
d5<-emmeans(b3,pairwise~stage, type="response")
cld(d5$emmeans,  Letters ='ABCDEFGHIJKLMNOPQRS')

ggplot(data=all.dat, aes(x=stage, y=prop_dw))+ 
	geom_point(position=position_jitter(width = 0.0), alpha=0.30, aes(color=chamber), size=2.5)+
	stat_summary(fun.data = "mean_se", colour="black", size=1)+
	theme_classic()+
	labs(x="", y="Total phenolics (prop. dw in GAE)")+
	theme(text = element_text(size=18), legend.position = "none")+
	scale_color_viridis(discrete = T, option = "D")+
	stat_summary(geom = 'text', label = c("A","B"),
				 fun = max, vjust = -1.5, size = 5.8)+
	scale_y_continuous(limits = c(0,.115))

###HERBIVORY---
#young leaves
all.dat2$prop_herb_Young1<-all.dat2$prop_herb_Young+0.00001
b5<-betareg(prop_herb_Young1~treat+growth+prop_dw_Young, dat=all.dat2)
summary(b5) #no sig, so no plots?

#mature leaves
all.dat2$prop_herb_Mature1<-all.dat2$prop_herb_Mature+0.00001
b8<-betareg(prop_herb_Mature1~treat+growth+prop_dw_Mature, dat=all.dat2)
summary(b8) #T+CO2 treat (pos), growth (neg), and phenolics (neg) all significant

#HERBIVORY + GROWTH (mature)
ggplot(all.dat2, aes(growth, prop_herb_Mature1))+
	geom_smooth(color="black",method = "glm")+
	geom_jitter(position=position_jitter(width = 0.0), alpha=0.30, aes(color=chamber), size=2.5)+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=19))+
	labs(x = "Proportion change in height", y = "Proportion leaf herbivory")+
	scale_color_viridis(discrete = T, option = "D")

#HERBIVORY + PHENOLICS (mature)
ggplot(all.dat2, aes(prop_dw_Mature, prop_herb_Mature1))+
	geom_smooth(color="black",method = "glm")+
	geom_jitter(position=position_jitter(width = 0.0), alpha=0.30, aes(color=chamber), size=2.5)+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=19))+
	labs(x = "Total phenolics (prop. dw in GAE)", y = "Proportion leaf herbivory")+
	scale_color_viridis(discrete = T, option = "D")

bmod1<-betareg(prop_herb_Mature1~prop_dw_Mature, dat=all.dat2)
plotPredy(data  = all.dat2,
		  y     = prop_herb_Mature1,
		  x     = prop_dw_Mature,
		  model = bmod1,
		  xlab  = "Phenolics",
		  ylab  = "Proportion herb")

#HERBIVORY + TREATMENT (mature)
#writing labels for plot

all.dat2 <- all.dat2[order(all.dat2$treat),]
lab1 <- c("Control", 
		  expression(CO["2"]),
		  "Temperature",
		  expression(CO["2"] + Temp))

#create cld
library(emmeans)
library(multcomp)
d2<-emmeans(b8,pairwise~treat, type="response")
cld(d2$emmeans,  Letters ='ABCDEFGHIJKLMNOPQRS')

ggplot(all.dat2, aes(treat, prop_herb_Mature1))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.40)+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=15))+
	labs(x = "Treatment", y = "Proportion leaf herbivory")+
	theme(text = element_text(size=18))+
	scale_x_discrete(labels=lab1)

ggplot(data=all.dat2, aes(x=treat, y=prop_herb_Mature1))+ 
	geom_point(aes(color=chamber),position=position_jitter(width = 0.04), alpha=0.40, size=2.5)+
	stat_summary(fun.data = "mean_se", colour="black", size=1)+
	theme_classic()+
	scale_color_viridis(discrete = T, option = "D")+
	theme(legend.position = "none",
		  text = element_text(size=20), axis.text.x = element_text(angle=20, hjust=1))+
	labs(x="", y="Proportion leaf herbivory")+
	stat_summary(geom = 'text', label = c("AB","AB","A","B"),
				 fun = max, vjust = -1.5, size = 5.8)+
	scale_y_continuous(limits = c(0, 0.25))+
	scale_x_discrete(labels=lab1)


##OTHER THINGS----

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
	geom_smooth(color="black",method = "glm")+
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
	geom_smooth(color="black",method = "glm")+
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

marginal <- emmeans(b8, ~ treat)
pairs(marginal, adjust="tukey")

Sum <- cld(marginal,
		  alpha = 0.05,
		  Letters = letters, ### Use lowercase letters for .group
		  adjust = "tukey") ### Tukey-adjusted comparisons
Sum

marginal2 <- emmeans(b11, ~ treatment)
pairs(marginal2, adjust="tukey")

Sum2 <- cld(marginal2,
		   alpha = 0.05,
		   Letters = letters, ### Use lowercase letters for .group
		   adjust = "tukey") ### Tukey-adjusted comparisons
Sum2

#The summary function in betareg produces a pseudo R-squared value for the model, and the
#recommended test for the p-value for the model is the lrtest function in the lmtest package.
library(lmtest)
lrtest(b8)
lrtest(b11)



