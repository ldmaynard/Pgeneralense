#Effects of climate change on allocation to growth and defense in a neotropical shrub

#LOADING LIBRARIES---
{library(lme4)
	library(ggplot2)
	library(car)
	library(dplyr)
	library(betareg)
	library(corrplot)
	library(viridis)
	library(emmeans)
	library(ggpubr)
	library(glmmTMB)
	library(RColorBrewer)
	library(multcomp)
	library(plyr)
}
#LOADING & WRANGLING DATA----------------------------------------------------
##GROWTH DATA----

#load growth data
grow <- read.csv(file="Piper_growth.csv",head=TRUE)

#renaming chamber column
colnames(grow)[1] <- "chamber"

#creating columns for: 
#total growth
grow$total_gro <- grow$ht.2018.09.cm-grow$ht.2018.04.cm
#proportion (relative) growth
grow$prop_gro<-grow$total_gro/grow$ht.2018.04.cm
#percent growth
grow$per_gro<-grow$prop_gro*100

#renaming dataset with all treatments, n=25
grow.all<-grow

#creating dataset without "No chamber" treatment, n=20 
grow <- grow[order(grow$Treatment),]
grow<-grow[-c(11:15),] #removing control (no chamber)


##CHEMISTRY DATA----

##Stanard curve math
{#Load standard curve data
	ga <- read.csv(file = "GA_StandardCurve.csv", head=T)
	
	ga.ab<-lm(ga$abs_avg~ga$ab_val_mg)
	summary(ga.ab)
	plot(ga$abs_avg~ga$ab_val_mg)
	#y=mx+b, y=3.511354x+0.056570, R^2=0.9974
	
	ga.conc<-lm(ga$abs_avg~ga$concen_mgml)
	summary(ga.conc)
	plot(ga$abs_avg~ga$concen_mgml)
	#y=mx+b, y=0.702271 +0.056570 , R^2=0.9974
}

##--
##Load total phenolics data
phen <- read.csv(file = "Piper_phenolics.csv", head=T)

#delete blanks/negative controls
phen <- phen[order(phen$sample),]
phen<-phen[-c(1:4),] 

#create new col for abs value in well, using standard curve results from above
phen$ab_val_mg<-(phen$abs_avg-0.056570)/3.511354

#defining data types
{phen$ab_val_mg<-as.numeric(phen$ab_val_mg)
	phen$abs_avg<-as.numeric(phen$abs_avg)
	phen$start_wt<-as.character(phen$start_wt)
	phen$start_wt<-as.numeric(phen$start_wt)
}

#create column for absolute value/%dw
#ab_val/0.2mL (vol in well) = y/1.1mL (total volume in tube)
phen$per_dw<-(((phen$ab_val_mg/0.2)*1.1)/(phen$start_wt))*100
phen$pdw<-(phen$per_dw)/100

#create new col for concentration
phen$concen<-(phen$abs_avg-0.056570)/0.702271

#change names for leaf stages
phen$stage<-as.character(phen$stage)
phen$stage[phen$stage=="y"]="Young"
phen$stage[phen$stage=="m"]="Mature"

#average triplicate readings for concentration
phen.all_conc<-aggregate(concen~treat+sample+stage+chamber, data=phen, FUN=mean)

#average triplicate readings for %dw, n=50
phen.all<-aggregate(pdw~treat+sample+stage+chamber, data=phen, FUN=mean)

#aggregate by chamber, n=25
phen.25<-aggregate(pdw~treat+chamber, data=phen.all, FUN=mean)

#create dataset without control (no chamber), n=40
phen.40<-phen.all[-c(1:10),]#removing control (no chamber)

#aggregating by chamber, n=20
phen.20<-aggregate(pdw~chamber+treat,data=phen.40,FUN=mean)


##HERBIVORY DATA----

#load herbivory data
herb.all <- read.csv(file="Piper_herbivory.csv",head=TRUE)

#creating col for proportion herbivory
herb.all$percent_herbivory<-as.numeric(herb.all$percent_herbivory)
herb.all$prop_herb<-(herb.all$percent_herbivory/100)
herb.all$prop_herb<-as.numeric(herb.all$prop_herb)

#creating dataset aggregated by leaf age  by age (n=40)
herb.50<-aggregate(prop_herb~chamber + age + treatment +sample,data=herb.all,FUN=mean)

#aggregating by chamber, n=25
herb.25<-aggregate(prop_herb~chamber+treatment,data=herb.all,FUN=mean)

#creating dataset without control (no chamber) = herb.80 (n=80)
herb.all <- herb.all[order(herb.all$treatment),]
herb.80<-herb.all[-c(41:60),]#removing control (no chamber)

#creating dataset with average herbivory  by age (n=40)
herb.40<-aggregate(prop_herb~chamber + age + treatment+sample,data=herb.80,FUN=mean)

#aggregating by chamber, n=20
herb.20<-aggregate(prop_herb~chamber+treatment,data=herb.80,FUN=mean)

##CREATING COMBINED DATASETS----

#four treatment data, not aggregated, n=80
all.dat80<-merge(herb.80, phen.40, by="sample", all = T)
head(all.dat80)
colnames(all.dat80)[4] <- "chamber"
all.dat80<-merge(all.dat80, grow, by="chamber", all=T)
#cleaning up, selecting columns
all.dat80<-subset(all.dat80, select = -c(11,13,15))

#adding small amount to herbivory to help beta regressions run with high number of zeros
all.dat80$prop_herb1<-all.dat80$prop_herb+0.00001

#four treatment data, aggregated by leaf age, n=40
all.dat40<-merge(herb.40, phen.40, by="sample", all = T)
head(all.dat40)
colnames(all.dat40)[2] <- "chamber"
all.dat40<-merge(all.dat40, grow, by="chamber", all=T)
#cleaning up, selecting columns
all.dat40<-subset(all.dat40, select = -c(8,10,15))

#four treatment data, aggregated by chamber, n=20
all.dat20<-merge(herb.20, phen.20, by="chamber", all = T)
all.dat20<-merge(all.dat20, grow, by="chamber", all=T)
#cleaning up, selecting columns
all.dat20<-subset(all.dat20, select = -c(4,6))

#for four treatment data, re-ordering factor levels so model will compare everything to control chamber
all.dat80$treatment<-as.factor(all.dat80$treatment)
levels(all.dat80$treatment)
all.dat80$treat <- factor(all.dat80$treat, levels=c("control chamber", "CO2", "TC", "TC + CO2" ))

all.dat40$treatment<-as.factor(all.dat40$treatment)
levels(all.dat40$treatment)
all.dat40$treat <- factor(all.dat40$treat, levels=c("control chamber", "CO2", "TC", "TC + CO2" ))

all.dat20$treatment<-as.factor(all.dat20$treatment)
levels(all.dat20$treatment)
all.dat20$treat <- factor(all.dat20$treat, levels=c("control chamber", "CO2", "TC", "TC + CO2" ))

#---

##CHECKING FOR COLINEARITY----
samp1<-all.dat80[,c(10,12,16)]#cat vars

L1 <- cor(samp1)#correlation matrix
corrplot(L1, method = "circle")
corrplot(L1, method = "number")
#low colinearity

check1 <- lm(prop_gro ~ prop_herb+pdw+stage+treatment, data=all.dat80)
summary(check1)
vif(check1)

##PRE-QUESTIONS----
##1) Is there a chamber effects between control treatments?
#Growth
#creating datasets with only two control treatments
grow.all <- grow.all[order(grow.all$Treatment),]
grow.control.c <-grow.all[c(6:10),]
grow.control.nc <-grow.all[c(11:15),]
t.test(grow.control.c$prop_gro, grow.control.nc$prop_gro, paired = F)
#t = -1.2361, df = 5.5473, p-value = 0.2662
#no chamber effect on growth

#Phenolics
#creating datasets with only two control treatments
phen.all <- phen.all[order(phen.all$treat),]
phen.all.cc <- phen.all[c(11:20),]
phen.all.nc <- phen.all[c(21:30),]
t.test(phen.all.cc$pdw, phen.all.nc$pdw, paired = F)
#t = -1.2182, df = 17.997, p-value = 0.2389

#Herbivory
#creating datasets with only two control treatments
herb.all <- herb.all[order(herb.all$treatment),]
herb.all.cc <- herb.all[c(21:40),]
herb.all.nc <- herb.all[c(41:60),]
t.test(herb.all.cc$prop_herb, herb.all.nc$prop_herb, paired = F)
#t = -1.089, df = 27.355, p-value = 0.2857

##2) Does leaf age affect phenolics and/or chemical defense?
#splitting phenolics data by leaf age
all.dat40<- all.dat40[order(all.dat40$stage),]
all.dat40_m <- all.dat40[c(1:20),]
all.dat40_y <- all.dat40[c(21:40),]

#splitting herbivory data by leaf age
all.dat80<- all.dat80[order(all.dat80$stage),]
all.dat80_m <- all.dat80[c(1:40),]
all.dat80_y <- all.dat80[c(41:80),]

t.test(all.dat40_m$pdw, all.dat40_y$pdw, paired = F)
#t = -4.8008, df = 37.647, p-value = 2.516e

t.test(all.dat80_m$prop_herb, all.dat80_y$prop_herb, paired = F)
#t = 3.2495, df = 62.412, p-value = 0.001864

#FIGURE A1a, PHENOLICS BY LEAF AGE
figA1a<-ggplot(data=all.dat40, aes(x=stage, y=pdw))+ 
	geom_boxplot()+
	geom_point(position=position_jitter(width = 0.025), alpha=0.6, size=2.5)+
	stat_summary(fun.data = "mean_se", colour="red", size=1)+
	theme_classic()+
	theme(legend.position = "none")+
	labs(y="Total phenolics (prop dw in GAE)", x="")+
	theme(text = element_text(size=16), axis.text.x = element_blank())
figA1a

#FIGURE A1b, HERBIVORY BY LEAF AGE
figA1b<-ggplot(data=all.dat80, aes(x=stage, y=prop_herb1))+ 
	geom_boxplot()+
	geom_point(position=position_jitter(width = 0.025), alpha=0.6, size=2.5)+
	stat_summary(fun.data = "mean_se", colour="red", size=1)+
	theme_classic()+
	theme(legend.position = "none")+
	labs(y="Proportion herbivory", x="Leaf age")+
	theme(text = element_text(size=16), axis.text.x = element_text(size=12))
figA1b

#Creating combined S1 figure
tiff('Maynard_etal_A1Fig.tiff', units="in", width=6, height=9, res=300)
ggarrange(figA1a, figA1b,
		  labels = c("a", "b"),heights = c(2, 2.2),
		  ncol = 1, nrow = 2)
dev.off()

#PHENOLICS SUMMARY STATS
phen.tab <- ddply(all.dat40, c("stage"), summarise,
				  N    = length(pdw),
				  mean = mean(pdw),
				  sd   = sd(pdw),
				  se   = sd / sqrt(N))
phen.tab
#young leaf avg pdw/old leaf avg pdw-old leaf avg pdw
(0.06859805-0.04970515)/0.04970515
#0.380099
#Young leaves had an average of 38% more total phenolics than mature leaves


#HERBIVORY SUMMARY STATS
herb.tab <- ddply(all.dat80, c("stage"), summarise,
				  N    = length(prop_herb),
				  mean = mean(prop_herb),
				  sd   = sd(prop_herb),
				  se   = sd / sqrt(N))
herb.tab
#Mature leaf avg pdw/young leaf avg pdw-young leaf avg pdw
(0.0528275-0.0103100)/0.0103100
#4.123909
#Mature leaves had an average of over four times more herbivory than young leaves


#QUESTION 1----
#Q1A GROWTH

gro.mod <- betareg(prop_gro ~ treatment, data = all.dat20)
Anova(gro.mod)
#treatment p=0.79

#Q1B PHENOLICS
#Ordering data so the models compare to the control treatment
#mature leaves data
all.dat40_m$treatment<-as.factor(all.dat40_m$treatment)
levels(all.dat40_m$treatment)
all.dat40_m$treatment <- factor(all.dat40_m$treatment, levels=c("control chamber", "CO2", "TC", "TC + CO2" ))

#young leaves data
all.dat40_y$treatment<-as.factor(all.dat40_y$treatment)
levels(all.dat40_y$treatment)
all.dat40_y$treatment <- factor(all.dat40_y$treatment, levels=c("control chamber", "CO2", "TC", "TC + CO2" ))

#analysis
#mature leaves
chem.mod.m <- betareg(pdw ~ treatment, data = all.dat40_m)
Anova(chem.mod.m)
#treatment p=0.70

#young leaves
chem.mod.y <- betareg(pdw ~ treatment, data = all.dat40_y)
Anova(chem.mod.y)
#treatment p=0.57


#Q1C HERBIVORY
#Ordering data so the models compare to the control treatment
{#mature leaves data
all.dat80_m$treatment<-as.factor(all.dat80_m$treatment)
levels(all.dat80_m$treatment)
all.dat80_m$treatment <- factor(all.dat80_m$treatment, levels=c("control chamber", "CO2", "TC", "TC + CO2" ))
#young leaves data
all.dat80_y$treatment<-as.factor(all.dat80_y$treatment)
levels(all.dat80_y$treatment)
all.dat80_y$treatment <- factor(all.dat80_y$treatment, levels=c("control chamber", "CO2", "TC", "TC + CO2" ))
}

#analysis, mixed model
#mature leaves
herb.mod.m <- glmmTMB(prop_herb1 ~ treatment + (1|chamber), data = all.dat80_m, family = "beta_family")
Anova(herb.mod.m)
#treatment p=0.56

#young leaves
herb.mod.y <- glmmTMB(prop_herb1 ~ treatment + (1|chamber), data = all.dat80_y, family = "beta_family")
Anova(herb.mod.y)
#treatment p=0.95

#analysis, regular betareg with averaged leaves by age
#adding small amount to herbivory to help beta regressions run with high number of zeros
all.dat40_y$prop_herb1<-all.dat40_y$prop_herb+0.0001
all.dat40_m$prop_herb1<-all.dat40_m$prop_herb+0.0001

#mature leaves
herb2.mod.m <- betareg(prop_herb1 ~ treatment, data = all.dat40_m)
Anova(herb2.mod.m)
#treatment p=0.80

#young leaves
herb2.mod.y <- betareg(prop_herb1 ~ treatment, data = all.dat40_y)
Anova(herb2.mod.y)
#treatment p=0.94

#QUESTION 2---- 
#Defense~growth * treatment

#Models
#mature leaves
mod.gro2m<-(betareg(pdw~treatment*prop_gro, dat=all.dat40_m))
Anova(mod.gro2m)

library(effects)
summary(allEffects(mod.gro2m))
plot(allEffects(mod.gro2m))

all.dat40_m <- all.dat40_m[order(all.dat40_m$treatment),]

CO2<-slice(all.dat40_m, 1:5) #CO2
control<-slice(all.dat40_m, 6:10) #control
TC<-slice(all.dat40_m, 11:15) #temperature
combo<-slice(all.dat40_m, 16:20) #combo treatment

modco2<-(betareg(pdw~prop_gro, dat=CO2))
modcon<-(betareg(pdw~prop_gro, dat=control))
modtc<-(betareg(pdw~prop_gro, dat=TC))
modcombo<-(betareg(pdw~prop_gro, dat=combo))

summary(modco2)#z=-0.44,p=0.663
summary(modcon)#z=-0.54, p=0.59
summary(modtc)#z=3.44, p=0.000591; estimate=2.0
#total phenolics increase 2% with every 1% increase in proportion growth in warmer environments
summary(modcombo)#z=0.001, p=0.999


#young leaves
mod.gro2y<-(betareg(pdw~treatment*prop_gro, dat=all.dat40_y))
Anova(mod.gro2y)
#interaction, chisq=8.31, p=0.040, significant

library(effects)
summary(allEffects(mod.gro2y))
plot(allEffects(mod.gro2y))

all.dat40_y <- all.dat40_y[order(all.dat40_y$treatment),]

CO2y<-slice(all.dat40_y, 1:5) #CO2
controly<-slice(all.dat40_y, 6:10) #control
TCy<-slice(all.dat40_y, 11:15) #temperature
comboy<-slice(all.dat40_y, 16:20) #combo treatment

modco2y<-(betareg(pdw~prop_gro, dat=CO2y))
modcony<-(betareg(pdw~prop_gro, dat=controly))
modtcy<-(betareg(pdw~prop_gro, dat=TCy))
modcomboy<-(betareg(pdw~prop_gro, dat=comboy))

summary(modco2y)#z=0.325,p=0.745
summary(modcony)#z=-1.217, p=0.224
summary(modtcy)#z=3.598, p=0.000321; estimate=0.48393
#total phenolics increase 0.5% with every 1% increase in proportion growth in warmer environments
summary(modcomboy)#z=5.183, p=2.19e-07; estimate=0.86345
##total phenolics increase 0.9% with every 1% increase in proportion growth in combo environments

#QUESTION 3----
#Climate change effects on herbivory and the effectiveness of defense

#Mature leaves
herb.chem.mod2_m<-glmmTMB(prop_herb1 ~ treatment * pdw + (1|chamber), data = all.dat80_m, family = "beta_family")
Anova(herb.chem.mod2_m)
#treatment  chisq=3.31, p=0.35
#chemistry  chisq=3.97, p=0.0464
#interaction chisq=1.01, p=0.80
#in mature leaves, increasing defenses had a negative effect on herbivory 

#Young leaves
herb.chem.mod2_y<-glmmTMB(prop_herb1 ~ treatment * pdw + (1|chamber), data = all.dat80_y, family = "beta_family")
Anova(herb.chem.mod2_y)
#nothing significant
#treatment, Chisq=0.53, p=0.91
#chemistry, chisq=0.24, p=0.63
#interaction, chisq=1.04, p=0.79
#in young leaves, increases chemical defense did not have an effect on herbivory

summary(herb.chem.mod2_m)
#Mature leaves  experienced ~20% less herbivory with every 1% increase in total phenolics

#non-mixed model (leaves averaged by age)
mod.herb2m<-(betareg(prop_herb1~treatment*pdw, dat=all.dat40_m))
Anova(mod.herb2m)
#treatment  chisq=3.09, p=0.38
#chemistry  chisq=3.72, p=0.054
#interaction chisq=4.70, p=0.20
#In mature leaves, increasing defences only has a marg sign neg effect on herbivory

mod.herb2y<-(betareg(prop_herb1~treatment*pdw, dat=all.dat40_y))
Anova(mod.herb2y)
#treatment  chisq=0.60, p=0.90
#chemistry  chisq=0.48, p=0.49
#interaction chisq=2.78, p=0.43

##GRAPHS Q1-3----

##Q1 GRAPH----
all.dat20$plant<-NA
all.dat20$plant<-"Whole plant"

gro.plot<-ggplot(data=all.dat20, aes(x=treatment, y=prop_gro, color=treatment))+ 
	geom_point(shape=18,aes(shape=plant),position=position_jitter(width = 0.025), alpha=0.6, 
			   size=5, show.legend = F)+
	stat_summary(aes(shape=plant),fun.data = "mean_se", size=1.1, color="black")+
	theme_classic()+
	labs(y="Proportion growth", x="")+
	theme(axis.text.x = element_blank(), 
		  text = element_text(size=16), legend.position = "right", legend.title = element_blank())+
	scale_x_discrete(limits=c("control chamber", "CO2", "TC", "TC + CO2"))+
	scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
	scale_color_brewer(palette=9, type = "div", direction = 1)
gro.plot

chem.plot<-ggplot(data=all.dat40, aes(x=treatment, y=pdw, color=treatment))+ 
	geom_point(aes(shape=stage),position=position_jitter(width = 0.025), alpha=0.6, 
			   size=4, show.legend = F)+
	stat_summary(aes(group=stage, shape=stage),fun.data = "mean_se", size=1.1, color="black")+
	theme_classic()+
	labs(y="Total phenolics (pdw GAE)", x="")+
	theme(axis.text.x = element_blank(),
		  text = element_text(size=15), legend.title = element_blank(),
		  legend.position = "right")+
	scale_color_brewer(palette=9, type = "div", direction = 1)+
	scale_x_discrete(limits=c("control chamber", "CO2", "TC", "TC + CO2"))
chem.plot

herb.plot<-ggplot(data=all.dat80, aes(x=treatment, y=prop_herb1, color=treatment))+ 
	geom_point(aes(shape=stage),position=position_jitter(width = 0.025), alpha=0.6, 
			   size=4, show.legend = F)+
	stat_summary(aes(group=stage, shape=stage, color=stage),fun.data = "mean_se", colour="black", size=1.1)+
	theme_classic()+
	labs(y="Proportion herbivory", x="")+
	theme(text = element_text(size=16), axis.text.x = element_text(angle=20, hjust=0.9, size=16),
		  legend.title = element_blank(), legend.position = "right")+
	scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
	scale_color_brewer(palette=9, type = "div", direction = 1)+
	scale_x_discrete(limits=c("control chamber", "CO2", "TC", "TC + CO2" ),
					 labels=c("control chamber"=expression(atop("Control")),
					 		 "TC"="Temperature",
					 		 CO2=expression(CO["2"]),
					 		 "TC + CO2"=expression(CO["2"] + Temp)))
herb.plot


tiff('Maynard_etal_Fig1.tiff', units="in", width=7, height=10, res=300)
ggarrange(gro.plot, chem.plot, herb.plot,
		  labels = c("a", "b", "c"),heights = c(1.5,1.5, 1.8),
		  ncol = 1, nrow = 3,
		  label.y = 1.025)
dev.off()

##Q2 GRAPH----

#legend text
leg.lab <- c("Control",
			 expression(CO["2"]),
			 "Temperature",
			 expression(CO["2"] + Temp))

#Growth*defense plot, mature leaves
gd_mature<-all.dat40_m %>%
	ggplot(aes(x=prop_gro, 
			   y=pdw,
			   color=treatment))+
	geom_smooth(aes(linetype=treatment),method="lm", se=F, size=1.8, show.legend=F)+
	geom_point(alpha=0.5, size=2.8)+
	labs(y="Total phenolics (prop dw GAE)", x="", title="Mature leaves")+
	theme_classic()+
	theme(legend.title = element_blank(), text = element_text(size=20), legend.position = c(0.88,0.88),
		  legend.text.align = 0,legend.spacing.y = unit(0, "mm"),
		  legend.box.background = element_rect(colour = "black"), plot.title = element_text(hjust = 0.5))+
	scale_linetype_manual(values = c("solid", "solid", "solid","solid"))+
	scale_color_brewer(palette=9, type = "div", direction = 1, labels=leg.lab)

gd_mature


gd_young<-all.dat40_y %>%
	ggplot(aes(x=prop_gro, 
			   y=pdw,
			   color=treatment))+
	geom_smooth(aes(linetype=treatment),method="lm", se=F, size=1.8, show.legend=F)+
	geom_point(alpha=0.5, size=3, shape=17)+
	labs(y="Total phenolics (prop dw GAE)", x="Proportion growth", title="Young leaves")+
	theme_classic()+
	theme(legend.title = element_blank(), text = element_text(size=20), legend.position = c(0.88,0.88),
		  legend.text.align = 0,legend.spacing.y = unit(0, "mm"),
		  legend.box.background = element_rect(colour = "black"), plot.title = element_text(hjust = 0.5))+
	scale_linetype_manual(values = c("solid", "solid", "solid","solid"))+
	scale_color_brewer(palette=9, type = "div", direction = 1, labels=leg.lab)
gd_young


#Export Fig2
tiff('Maynard_etal_Fig2.tiff', units="in", width=8, height=10, res=300)
ggarrange(gd_mature, gd_young,
		  labels = c("a", "b"),heights = c(2, 2),
		  ncol = 1, nrow = 2,
		  font.label=list(color="black",size=20))
dev.off()


##Q3 GRAPH----
fig3<-all.dat80_m %>%
	ggplot(aes(x=pdw, 
			   y=prop_herb1))+ 
	geom_point(alpha=0.6, size=2.5, aes(color=treatment))+
	geom_smooth(method = 'lm', fill="light grey", linetype="solid", color="#440154FF")+
	theme_classic()+
	labs(y="Proportion herbivory", x="Total phenolics (prop dw GAE)")+
	theme(text = element_text(size=16))+
	scale_color_brewer(palette=9, type = "div", direction = 1, labels=leg.lab)+
	theme(legend.title = element_blank(), text = element_text(size=18), legend.position = c(0.8,0.88),
		  legend.text.align = 0,legend.spacing.y = unit(0, "mm"),
		  legend.box.background = element_rect(colour = "black"))
fig3

#EXPORT PLOT
tiff('Maynard_etal_Fig3.tiff', units="in", width=7, height=5, res=300)
fig3
dev.off()


##OTHER/OLDER STUFF----
#OLD QUESTION 1
#1A. GROWTH----
mod.gro2<-(betareg(prop_gro~treatment, dat=all.dat20))
Anova(mod.gro2)
#df=3, chisq=1.0444, p=0.79
#There is no effect of treatment on plant growth

#1B. TOTAL PHENOLICS----
#both leaf ages
chem.mod2 <- glmmTMB(pdw ~ treatment + (1|chamber) + (1|stage), data = all.dat40, family = "beta_family")
Anova(chem.mod2)
#df=3, chisq=1.73,  p=0.63
#There is no effect of treatment on leaf chemical defense

#young leaves
chem.mod2y <- glmmTMB(pdw ~ treatment + (1|chamber), data = all.dat40_y, family = "beta_family")
Anova(chem.mod2y)
#df=3, chisq=1.95,  p=0.58
#There is no effect of treatment on young leaf chemical defense

#mature leaves
chem.mod2m <- glmmTMB(pdw ~ treatment + (1|chamber), data = all.dat40_m, family = "beta_family")
Anova(chem.mod2m)
#df=3, chisq=1.24,  p=0.74
#There is no effect of treatment on mature leaf chemical defense

#1C. HERBIVORY----
herb.mod2 <- glmmTMB(prop_herb1 ~ treatment + (1|chamber) + (1|stage), data = all.dat80, family = "beta_family")
Anova(herb.mod2)
#df=3, chisq=0.94,  p=0.82
#There is no effect of treatment on leaf herbivory

#young leaves
herb.mod2y <- glmmTMB(prop_herb1 ~ treatment + (1|chamber), data = all.dat80_y, family = "beta_family")
Anova(herb.mod2y)
#df=3, chisq=0.37,  p=0.95
#There is no effect of treatment on young leaf herbivory

#mature leaves
herb.mod2m <- glmmTMB(prop_herb1 ~ treatment + (1|chamber), data = all.dat80_m, family = "beta_family")
Anova(herb.mod2m)
#df=3, chisq=2.08,  p=0.56
#There is no effect of treatment on mature leaf herbivory



##Old stats----
#SUMMARY STATS
sum.tab.m <- ddply(all.dat80_m, c("treatment"), summarise,
				   N    = length(prop_herb),
				   mean = mean(prop_herb),
				   sd   = sd(prop_herb),
				   se   = sd / sqrt(N))
sum.tab.m
#combo herbivory/temp herb - temp herb
(0.10493750/0.02531667)-0.02531667
#mature leaves in combo treatments had 4.12 times more herbivory than leaves in temperature

(0.10493750/0.03472000)-0.03472000
#mature leaves in combo treatment had 2.99 times more herbivory than leaves in control

(0.10493750/0.06226000)-0.06226000
#mature leaves in combo treatment had 1.62 times more herbivory relative to leaves in CO2
#herbivory plot mature
ggplot(data=all.dat80_m, aes(x=treatment, y=prop_herb1))+ 
	geom_point(position=position_jitter(width = 0.025), alpha=0.6, aes(color=treatment), size=2.5)+
	stat_summary(fun.data = "mean_se", colour="black", size=1)+
	theme_classic()+
	theme(legend.position = "none")+
	labs(y="Proportion herbivory", x="")+
	theme(text = element_text(size=16), axis.text.x = element_text(angle=20, hjust=0.9, size=12),
		  legend.position = "none")+
	scale_x_discrete(limits=c("control chamber","CO2", "TC", "TC + CO2"),
					 labels=c("control chamber"=expression(atop("Control")),
					 		 "TC"="Temperature",
					 		 CO2=expression(CO["2"]),
					 		 "TC + CO2"=expression(CO["2"] + Temp)))+
	scale_color_brewer(palette=9, type = "div", direction = 1)

#R^2 data
mod1<-lm(data = all.dat80_m, prop_herb1~pdw*treatment)
summary(mod1)#Multiple R-squared:  0.03
#y=mx+b, y=-1.03027 + 0.10405
plot(all.dat80_m$prop_herb1~all.dat80_m$pdw)
#Mature leaves in experienced 1.03% less herbivory with every 1% increase in total phenolics


#Maure leaves, R^2 data
all.dat40_m <- all.dat40_m[order(all.dat40_m$treatment),]
q2_temp_m<-slice(all.dat40_m, 11:16)
sum.q2<-lm(data = q2_temp_m, pdw~prop_gro)
summary(sum.q2)#Multiple R-squared:  0.7486, p=0.0260
#y=mx+b, y=0.078748x + 0.024343
plot(q2_temp_m$pdw~q2_temp_m$prop_gro)
#Phenolics of Mature leaves in elevated temperatures inc 0.8% with every 1% increase in growth?


#Young leaves, R^2 data
all.dat40_y <- all.dat40_y[order(all.dat40_y$treatment),]
ydat_combo<-slice(all.dat40_y, 17:20)
sum.q2S<-lm(data = ydat_combo, pdw~prop_gro)
summary(sum.q2S)#Multiple R-squared:  0.8526, p=0.0766
#y=mx+b, y=0.053937x + 0.049805
plot(ydat_combo$pdw~ydat_combo$prop_gro)
#Chemical defense of Young leaves in combo treatments inc 0.05% with every 1% increase in growth

#defense~growth*treatment
#mature leaves, interactive mixed model
grow.mod2m_i<-glmmTMB(pdw ~ treatment * prop_gro + (1|chamber), data = all.dat40_m, family = "beta_family")
Anova(grow.mod2m_i)
#treatment, Chisq=0.73, p=.87
#growth, chisq=0.50, p=0.48
#interaction, chisq=8.44, p=0.038
joint_tests((grow.mod2m_i), by = "treatment")
#temperature, F=8.32, p=0.016
#significant effect of growth on defense in  temperature treatment
#while usually there is a negative relationship between growth and defense, 
#for mature leaves in temperature experiment, plants that grew more also had higher defenses in mature leaves
#young leaves, interactive model mixed model
grow.mod2y_i<-glmmTMB(pdw ~ treatment * prop_gro + (1|chamber), data = all.dat40_y, family = "beta_family")
Anova(grow.mod2y_i)
#treatment, Chisq=2.93, p=0.40
#growth, chisq=1.37, p=0.24
#interaction, chisq=7.73, p=0.0519, marginally significant
joint_tests((grow.mod2y_i), by = "treatment")
#combo treatment marginally significant
#F=4.88, p=0.0516

#Extrapolating slopes from the mixed models
coef(grow.mod2m_i)
#Phenolics of Mature leaves in elevated temperatures inc 2.2% with every 1% increase in growth?
coef(grow.mod2y_i)
#Chemical defense of Young leaves in combo treatments inc 0.72% with every 1% increase in growth
#or relative to control??

gd_mature_r<-gd_mature+
	annotate("text", x = 0.47, y = 0.055,
			 label = "paste(italic(R) ^ 2, \" = 0.75\")", parse = TRUE, size =5)
gd_mature_r
gd_young_r<-gd_young+
	annotate("text", x = 0.53, y = 0.083,
			 label = "paste(italic(R) ^ 2, \" = 0.86\")", parse = TRUE, size =5)
gd_young_r