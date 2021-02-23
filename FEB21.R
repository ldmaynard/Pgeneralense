#Interactive effects of climate change, leaf age, and secondary metabolites 
#on plant growth, defense, and herbivory.

library(lme4)
library(ggplot2)
library(car)
library(multcomp)
library(plyr)
library(AICcmodavg)
library(betareg)

#GROWTH ANALYSIS----
grow <- read.csv(file="Piper_growth.csv",head=TRUE)
colnames(grow)[1] <- "Casa"

grow$rel_gro<-((grow$ht.2018.09.cm-grow$ht.2018.04.cm)/grow$ht.2018.04.cm)
grow$total_gro <- grow$ht.2018.09.cm-grow$ht.2018.04.cm

grow <- grow[order(grow$Treatment),]
grow<-grow[-c(11:15),]#removing control (no chamber)

hist(grow$rel_gro)  #SW: not normal, this is a typical proportion
shapiro.test(grow$rel_gro)#LM: it is technically normal
hist(grow$total_gro)
shapiro.test(grow$total_gro)#LM: but so is this one. Both yield the same results, no treatment effect

#re-ordering factor levels so lm will compare everything to control
grow$Treatment <- factor(grow$Treatment, levels=c("control chamber", "CO2", "T°C", "T°C + CO2" ))

gro3 <- lm(rel_gro ~ Treatment, data=grow)
summary(gro3)
shapiro.test(resid(gro3))#normally dist resids

gro4 <- lm(total_gro ~ Treatment, data=grow)
summary(gro4)
shapiro.test(resid(gro4))#also normal

#Growth plot----
#relative growth
ggplot(grow, aes(Treatment, rel_gro))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	theme(legend.title = element_blank(),
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "Relative growth (cm)")

#total growth
ggplot(grow, aes(Treatment, total_gro))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	theme(legend.title = element_blank(),
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "Total growth (cm)")

#PHENOLICS ANALYSIS----

#Script that can be used to quantify compounds as concentration or absolute value
##Load standard curve data
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

#models for pdw
phen1<-lmer(pdw ~ treat * stage + (1|chamber), data=phen_ag2, na.action = "na.fail")
phen2<-lmer(pdw ~ treat + stage + (1|chamber), data=phen_ag2, na.action = "na.fail")
phen3<-lmer(pdw ~ treat + (1|chamber), data=phen_ag2, na.action = "na.fail")
phen4<-lmer(pdw ~ stage + (1|chamber), data=phen_ag2, na.action = "na.fail")
phen.null<-lmer(pdw ~ 1 + (1|chamber), data=phen_ag2, na.action = "na.fail")

modcomp.phen<-aictab(cand.set=list(phen1, phen2, phen3, phen4, phen.null),
					 modnames=c("interaxn","add","treat", "stage", "null"), REML=F)#AIC table
modcomp.phen #best fit model=only stage, next best is additive, dAIC=4.04

d.phen<-dredge(phen1)
d.phen
d.phen.avg<-model.avg(d.phen, subset=delta<4)#only consists of one model
summary(d.phen.avg)#thus, can't model average

d.phen.avg1<-model.avg(d.phen)
summary(d.phen.avg1)#stage is only significant p val

##Phenolics plots----
#leaf age
ggplot(phen_ag2, aes(stage, pdw))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.40)+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=15))+
	labs(x = "Leaf age", y = "Total phenolics (%dw in gallic acid equivalents)")

#treatment
ggplot(phen_ag2, aes(treat, pdw))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.position = "none",
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "Total phenolics (%dw in gallic acid equivalents)")



#HERBIVORY ANALYSIS----
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

library(MuMIn)
ph6<-lmer(percent_herbivory ~ stage * pdw * treat + (1|chamber), data=ph, na.action = "na.fail")
d1<-dredge(ph6)
d1
# Model average models with delta AICc < 4
d1.avg<-model.avg(d1, subset=delta<4)
summary(d1.avg)
#sig interaction b/w phenolics and temp+Co2 treatment

#all models
d2.avg<-model.avg(d1)
summary(d2.avg)
#moderately sig interaxn bw chem and T+CO2 treatment
#mod sig interaxn among chem, young leaves, and T+CO2 treatment

#HERBIVORY PLOT----
lab1 <- c(expression(CO["2"]),
		  "Control chamber", 
		  expression(Temp + CO["2"]),
		  "Temperature")

ph$Treatment<-ph$treat

herbplot<-ggplot(ph)+
	geom_smooth(aes(pdw, percent_herbivory, color=Treatment, linetype=Treatment),method = "lm", se=F, formula = "y~x")+
	geom_jitter(aes(pdw, percent_herbivory, color=Treatment),position=position_jitter(width = 0.04), alpha=0.30)+
	theme_classic()+
	theme(legend.position = c(0.85, 0.8),
		  legend.direction = "vertical",
		  legend.text = element_text(size = 10, hjust = 0),
		  legend.title = element_text(size = 10),
		  text = element_text(size=15))+
	labs(x = "%dw total phenolics", y = "% herbivory")+
	scale_linetype_manual(values=c("twodash", "twodash", "solid", "twodash"), labels = lab1)+
	scale_color_viridis(discrete = T, option = "D", labels = lab1)
herbplot

herbplot+annotate("text", x = 3.5, y = 14,
			 label = "paste(italic(R) ^ 2, \" = 0.42\")", parse = TRUE, size = 4)
herbplot+facet_grid(.~stage)

ph <- ph[order(ph$treat),]
inter <- ph[21:28,]
ph.inter.mod<-lm(inter$percent_herbivory~inter$pdw)
summary(ph.inter.mod)
plot(inter$prop_herb~inter$pdw)
#y=mx+b, y=-4.708x+34.151, R^2=0.419
#When exposed to increased temperature and CO2,
#Herbivory decreases 1% with every 4.7% increase in total phenolics in increased CO2 and temp conditions

inter <- inter[order(inter$stage),]
inter.mature <- inter[1:4,]
inter.yg <- inter[5:8,]
ph.inter.mod.m<-lm(inter.mature$percent_herbivory~inter.mature$pdw)
summary(ph.inter.mod.m)
plot(inter.mature$percent_herbivory~inter.mature$pdw)
#In mature leaves, herbivory decerases 7.25% with every  1% increase in total phenolics in inc CO2 and temp environment
#R2=0.2827

ph.inter.mod.y<-lm(inter.yg$percent_herbivory~inter.yg$pdw)
summary(ph.inter.mod.y)
plot(inter.mature$percent_herbivory~inter.mature$pdw)
#In young leaves, herbivory decerases 3.34% with every  1% increase in total phenolics in inc CO2 and temp environment
#R2=0.439
