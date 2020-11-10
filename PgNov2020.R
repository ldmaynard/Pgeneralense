#Climate change & Piper generalense

library(lme4)
library(ggplot2)
library(car)

#GROWTH----
grow <- read.csv(file="Piper_growth.csv",head=TRUE)

grow$growth<-grow$ht.2018.09.cm-grow$ht.2018.04.cm
grow$rel_gro<-((grow$ht.2018.09.cm-grow$ht.2018.04.cm)/grow$ht.2018.04.cm)

grow <- grow[order(grow$Treatment),]
grow<-grow[-c(11:15),]#removing control (no chamber)

library(lme4)
gro1<-lmer(rel_gro ~ Treatment + (1|Casa), data=grow, na.action = "na.omit")
#Error: number of levels of each grouping factor must be < number of observations
#not enough data to run lmer

gro2<-aov(grow$rel_gro~grow$Treatment)
summary.aov(gro2)
#p=0.927

#plot
ggplot(grow, aes(Treatment, rel_gro))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	theme(legend.title = element_blank(),
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "Relative growth (cm)")

#CHEMISTRY----

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

#Load total phenolics data
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
phen$stage[phen$stage=="m"]="Old"

#change names for treatments
phen$treat<-as.character(phen$treat)
phen$treat[phen$treat=="control_chamber"]="Control (chamber)"
phen$treat[phen$treat=="control_nat"]="Control (no chamber)"
phen$treat[phen$treat=="TC"]="Temperature"
phen$treat[phen$treat=="TC+CO2"]="Temp + CO2"


ggplot(phen, aes(x=treat, y=concen))+geom_boxplot()+geom_point()

phen<-phen[-c(1:30),]#removing control (no chamber)
phen <- phen[order(phen$treat),]#check

#combine triplicate readings for concentration
phen_ag<-aggregate(concen~treat+sample+stage+chamber, data=phen, FUN=mean)

#combine triplicate readings for %dw
phen_ag2<-aggregate(pdw~treat+sample+stage+chamber, data=phen, FUN=mean)

phen1<-lmer(pdw ~ treat * stage + (1|chamber), data=phen_ag2, na.action = "na.omit")
summary(phen1)
Anova(phen1)

phen2<-lmer(pdw ~ treat + stage + (1|chamber), data=phen_ag2, na.action = "na.omit")
summary(phen2)
Anova(phen2)
shapiro.test(resid(phen2))

phen3<-lmer(pdw ~ treat + (1|chamber), data=phen_ag2, na.action = "na.omit")
phen4<-lmer(pdw ~ stage + (1|chamber), data=phen_ag2, na.action = "na.omit")
phen.null<-lmer(pdw ~ 1 + (1|chamber), data=phen_ag2, na.action = "na.omit")

library(AICcmodavg)
modcomp.phen<-aictab(cand.set=list(phen1, phen2, phen3, phen4, phen.null),
				modnames=c("interaxn","add","treat", "stage", "null"), REML=F)#AIC table
#don't know why it's not accepting REML=F?
modcomp.phen
#best fit model=only stage

##PLOTS
#Stage 
ggplot(phen_ag2, aes(stage, pdw))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.40, aes(color=stage))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.position = "none",
		  text = element_text(size=15))+
	labs(x = "", y = "%dw in gallic acid equivalents")

#interaction
ggplot(phen_ag2, aes(treat, pdw, color=stage))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.position = "none",
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "%dw in gallic acid equivalents")#+
	#facet_wrap(~stage)

#treatment
ggplot(phen_ag2, aes(treat, pdw))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.position = "none",
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "%dw in gallic acid equivalents")

#HERBIVORY----
pg <- read.csv(file="Piper_herbivory.csv",head=TRUE)

#creating col for proportion herbivory
pg$percent_herbivory<-as.numeric(pg$percent_herbivory)
pg$prop_herb<-(pg$percent_herbivory/100)
pg$prop_herb<-as.numeric(pg$prop_herb)

herb1<-lmer(prop_herb ~ treatment * age + (1|chamber), data=pg, na.action = "na.omit")
summary(herb1)#warning messages
Anova(herb1)#age signif
shapiro.test(resid(herb1))#residuals not normally distributed
hist(pg$prop_herb)#skewed, zero-inflated

ggplot(pg, aes(x=age, y=percent_herbivory))+geom_boxplot()+geom_point()

pg <- pg[order(pg$treatment),]
pg1<-pg[-c(41:60),]#removing control (no chamber)

herb2<-lmer(prop_herb ~ treatment * age + (1|chamber), data=pg1, na.action = "na.omit")
summary(herb2)#warnings
Anova(herb2)#age signif

shapiro.test(resid(herb2))#residuals not normally distributed
hist(pg$prop_herb)#skewed, zero-inflated

ggplot(pg1, aes(x=treatment, y=percent_herbivory))+geom_boxplot()+geom_point()
ggplot(pg1, aes(x=age, y=percent_herbivory))+geom_boxplot()+geom_point()

#create col for herbivory pres/abs
pg1$pa_herb <- NA
for(i in 1:length(pg1$percent_herbivory)){
	if(pg1$percent_herbivory[i]==0){pg1$pa_herb[i]=0}
	if(pg1$percent_herbivory[i]>0){pg1$pa_herb[i]=1}
}
hist(pg1$pa_herb)

herb.pa.mod <- glmer(pa_herb ~ treatment * age + (1|chamber), data=pg1, family=binomial, na.action="na.fail")
#singular, overfitting
?isSingular
summary(herb.pa.mod)
Anova(herb.pa.mod)
shapiro.test(resid(herb.pa.mod))#not normal

herb.pa.mod2 <- glmer(pa_herb ~ treatment + age + (1|chamber), data=pg1, family=binomial, na.action="na.fail")
#still singular, overfitting
summary(herb.pa.mod2)
Anova(herb.pa.mod2)
shapiro.test(resid(herb.pa.mod2))#not normal

#all data, including non herb
herb5<-lmer(prop_herb ~ treatment * age + (1|chamber), data=pg1, na.action = "na.omit")
summary(herb5)
Anova(herb5)#age signif
shapiro.test(resid(herb5))#not normal

#zero skewed wont run
library(betareg)
betaherb<-betareg(prop_herb ~ treatment + age, dat=pg1)
summary(betaherb)

pg1 <- pg1[order(pg1$pa_herb),]
pg.herb.pres<-pg1[-c(1:39),]#removing no herb

herb6<-lmer(prop_herb ~ treatment * age + (1|chamber), data=pg.herb.pres, na.action = "na.omit")
summary(herb6)
Anova(herb6)#no signif, interaction marginally
shapiro.test(resid(herb6))#normal

betaherb1<-betareg(prop_herb ~ treatment * age, dat=pg.herb.pres)
summary(betaherb1)

betaherb2<-betareg(prop_herb ~ treatment + age, dat=pg.herb.pres)
summary(betaherb2)

betaherb3<-betareg(prop_herb ~ treatment, dat=pg.herb.pres)
betaherb4<-betareg(prop_herb ~ age, dat=pg.herb.pres)
betaherb.null<-betareg(prop_herb ~ 1, dat=pg.herb.pres)

modcomp.herb<-aictab(cand.set=list(betaherb1, betaherb2, betaherb3, betaherb4, betaherb.null),
					 modnames=c("interaxn","add","treat", "stage", "null"), REML=F)#AIC table
modcomp.herb

pg.herb.pres$yhat<-predict(betaherb3)
predplot<-ggplot(pg.herb.pres)+
	geom_point(aes(x=treatment, y=prop_herb))+
	geom_point(aes(x=treatment, y=yhat), color="red", size=2)+
	geom_line(aes(x=treatment, y=yhat) ,color="red", size=1)
predplot

pg1$prop_herb2<-(pg1$prop_herb)+0.1

betaherb10<-betareg(prop_herb2 ~ treatment * age, dat=pg1)
summary(betaherb10)
betaherb11<-betareg(prop_herb2 ~ treatment + age, dat=pg1)
betaherb12<-betareg(prop_herb2 ~ treatment, dat=pg1)
betaherb13<-betareg(prop_herb2 ~ age, dat=pg1)
betaherb14<-betareg(prop_herb2 ~ 1, dat=pg1)

modcomp.herb2<-aictab(cand.set=list(betaherb10, betaherb11, betaherb12, betaherb13, betaherb14),
					 modnames=c("interaxn","add","treat", "stage", "null"), REML=F)#AIC table
modcomp.herb2

summary(betaherb11)
shapiro.test(resid(betaherb11))#not normal though

#changing names for plot
pg1$age<-as.character(pg1$age)
pg1$age[pg1$age=="Y"]="Young"
pg1$age[pg1$age=="M"]="Old"

pg1$treatment<-as.character(pg1$treatment)
pg1$treatment[pg1$treatment=="control chamber"]="Control chamber"
pg1$treatment[pg1$treatment=="T°C"]="Temperature"
pg1$treatment[pg1$treatment=="T°C + CO2"]="Temp + CO2"

perherb_all<-ggplot(pg1, aes(age, percent_herbivory))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.30, aes(color=age))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.position = "none",
		  text = element_text(size=15))+
	labs(x = "", y = "% herbivory")
perherb_all

perherb_treat<-ggplot(pg1, aes(treatment, percent_herbivory))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.25, aes(color=treatment))+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "% herbivory")+
	scale_x_discrete(limits=c("Control chamber", "CO2", "Temperature", "Temp + CO2"))+
	scale_color_manual(values = c("#6baed6", "#969696", "#810f7c", "#fb6a4a"))
perherb_treat


