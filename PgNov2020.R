#Climate change & Piper generalense

library(lme4)
library(ggplot2)
library(car)
library(multcomp)
library(plyr)

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
phen <- phen[order(phen$treat),]#check removed correct rows

#combine triplicate readings for concentration
phen_ag<-aggregate(concen~treat+sample+stage+chamber, data=phen, FUN=mean)

#combine triplicate readings for %dw
phen_ag2<-aggregate(pdw~treat+sample+stage+chamber, data=phen, FUN=mean)

phen1<-lmer(pdw ~ treat * stage + (1|chamber), data=phen_ag2, na.action = "na.omit")

phen2<-lmer(pdw ~ treat + stage + (1|chamber), data=phen_ag2, na.action = "na.omit")

phen3<-lmer(pdw ~ treat + (1|chamber), data=phen_ag2, na.action = "na.omit")
phen4<-lmer(pdw ~ stage + (1|chamber), data=phen_ag2, na.action = "na.omit")
phen.null<-lmer(pdw ~ 1 + (1|chamber), data=phen_ag2, na.action = "na.omit")

library(AICcmodavg)
modcomp.phen<-aictab(cand.set=list(phen1, phen2, phen3, phen4, phen.null),
				modnames=c("interaxn","add","treat", "stage", "null"), REML=F)#AIC table
#don't know why it's not accepting REML=F?
modcomp.phen
#best fit model=only stage
summary(phen4)
Anova(phen4)
shapiro.test(resid(phen4))#normal!!
#but silly to exclude treatment bc that's one of main questions of study

modcomp.phen2<-aictab(cand.set = list(phen1, phen2, phen.null),
					 modnames = c("interaxn", "add", "null"), REML=F)
#don't know why it's not accepting REML=F?
modcomp.phen2
#additive=model of better fit

summary(phen2)
Anova(phen2)
shapiro.test(resid(phen2))#normal dist!!

#SUMMARY STATS
library(plyr)
phen.tab <- ddply(phen_ag2, c("stage"), summarise,
					N    = length(pdw),
					mean = mean(pdw),
					sd   = sd(pdw),
					se   = sd / sqrt(N))
phen.tab

#young leaf avg pdw/old leaf avg pdw
6.859805/4.970515
#1.380099
#Young leaves had an average of 1.4 times more total phenolics

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

##couldn't run mixed model with chamber as random effect, so aggregated dataset by chamber to avoid pseudorep
##run LMs and betaregs, but model comparison identified the null model as mod of best fit
##then I sep data, running model for  presence/absense herbivory and  another model for only leaves that had herbivory
##for proportion of leaves that received herbivory, leaf age was a clear predictor 
##and the additive model was mod of best fit.
##for proportion of herbivory on leaves, the null model was always model of better fit

pg <- read.csv(file="Piper_herbivory.csv",head=TRUE)

#creating col for proportion herbivory
pg$percent_herbivory<-as.numeric(pg$percent_herbivory)
pg$prop_herb<-(pg$percent_herbivory/100)
pg$prop_herb<-as.numeric(pg$prop_herb)

pg <- pg[order(pg$treatment),]
pg1<-pg[-c(41:60),]#removing control (no chamber)
#pg1=data without no chamber control group

#combine chambers to avoid pseudorep and take out of models, averaging prop. herb, N=8 
pg2<-aggregate(prop_herb~treatment + age,data=pg1,FUN=mean)

herb.3<-lm(prop_herb ~ treatment + age, data=pg2)
summary(herb.3)#effect of age
shapiro.test(resid(herb.3))#normal!

herb.4<-lm(prop_herb ~ treatment * age, data=pg2)
summary(herb.4)#all NAs, not enough data?

herb.5<-lm(prop_herb ~ 1, data=pg2)

modcomp.herb<-aictab(cand.set=list(herb.3, herb.5),
					 modnames=c("add", "null"), REML=F)#AIC table
modcomp.herb
#null is model of better fit...

library(betareg)
betaherb<-betareg(prop_herb ~ treatment + age, dat=pg2)
summary(betaherb)
shapiro.test(resid(betaherb))#normal! 
Anova(betaherb)#both age and treatment significant...

betaherb2<-betareg(prop_herb ~ treatment * age, dat=pg2)
summary(betaherb2)#still can't do that interaction

betaherb5<-betareg(prop_herb~1, dat=pg2)

modcomp.herb.b<-aictab(cand.set=list(betaherb, betaherb5),
					   modnames=c("add", "null"), REML=F)#AIC table
modcomp.herb.b
#null is model of better fit

#plots
ggplot(pg2, aes(x=treatment, y=prop_herb))+geom_boxplot()+geom_point()
ggplot(pg2, aes(x=age, y=prop_herb))+geom_boxplot()+geom_point()

#separating datasets
#create col for herbivory pres/abs 
pg1$pa_herb <- NA
for(i in 1:length(pg1$percent_herbivory)){
	if(pg1$percent_herbivory[i]==0){pg1$pa_herb[i]=0}
	if(pg1$percent_herbivory[i]>0){pg1$pa_herb[i]=1}
}
hist(pg1$pa_herb)

herb.pa.mod2 <-glm (pa_herb~treatment+age, data=pg1,  family = binomial)
summary(herb.pa.mod2)#age sig
herb.pa.mod3 <-glm (pa_herb~treatment*age, data=pg1,  family = binomial)
herb.pa.mod4 <-glm (pa_herb~1, data=pg1,  family = binomial)

modcomp.herb.pa<-aictab(cand.set=list(herb.pa.mod2, herb.pa.mod3, herb.pa.mod4),
					 modnames=c("add","interaxn", "null"), REML=F)#AIC table
modcomp.herb.pa
#add=model of better fit

summary(herb.pa.mod2)#age sig, treat not
Anova(herb.pa.mod2)#treamtnet p=0.4244, age p < 0.0001
shapiro.test(resid(herb.pa.mod2))#not normal :(

#PLOT
#changing names for plot
pg1$age<-as.character(pg1$age)
pg1$age[pg1$age=="Y"]="Young"
pg1$age[pg1$age=="M"]="Old"

ggplot(pg1, aes(age, pa_herb,fill=age,color))+
	geom_bar(stat = "summary")+
	theme_classic()+
	scale_fill_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.position = "none",
		  text = element_text(size=15))+
	labs(x = "", y = "Proportion of leaves with herbivory present")

herb.pa.tab <- ddply(pg1, c("age"), summarise,
				  N    = length(pa_herb),
				  mean = mean(pa_herb),
				  sd   = sd(pa_herb),
				  se   = sd / sqrt(N))
herb.pa.tab
0.825/0.200

#prop herbivory with leaves that have herbivory
pg1 <- pg1[order(pg1$pa_herb),]
pg.herb.pres<-pg1[-c(1:39),]#removing no herb

herb20<-glm(prop_herb~treatment+age, data=pg.herb.pres)
herb21<-glm(prop_herb~treatment*age, data=pg.herb.pres)
herb22<-glm(prop_herb~1, data=pg.herb.pres)

modcomp.herb.20<-aictab(cand.set=list(herb20, herb21, herb22),
						modnames=c("add","interaxn", "null"), REML=F)#AIC table
modcomp.herb.20
#null=model of better fit

summary(herb20)#nothing sig
shapiro.test(resid(herb20))#not normal

#betaregressions
betaherb1<-betareg(prop_herb ~ treatment * age, dat=pg.herb.pres)
summary(betaherb1)

betaherb2<-betareg(prop_herb ~ treatment + age, dat=pg.herb.pres)
summary(betaherb2)

betaherb.null<-betareg(prop_herb ~ 1, dat=pg.herb.pres)

modcomp.herb.beta<-aictab(cand.set=list(betaherb1, betaherb2, betaherb.null),
						  modnames=c("interaxn","add", "null"), REML=F)#AIC table
modcomp.herb.beta
#null model is best fit...


herb1<-lmer(prop_herb ~ treatment * age + (1|chamber), data=pg1, na.action = "na.omit")
summary(herb1)#warnings, prob overfitting
Anova(herb1)#age signif

shapiro.test(resid(herb1))#residuals not normally distributed
hist(pg$prop_herb)#skewed, zero-inflated

herb.1<-lm(prop_herb ~ treatment + age + chamber, data=pg1)
summary(herb.1)#effect of age and chamber(H & I)
shapiro.test(resid(herb.1))#not normal


herb.2<-lm(prop_herb ~ treatment + age, data=pg1)
summary(herb.2)#effect of age and marginal effect of temp+CO2 treatment
shapiro.test(resid(herb.2))#not normal

#combine chambers to avoid pseudorep and take out of models, averaging prop. herb, N=8 
pg2<-aggregate(prop_herb~treatment + age,data=pg1,FUN=mean)

herb.3<-lm(prop_herb ~ treatment + age, data=pg2)
summary(herb.3)#effect of age, p=0.0394, t=-3.502
shapiro.test(resid(herb.3))#normal!

herb.4<-lm(prop_herb ~ treatment * age, data=pg2)
summary(herb.4)#all NAs, not enough data?



herb.pa.mod <- glmer(pa_herb ~ treatment + age + (1|chamber), data=pg1, family=binomial, na.action="na.fail")
#singular, overfitting
?isSingular
summary(herb.pa.mod)
Anova(herb.pa.mod)
shapiro.test(resid(herb.pa.mod))#not normal

herb6<-lmer(prop_herb ~ treatment * age + (1|chamber), data=pg.herb.pres, na.action = "na.omit")
summary(herb6)#warnings
Anova(herb6)#no signif, interaction marginally

#changing names for plots
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

herb_plant<-aggregate(pa_herb~treatment + age ,data=pg1,FUN=mean)
herb_plant$pa_herb<-as.numeric(herb_plant$pa_herb)
herb_plant$pa_herb[herb_plant$pa_herb==1] <- 0.99999
herb_plant$pa_herb[herb_plant$pa_herb==0] <- 0.00001
b1<-betareg(pa_herb ~ treatment + age, dat=herb_plant)
summary(b1)

m1 <- glmer(pa_herb ~ treatment * age + (1|chamber), data=herb_plant, family=binomial, na.action="na.fail")
summary(m1)
Anova(m1)

ggplot(herb_plant, aes(age, pa_herb))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.30, aes(color=age))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.position = "none",
		  text = element_text(size=15))+
	labs(x = "", y = "% herbivory")#+
	#scale_y_continuous(limits =  c(0,1))

ggplot(herb_plant, aes(treatment, pa_herb))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.25, aes(color=treatment))+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "% herbivory")

