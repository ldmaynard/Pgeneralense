#Climate change & Piper generalense

library(lme4)
library(ggplot2)
library(car)
library(multcomp)
library(plyr)
library(AICcmodavg)
library(betareg)

#GROWTH----
grow <- read.csv(file="Piper_growth.csv",head=TRUE)

grow$rel_gro<-((grow$ht.2018.09.cm-grow$ht.2018.04.cm)/grow$ht.2018.04.cm)

grow <- grow[order(grow$Treatment),]
grow<-grow[-c(11:15),]#removing control (no chamber)

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

#Script that can be used to quantify compounds as concentration or absolute value
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
ggplot(phen, aes(x=treat, y=concen))+geom_boxplot()+geom_point()

#average triplicate readings for concentration
phen_ag<-aggregate(concen~treat+sample+stage+chamber, data=phen, FUN=mean)

#combine triplicate readings for %dw
phen_ag2<-aggregate(pdw~treat+sample+stage+chamber, data=phen, FUN=mean)

#models for pdw
phen1<-lmer(pdw ~ treat * stage + (1|chamber), data=phen_ag2, na.action = "na.omit")
phen2<-lmer(pdw ~ treat + stage + (1|chamber), data=phen_ag2, na.action = "na.omit")
phen3<-lmer(pdw ~ treat + (1|chamber), data=phen_ag2, na.action = "na.omit")
phen4<-lmer(pdw ~ stage + (1|chamber), data=phen_ag2, na.action = "na.omit")
phen.null<-lmer(pdw ~ 1 + (1|chamber), data=phen_ag2, na.action = "na.omit")

modcomp.phen<-aictab(cand.set=list(phen1, phen2, phen3, phen4, phen.null),
				modnames=c("interaxn","add","treat", "stage", "null"), REML=F)#AIC table
#don't know why it's not accepting REML=F?
modcomp.phen #best fit model=only stage, next best is additive, dAIC=4.04

write.table(modcomp.phen, file = "aic_totalphenolics.csv", sep = ",", quote = FALSE, row.names = F)

summary(phen4)
shapiro.test(resid(phen4))#normal!!
Anova(phen4)#stage p=3.9 e-15
Anova(phen3)#treatment=0.72
#don't need to do mod avg bc next model has deltaAIC>4 for add and treat univar mod was last in mod comp

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

ggplot(phen_ag2, aes(stage, pdw))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.40)+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=15))+
	labs(x = "", y = "Total phenolics (%dw in gallic acid equivalents)")

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

##ran LMM with chamber as random effect, model comp showed additive was model of best fit with null behind and dAIC>3
##but residuals for additive model don't have normal distribution
##so I aggregated the dataset (pg2) by chamber (to avoid pseudorep) and ran LMs and betaregs
##null model was always the clear top model in LM. Can't run interactive model
##additive was top model in beta reg, null close behind. Again, can't run interactive model
##Good news: residuals of age betareg are normal. So.
##Next, I sep data, running model for  presence/absense herbivory (binomial GLM, not enough data to fit GLMM, so not controling for chamber)
##and another set of models for only leaves that had herbivory (LMMs)
##for proportion of leaves that received herbivory, leaf age was mod of best fit again (dAIC>3)
##for proportion of herbivory on leaves, the null model was model of best fit (dAIC>8)
##The best options seem to be the first one (LMM) and how I did the phenolics, but still have non-normal dist residuals for the top model
##Or the betareg with the aggregated data (norm dist resids)
##Either options gives the same outcome. Top models=age univariate
##Hurdle method isn't working with these data, so I could transform the data, logit transformation?

pg <- read.csv(file="Piper_herbivory.csv",head=TRUE)

#creating col for proportion herbivory
pg$percent_herbivory<-as.numeric(pg$percent_herbivory)
pg$prop_herb<-(pg$percent_herbivory/100)
pg$prop_herb<-as.numeric(pg$prop_herb)

pg <- pg[order(pg$treatment),]
pg1<-pg[-c(41:60),]#removing control (no chamber)
#pg1=data without no chamber control group

#logit transformation
logitTransform <- function(p) { log(p/(1-p)) }
pg1$prop_herb_add <- (pg1$prop_herb)+0.001 #adding small amount to get rid of zeros 
pg1$pLogit <- logitTransform(pg1$prop_herb_add)
pg1$pLogit <- as.numeric(pg1$pLogit)

#Arcsine transformation
asinTransform <- function(p) { asin(sqrt(p)) }
pg1$pAsin <- asinTransform(pg1$prop_herb_add)
pg1$pAsin <- as.numeric(pg1$pAsin)

#LMM with chamber as random effect
lmm1<-lmer(prop_herb ~ age + (1|chamber), data=pg1, na.action = "na.omit")
lmm2<-lmer(prop_herb ~ treatment + (1|chamber), data=pg1, na.action = "na.omit")
lmm3<-lmer(prop_herb ~ age + treatment + (1|chamber), data=pg1, na.action = "na.omit")
lmm4<-lmer(prop_herb ~ age * treatment + (1|chamber), data=pg1, na.action = "na.omit")#warning
lmm5<-lmer(prop_herb ~ (1|chamber), data=pg1, na.action = "na.omit")

modcomp.lmm<-aictab(cand.set=list(lmm1, lmm2, lmm3, lmm4, lmm5),
				 modnames=c("age", "treat", "add", "interactive", "null"), REML=F)#AIC table
modcomp.lmm#error about fixed effects bening different?
#best model is age, next model = null and dAIC>3

summary(lmm1)
shapiro.test(resid(lmm1))#not normal 
hist(resid(lmm1))
qqnorm(resid(lmm1))
qqline(resid(lmm1))
#transform data? or try hurdle



#summary stats
herb.sum <- ddply(pg1, c("age"), summarise,
					 N    = length(prop_herb),
					 mean = mean(prop_herb),
					 sd   = sd(prop_herb),
					 se   = sd / sqrt(N))
herb.sum
0.0528275/0.0103100
#mature leaves had 5.1 times more herbivory than younger leaves

#PLOT
ggplot(pg1, aes(age, prop_herb))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.30)+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=15))+
	labs(x = "", y = "Proportion herbivory")

#plots
ggplot(pg1, aes(x=treatment, y=prop_herb))+geom_boxplot()+geom_point()
ggplot(pg1, aes(x=age, y=prop_herb))+geom_boxplot()+geom_point()

#LMM with chamber as random effect with logit transformed data
lmm1t<-lmer(pLogit ~ age + (1|chamber), data=pg1, na.action = "na.omit")
lmm2t<-lmer(pLogit ~ treatment + (1|chamber), data=pg1, na.action = "na.omit")
lmm3t<-lmer(pLogit ~ age + treatment + (1|chamber), data=pg1, na.action = "na.omit")
lmm4t<-lmer(pLogit ~ age * treatment + (1|chamber), data=pg1, na.action = "na.omit")#warning
lmm5t<-lmer(pLogit ~ (1|chamber), data=pg1, na.action = "na.omit")

modcomp.lmmt<-aictab(cand.set=list(lmm1t, lmm2t, lmm3t, lmm4t, lmm5t),
					modnames=c("age", "treat", "add", "interactive", "null"), REML=F)#AIC table
modcomp.lmmt#error about fixed effects bening different?
#best model is age, next model = null and dAIC>3

summary(lmm1t)
shapiro.test(resid(lmm1t))#not normal, but more normal than untransformed data 

#LMM with chamber as random effect with arcsine transformed data
lmm1a<-lmer(pAsin ~ age + (1|chamber), data=pg1, na.action = "na.omit")
lmm2a<-lmer(pAsin ~ treatment + (1|chamber), data=pg1, na.action = "na.omit")
lmm3a<-lmer(pAsin ~ age + treatment + (1|chamber), data=pg1, na.action = "na.omit")
lmm4a<-lmer(pAsin ~ age * treatment + (1|chamber), data=pg1, na.action = "na.omit")#warning
lmm5a<-lmer(pAsin ~ (1|chamber), data=pg1, na.action = "na.omit")

modcomp.lmma<-aictab(cand.set=list(lmm1a, lmm2a, lmm3a, lmm4a, lmm5a),
					 modnames=c("age", "treat", "add", "interactive", "null"), REML=F)#AIC table
modcomp.lmma#error about fixed effects bening different?
#best model is age, next model = null and dAIC>3

summary(lmm1a)
shapiro.test(resid(lmm1a))#not normal

#create dataset that combines chambers to avoid pseudorep and remove random effect, averaging prop. herb, N=8 
pg2<-aggregate(prop_herb~treatment + age,data=pg1,FUN=mean)

herb.1<-lm(prop_herb ~ treatment, data=pg2)
herb.2<-lm(prop_herb ~ age, data=pg2)
herb.5<-lm(prop_herb ~ 1, data=pg2)
herb.3<-lm(prop_herb ~ treatment + age, data=pg2)
herb.4<-lm(prop_herb ~ treatment * age, data=pg2)
summary(herb.4)#all NAs, not enough data?, breaks AIC mod comp

modcomp.herb<-aictab(cand.set=list(herb.1, herb.2, herb.3, herb.5),
					 modnames=c("treatment", "age", "add", "null"), REML=F)#AIC table
modcomp.herb
#null is model of better fit, age close behind, dAICc=-0.4

#plots
ggplot(pg2, aes(x=treatment, y=prop_herb))+geom_boxplot()+geom_point()
ggplot(pg2, aes(x=age, y=prop_herb))+geom_boxplot()+geom_point()

#hurdle
#separating datasets
#create col for herbivory pres/abs 
pg1$pa_herb <- NA
for(i in 1:length(pg1$percent_herbivory)){
	if(pg1$percent_herbivory[i]==0){pg1$pa_herb[i]=0}
	if(pg1$percent_herbivory[i]>0){pg1$pa_herb[i]=1}
}
hist(pg1$pa_herb)

herb.pa.mod <-glmer (pa_herb ~ treatment + (1|chamber), data=pg1,  family = binomial)
#singular fit, not enough data

herb.pa.mod1 <-glm (pa_herb~treatment, data=pg1,  family = binomial)
herb.pa.mod2 <-glm (pa_herb~age, data=pg1,  family = binomial)
herb.pa.mod3 <-glm (pa_herb~treatment+age, data=pg1,  family = binomial)
herb.pa.mod4 <-glm (pa_herb~treatment*age, data=pg1,  family = binomial)
herb.pa.mod5 <-glm (pa_herb~1, data=pg1,  family = binomial)

modcomp.herb.pa<-aictab(cand.set=list(herb.pa.mod1, herb.pa.mod2, herb.pa.mod3, herb.pa.mod4, herb.pa.mod5),
					 modnames=c("treat", "age","add","interaxn", "null"), REML=F)#AIC table
modcomp.herb.pa
#age=model of best fit, add dAIC>3

summary(herb.pa.mod2)
shapiro.test(resid(herb.pa.mod2))#not normal :(
hist(resid(herb.pa.mod2))
qqnorm(resid(herb.pa.mod2))
qqline(resid(herb.pa.mod2))


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
#old leaves were 4.1 times more likely to have herbivore damage compared to young leaves

#prop herbivory with leaves that have herbivory
pg1 <- pg1[order(pg1$pa_herb),]
pg.herb.pres<-pg1[-c(1:39),]#removing no herb

lmm20<-lmer(prop_herb ~ age + (1|chamber), data=pg.herb.pres, na.action = "na.omit")
lmm21<-lmer(prop_herb ~ treatment + (1|chamber), data=pg.herb.pres, na.action = "na.omit")
lmm22<-lmer(prop_herb ~ age + treatment + (1|chamber), data=pg.herb.pres, na.action = "na.omit")
lmm23<-lmer(prop_herb ~ age * treatment + (1|chamber), data=pg.herb.pres, na.action = "na.omit")
lmm24<-lmer(prop_herb ~ (1|chamber), data=pg.herb.pres, na.action = "na.omit")

modcomp.lmm2<-aictab(cand.set=list(lmm20, lmm21, lmm22, lmm23, lmm24),
					modnames=c("age", "treat", "add", "interactive", "null"), REML=F)#AIC table
modcomp.lmm2
#best model is null, next model = age and dAIC>8

pg.herb.pres2<-aggregate(prop_herb~treatment + age,data=pg.herb.pres,FUN=mean)
lm25<-lm(prop_herb ~ age, data=pg.herb.pres2)
lm26<-lm(prop_herb ~ treatment, data=pg.herb.pres2)
lm27<-lm(prop_herb ~ age + treatment, data=pg.herb.pres2)
lm28<-lm(prop_herb ~ age * treatment, data=pg.herb.pres2)
lm29<-lm(prop_herb ~ 1, data=pg.herb.pres2)

summary(lm28)#NAs again, removing from AIC mod comp

modcomp.lmm3<-aictab(cand.set=list(lm25, lm26, lm27, lm29),
					 modnames=c("age", "treat", "add", "null"), REML=F)#AIC table
modcomp.lmm3
#again null is best model, age is next but dAIC>5


#plot
ggplot(pg.herb.pres, aes(treatment, prop_herb))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.30, aes(color=treatment))+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=15))+
	labs(x = "", y = "Proportion herbivory")

#OTHER THINGS I TRIED####
#betaregressions
#using data with herb present and combined chambers (bc can't add random effect in beta regs)
#if you run it with full dataset, null is still best model, but add model has lower dAIC
library(betareg)
betaherb1<-betareg(prop_herb ~ treatment, dat=pg.herb.pres2)
betaherb2<-betareg(prop_herb ~ age, dat=pg.herb.pres2)
betaherb3<-betareg(prop_herb ~ treatment + age, dat=pg.herb.pres2)
betaherb4<-betareg(prop_herb ~ treatment * age, dat=pg.herb.pres2)#doesn't like interactive model
betaherb.null<-betareg(prop_herb ~ 1, dat=pg.herb.pres2)

modcomp.herb.beta<-aictab(cand.set=list(betaherb1, betaherb2, betaherb3, betaherb.null),
						  modnames=c("treat", "age", "add", "null"), REML=F)#AIC table
modcomp.herb.beta
#null model is best fit, age dAIC>5

betaherb10<-betareg(prop_herb ~ treatment, dat=pg1)#can't run with full dataset bc zero-inflated. also doesn't account for chamber

#running with aggregated data
betaherb11<-betareg(prop_herb ~ treatment, dat=pg2)
betaherb12<-betareg(prop_herb ~ age, dat=pg2)
betaherb13<-betareg(prop_herb ~ treatment + age, dat=pg2)
betaherb14<-betareg(prop_herb ~ treatment * age, dat=pg2)#doesn't like interactive model
betaherb.null1<-betareg(prop_herb ~ 1, dat=pg2)

modcomp.herb.beta1<-aictab(cand.set=list(betaherb11, betaherb12, betaherb13, betaherb.null1),
						  modnames=c("treat", "age", "add", "null"), REML=F)#AIC table
modcomp.herb.beta1
#age is top model, null close behind dAIC=1.23

shapiro.test(resid(betaherb12))#normal!


#dataset not accounting for chamber, treatment only model is best model
herb20<-glm(prop_herb~treatment, data=pg.herb.pres)
herb21<-glm(prop_herb~age, data=pg.herb.pres)
herb22<-glm(prop_herb~treatment+age, data=pg.herb.pres)
herb23<-glm(prop_herb~treatment*age, data=pg.herb.pres)
herb24<-glm(prop_herb~1, data=pg.herb.pres)

modcomp.herb.20<-aictab(cand.set=list(herb20, herb21, herb22, herb23, herb24),
						modnames=c("treat", "age","add","interaxn", "null"), REML=F)#AIC table
modcomp.herb.20
#treat=model of better fit, null close behind dAIC=2.11

summary(herb20)#nothing sig, temp+CO2 p=0.09
Anova(herb20)#but when you run the ANOVA, treatment is signficant
summary(glht(herb20, linfct=mcp(treatment="Tukey")))
#tukey test shoes sig diff bw temp + CO2 chamber and control chamber, and diff bw T+C chamber + temp chamber

shapiro.test(resid(herb20))#not normal, p=0.02
hist(resid(herb20))
qqnorm(resid(herb20))
qqline(resid(herb20))


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

