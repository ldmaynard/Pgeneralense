pg <- read.csv(file="Piper_herbivory.csv",head=TRUE)

library(ggplot2)
library(viridis)
library(ggsignif)

#creating col for proportion herbivory
pg$percent_herbivory<-as.numeric(pg$percent_herbivory)
pg$prop_herb<-(pg$percent_herbivory/100)
pg$prop_herb<-as.numeric(pg$prop_herb)

library(lme4)
herb1<-lmer(prop_herb ~ treatment * age + (1|chamber), data=pg, na.action = "na.omit")
summary(herb1)
library(car)
Anova(herb1)

ggplot(pg, aes(x=treatment, y=percent_herbivory, fill=age))+geom_boxplot()+geom_point()

shapiro.test(resid(herb1))#residuals not normally distributed
hist(pg$prop_herb)#skewed, zero-inflated

pg <- pg[order(pg$treatment),]
pg1<-pg[-c(41:60),]#removing control (no chamber)

herb2<-lmer(prop_herb ~ treatment * age + (1|chamber), data=pg1, na.action = "na.omit")
summary(herb2)
Anova(herb2)

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

#nonparametric test of correlation
tbl = table(pgag$percent_herbivory, phen_ag$pdw) 
tbl
chisq.test(tbl)
#not working

pg1 <- pg1[order(pg1$pa_herb),]
pg.herb.pres<-pg1[-c(1:39),]#removing no herb

herb5<-lmer(prop_herb ~ treatment * age + (1|chamber), data=pg.herb.pres, na.action = "na.omit")
summary(herb5)
Anova(herb5)

shapiro.test(resid(herb5))#normal dist
hist(pg.herb.pres$prop_herb)#skewed, zero-inflated

library(betareg)
betaherb<-betareg(prop_herb ~ treatment + age, dat=pg.herb.pres)
summary(betaherb)

ggplot(pg.herb.pres, aes(x=treatment, y=percent_herbivory))+geom_boxplot()+geom_point()
ggplot(pg.herb.pres, aes(x=age, y=percent_herbivory))+geom_boxplot()+geom_point()
ggplot(pg.herb.pres, aes(x=treatment, y=percent_herbivory, class=age))+geom_boxplot()+geom_point()

#summarizing data by age and not by individual leaf, n=100->n=50
pgag<-aggregate(prop_herb~chamber+treatment+age, data=pg, FUN=mean)
ggplot(pgag, aes(x=treatment, y=prop_herb, fill=age))+geom_boxplot()+geom_point()

pgag2<-aggregate(percent_herbivory~chamber+treatment+age, data=pg, FUN=mean)


pg1<-aov(pgag$prop_herb~pgag$treatment)
summary.aov(pg1)

pg2<-aov(pg$prop_herb~pg$treatment)
summary.aov(pg2)

#THESE TWO
pg3<-aov(pgag$prop_herb~pgag$treatment+pgag$age)
summary.aov(pg3)

pg4<-aov(pg$prop_herb~pg$treatment+pg$age)
summary.aov(pg4)
##

pg10<-aov(percent_herbivory~treatment+age+chamber, data = pgag2)
summary.aov(pg10)

pg11<-aov(percent_herbivory~treatment+age, data = pgag2)

pg12<-aov(percent_herbivory~treatment, data = pgag2)

pg13<-aov(percent_herbivory~age, data = pgag2)

pg14<-aov(percent_herbivory~chamber, data = pgag2)

pg15<-aov(percent_herbivory~1, data = pgag2)

library(AICcmodavg)
modcomp<-aictab(cand.set=list(pg10,pg11,pg12,pg13, pg14, pg15),
				modnames=c("all","t+a","treat", "age", "chamber", "null"))#AIC table

modcomp

summary.aov(pg13)
summary.aov(pg11)

#summarizing % herb by age to add row to data
library(plyr)
sumtab <- ddply(pg, c("age"), summarise,
				N    = length(percent_herbivory),
				mean = mean(percent_herbivory),
				sd   = sd(percent_herbivory),
				se   = sd / sqrt(N))
sumtab

#adding summary rows to data
library(tidyverse)
#pg<-pg %>% add_row(chamber = "ALL", treatment = "All", age = "Y", percent_herbivory = 1.4608)
#pg<-pg %>% add_row(chamber = "ALL", treatment = "All", age = "M", percent_herbivory = 5.3598)

#changing names for plot
pg$age<-as.character(pg$age)
pg$age[pg$age=="Y"]="Young"
pg$age[pg$age=="M"]="Old"

pg$treatment<-as.character(pg$treatment)
pg$treatment[pg$treatment=="control chamber"]="Control (chamber)"
pg$treatment[pg$treatment=="natural control"]="Control (no chamber)"
pg$treatment[pg$treatment=="T°C"]="Temperature"
pg$treatment[pg$treatment=="T°C + CO2"]="Temp + CO2"

ggplot(pgag2, aes(x=treatment, y=percent_herbivory, fill=age))+geom_boxplot()+geom_point()

ggplot(pg, aes(x=treatment, y=prop_herb, fill=age))+geom_boxplot()+geom_point()

perherb_treat<-ggplot(pg, aes(treatment, percent_herbivory, color=age))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.position = "none",
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "% herbivorized")+ 
	facet_wrap(~ age)
perherb_treat

#EXPORT PLOT
tiff('perherb_treat.tiff', units="in", width=6, height=5, res=400)
perherb_treat
dev.off()

perherb_all<-ggplot(pg, aes(age, percent_herbivory, color=age))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.position = "none",
		  text = element_text(size=15), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "% herbivorized")
perherb_all

tiff('perherb_all.tiff', units="in", width=6, height=5, res=400)
perherb_all
dev.off()


library(lme4)
pg2<-lmer(prop_herb ~ treatment + age + (1|chamber), data=pgag, na.action = "na.omit")
summary(pg2)
library(car)
Anova(pg2)

lm1<-lm(prop_herb~treatment+age, data=pg)
summary(lm1)

library(multcomp)
summary(glht(lm1, linfct=mcp(treatment="Tukey")))


library(ggplot2)
library(Rmisc)
sebars <- summarySE(pgag, measurevar="prop_herb", groupvars=c("treatment", "age"))

pgplot_bar<-ggplot(sebars, aes(x=treatment, y=prop_herb, color=age))+geom_bar(stat = "identity")+
	geom_errorbar(aes(ymin=prop_herb-se, ymax=prop_herb+se), width=.2,
				  position=position_dodge(.9))
pgplot_bar

pgplot_point<-ggplot(sebars, aes(x=treatment, y=prop_herb, color=age))+geom_point(stat = "identity")+
	geom_errorbar(aes(ymin=prop_herb-se, ymax=prop_herb+se), width=.2,
				  position=position_dodge(.9))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.title = element_blank(),
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "Proportion herbivorized")
pgplot_point

ggplot(sebars, aes(x=treatment, y=prop_herb, color=age))+geom_point(stat = "identity")+
	geom_errorbar(aes(ymin=prop_herb-se, ymax=prop_herb+se), width=.2,
				  position=position_dodge(.9))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.position = "none", text = element_text(size=12), 
		  axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "Proportion herbivorized")+
	facet_wrap(~age)


pg$herb_pa[pg$percent_herbivory==0] <- 0
pg$herb_pa[pg$percent_herbivory>0] <- 1



pg$chamber<-as.character(pg$chamber)
modz<- glmer(herb_pa ~ treatment + treatment|chamber, family = binomial, data=pg, na.action = "na.fail")
#don't have sufficient no. of obs to support this model
summary(modz)
modq<- glm(herb_pa ~ treatment, family = binomial, data=pg)
summary(modq)
ggplot(pg, aes(x=treatment, y=herb_pa))+geom_bar(stat = "identity")

pa.ag<-aggregate(herb_pa~chamber+treatment, data=pg, FUN=mean)
pa.ag$herb_pa<-as.numeric(pa.ag$herb_pa)
ggplot(pa.ag, aes(treatment,herb_pa))+geom_boxplot()+geom_point()+theme_classic()
modx<- betareg(herb_pa ~ treatment, data=pa.ag, na.action = "na.fail")
summary(modx)
pg3<-aov(pa.ag$herb_pa~pa.ag$treatment)
summary.aov(pg3)
sebars3 <- summarySE(pa.ag, measurevar="herb_pa", groupvars=c("treatment"))

ggplot(sebars3, aes(x=treatment, y=herb_pa))+geom_point(stat = "identity")+
	geom_errorbar(aes(ymin=herb_pa-se, ymax=herb_pa+se), width=.2,
				  position=position_dodge(.9))
pgplot_bar

pg$prop_herb<-(pg$percent_herbivory/100)
pg$prop_herb<-as.numeric(pg$prop_herb)
pg$prop_herb[pg$prop_herb==0] <- NA

New <- data.frame("prop_herb" = pg$prop_herb, 
				  "treatment" = pg$treatment, 
				  "chamber" = pg$chamber) 

New<-na.omit(New)

modr<- betareg(prop_herb ~ treatment, data=New, na.action = "na.fail")
summary(modr)
emmeans(modr,pairwise~treatment,type="response")

sebars1 <- summarySE(New, measurevar="prop_herb", groupvars=c("treatment"))

plot1<-ggplot(sebars1, aes(x=treatment, y=prop_herb))+geom_point(stat = "identity")+
	geom_errorbar(aes(ymin=prop_herb-se, ymax=prop_herb+se), width=.2,
				  position=position_dodge(.9))+
	theme_classic()
plot1

ggplot(New, aes(x=treatment, y=prop_herb))+geom_boxplot()

a1<-aov(prop_herb~treatment,data=New)
summary.aov(a1)
summary(glht(a1, linfct=mcp(treatment="Tukey")))

#Phenolics

ga <- read.csv(file = "GA_StandardCurve.csv", head=T)

ga.ab<-lm(ga$abs_avg~ga$ab_val_mg)
summary(ga.ab)
plot(ga$abs_avg~ga$ab_val_mg)
#y=mx+b, y=3.511354x+0.056570, R^2=0.9974

ga.conc<-lm(ga$abs_avg~ga$concen_mgml)
summary(ga.conc)
plot(ga$abs_avg~ga$concen_mgml)
#y=mx+b, y=0.702271 +0.056570 , R^2=0.9974

#x=(y-b)/m

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

#ab_val/0.2mL (vol in well) = y/1.1mL (total volume in tube)

phen$pdw<-(((phen$ab_val_mg/0.2)*1.1)/(phen$start_wt))*100

#create new col for concentration
phen$concen<-(phen$abs_avg-0.056570)/0.702271

phen$stage<-as.character(phen$stage)
phen$stage[phen$stage=="y"]="Young"
phen$stage[phen$stage=="m"]="Old"

phen$treat<-as.character(phen$treat)
phen$treat[phen$treat=="control_chamber"]="Control (chamber)"
phen$treat[phen$treat=="control_nat"]="Control (no chamber)"
phen$treat[phen$treat=="TC"]="Temperature"
phen$treat[phen$treat=="TC+CO2"]="Temp + CO2"


ggplot(phen, aes(x=treat, y=concen))+geom_boxplot()+geom_point()

phen_ag<-aggregate(concen~treat+sample+stage, data=phen, FUN=mean)

phen_ag2<-aggregate(pdw~treat+sample+stage, data=phen, FUN=mean)



pg5<-aov(phen_ag$concen~phen_ag$treat+phen_ag$stage)
summary.aov(pg5)

pg6<-aov(pdw~treat+stage, data=phen_ag2)
summary.aov(pg6)

##

ggplot(phen_ag2, aes(treat, pdw, color=stage))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.position = "none",
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "%dw in gallic acid equivalents")+
	facet_wrap(~stage)

ggplot(phen_ag2, aes(stage, pdw, color=stage))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.position = "none",
		  text = element_text(size=15))+
	labs(x = "", y = "%dw in gallic acid equivalents")


ggplot(phen_ag, aes(treat, concen, color=stage))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.title = element_blank(),
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "Concentration (mg/mL)")+ 
	facet_wrap(~ stage)

ggplot(phen_ag, aes(treat, concen, color=stage))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.title = element_blank(),
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "Concentration (mg/mL)")

ggplot(phen_ag, aes(treat, concen))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	scale_color_manual(values = c("#006d2c", "#66c2a4"))+
	theme(legend.title = element_blank(),
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "Concentration (mg/mL)")

#growth
grow <- read.csv(file="Piper_growth.csv",head=TRUE)

grow$growth<-grow$ht.2018.09.cm-grow$ht.2018.04.cm
grow$rel_gro<-((grow$ht.2018.09.cm-grow$ht.2018.04.cm)/grow$ht.2018.04.cm)

grow <- grow[order(grow$Treatment),]
grow<-grow[-c(11:15),]#removing control (no chamber)

gro1<-aov(grow$growth~grow$Treatment)
summary.aov(gro1)

gro2<-aov(grow$rel_gro~grow$Treatment)
summary.aov(gro2)

ggplot(grow, aes(Treatment, growth))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	theme(legend.title = element_blank(),
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "Total growth (cm)")

ggplot(grow, aes(Treatment, rel_gro))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
	theme_classic()+
	theme(legend.title = element_blank(),
		  text = element_text(size=12), axis.text.x = element_text(angle=45, hjust=1))+
	labs(x = "", y = "Relative growth (cm)")
