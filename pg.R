pg <- read.csv(file="Piper_herbivory.csv",head=TRUE)

pg$percent_herbivory<-as.numeric(pg$percent_herbivory)
pg$treatment<-as.character(pg$treatment)
pg <- pg[order(pg$percent_herbivory),]

#pg <- pg[1:98,] #removing 2 rows with NAs

ggplot(pg, aes(x=treatment, y=percent_herbivory))+geom_boxplot()+geom_point()

hist(pg$percent_herbivory)#skewed, zero-inflated
shapiro.test(pg$percent_herbivory)#very not normal

pg$prop_herb<-(pg$percent_herbivory/100)
pg$prop_herb<-as.numeric(pg$prop_herb)

library(betareg)
beta.pg<-betareg(prop_herb~treatment, data = pg)
#won't run bc zero-inflated

pgag<-aggregate(prop_herb~chamber+treatment, data=pg, FUN=mean)
ggplot(pgag, aes(x=treatment, y=prop_herb))+geom_boxplot()+geom_point()

betaa<-betareg(prop_herb~treatment, data=pgag)
summary(betaa)
library(car)
Anova(betaa)
library(emmeans)
emmeans(betaa,pairwise~treatment,type="response")


#figure out how to load these to library
#library(rjags)
#library(zoib)
#zo<-zoib(prop_herb~treatment, data = pg)
#zero-inflated models

	
pg1<-aov(pgag$prop_herb~pgag$treatment)
summary.aov(pg1)

pg2<-aov(pg$prop_herb~pg$treatment)
summary.aov(pg2)

library(lme4)
pg2<-lmer(prop_herb ~ treatment + (1|chamber), data=pgag, na.action = "na.omit")
summary(pg2)
library(car)
Anova(pg2)

lm1<-lm(prop_herb~treatment, data=pg)
summary(lm1)

library(multcomp)
summary(glht(lm1, linfct=mcp(treatment="Tukey")))


library(ggplot2)
library(Rmisc)
sebars <- summarySE(pgag, measurevar="prop_herb", groupvars=c("treatment"))

pgplot_bar<-ggplot(sebars, aes(x=treatment, y=prop_herb))+geom_bar(stat = "identity")+
	geom_errorbar(aes(ymin=prop_herb-se, ymax=prop_herb+se), width=.2,
				  position=position_dodge(.9))
pgplot_bar

pgplot_point<-ggplot(sebars, aes(x=treatment, y=prop_herb))+geom_point(stat = "identity")+
	geom_errorbar(aes(ymin=prop_herb-se, ymax=prop_herb+se), width=.2,
				  position=position_dodge(.9))+
	theme_classic()
pgplot_point


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
