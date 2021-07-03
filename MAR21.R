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
library(ggsignif)

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
#SRW: removed row below, I think you did this above?? I was getting an error
#phen_ag$pdw<-(phen_ag$pdw)*100

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
b10<-betareg(prop_gro~treatment+prop_herb+prop_dw, dat=all.dat3)
shapiro.test(resid(b10)) #normal!!

summary(b10)
#Pseudo R-squared: 0.3165
#y=mx+b, y= -4.774x - 0.099
#growth decreased 4.8% with every 1% increase in leaf herbivory 
Anova(b10)#herb sig (neg)

#USE THIS ONE
b11<-betareg(prop_gro~treatment+prop_herb, dat=all.dat3)
shapiro.test(resid(b11)) #normal!!
Anova(b11)#herb sig (neg)

summary(b11)
#Pseudo R-squared: 0.3165
#y=mx+b, y= -11.66x - 0.424
#growth decreased ~12% with every 1% increase in leaf herbivory 

summary(betareg(prop_gro~prop_herb, dat=all.dat3))
#y= -8.95 - 0.65
#growth decreased ~9% with every 1% increase in leaf herbivory

#GROWTH + HERBIVORY PLOT
ggplot(all.dat3, aes(prop_herb, prop_gro))+
	geom_smooth(color="black",method = "lm")+
	geom_jitter(position=position_jitter(width = 0.0), alpha=0.45, size=2.5, aes(color=treatment))+
	theme_classic()+
	theme(legend.position = "right",
		  text = element_text(size=19))+
	labs(x = "Proportion leaf herbivory", y = "Proportion change in height")+
	scale_color_viridis(discrete = T, option = "D")+
	scale_x_continuous(limits = c(0,.15))
#pseudo R^2=0.178

bmod<-betareg(prop_gro~prop_herb, dat=all.dat3)
library(rcompanion)
plotPredy(data  = all.dat3,
		  y     = prop_gro,
		  x     = prop_herb,
		  model = bmod,
		  xlab  = "Phenolics",
		  ylab  = "Proportion herb")

library(lmtest)
lrtest(b10)
lrtest(bmod)

summary(bmod)
#y=mx+b, y= -8.95 - 0.65
#Change in height decreased 8.95% with every 1% increase in herbivory

###PHENOLICS----
b3<-betareg(prop_dw~treat+stage, dat=all.dat)
shapiro.test(resid(b3)) #normal

Anova(b3, test.statistic="F")#stage sig, p=4.32e-07 ***
Anova(b3)

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
	geom_point(position=position_jitter(width = 0.025), alpha=0.4, aes(color=treat), size=2.5)+
	stat_summary(fun.data = "mean_se", colour="black", size=1)+
	theme_classic()+
	labs(x="", y="Total phenolics (prop. dw in GAE)")+
	theme(text = element_text(size=18), legend.position = "none")+
	scale_color_viridis(discrete = T, option = "D")+
	stat_summary(geom = 'text', label = c("A","B"),
				 fun = max, vjust = -1.5, size = 5.8)+
	scale_y_continuous(limits = c(0,.115))

#PHENOLICS SUMMARY STATS
library(plyr)
phen.tab <- ddply(all.dat, c("stage"), summarise,
				  N    = length(prop_dw),
				  mean = mean(prop_dw),
				  sd   = sd(prop_dw),
				  se   = sd / sqrt(N))
phen.tab

#young leaf avg pdw/old leaf avg pdw
(0.06859805-0.04970515)/0.04970515
#0.380099
#Young leaves had an average of 38% times more total phenolics

###HERBIVORY---
#young leaves
all.dat2$prop_herb_Young1<-all.dat2$prop_herb_Young+0.00001
b5<-betareg(prop_herb_Young1~treat+growth+prop_dw_Young, dat=all.dat2)
summary(b5) #no sig, so no plots?
Anova(b5)

#mature leaves
all.dat2$prop_herb_Mature1<-all.dat2$prop_herb_Mature+0.00001
b8<-betareg(prop_herb_Mature1~treat+prop_dw_Mature+growth, dat=all.dat2)
summary(b8) #T+CO2 treat (pos), growth (neg), and phenolics (neg) all significant
##herbivory decreased ~30% with every 1% increase in phenolic concentration
Anova(b8)
Anova(b8, test.statistic="F")

b80<-betareg(prop_herb_Mature1~treat+prop_dw_Mature, dat=all.dat2)
summary(b80)
#if remove growth from model, nothing is significant....

summary(betareg(prop_herb_Mature1~prop_dw_Mature, dat=all.dat2))
#herb~chem
#y= -17.14 - 2.05
#herbivory decreased ~17% with every 1% increase in phenolic concentration
summary(betareg(prop_herb_Mature1~growth, dat=all.dat2))
#herb~growth
#y = -2.36 -1.96
#herbivory decreased 2.4% with every 1% increase in growth


#HERBIVORY + GROWTH (mature)
ggplot(all.dat2, aes(growth, prop_herb_Mature1))+
	geom_smooth(color="black",method = "glm")+
	geom_jitter(position=position_jitter(width = 0.0), alpha=0.4, aes(color=treat), size=2.5)+
	theme_classic()+
	theme(legend.position = "none",
		  text = element_text(size=19))+
	labs(x = "Proportion change in height", y = "Proportion leaf herbivory")+
	scale_color_viridis(discrete = T, option = "D")


#HERBIVORY + PHENOLICS (mature)
ggplot(all.dat2, aes(prop_dw_Mature, prop_herb_Mature1))+
	geom_smooth(color="black",method = "glm")+
	geom_jitter(position=position_jitter(width = 0.0), alpha=0.4, aes(color=treat), size=2.5)+
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
	geom_point(aes(color=treat),position=position_jitter(width = 0.04), alpha=0.45, size=2.5)+
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

treat.tab <- ddply(all.dat2, c("treat"), summarise,
				  N    = length(prop_herb_Mature1),
				  mean = mean(prop_herb_Mature1),
				  sd   = sd(prop_herb_Mature1),
				  se   = sd / sqrt(N))
treat.tab

(0.10494750-0.02532667)/0.02532667
#3.14
#Mature leaves in CO2+temp experienced 3.14 times more herbivory than mature leaves in 
#environments with increased temp

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





####SRW trying some things

hist(all.dat$percent_herbivory)
hist(all.dat$growth)
hist(all.dat$pdw)

hist(all.dat3$prop_herb)
hist(all.dat3$prop_gro)
hist(all.dat3$pdw)


m1<-betareg(prop_gro~pdw*treatment, dat=all.dat3)
shapiro.test(resid(m1)) #normal!!
Anova(m1)  #significant interaction

library(dplyr)
all.dat3 %>%
	ggplot(aes(x=pdw, 
			   y=prop_gro, 
			   color=treatment))+
	geom_point()+
	geom_smooth(method="lm")



m2<-betareg(prop_herb~prop_gro*treatment, dat=all.dat3)
shapiro.test(resid(m2)) #normal!!
Anova(m2)  #significant interaction

all.dat3 %>%
	ggplot(aes(x=prop_gro, 
			   y=prop_herb, 
			   color=treatment))+
	geom_point()+
	geom_smooth(method="lm")



m3<-betareg(prop_herb~pdw*treatment, dat=all.dat3)
shapiro.test(resid(m3)) #normal!!
Anova(m3)  #significant interaction

all.dat3 %>%
	ggplot(aes(x=pdw, 
			   y=prop_herb, 
			   color=treatment))+
	geom_point()+
	geom_smooth(method="lm")


###Okay, so there are some interesting interactions going on here...
#In CO2+temp treatment, the plants that are doing the best in terms of 
#growth are getting low levels of herbivory and also producing high
#phenolics. Also we see the predicted relationship between high phenolics-low herbivory

#In control treatments, and most other treatments, 
#basically we don't see these trade-offs

#I am thinking maybe this is just driven by a general increase in herbivory in 
#the CO2+temp chamber. So herbivory is becoming more apparent and the 
#role of defense is becoming more important

#We never see any growth-defense trade-offs, where plants that are 
#growing more producing fewer defenses, though it is possible
#this would have been more apparent in the absence of any herbivory


#I also think it would be worth repeating this, but looking separately at
#young leaves and old leaves rather than




#Another thing to look at...are changes in herbivory mediated by shifts in 
#ontogenetic patterns?? It could be that in rapidly growing plants
#under climate change, more of the herbivory is shifted to young leaves??

all.dat2$prop_dw_YminusM <- all.dat2$prop_dw_Young - all.dat2$prop_dw_Mature
all.dat2$prop_herb_YminusM <- all.dat2$prop_herb_Young - all.dat2$prop_herb_Mature

hist(all.dat2$prop_dw_YminusM)
hist(all.dat2$prop_herb_YminusM)
plot(prop_dw_YminusM ~ treat, data=all.dat2)
plot(prop_herb_YminusM ~ treat, data=all.dat2)

plot(prop_herb_Young ~ treat, data=all.dat2)
plot(prop_herb_Mature ~ treat, data=all.dat2)


m1 <- lm(prop_dw_YminusM ~ treat, data=all.dat2)
summary(m1)
anova(m1)

m2 <- lm(prop_herb_YminusM ~ treat, data=all.dat2)
summary(m2)
anova(m2)


all.dat %>%
	ggplot(aes(stage,pdw, color=treat)) +
	geom_point(aes(fill=treat),size=3) +
	geom_line(aes(group = chamber))


all.dat %>%
	ggplot(aes(stage,percent_herbivory, color=treat)) +
	geom_point(aes(fill=treat),size=3) +
	geom_line(aes(group = chamber))


