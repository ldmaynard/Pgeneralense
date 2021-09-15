#Interactive effects of climate change, leaf age, and secondary metabolites 
#on plant growth, defense, and herbivory.

#LOADING LIBRARIES----------------------------------------------------
{library(lme4)
library(ggplot2)
library(car)
library(multcomp)
library(dplyr)
library(AICcmodavg)
library(betareg)
library(MuMIn)
library(car)
library(corrplot)
library(viridis)
library(ggsignif)
library(emmeans)
library(multcomp)
library(plyr)
	}

#LOADING & WRANGLING DATA----------------------------------------------------

##GROWTH

#load growth data
grow <- read.csv(file="Piper_growth.csv",head=TRUE)
table(grow$Treatment)   
#SRW: Why are there 6 T chambers and 4 T+CO2???
#LDM: Unsure, but its consistent across all data. Assuming the CO2 failed in one of the T+CO2 treatments

#renaming column
colnames(grow)[1] <- "chamber"

grow$total_gro <- grow$ht.2018.09.cm-grow$ht.2018.04.cm
grow$prop_gro<-grow$total_gro/grow$ht.2018.04.cm
grow$per_gro<-grow$prop_gro*100

#creating dataset with all treatments
grow.all<-grow

#creating dataset without "No chamber" treatment
grow <- grow[order(grow$Treatment),]
grow<-grow[-c(11:15),] #removing control (no chamber)

hist(grow$prop_gro) 
shapiro.test(grow$prop_gro)#normal
hist(grow$total_gro)
shapiro.test(grow$total_gro)#so is this one. Both yield the same results


##CHEMISTRY

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


##HERBIVORY

#load herbivory data
herb.all <- read.csv(file="Piper_herbivory.csv",head=TRUE)

#creating col for proportion herbivory
herb.all$percent_herbivory<-as.numeric(herb.all$percent_herbivory)
herb.all$prop_herb<-(herb.all$percent_herbivory/100)
herb.all$prop_herb<-as.numeric(herb.all$prop_herb)

#creating dataset aggregated by leaf age  by age (n=40)
herb.50<-aggregate(prop_herb~chamber + age + treatment +sample,data=herb.all,FUN=mean)

#aggregating by chamber, n=20
herb.25<-aggregate(prop_herb~chamber+treatment,data=herb.all,FUN=mean)

#creating dataset without control (no chamber) = herb.80 (n=80)
herb.all <- herb.all[order(herb.all$treatment),]
herb.80<-herb.all[-c(41:60),]#removing control (no chamber)

#creating dataset with average herbivory  by age (n=40)
herb.40<-aggregate(prop_herb~chamber + age + treatment+sample,data=herb.80,FUN=mean)

#aggregating by chamber, n=20
herb.20<-aggregate(prop_herb~chamber+treatment,data=herb.80,FUN=mean)

##CREATING COMBINED DATASETS----

#all data, n=100
all.dat100<-merge(phen.all, herb.all, by="sample", all = F)
head(all.dat100)
#this creates two chamber columns=chamber.x and chamber.y
colnames(all.dat100)[4] <- "chamber"
all.dat100<-merge(all.dat100, grow.all, by="chamber", all=T)
#cleaning up, selecting columns
all.dat100<-select(all.dat100, chamber, stage, pdw, treatment, sample, ID, total.area.cm2,
				   real.area.cm2, percent_herbivory, prop_herb, ht.2018.04.cm, ht.2018.09.cm,
				   total_gro, prop_gro)

#all treatment data, aggregated by leaf age, n=50
##adding prop dry weight column to herbivory data
all.dat50<-merge(phen.all, herb.50, by="sample", all = T)
head(all.dat50)
colnames(all.dat50)[4] <- "chamber"
all.dat50<-merge(all.dat50, grow.all, by="chamber", all=T)
#cleaning up, selecting columns
all.dat50<-select(all.dat50, chamber, stage, pdw, treatment, sample, 
				  prop_herb, ht.2018.04.cm, ht.2018.09.cm,total_gro, prop_gro)

#all treatment data, aggregated by chamber, n=25
all.dat25<-merge(grow.all, herb.25,by="chamber", all = T)
all.dat25<-merge(all.dat25, phen.25,by="chamber", all = T)
#cleaning up, selecting columns
all.dat25<-select(all.dat25, chamber, pdw, treatment, 
				  prop_herb, ht.2018.04.cm, ht.2018.09.cm,total_gro, prop_gro)

#four treatment data, not aggregated, n=80
all.dat80<-merge(herb.80, phen.40, by="sample", all = T)
head(all.dat80)
colnames(all.dat80)[4] <- "chamber"
all.dat80<-merge(all.dat80, grow, by="chamber", all=T)
#cleaning up, selecting columns
all.dat80<-select(all.dat80, chamber, stage, pdw, treatment, sample, ID, total.area.cm2,
				   real.area.cm2, percent_herbivory, prop_herb, ht.2018.04.cm, ht.2018.09.cm,
				   total_gro, prop_gro)

#four treatment data, aggregated by leaf age, n=40
all.dat40<-merge(herb.40, phen.40, by="sample", all = T)
head(all.dat40)
colnames(all.dat40)[2] <- "chamber"
all.dat40<-merge(all.dat40, grow, by="chamber", all=T)
#cleaning up, selecting columns
all.dat40<-select(all.dat40, chamber, stage, pdw, treatment, sample, prop_herb, ht.2018.04.cm, 
				  ht.2018.09.cm, total_gro, prop_gro)

#four treatment data, aggregated by chamber, n=20
all.dat20<-merge(herb.20, phen.20, by="chamber", all = T)
all.dat20<-merge(all.dat20, grow, by="chamber", all=T)
#cleaning up, selecting columns
all.dat20<-select(all.dat20, chamber, pdw, treatment, 
				  prop_herb, ht.2018.04.cm, ht.2018.09.cm,total_gro, prop_gro)

#for four treatment data, re-ordering factor levels so model will compare everything to control chamber
all.dat80$treatment<-as.factor(all.dat80$treatment)
levels(all.dat80$treatment)
all.dat80$treat <- factor(all.dat80$treat, levels=c("control chamber", "CO2", "T°C", "T°C + CO2" ))

all.dat40$treatment<-as.factor(all.dat40$treatment)
levels(all.dat40$treatment)
all.dat40$treat <- factor(all.dat40$treat, levels=c("control chamber", "CO2", "T°C", "T°C + CO2" ))

all.dat20$treatment<-as.factor(all.dat20$treatment)
levels(all.dat20$treatment)
all.dat20$treat <- factor(all.dat20$treat, levels=c("control chamber", "CO2", "T°C", "T°C + CO2" ))

#---

##CHECKING FOR COLINEARITY----
samp1<-all.dat100[,-c(1,2,4:9,11:13)]#cat vars

L1 <- cor(samp1)#correlation matrix
corrplot(L1, method = "circle")
corrplot(L1, method = "number")
#low colinearity

check1 <- lm(prop_gro ~ prop_herb+pdw+stage+treatment, data=all.dat100)
summary(check1)
vif(check1)


##JULY ANALYSIS----
##QUESTION 1----

#GROWTH
#all treatments, n=25
mod.gro1<-(betareg(prop_gro~treatment, dat=all.dat25))
Anova(mod.gro1)
#p=0.5431, no effect of treatment on growth, no effect of chamber
shapiro.test(resid(mod.gro1)) #residuals  normally distributed

#growth plot, n=25
ggplot(data=all.dat25, aes(x=treatment, y=prop_gro))+ 
	geom_point(position=position_jitter(width = 0.025), alpha=0.4, size=2.5)+
	stat_summary(fun.data = "mean_se", colour="black", size=1)+
	theme_classic()

#four treatment data
mod.gro2<-(betareg(prop_gro~treatment, dat=all.dat20))
Anova(mod.gro2)
#p=0.82, no effect of treatment 
shapiro.test(resid(mod.gro2)) #residuals  normally distributed

#growth plot, n=20
ggplot(data=all.dat20, aes(x=treatment, y=prop_gro))+ 
	geom_point(position=position_jitter(width = 0.025), alpha=0.4, size=2.5)+
	stat_summary(fun.data = "mean_se", colour="black", size=1)+
	theme_classic()

#CHEMISTRY
library(glmmTMB)

#all treatments, n=50
chem.mod1 <- glmmTMB(pdw ~ treatment + (1|chamber), data = all.dat50, family = "beta_family")
summary(chem.mod1)
Anova(chem.mod1)
#treatment p=0.54
drop1(chem.mod1, test="Chisq")
#treatment p=0.57

#chemistry plot, n=50
ggplot(data=all.dat50, aes(x=treatment, y=pdw))+ 
	geom_point(position=position_jitter(width = 0.025), alpha=0.4, size=2.5)+
	stat_summary(fun.data = "mean_se", colour="black", size=1)+
	theme_classic()+
	labs(y="Phenolics (proportion dw)")

#four treatments, n=40
chem.mod2 <- glmmTMB(pdw ~ treatment + (1|chamber), data = all.dat40, family = "beta_family")
summary(chem.mod2)
Anova(chem.mod2)
#treatment p=0.76
drop1(chem.mod2, test="Chisq")
#treatment p=0.77

#chemistry plot, n=40
ggplot(data=all.dat40, aes(x=treatment, y=pdw))+ 
	geom_point(position=position_jitter(width = 0.025), alpha=0.4, size=2.5)+
	stat_summary(fun.data = "mean_se", colour="black", size=1)+
	theme_classic()+
	labs(y="Phenolics (proportion dw)")



#HERBIVORY

all.dat100$prop_herb1<-all.dat100$prop_herb+0.0001
all.dat80$prop_herb1<-all.dat80$prop_herb+0.0001

#all treatments, n=100
herb.mod1 <- glmmTMB(prop_herb1 ~ treatment + (1|chamber), data = all.dat100, family = "beta_family")
summary(herb.mod1)
Anova(herb.mod1)
#treatment p=0.62
drop1(herb.mod1, test="Chisq")
#treatment p=0.61

#herbivory plot, n=100
ggplot(data=all.dat100, aes(x=treatment, y=prop_herb1))+ 
	geom_point(position=position_jitter(width = 0.025), alpha=0.4, size=2.5)+
	stat_summary(fun.data = "mean_se", colour="black", size=1)+
	theme_classic()+
	labs(y="Proportion herbivory")

#four treatments, n=80
herb.mod2 <- glmmTMB(prop_herb1 ~ treatment + (1|chamber), data = all.dat80, family = "beta_family")
summary(herb.mod2)
Anova(herb.mod2)
#treatment p=0.68
drop1(herb.mod2, test="Chisq")
#treatment p=0.68
#kind of strange the temp+co2 treatment isn't coming out significant

#herbivory plot, n=80
ggplot(data=all.dat80, aes(x=treatment, y=prop_herb1))+ 
	geom_point(position=position_jitter(width = 0.025), alpha=0.4, size=2.5)+
	stat_summary(fun.data = "mean_se", colour="black", size=1)+
	theme_classic()+
	labs(y="Proportion herbivory")


##QUESTION 2A----

#Data not including natural control, n=40
chem.mod2a <- glmmTMB(pdw ~ treatment + stage + (1|chamber), data = all.dat40, family = "beta_family")
summary(chem.mod2a)
#stage p<0.00001
Anova(chem.mod2a)
#stage p<0.00001, treatment p=0.71
drop1(chem.mod2a, test="Chisq")
#stage p<0.00001, treatment p=0.72

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
#Young leaves had an average of 38% times more total phenolics than mature leaves

#plot
all.dat40 %>%
	ggplot(aes(stage,pdw, color=treatment)) +
	geom_point(aes(fill=treatment),size=3) +
	geom_line(aes(group = chamber))+
	theme_classic()+
	labs(y="Total phenolics (prop. dw)")

##QUESTION 2B----
#Growth-defense trade-off

#Opcion uno
grow.mod2b<-glmmTMB(pdw ~ treatment * prop_gro + (1|chamber), data = all.dat40, family = "beta_family")
summary(grow.mod2b)
#based on summary, temp treatment (p=0.01) and it's interaction with growth (p=0.02) are sig
Anova(grow.mod2b)
#However, not based on an ANOVA, growth p=0.33, treatment p=0.82, interaction p=0.12
drop1(grow.mod2b, test="Chisq")
#Or drop1, growth p=0.34, treatment p=0.83, interaction p=0.14

joint_tests((grow.mod2b), by = "treatment")
#significant effect of growth on defense in  temperature treatment
#while usually there is a negative relationship between growth and defense, 
#plants in the temp treatment that grew more also had higher defenses in mature leaves

#Growth*defense plot, all leaves
all.dat40 %>%
	ggplot(aes(x=prop_gro, 
			   y=pdw,
			   color=treatment))+
	geom_point()+
	geom_smooth(method="lm")

#Opcion dos
#Two models, where leaf age is split
all.dat40<- all.dat40[order(all.dat40$stage),]
all.dat40_m <- all.dat40[c(1:20),]
all.dat40_y <- all.dat40[c(21:40),]

#mature leaves, interactive model
grow.mod2m_i<-glmmTMB(pdw ~ treatment * prop_gro + (1|chamber), data = all.dat40_m, family = "beta_family")
summary(grow.mod2m_i)
Anova(grow.mod2m_i)
#treatment p=.93, growth=0.45, interaction p=0.034

joint_tests((grow.mod2m_i), by = "treatment")
#significant effect of growth on defense in  temperature treatment
#while usually there is a negative relationship between growth and defense, 
#for mature leaves in temperature experiment, plants that grew more also had higher defenses in mature leaves

#Growth*defense plot, mature leaves
all.dat40_m %>%
	ggplot(aes(x=prop_gro, 
			   y=pdw,
			   color=treatment))+
	geom_point()+
	geom_smooth(method="lm")+
	labs(title = "Mature leaves", y="Total phenolics (prop. dw)", x="Proportion growth")+
	theme_classic()

#young leaves, interactive model
grow.mod2y_i<-glmmTMB(pdw ~ treatment * prop_gro + (1|chamber), data = all.dat40_y, family = "beta_family")
summary(grow.mod2y_i)
Anova(grow.mod2y_i)
#treatment p=.42, chemistry p=0.23, interaction p=0.055

joint_tests((grow.mod2y_i), by = "treatment")
#marginally significant effect of growth on defense in combo treatment
#while usually there is a negative relationship between growth and defense, 
#young leaves in combo treatment, plants that grew more also had higher defenses in young leaves

#Growth*defense plot, young leaves
all.dat40_y %>%
	ggplot(aes(x=prop_gro, 
			   y=pdw,
			   color=treatment))+
	geom_point()+
	geom_smooth(method="lm")+
	labs(title = "Young leaves", y="Total phenolics (prop. dw)", 
		 	x="Proportion growth")+
	theme_classic()


##SRW: I guess I like the opcion uno better for this one?? But we could say it is mostly driven by the 
#young leaves and could even show those separately in a supplement. 
#I am still struggling with how to interpret this though...basically when you have increased temp, 
#plants that do better just do better in terms of everything...growth and phenolics, but there is no "trade-off"
#i.e. a negative relationship for any of the plants or treatments.
#Maybe this just has to do with physiological differences across individuals in heat tolerance?? Where those
#that are tolerant of higher temps just do better generally

#I would also be curious to see the data for the no chamber control here...though it would probably just make things
#more confusing

#WITH ALL TREATMENTS
#Opcion uno
grow.mod2b_all<-glmmTMB(pdw ~ treatment * prop_gro + (1|chamber), data = all.dat50, family = "beta_family")
summary(grow.mod2b_all)
#again, based on summary, temp treatment and it's interaction with growth are sig
Anova(grow.mod2b_all)
#However, not based on an ANOVA

joint_tests((grow.mod2b_all), by = "treatment")
#same result, significant effect of growth on defense in  temperature treatment
#while usually there is a negative relationship between growth and defense, 
#plants in the temp treatment that grew more also had higher defenses in mature leaves

#Growth*defense plot, all leaves
all.dat50 %>%
	ggplot(aes(x=prop_gro, 
			   y=pdw,
			   color=treatment))+
	geom_point()+
	geom_smooth(method="lm")

##QUESTION 2C----
#Relative change in the effectiveness of defense

herb.chem.mod2c<-glmmTMB(prop_herb1 ~ treatment * pdw + (1|chamber), data = all.dat80, family = "beta_family")
summary(herb.chem.mod2c)
Anova(herb.chem.mod2c)
#chemistry p=0.0003, treatment p=0.22, interaction p=0.28


#Chem plot
all.dat80 %>%
	ggplot(aes(x=pdw, 
			   y=prop_herb1))+
	geom_point(aes(color=stage))+
	geom_smooth(method = 'lm')+
	labs(x="Total phenolics (prop. dw)", y="Proportion herbivory")+
	theme_classic()


#Option 2, split leave ages
all.dat80<- all.dat80[order(all.dat80$stage),]
all.dat80_m <- all.dat80[c(1:40),]
all.dat80_y <- all.dat80[c(41:80),]

#Young leaves
herb.chem.mod2c_y<-glmmTMB(prop_herb1 ~ treatment * pdw + (1|chamber), data = all.dat80_y, family = "beta_family")
summary(herb.chem.mod2c_y)
Anova(herb.chem.mod2c_y)
#nothing significant

#Mature leaves
herb.chem.mod2c_m<-glmmTMB(prop_herb1 ~ treatment * pdw + (1|chamber), data = all.dat80_m, family = "beta_family")
summary(herb.chem.mod2c_m)
Anova(herb.chem.mod2c_m)
#treatment significant p=0.047, chemistry marg sig p=0.053, interaction not sig
#in mature leaves, increasing defenses had a marginally sign neg effect on herbivory
#in mature leaves, 

#Chem plot
all.dat80_m %>%
	ggplot(aes(x=pdw, 
			   y=prop_herb1))+
	geom_point(aes(color=stage))+
	geom_smooth(method = 'lm')

#herbivory plot, n=80
ggplot(data=all.dat80_m, aes(x=treatment, y=prop_herb1))+ 
	geom_point(position=position_jitter(width = 0.025), alpha=0.4, aes(color=treatment), size=2.5)+
	stat_summary(fun.data = "mean_se", colour="black", size=1)+
	theme_classic()+
	theme(legend.position = "none")+
	labs(y="Proportion herbivory", title="Mature leaves")

#clds
m2c<-emmeans(herb.chem.mod2c_m,pairwise~treatment, type="response")
cld(m2c$emmeans,  Letters ='abcde')
#all in same group

#additive model
a1<-glmmTMB(prop_herb1 ~ treatment + pdw + (1|chamber), data = all.dat80_m, family = "beta_family")
a2<-emmeans(a1,pairwise~treatment, type="response")
cld(a2$emmeans,  Letters ='abcde')
#same result

#SUMMARY STATS
sum.tab <- ddply(all.dat80_m, c("treatment"), summarise,
				  N    = length(prop_herb1),
				  mean = mean(prop_herb1),
				  sd   = sd(prop_herb1),
				  se   = sd / sqrt(N))
sum.tab

#trying tukey
summary(glht(herb.chem.mod2c_m, linfct=mcp(treatment="Tukey")))
summary(glht(a1, linfct=mcp(treatment="Tukey")))
#in additive model, combo treatment is sig higher than temp treatment


##OLD ANALYSES/BRAIN DUMPS----
#SRW: there is one outlier here, which I would remove--I feel something
#definitely went wrong with that one...prob should remove from aggregate values also
d.temp <- phen.all[which(phen.all$pdw>0),]
hist(d.temp$pdw)  #also this is beautifully normal, so we I don't think we need betareg anyway
#SRW: another approach for these initial Q1 analyses--don't aggregate any given
#dataset, i.e. use all the phenolics samples you have when that is the response
#I think we discussed this and a major issue was that you can't use the random effects in the betareg. 
#However, I was just looking at this very helpful paper
#Douma and Weedon 2019: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13234
#they give some examples of beta regression with random effects using the package glmmTMB

#note in the paper I mention above, they suggest rescaling the dataset instead of just adding a constant
#From the paper, Appendix S3--
#A suggested rescaling equation is:
#	x∗i=xi(n−1)+0.5n
#Where x∗i is the transformation of xi and n is the total number of observations in the dataset.
#For convenience we define this as a custom function tranform01 and apply it to the dataset:
#	transform01 <- function(x) {
#		(x * (length(x) - 1) + 0.5) / (length(x))
#	}
#andrew2$ALGAE.scaled <- transform01(andrew2$ALGAE.mean)

#I did not mess with this...

#SRW: I like using all the data, but it really bugs me that we can't use a random
#effect with the betareg because without it the data are pseudoreplicated. Maybe we could try running this with
#the glmmTMB package??

#another option with this distribution could be a hurdle model, where we 
#first assess whether leaves had herbivory or not 0/1 with a binomial model, then do 
#another model just for those that received herbivory. I think that could work with lmer or with glmmTMB
#I started this below with lmer, but later found the glmmTMB and I have not tried that

hist(pg1$percent_herbivory)
pg1$herb_pa <- ifelse(pg1$percent_herbivory==0, 0, 1)

m1  <- glmer(herb_pa ~ treatment * age + (1|chamber), data=pg1, family=binomial)
m2  <- glmer(herb_pa ~ treatment + age + (1|chamber), data=pg1, family=binomial)
#singular fit error
drop1(m1, test="Chisq") #but no effect
drop1(m2, test="Chisq")

d.temp <- pg1[which(pg1$herb_pa==1),]
hist(d.temp$percent_herbivory)  #still very not normal...would have to transform for lmer
d.temp$ph_tr <- asin(sqrt(d.temp$prop_herb))
hist(d.temp$ph_tr) #better
d.temp$ph_tr2 <- logit(d.temp$prop_herb)
hist(d.temp$ph_tr2) #much better

m1 <- lmer(ph_tr2 ~ treatment * age + (1|chamber), data=d.temp)
m2 <- lmer(ph_tr2 ~ treatment + age + (1|chamber), data=d.temp)
summary(m1)
drop1(m1, test="Chisq") #no effect
drop1(m2, test="Chisq") #no effect

#Honestly I don't think any of these things will change the outcome, but I still think dialing in the 
#statistical approach for proportions would be good for this paper. It seems to keep coming
#up in our lives!


#I also have the feeling that age should always be in the models
#since we know it has such a huge effect on everything

#Playing around with these things below...

hist(phen$pdw)

#there is one outlier here, which I would remove--I feel something
#definitely went wrong with that one...prob should remove from aggregate values also
d.temp <- phen[which(phen$pdw>0),]
hist(d.temp$pdw)  #also this is beautifully normal, so we I don't think we need betareg anyway


library(glmmTMB)

d.temp$prop_dw <- d.temp$pdw/100

m1 <- glmmTMB(prop_dw ~ treat + stage + (1|chamber), data = d.temp, family = "beta_family")
#I am getting some warnings here about the way family is specified, but it is running

summary(m1)
drop1(m1, test="Chisq")
library(car)
Anova(m1)


#Also could try with lmer, since the data are so normal
m1 <- lmer(pdw ~ treat + stage + (1|chamber), data=d.temp)
summary(m1)
drop1(m1, test="Chisq")
boxplot(pdw ~ treat, data=d.temp)


#Option uno
#N=20, all leaves in each chamber are combined
Anova((betareg(prop_gro~treatment*pdw, data = all.dat3)))
#interaction significant p=0.02; treatment and chemistry alone are not
shapiro.test(resid(betareg(prop_gro~treatment*pdw, data = all.dat3)))
#normal

#joint_tests() function that obtains and tests the interaction contrasts 
#for all effects in the model and compiles them in one Type-III-ANOVA-like table
library (emmeans)
joint_tests((betareg(prop_gro~treatment*pdw, data = all.dat3)), by = "treatment")
#significant effect of  pdw on growth in both temperature treatments

#Growth*defense plot, all leaves combined
all.dat3 %>%
	ggplot(aes(x=prop_gro, 
			   y=pdw,
			   color=treatment))+
	geom_point()+
	geom_smooth(method="lm")



###GROWTH (using data averaged across all four leaves)
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

###PHENOLICS
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


m1<-betareg(prop_gro~pdw*treatment, dat=all.dat3,na.action = "na.fail")
shapiro.test(resid(m1)) #normal!!
Anova(m1)  #significant interaction

m1.add<-betareg(prop_gro~pdw+treatment, dat=all.dat3)
m1.t<-betareg(prop_gro~treatment, dat=all.dat3)
m1.c<-betareg(prop_gro~pdw, dat=all.dat3)
m1.null<-betareg(prop_gro~1, dat=all.dat3)

modcomp.grow<-aictab(cand.set=list(m1, m1.add, m1.t, m1.c, m1.null),
					 modnames=c("interact", "add","treat", "chem", "null"), REML=F)#AIC table
modcomp.grow
#better model of fit = null model, next is chem with dAICc=2.41

growd<-dredge(m1)
grow.avg<-model.avg(growd, subset=delta<4)#null and chem models
summary(grow.avg)#no sig

q1<-emmip(m1,pairwise~treatment|pdw, type="response")
q1



m1.emm<-emmeans(m1, ~treatment*pdw)
pairs(m1.emm, simple = "treatment")

#joint_tests() function that obtains and tests the interaction contrasts 
#for all effects in the model and compiles them in one Type-III-ANOVA-like table
joint_tests(m1, by = "treatment")
#significant effect of  pdw on growth in both co2 treatments



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

joint_tests(m2, by = "treatment")
#sig effect of growth on herbivory in co2+temp treatment

m2.add<-betareg(prop_herb~prop_gro+treatment, dat=all.dat3)
m2.t<-betareg(prop_herb~treatment, dat=all.dat3)
m2.g<-betareg(prop_herb~prop_gro, dat=all.dat3)
m2.null<-betareg(prop_herb~1, dat=all.dat3)

modcomp.herb<-aictab(cand.set=list(m2, m2.add, m2.t, m2.g, m2.null),
					 modnames=c("interact", "add","treat", "grow", "null"), REML=F)#AIC table
modcomp.herb
#top model is growth, followed narrowly by null (.97) then add (5.12)


all.dat3 %>%
	ggplot(aes(x=prop_gro, 
			   y=prop_herb, 
			   color=treatment))+
	geom_point()+
	geom_smooth(method="lm")



m3<-betareg(prop_herb~pdw*treatment, dat=all.dat3, na.action = "na.fail")
shapiro.test(resid(m3)) #normal!!
Anova(m3)  #significant interaction

joint_tests(m3, by = "treatment")
#sig effect of chemistry on herbivory in the additive cc treatment

q2<-emmeans(m3,pairwise~treatment, type="response")

m3.add<-betareg(prop_herb~pdw+treatment, dat=all.dat3)
m3.t<-betareg(prop_herb~treatment, dat=all.dat3)
m3.c<-betareg(prop_herb~pdw, dat=all.dat3)
m3.null<-betareg(prop_herb~1, dat=all.dat3)

modcomp.herb2<-aictab(cand.set=list(m3, m3.add, m3.t, m3.c, m3.null),
					 modnames=c("interact", "add","treat", "chem", "null"), REML=F)#AIC table
modcomp.herb2

herb2d<-dredge(m3)
herb2.avg<-model.avg(herb2d, subset=delta<4)#null and chem models
summary(herb2.avg)#no sig

all.dat3 %>%
	ggplot(aes(x=pdw, 
			   y=prop_herb, 
			   color=treatment))+
	geom_point()+
	geom_smooth(method="lm")

all.dat3 %>%
	ggplot(aes(x=pdw, 
			   y=prop_herb))+
	geom_point()+
	geom_smooth(method="lm")


all.dat3 %>%
	ggplot(aes(x=treatment, 
			   y=prop_herb))+
	geom_boxplot()+
	geom_point()

m4<-betareg(prop_herb_Mature1~prop_dw_Mature*treat, dat=all.dat2, na.action = "na.fail")
shapiro.test(resid(m4)) #normal!!
Anova(m4)  #no sig

m5<-betareg(prop_herb_Young1~prop_dw_Young*treat, dat=all.dat2, na.action = "na.fail")
#won't run

m6<-betareg(prop_herb_Young1~prop_dw_Young+treat, dat=all.dat2, na.action = "na.fail")
shapiro.test(resid(m6)) #normal!!
Anova(m6)  #no sig

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
plot(prop_dw ~ treat, data=all.dat3)
plot(prop_herb_YminusM ~ treat, data=all.dat2)

plot(prop_herb_Young ~ treat, data=all.dat2)
plot(prop_herb_Mature ~ treat, data=all.dat2)

plot(prop_herb~treatment, data=all.dat3)

m1z <- lm(prop_dw_YminusM ~ treat, data=all.dat2)
summary(m1z)
anova(m1z)

m2z <- lm(prop_herb_YminusM ~ treat, data=all.dat2)
summary(m2z)
anova(m2)z


all.dat %>%
	ggplot(aes(stage,pdw, color=treat)) +
	geom_point(aes(fill=treat),size=3) +
	geom_line(aes(group = chamber))


all.dat %>%
	ggplot(aes(stage,percent_herbivory, color=treat)) +
	geom_point(aes(fill=treat),size=3) +
	geom_line(aes(group = chamber))


#re-ordering factor levels so lm will compare everything to control
all.dat3$treatment <- factor(all.dat3$treatment, levels=c("control chamber", "CO2", "T°C", "T°C + CO2" ))

#Phenolics
summary(betareg(prop_dw~treatment, data = all.dat3))
Anova(betareg(prop_dw~treatment, data = all.dat3))
#No effect of treatment on leaf chemistry 

summary(betareg(prop_dw~treat*stage, data = all.dat))
Anova(betareg(prop_dw~treat*stage, data = all.dat))
#significant effect of leaf age, no effected by climate change

summary(betareg(prop_dw_Mature~treat, data = all.dat2))
Anova(betareg(prop_dw_Mature~treat, data = all.dat2))
#No effect of treatment on mature leaf chemistry 

summary(betareg(prop_dw_Young~treat, data = all.dat2))
Anova(betareg(prop_dw_Young~treat, data = all.dat2))
#No effect of treatment on mature young chemistry 

#Growth
summary(betareg(prop_gro~treatment, data = all.dat3))
Anova(betareg(prop_gro~treatment, data = all.dat3))

#Herbivory
summary(betareg(prop_herb~treatment, data = all.dat3))
Anova(betareg(prop_herb~treatment, data = all.dat3))

summary(betareg(prop_herb_Young1~treat, data = all.dat2))
Anova(betareg(prop_herb_Young1~treat, data = all.dat2))
#No effect of treatment on young leaf herbivory

summary(betareg(prop_herb_Mature~treat, data = all.dat2))
Anova(betareg(prop_herb_Mature~treat, data = all.dat2))
#No effect of treatment on mature leaf herbivory

#Growth+defense
summary(betareg(prop_gro~treatment+pdw, data = all.dat3))
Anova((betareg(prop_gro~treatment+pdw, data = all.dat3)))

#Growth*defense
summary(betareg(prop_gro~treatment*pdw, data = all.dat3))
Anova((betareg(prop_gro~treatment*pdw, data = all.dat3)))
#sig interaction

#Herbivory+defense
summary(betareg(prop_herb~treatment+pdw, data = all.dat3))
Anova((betareg(prop_herb~treatment+pdw, data = all.dat3)))

#Herbivory*defense
summary(betareg(prop_herb~treatment*pdw, data = all.dat3))
Anova((betareg(prop_herb~treatment*pdw, data = all.dat3)))
#sig interaction

all.dat$prop_herb1<-all.dat$prop_herb + 0.00001
#Herbivory*defense, all data
summary(betareg(prop_herb1~treat*pdw, data = all.dat))
Anova((betareg(prop_herb1~treat*pdw, data = all.dat)))
#sig interaction

all.dat %>%
	ggplot(aes(x=pdw, 
			   y=prop_herb1,
			   color=treat))+
	geom_point()+
	geom_smooth(method="lm")
