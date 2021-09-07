#Interactive effects of climate change, leaf age, and secondary metabolites 
#on plant growth, defense, and herbivory.

#LOADING LIBRARIES----------------------------------------------------
{library(lme4)
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
}

#LOADING & WRANGLING DATA----------------------------------------------------

##GROWTH

#load growth data
grow <- read.csv(file="Piper_growth.csv",head=TRUE)
table(grow$Treatment)   #SRW: Why are there 6 T chambers and 4 T+CO2???
colnames(grow)[1] <- "Casa"

grow$Treatment[which(grow$Treatment=="TÂ°C")] <- "T°C"
grow$Treatment[which(grow$Treatment=="TÂ°C + CO2")] <- "T°C + CO2"

grow$total_gro <- grow$ht.2018.09.cm-grow$ht.2018.04.cm
grow$rel_gro<-grow$total_gro/grow$ht.2018.04.cm
grow$per_gro<-grow$rel_gro*100

grow <- grow[order(grow$Treatment),]
grow<-grow[-c(11:15),] #removing control (no chamber)

hist(grow$rel_gro)  #SW: not normal, this is a typical proportion
shapiro.test(grow$rel_gro)#LM: it is technically normal
hist(grow$total_gro)
shapiro.test(grow$total_gro)#LM: but so is this one. Both yield the same results, no treatment effect

#re-ordering factor levels so lm will compare everything to control
grow$Treatment[]
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
pg$treatment[which(pg$treatment=="TÂ°C")] <- "T°C"
pg$treatment[which(pg$treatment=="TÂ°C + CO2")] <- "T°C + CO2"

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
samp1<-samp1[,-4]#remove second herb var   #SRW: getting an error here, there is no col 4

L1 <- cor(samp1)#correlation matrix
corrplot(L1, method = "circle")
corrplot(L1, method = "number")
#low colinearity

check1 <- lm(prop_gro ~ prop_herb+prop_dw+stage+treat, data=all.dat)
summary(check1)
vif(check1)



##JULY ANALYSIS----
##QUESTION 1----

#GROWTH
#all.dat, N=20, only one measurement per plant/per chamber

#levels(all.dat3$treatment) <- factor(all.dat3$treatment, levels=c("control chamber", "CO2", "T°C", "T°C + CO2"))

Anova(betareg(prop_gro~treatment, dat=all.dat3))
#p=0.82, no effect of treatment 


shapiro.test(resid(betareg(prop_gro~treatment, dat=all.dat3)))
#residuals normally distributed

#growth plot
plot(all.dat3$prop_gro~all.dat3$treatment)


##SRW: To send to Maaike, it might be good to also include some alternative plots/analyses with the natural 
#control as well. I did a quick check on this just not running that line above (L35) and using this code

grow$prop_gro <- grow$per_gro/100
Anova(betareg(prop_gro~Treatment, dat=grow))
#p=0.54, no effect of treatment 


shapiro.test(resid(betareg(prop_gro~Treatment, dat=grow)))
#residuals normally distributed

#growth plot
boxplot(grow$prop_gro~grow$Treatment)




#CHEMISTRY
#all.dat, N=40, leaves of same stage on each plant were combined for chem analysis
all.dat$treat <- factor(all.dat$treat, levels=c("control_chamber", "CO2", "TC", "TC+CO2" ))
Anova(betareg(prop_dw~treat, dat=all.dat))
#p=0.76, no effect of treatment

shapiro.test(resid(betareg(prop_dw~treat, dat=all.dat)))
#residuals normally distributed

#chemistry plot
plot(all.dat$prop_dw~all.dat$treat)


#SRW: another approach for these initial Q1 analyses--don't aggregate any given
#dataset, i.e. use all the phenolics samples you have when that is the response
#I think we discussed this and a major issue was that you can't use the random effects in the betareg. 
#However, I was just looking at this very helpful paper
#Douma and Weedon 2019: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13234
#they give some examples of beta regression with random effects using the package glmmTMB

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





#HERBIVORY
#pg1 data, N=80, all herbivory data
#pg1$treatment <- factor(pg1$treatment, levels=c("control chamber", "CO2", "T°C", "T°C + CO2" ))
pg1$treatment<-as.factor(pg1$treatment)
#adding small number to avoid error (must be between 0,1)
pg1$prop_herb1<-pg1$prop_herb+0.0001

#note in the paper I mention above, they suggest rescaling the dataset instead of just adding a constant
#From the paper, Appendix S3-------------------
#A suggested rescaling equation is:
#	x∗i=xi(n−1)+0.5n
#Where x∗i is the transformation of xi and n is the total number of observations in the dataset.
#For convenience we define this as a custom function tranform01 and apply it to the dataset:
#	transform01 <- function(x) {
#		(x * (length(x) - 1) + 0.5) / (length(x))
#	}
#andrew2$ALGAE.scaled <- transform01(andrew2$ALGAE.mean)

#I did not mess with this...

Anova(betareg(prop_herb1~treatment, dat=pg1))#p=0.6813, no effect of treatment
shapiro.test(resid(betareg(prop_herb1~treatment, dat=pg1)))
#residuals NOT normally distributed

#herbivory plot N=80
plot(pg1$prop_herb1~pg1$treatment)

#running model where leaves of each stage of combined, N=40, all.dat
all.dat$prop_herb1<-all.dat$prop_herb+0.0001
Anova(betareg(prop_herb1~treat, dat=all.dat))#p=0.772
shapiro.test(resid(betareg(prop_herb1~treat, dat=all.dat)))
#residuals still not normally distributed, but less non normal...

#herbivory plot N=40
plot(all.dat$prop_herb1~all.dat$treat)

all.dat3$prop_herb1<-all.dat3$prop_herb+0.0001
Anova(betareg(prop_herb1~treatment, dat=all.dat3))#p=0.6528
shapiro.test(resid(betareg(prop_herb1~treatment, dat=all.dat3)))
#residuals normally distributed

#herbivory plot N=20
plot(all.dat3$prop_herb1~all.dat3$treat)


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





##QUESTION 2A----

#N=40, only data choice for this analysis
Anova(betareg(prop_dw~stage*treat, dat=all.dat))
#stage signficant p<0.0001
shapiro.test(resid(betareg(prop_dw~stage*treat, dat=all.dat)))#normal

#plot
library(dplyr)
all.dat %>%
	ggplot(aes(stage,pdw, color=treat)) +
	geom_point(aes(fill=treat),size=3) +
	geom_line(aes(group = chamber))


#SRW: see above, would be great to do this with the random effect included.
#also the chem data are so normal, could also do this with lmer
hist(all.dat$prop_dw)
Anova(lmer(pdw~stage*treat + (1|chamber), dat=all.dat))
#not that it matters, I just hate ignoring the random effects



##QUESTION 2B----
#Growth-defense trade-off
#Because only one measurement for growth, N must be 20 since we aren't including random effects
#to account for repeated measurements

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

#Option dos
#Two models, where leaf age is split

#Mature leaves
Anova((betareg(growth~treat*prop_dw_Mature, data = all.dat2)))
#no significance
shapiro.test(resid(betareg(growth~treat*prop_dw_Mature, data = all.dat2)))#normal

#Young leaves
Anova((betareg(growth~treat*prop_dw_Young, data = all.dat2)))
#interaction sig, p=0.0065
shapiro.test(resid(betareg(growth~treat*prop_dw_Young, data = all.dat2)))#normal

joint_tests((betareg(growth~treat*prop_dw_Young, data = all.dat2)), by = "treat")
#significant effect of phenolics on growth in both temperature treatments
#defense decreased faster as growth increased in temp treatments

#Growth*defense plot, young leaves
all.dat2 %>%
	ggplot(aes(x=growth, 
			   y=prop_dw_Young,
			   color=treat))+
	geom_point()+
	geom_smooth(method="lm")


##SRW: I guess I like the opcion uno better for this one?? But we could say it is mostly driven by the 
#young leaves and could even show those separately in a supplement. 
#I am still struggling with how to interpret this though...basically when you have increased temp, 
#plants that do better just do better in terms of everything...growth and phenolics, but there is no "trade-off"
#i.e. a negative relationship for any of the plants or treatments.
#Maybe this just has to do with physiological differences across individuals in heat tolerance?? Where those
#that are tolerant of higher temps just do better generally

#I would also be curious to see the data for the no chamber control here...though it would probably just make things
#more confusing





##QUESTION 2C----
#Defense-herbivory tradeoff

#SRW: I guess we would not really think of this as a "trade-off" but as  relative change in the effectiveness
#of defense

#Option uno, N=20
Anova((betareg(prop_herb~treatment*pdw, data = all.dat3)))
#all are significant

shapiro.test(resid(((betareg(prop_herb~treatment*pdw, data = all.dat3)))))#normal

#Defense plot
all.dat3 %>%
	ggplot(aes(x=pdw, 
			   y=prop_herb))+
	geom_point()+
	geom_smooth(method = 'lm')


#Treatment plot
all.dat3 %>%
	ggplot(aes(x=treatment, 
			   y=prop_herb))+
	geom_boxplot()+
	geom_point()

#Defense*treatment plot
all.dat3 %>%
	ggplot(aes(x=pdw, 
			   y=prop_herb,
			   color=treatment))+
	geom_point()+
	geom_smooth(method="lm")

#Option dos, N=40
Anova((betareg(prop_herb1~treat*pdw, data = all.dat)))
#only chemistry signficant, p<0.0001
#this doesn't account/control for leaf age
#Chem plot
all.dat %>%
	ggplot(aes(x=pdw, 
			   y=prop_herb1))+
	geom_point(aes(color=stage))+
	geom_smooth(method = 'lm')


#SRW: I think we need to account for leaf age in whatever we do, so this one doesn't really work. The pattern
#here could just be driven by age (old leaves have lower pdw and also higher herbivory)--I added color to the
#graph above to see this and it does seem to be the case


#Option 3, split leave ages
#Young leaves
all.dat2$prop_herb_Young1<-all.dat2$prop_herb_Young+0.0001
Anova((betareg(prop_herb_Young1~treat*prop_dw_Young, data = all.dat2)))
#nothing signif
shapiro.test(resid(betareg(prop_herb_Young1~treat*prop_dw_Young, data = all.dat2)))#normal

#Mature leaves
all.dat2$prop_herb_Mature1<-all.dat2$prop_herb_Mature+0.0001
Anova((betareg(prop_herb_Mature1~treat*prop_dw_Mature, data = all.dat2)))
#chemistry marginally signif
shapiro.test(resid(betareg(prop_herb_Mature1~treat*prop_dw_Mature, data = all.dat2)))#normal

#Plot, mature leaves defense-herbivory
all.dat2 %>%
	ggplot(aes(x=prop_dw_Mature, 
			   y=prop_herb_Mature1,
			   color=treat))+
	geom_point()+
	geom_smooth(method="lm")
#purple line marginally signficant



#SRW: Not sure what to make of this--we would basically have to conclude that phenolics are effective as a defense 
#against herbivores only in a climate change scenario with T+CO2??? This does not really make sense to me and I am 
#struggling because I think there are just not enough data to really conclude much...only 4 points for the 
#T+CO2 treatment where we see that negative relationship




##OLD ANALYSES/BRAIN DUMPS----


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
