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
ggplot(phen_ag2, aes(stage, pdw, color=stage))+
	geom_boxplot(outlier.shape = NA)+
	geom_jitter(position=position_jitter(width =0.04))+
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


