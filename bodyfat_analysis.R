library(lme4)
library(car)
library(nortest)
library(lsmeans)

#read the file containing  fat content data
lipid=read.table("bodyFatContent.txt",sep="\t",h=TRUE,dec=",")

head(lipid)
str(lipid)

#change Replicate of Base to 0 so not be confused with replicate 1 of hot-evolved Pop
lipid$rep_id <- ifelse(lipid$Pop == 'B', 0, lipid$Rep)

lipid$Rep <- as.factor(lipid$Rep)
lipid$Wet_Weight <- as.numeric(as.character(lipid$Wet_Weight))
lipid$Lipid <- as.numeric(as.character(lipid$Lipid))
lipid$LperFly <- as.numeric(as.character(lipid$LperFly))
lipid$rep_id <- as.factor(lipid$rep_id)

#######################################################################################
#A. test for significance of fat content  between ancestral and hot-evolved populations 
#######################################################################################

#####################################
#set contrast matrix of fixed factors
#####################################

lipid$Pop <- factor(lipid$Pop , levels = c( "H","B")) 
contrasts(lipid$Pop) = contr.sum(2)

lipid$Sex <- factor(lipid$Sex , levels = c( "F","M")) 
contrasts(lipid$Sex) = contr.sum(2)

lipid$rep_id <- factor(lipid$rep_id, levels = c( "1","2","3","4","5","6","7","8","9","10","0"))
contrasts(lipid$rep_id) = contr.sum(11) 

#####################
##fit the mixed model
#####################

lipid_lmer <- lmer(formula = log10(LperFly) ~ Pop + (1|rep_id) + Sex + Pop:Sex + Wet_Weight,data = lipid, REML = FALSE,contrasts = 'contr.sum')
summary(lipid_lmer) # (1|rep_id) doesn't explain much os the variance

#to check whether (1|rep_id) is significant
lipid_lmer1 <- lmer(formula = log10(LperFly) ~ Pop + Sex + Pop:Sex + (1|rep_id) + Wet_Weight,data = lipid, REML = FALSE,contrasts = 'contr.sum')
lipid_lmer2 <- lm(formula = log10(LperFly) ~ Pop + Sex + Pop:Sex+ Wet_Weight,data = lipid ,contrasts = 'contr.sum')
anova(lipid_lmer1,lipid_lmer2) #p-value=0.4853 so it's better to use lm

lipid_lmer <- lm(formula = log10(LperFly) ~ Pop + Sex + Pop:Sex + Wet_Weight, data = lipid, contrasts = 'contr.sum')
summary(lipid_lmer) #Wet_Weight is not significant so will be dropped from the model

#############
#final model
#############

lipid_lmer <- lm(formula = log10(LperFly) ~ Pop + Sex + Pop:Sex  ,data = lipid ,contrasts = 'contr.sum')
model.matrix(lipid_lmer)
plot(lipid_lmer,which=1)#variance is homogeneous
plot(lipid_lmer,which=2)#the residuals have normal distribution
summary(lipid_lmer) 
Anova(lipid_lmer, type = 3)

#check for the homogeneity of variance
var.test(lipid$LperFly[lipid$Pop == 'B'],lipid$LperFly[lipid$Pop == 'H'])
#variance is homogeneous

#test for normality of residuals
ad.test(residuals(lipid_lmer))#the residuals have normal distribution

##############
#get p-values
##############
(lsmeans <- lsmeans(lipid_lmer,"Pop", transform = 'response', by = 'Sex'))
contrast(lsmeans, "pairwise", transform = 'response')
cld(lsmeans) 

(lsmeans <- lsmeans(lipid_lmer, pairwise ~ Pop | Sex, transform = 'response'))
contrast(lsmeans, "pairwise", transform = 'response')
cld(lsmeans) 

##########################################################################
#B. test for significance of fat content between hot-evolved replicates 
##########################################################################

#####################################
#set contrast matrix of fixed factors
#####################################

lipid_H <- subset(lipid, lipid$Pop == 'H')

lipid_H$Rep <- factor(lipid_H$Rep, levels = c( "1","2","3","4","5","6","7","8","9","10"))
lipid_H$Sex <- factor(lipid_H$Sex , levels = c( "F","M")) 

contrasts(lipid_H$Rep) = contr.sum(10) 
contrasts(lipid_H$Sex) = contr.sum(2) 

#####################
##fit the mixed model
#####################

lipid_lmer <- lmer(formula = log10(LperFly) ~ Rep + (1|Rep) + Sex + Rep:Sex + Wet_Weight,data = lipid_H, REML = FALSE,contrasts = 'contr.sum')
summary(lipid_lmer) # (1|rep_id) doesn't explain much of the variance

#to check whether (1|rep_id) is significant
lipid_lmer1 <- lmer(formula = log10(LperFly) ~ Rep + Sex + Rep:Sex + (1|Rep) + Wet_Weight,data = lipid_H, REML = FALSE,contrasts = 'contr.sum')
lipid_lmer2 <- lm(formula = log10(LperFly) ~ Rep + Sex + Rep:Sex+ Wet_Weight,data = lipid_H ,contrasts = 'contr.sum')
anova(lipid_lmer1,lipid_lmer2) #p-value=1 so it's better to use lm

lipid_lmer <- lm(formula = log10(LperFly) ~ Rep + Sex + Rep:Sex + Wet_Weight, data = lipid_H, contrasts = 'contr.sum')
summary(lipid_lmer) 
Anova(lipid_lmer, type = 3) #Wet_Weight and Rep:Sex are not significant so will be dropped from the model

#############
#final model
#############

lipidH_lmer <- lm(formula = log10(LperFly) ~ Rep + Sex ,data = lipid_H,contrasts = 'contr.sum')
summary(lipidH_lmer)
Anova(lipidH_lmer, type = 3) 
model.matrix(lipidH_lmer)
plot(lipidH_lmer,which=1)#variance is homogeneous
plot(lipidH_lmer,which=2)#the residuals have normal distribution

#test for normality of residuals
ad.test(residuals(lipidH_lmer))#the residuals have normal distribution

#############################################
#Tukey’s HSD to correct for multiple testing
#############################################

(lsmeans <- lsmeans(lipidH_lmer,"Rep", transform = 'response',by = 'Sex'))
contrast(lsmeans, "pairwise")
cld(lsmeans)

(lsmeans <- lsmeans(lipidH_lmer, pairwise ~ Rep | Sex, transform = 'response'))
cld(lsmeans)
