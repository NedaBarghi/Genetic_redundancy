library(lme4)
library(car)
library(nortest)
library(lsmeans)

#read the file containing metabolic rate data
metabol <- read.table('metabolism.txt',header=TRUE, sep='\t')

head(metabol)
str(metabol)

#change Replicate of Base to 0 so not be confused with replicate 1 of hot-evolved population
metabol$rep_id <- ifelse(metabol$Pop == 'B', 0, metabol$Rep)

metabol$Rep <- as.factor(metabol$Rep)
metabol$rep_id <- as.factor(metabol$rep_id)

##########################################################################################
#A. test for significance of metabolic rate  between ancestral and hot-evolved populations 
##########################################################################################

#####################################
#set contrast matrix of fixed factors
#####################################

metabol$Pop <- factor(metabol$Pop, levels = c( "H","B")) 
contrasts(metabol$Pop) = contr.sum(2)

metabol$Sex <- factor(metabol$Sex, levels = c( "F","M"))
contrasts(metabol$Sex) = contr.sum(2)

metabol$rep_id <- factor(metabol$rep_id, levels = c( "1","2","3","4","5","6","7","8","9","10","0"))
contrasts(metabol$rep_id) = contr.sum(11) 

############################
#fit the mixed effect model 
###########################

met_lmer <- lmer(formula = CO2PerMg ~ Pop + (1|rep_id) + Sex + Sex:Pop, data = metabol, REML = FALSE,contrasts = 'contr.sum')
summary(met_lmer)# (1|rep_id) doesn't explain much of the variance
plot(met_lmer,which=1)

#to check whether (1|rep_id) is significant
met_lmer1 <- lmer(formula = CO2PerMg ~ Pop + (1|rep_id) + Sex + Sex:Pop, data = metabol, REML = FALSE,contrasts = 'contr.sum')
met_lmer2 <- lm(formula = CO2PerMg ~ Pop + Sex + Sex:Pop, data = metabol,contrasts = 'contr.sum')
anova(met_lmer1,met_lmer2) #p-value= 1 so it's better to use lm

#############
#final model
#############

met_lmer <- lm(formula = CO2PerMg ~ Pop + Sex + Sex:Pop, data = metabol,contrasts = 'contr.sum')
model.matrix(met_lmer)
summary(met_lmer)
Anova(met_lmer, type = 3)
plot(met_lmer,which=1)#variance is homogeneous
plot(met_lmer,which=2)#the residuals have normal distribution

#test for normality of residuals
ad.test(residuals(met_lmer)) #the residuals have normal distribution

##############
#get p-values
##############

(lsmeans <- lsmeans(met_lmer,"Pop", by = 'Sex'))
contrast(lsmeans, "pairwise")
cld(lsmeans)

(lsmeans <- lsmeans(met_lmer, pairwise ~ Pop | Sex)) 
contrast(lsmeans, "pairwise")
cld(lsmeans)

##########################################################################
#B. test for significance of metabolic rate between hot-evolved replicates 
##########################################################################

#####################################
#set contrast matrix of fixed factors
#####################################

metabol_H <- subset(metabol, metabol$Pop == 'H')

metabol_H$Rep <- factor(metabol_H$Rep, levels = c( "1","2","3","4","5","6","7","8","9","10"))
contrasts(metabol_H$Rep) = contr.sum(10) 

metabol_H$Sex <- factor(metabol_H$Sex, levels = c( "F","M"))
contrasts(metabol_H$Sex) = contr.sum(2) 

#####################
##fit the mixed model
#####################

#to check whether  (1|Rep) is significant
metH_lmer1 <- lmer(formula = CO2PerMg ~ Rep + (1|Rep) + Sex + Sex:Rep, data = metabol_H, REML = FALSE,contrasts = 'contr.sum')
metH_lmer2 <- lm(formula = CO2PerMg ~ Rep + Sex + Sex:Rep, data = metabol_H,contrasts = 'contr.sum')
anova(metH_lmer1,metH_lmer2) #p-value= 1 so it's better to use lm

#############
#final model
#############

metH_lmer <- lm(formula = CO2PerMg ~ Rep + Sex + Sex:Rep, data = metabol_H,contrasts = 'contr.sum')
plot(metH_lmer,which=1)#variance is homogeneous
plot(metH_lmer,which=2)#the residuals have normal distribution
summary(metH_lmer)
Anova(metH_lmer, type = 3) 
model.matrix(metH_lmer)

#test for normality of residuals
ad.test(residuals(metH_lmer)) #the residuals have normal distribution

#############################################
#Tukey’s HSD to correct for multiple testing
#############################################

(lsmeans <- lsmeans(metH_lmer,"Rep", by = 'Sex'))
contrast(lsmeans, "pairwise")
cld(lsmeans)


