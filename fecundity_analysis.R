library(lme4)
library(car)
library(nortest)
library(lsmeans)

#read the file containing fecundity data
HOT_Base_pool <- read.table('fecundity.txt',header=TRUE, sep='\t')

head(HOT_Base_pool)
str(HOT_Base_pool)

#change Replicate of Base to 0 so not be confused with replicate 1 of hot-evolved population
HOT_Base_pool$rep_id <- ifelse(HOT_Base_pool$Pop == 'B', 0, HOT_Base_pool$Rep)

#Dry_Weight is the total weight of all Female_Num flies, convert it to the average weight of each female
HOT_Base_pool$wt_perfemale <- HOT_Base_pool$Dry_Weight/HOT_Base_pool$Female_Num

HOT_Base_pool$Rep <- as.factor(HOT_Base_pool$Rep)
HOT_Base_pool$Subrep <- as.factor(HOT_Base_pool$Subrep)
HOT_Base_pool$rep_id <- as.factor(HOT_Base_pool$rep_id)

####################################################################################
#A. test for significance of fecundity between ancestral and hot-evolved populations 
####################################################################################

#####################################
#set contrast matrix of fixed factors
#####################################

HOT_Base_pool$Pop <- factor(HOT_Base_pool$Pop, levels = c("H","B")) 
contrasts(HOT_Base_pool$Pop) = contr.sum(2)

HOT_Base_pool$rep_id <- factor(HOT_Base_pool$rep_id, levels = c( "1","2","3","4","5","6","7","8","9","10",'0')) 
contrasts(HOT_Base_pool$rep_id) = contr.sum(11)

############################
#fit the mixed effect model 
###########################

fec_lmer <- lmer(log10(Egg_Num/Female_Num) ~ Pop + log10(wt_perfemale) + (1|rep_id) , REML = FALSE,data = HOT_Base_pool,contrasts = 'contr.sum')
summary(fec_lmer) # (1|rep_id) doesn't explain much of the variance
Anova(fec_lmer, type = 3)

#to check whether (1|rep_id) is significant
fec_lmer1 <- lmer(log10(Egg_Num/Female_Num) ~ Pop + log10(wt_perfemale) + (1|rep_id) , REML = FALSE,data = HOT_Base_pool,contrasts = 'contr.sum')
fec_lmer2 <- lm(log10(Egg_Num/Female_Num) ~ Pop + log10(wt_perfemale) ,data = HOT_Base_pool,contrasts = 'contr.sum')
anova(fec_lmer1,fec_lmer2) #p-value= 1 so it's better to use lm

#############
#final model
#############

fec_lmer <- lm(log10(Egg_Num/Female_Num) ~ Pop + log10(wt_perfemale) ,data = HOT_Base_pool,contrasts = 'contr.sum')

summary(fec_lmer)
Anova(fec_lmer, type = 3) 
model.matrix(fec_lmer)

#test for normality of residuals
ad.test(residuals(fec_lmer)) #the residuals have normal distribution

#check for the homogeneity of variance
var.test(HOT_Base_pool$Egg_Num[HOT_Base_pool$Pop == 'B']/HOT_Base_pool$Female_Num[HOT_Base_pool$Pop == 'B'],HOT_Base_pool$Egg_Num[HOT_Base_pool$Pop == 'H']/HOT_Base_pool$Female_Num[HOT_Base_pool$Pop == 'H'])
#variance is homogeneous

##############
#get p-values
##############

(lsmeans <- lsmeans(fec_lmer,"Pop", transform = 'response'))
contrast(lsmeans, "pairwise")
cld(lsmeans)

######################################################################
#B. test for significance of fecundity between hot-evolved replicates 
######################################################################

#####################################
#set contrast matrix of fixed factors
#####################################

HOT_pool <- subset(HOT_Base_pool,HOT_Base_pool$Pop == "H")

HOT_pool$Rep <- factor(HOT_pool$Rep, levels = c( "1","2","3","4","5","6","7","8","9","10"))
contrasts(HOT_pool$Rep) = contr.sum(10)

#####################
##fit the mixed model
#####################

fec_lmer <- lmer(log10(Egg_Num/Female_Num) ~ Rep + log10(wt_perfemale) + (1|Rep) ,data = HOT_pool, REML = FALSE, contrasts = 'contr.sum') 
summary(fec_lmer) # (1|Rep) doesn't explain much os the variance
Anova(fec_lmer, type = 3) 

#to check whether  (1|Rep) is significant
fec_lmer1 <- lmer(log10(Egg_Num/Female_Num) ~ Rep + log10(wt_perfemale) + (1|Rep) ,data = HOT_pool, REML = FALSE, contrasts = 'contr.sum') 
fec_lmer2 <- lm(log10(Egg_Num/Female_Num) ~ Rep + log10(wt_perfemale) ,data = HOT_pool, contrasts = 'contr.sum') 
anova(fec_lmer1,fec_lmer2) #p-value= 1 so it's better to use lm

#############
#final model
#############

fec_lmer <- lm(log10(Egg_Num/Female_Num) ~ Rep + log10(wt_perfemale) ,data = HOT_pool, contrasts = 'contr.sum') 
plot(fec_lmer, which = 1) #variance is homogeneous
plot(fec_lmer, which = 2) #residuals are normal
summary(fec_lmer)
Anova(fec_lmer, type = 3) 
model.matrix(fec_lmer)

#############################################
#Tukey’s HSD to correct for multiple testing
#############################################

(lsmeans <- lsmeans(fec_lmer,"Rep", transform = 'response'))
contrast(lsmeans, "pairwise")
cld(lsmeans)