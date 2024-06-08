
library(GLMMadaptive)
require(lme4)
require(lmerTest)

#set working directory to "analysis_data" folder
#setwd("~/analysis_data")


#format data
d= read.csv( "analysis_ready_7.csv")
d$name = as.factor(d$name)
d$pressure = d$pressure_min 
d$locGrp = as.factor(d$locGrp)
d$per = NA
for(i in 1:length(d$per)){
  if(d$yr[i] < 1980){
    d$per[i] <-1
  } else if (d$yr[i] < 1985){
    d$per[i] <-2
  } else if (d$yr[i] < 1990){
    d$per[i] <-3
  } else if (d$yr[i] < 1995){
    d$per[i] <-4
  } else if (d$yr[i] < 2000){
    d$per[i] <-5
  } else if (d$yr[i] < 2005){
    d$per[i] <-6
  } else if (d$yr[i] < 2010){
    d$per[i] <-7
  } else if (d$yr[i] < 2015){
    d$per[i] <-8
  } else {
    d$per[i] <-9
  }
}    

d$yr = as.factor(d$yr)

vars = gsub('[[:digit:]]+', '', d$name)

d_fs<-d

#Calculate fs' and add to original d
for(i in unique(d$locGrp)){
  for(q in unique(d$per)){
    
    d$failTF_sum[which(d$locGrp==i & d$per==q)] = round(d$failTF_sum[which(d$locGrp==i & d$per==q)] - mean(d$failTF_sum[which(d$locGrp==i & d$per==q & vars=="NoStorm")], na.rm=FALSE)/3, digits = 0)
  }}

d$failTF_sum[d$failTF_sum<0]<-0
d <- d[!(vars=="NoStorm"),]
d <- d[complete.cases(d),]


############################################################################################################################
#
#                         Principle Component Analysis
#
############################################################################################################################

#Normalize hurricane variables
#dn is a dataframe of hurricane data
dn <- scale(d[,c("wind_max","pressure_min","minlat_min")])

#name rownames by storm ID
rownames(dn)<-paste(d$name, d$locGrp, sep="_")

#Calculate principal components of three hurricane variables
pca <- princomp(dn)

#add first PC to analysis ready dataset, variable name HurComp
d$HurComp = pca$scores[,1]
d<-d[complete.cases(d),]


############################################################################################################################
#
#                     Binomial Mixed-Effects Model
#
############################################################################################################################
#initiate dataframe for binomial (no failures, one or more failures) dataset
d0 = d  
d0$failTF_sum[d0$failTF_sum > 0]<-1
d0$failTF_sum<-as.factor(d0$failTF_sum)

#Train models
mb1 <- glmer(failTF_sum ~ 1 + (1 | yr) + (1|locGrp) , data = d0, family = binomial(link = "logit"))
mb2 <- glmer(failTF_sum ~  HurComp + (1 | yr) + (1|locGrp), data = d0, family = binomial(link = "logit"))
#I have a feeling this last model will help us on review
#mb3 <- glmer(failTF_sum ~  HurComp + (1 | yr) + (HurComp|locGrp), data = d0, family = binomial(link = "logit"))

#Summarize models
summary(mb1)
summary(mb2)
#See comment above
#summary(mb3)

#Intepreting the log odds ratio
plogis(fixef(mb2)[2])
#a 1 sd increase in HurComp per TC is associated with a 0.5880979*100 percent increase in likelihood of at least one failure 

############################################################################################################################
#
#                     Poisson Mixed-Effects Model
#
############################################################################################################################
m1 <- glmer(failTF_sum ~ 1 + (1 | yr) + (1|locGrp) , data = d, family = poisson(link = "log"))
m2 <- glmer(failTF_sum ~  HurComp + (1 | yr) + (1|locGrp), data = d, family = poisson(link = "log"))

#look at model residuals
plot(exp(predict(m2))~d$failTF_sum, xlim=c(0, 25), ylim=c(0, 25), xlab="observed fs'", ylab="Poisson predicted fs'")
abline(0,1)

summary(m1)
summary(m2)

##############################################################################################################################
#
#                   Save models and data derivatives to model_outputs folder
#
##############################################################################################################################
save(d_fs,d,d0,dn,m1,m2,mb1,mb2,pca,vars, file="../model_outputs/model_outputs.RData")
