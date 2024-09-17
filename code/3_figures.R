############################################################################################################################
#           TODO
#           *make this pretty (maybe an R Markdown?)
#           *much finishing of figures, see "TODOs" throughout
#           *generate Figure 2
#           *description of each dataset loaded in .RData
#           *script save all formatted files to figures folder
#############################################################################################################################
library(GLMMadaptive)
require(lme4)
require(lmerTest)
library(lattice)
library(sjPlot)
library(factoextra)
library('corrr')
library(ggcorrplot)
library("FactoMineR")
library(gdata)
library(data.table)
library(multcompView)
library(RColorBrewer)

#Set working directory to model_outputs
#setwd("./model_outputs")

#load data generated in step 2 (./code/2_modelTrain.R)
load("model_outputs.RData")


############################################################################################################################
#
#                                       Figure S1
#
############################################################################################################################
#Figure S1a
custom_breaks <- c(0,1, 6, 11, 16,21,26,31,36,41,46,51,56,61,66,71,76,81,86,91,96,101,106,111,116,121,126)
hist(d_fs$failTF_sum, breaks=custom_breaks, freq=TRUE, xlab="fs", main="", ylim=c(0,300))

#Figure S1b
hist(d$failTF_sum, breaks=custom_breaks, freq=TRUE, xlab="fs'", main="", ylim=c(0,300))

############################################################################################################################
#
#                                       Figure S2
#
############################################################################################################################

corr_matrix <- cor(dn)
ggcorrplot(corr_matrix)


############################################################################################################################
#
#                                       Figure S3
#
############################################################################################################################

#Use the "label" and "select.ind" arguments to select only larger storms or high impact storms to label
fviz_pca_biplot(pca, 
                fill.ind=d$failTF_sum,
                gradient.cols=c("white",brewer.pal(9, "YlOrRd")),
                ggtheme = theme_minimal(),
                pointsize=3,
                pointshape=21,
                label=FALSE,
                repel = TRUE)


############################################################################################################################
#
#                                       Figure S8
#
############################################################################################################################


#TODO change axis labels, add connecting letter report to figure if there's time
ggplot(d0,aes(locGrp,HurComp))+geom_boxplot(aes(fill=failTF_sum))+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Figure S8 Tukey's HSD test
anova <- aov(HurComp~failTF_sum*locGrp, data = d0)
tukey <- TukeyHSD(anova)
(cld <- multcompLetters4(anova, tukey))


############################################################################################################################
#
#                                       Figure 1
#
############################################################################################################################

#Figure 1, table component (TODO: format these model summaries into a table, add them as a panel to figure 2, you don't need to use R)
summary(mb1)
summary(mb2)


#Figure 1
#year caterpillar plot
x <- ranef(mb1,condVar=TRUE)[[1]]
pv   <- attr(x, "postVar")
cols <- 1:(dim(pv)[1])
se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
ord  <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))

pDf  <- data.frame(y=unlist(x)[ord],
                   ci=1.96*se[ord],
                   nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                   ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                   ind=gl(ncol(x), nrow(x), labels=names(x)), 
                   mod=rep("Null", nrow(x)))


x <- ranef(mb2,condVar=TRUE)[[1]]
pv   <- attr(x, "postVar")
cols <- 1:(dim(pv)[1])
se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))

pDf2  <- data.frame(y=unlist(x)[ord],
                    ci=1.96*se[ord],
                    nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                    ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                    ind=gl(ncol(x), nrow(x), labels=names(x)),
                    mod=rep("Poisson", nrow(x)))


test<-rbindlist(list(pDf, pDf2))[order(ID)]
test$grp <-paste(test$ID, test$mod, sep="_")
test$grp <- factor(test$grp, levels=test$grp)



shapel = rep(c(1, 16), 46)
sizel = rep(c(4.5, 2.5),46)
p <- ggplot(test, aes(ID,y, group=grp)) + coord_flip()

p <- p + geom_hline(yintercept=0)
p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0.5, colour=rep(c("darkgray", "black"), 46))
p <- p + geom_point(aes(colour = ID), shape=shapel, colour=rep(topo.colors(46)[ord], each=2), size=sizel)
p <- p +  scale_fill_continuous(guide = guide_colourbar())
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "Year",
       y = "Intercept",
       color = "Legend") +
  scale_color_manual(values = terrain.colors(11), breaks= seq(1975,2021, by=5))


#Figure 1, location caterpillar plot
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")


x <- ranef(mb1,condVar=TRUE)[[2]]
pv   <- attr(x, "postVar")
cols <- 1:(dim(pv)[1])
se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
ord  <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))

pDf  <- data.frame(y=unlist(x)[ord],
                   ci=1.96*se[ord],
                   nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                   ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                   ind=gl(ncol(x), nrow(x), labels=names(x)), 
                   mod=rep("Null", nrow(x)))


x <- ranef(mb2,condVar=TRUE)[[2]]
pv   <- attr(x, "postVar")
cols <- 1:(dim(pv)[1])
se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))

pDf2  <- data.frame(y=unlist(x)[ord],
                    ci=1.96*se[ord],
                    nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                    ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                    ind=gl(ncol(x), nrow(x), labels=names(x)),
                    mod=rep("Poisson", nrow(x)))


test<-rbindlist(list(pDf, pDf2))[order(ID)]
test$grp <-paste(test$ID, test$mod, sep="_")
test$grp <- factor(test$grp, levels=test$grp)



shapel = rep(c(1, 16), 6)
sizel = rep(c(4.5, 2.5),6)

p <- ggplot(test, aes(ID,y)) + coord_flip()
p <- p + theme(legend.position="none")
p <- p + geom_hline(yintercept=0)
p <- p + geom_point(size=sizel, shape=shapel, 
                    colour=rep(cmap[ord], each=2)) 
p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0.15, colour=rep(c("gray","black"), 6))#rep(terrain.colors(46)[rep(ord, each=2)], each=2))

p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))                      #  scale_fill_gradientn(



############################################################################################################################
#
#                                       Figure 2
#
############################################################################################################################

#Figure 2, table component (TODO: format these model summaries into a table, add them as a panel to figure 2, you don't need to use R)
summary(m1)
summary(m2)

#Figure 2 year caterpillar plot
x <- ranef(m1,condVar=TRUE)[[1]]
pv   <- attr(x, "postVar")
cols <- 1:(dim(pv)[1])
se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
ord  <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))

pDf  <- data.frame(y=unlist(x)[ord],
                   ci=1.96*se[ord],
                   nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                   ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                   ind=gl(ncol(x), nrow(x), labels=names(x)), 
                   mod=rep("Null", nrow(x)))


x <- ranef(m2,condVar=TRUE)[[1]]
pv   <- attr(x, "postVar")
cols <- 1:(dim(pv)[1])
se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))

pDf2  <- data.frame(y=unlist(x)[ord],
                    ci=1.96*se[ord],
                    nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                    ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                    ind=gl(ncol(x), nrow(x), labels=names(x)),
                    mod=rep("Poisson", nrow(x)))


test<-rbindlist(list(pDf, pDf2))[order(ID)]
test$grp <-paste(test$ID, test$mod, sep="_")
test$grp <- factor(test$grp, levels=test$grp)



shapel = rep(c(1, 16), 46)
sizel = rep(c(4.5, 2.5),46)
p <- ggplot(test, aes(ID,y, group=grp)) + coord_flip()

p <- p + geom_hline(yintercept=0)
p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0.5, colour=rep(c("darkgray", "black"), 46))
p <- p + geom_point(aes(colour = ID), shape=shapel, colour=rep(topo.colors(46)[ord], each=2), size=sizel)
p <- p +  scale_fill_continuous(guide = guide_colourbar())
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "Year",
       y = "Intercept",
       color = "Legend") +
  scale_color_manual(values = terrain.colors(11), breaks= seq(1975,2021, by=5))



#Figure 2, location catepillar plot
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")


x <- ranef(m1,condVar=TRUE)[[2]]
pv   <- attr(x, "postVar")
cols <- 1:(dim(pv)[1])
se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
ord  <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))

pDf  <- data.frame(y=unlist(x)[ord],
                   ci=1.96*se[ord],
                   nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                   ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                   ind=gl(ncol(x), nrow(x), labels=names(x)), 
                   mod=rep("Null", nrow(x)))


x <- ranef(m2,condVar=TRUE)[[2]]
pv   <- attr(x, "postVar")
cols <- 1:(dim(pv)[1])
se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))

pDf2  <- data.frame(y=unlist(x)[ord],
                    ci=1.96*se[ord],
                    nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                    ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                    ind=gl(ncol(x), nrow(x), labels=names(x)),
                    mod=rep("Poisson", nrow(x)))


test<-rbindlist(list(pDf, pDf2))[order(ID)]
test$grp <-paste(test$ID, test$mod, sep="_")
test$grp <- factor(test$grp, levels=test$grp)



shapel = rep(c(1, 16), 6)
sizel = rep(c(4.5, 2.5),6)

p <- ggplot(test, aes(ID,y)) + coord_flip()
p <- p + theme(legend.position="none")
p <- p + geom_hline(yintercept=0)
p <- p + geom_point(size=sizel, shape=shapel, 
                    colour=rep(cmap[ord], each=2)) 
p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0.15, colour=rep(c("gray","black"), 6))#rep(terrain.colors(46)[rep(ord, each=2)], each=2))

p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))                      #  scale_fill_gradientn(


############################################################################################################################
#
#                                       Figure 3
#
############################################################################################################################

#Figure 3 component
fviz_pca_var(pca, col.var = "cos2",
             #gradient.cols = c("orange", "green", "blue"),
             repel = FALSE)



#Figure 3 component
hc <- as.data.frame(list(yr=unique(d$yr), hc=NA))
hc$yr <-as.numeric(as.character(hc$yr))
hc <-hc[order(hc$yr),]
for(i in 1:length(hc$yr)){
  hc$hc[i]<-max(d$HurComp[d$yr %in% hc$yr[i]])
}
p <- ggplot(hc, aes(yr, hc)) + geom_line()
p + geom_smooth(method = "lm")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))# + 

#Figure 3 component TODO: add these model stats to plot just like previous figures
summary(lm(hc~yr, data=hc))                         


############################################################################################################################
#
#                                       Figure 4
#
############################################################################################################################

library(GLMMadaptive)
require(lme4)
require(lmerTest)

#set working directory to "analysis_data" folder
setwd("~/analysis_data")


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

############################################################################################################################
#
#                     Predictions 
#
############################################################################################################################

# Historical statistics 
historical_mean_windspeed <- 53.704824202780046
historical_sd_windspeed <- 26.272791647341982
historical_mean_pressure <- 991.9623875715454
historical_sd_pressure <- 19.597726204224276
historical_mean_latitude <- 24.782739165985284
historical_sd_latitude <- 8.547811621083

# Future predictions for 2050 
future_windspeed_2050 <- 57.879812176803924
future_pressure_2050 <- 987.3967669704609
future_latitude_2050 <- 22.28750042690649

# Standardize future predictions 
wind_max <- (future_windspeed_2050 - historical_mean_windspeed) / historical_sd_windspeed
pressure_min <- (future_pressure_2050 - historical_mean_pressure) / historical_sd_pressure
minlat_min <- (future_latitude_2050 - historical_mean_latitude) / historical_sd_latitude

# Transform future predictions into pc space
# https://cran.r-project.org/web/packages/LearnPCA/vignettes/Vig_04_Scores_Loadings.pdf
# https://cran.r-project.org/web/packages/pls/vignettes/pls-manual.pdf
# Section 8
# Combine into a single HurComp value 
future_hurcomp_2050 <- predict(pca, newdata=data.frame(wind_max, pressure_min, minlat_min))[1]

# Create new data frame for prediction
newdata <- data.frame(
  HurComp = future_hurcomp_2050,
  yr = factor(rep(NA, 6)),  
  locGrp = factor(0:5)    
)



# Ensure locGrp in newdata matches levels in original model
levels(newdata$locGrp) <- levels(d$locGrp)  

# Make predictions using the binomial model
binom_predictions <- predict(mb2, newdata = newdata,  type = "response", allow.new.levels = TRUE)

# Make predictions using the Poisson model
poisson_predictions <- exp(predict(m2, newdata = newdata,  type = "response", allow.new.levels = TRUE))

# Print the predictions
cat("Binomial model predictions (probability of failure):\n")
print(binom_predictions)

cat("\nPoisson model predictions (expected number of failures):\n")
print(poisson_predictions)

# Plotting the predictions
plot(exp(predict(m2, newdata = newdata, allow.new.levels = TRUE)), type = "p", col = "blue", pch = 19,
     xlab = "Location Group", ylab = "Predicted Failures",
     main = "Predicted Pipeline Failures by Location Group (2050)")


#explore predictions
plot(exp(predict(m2, data.frame(HurComp=d$HurComp, locGrp=0, yr=2014),
                 na.action=na.exclude))~((d$wind_max - historical_mean_windspeed) / historical_sd_windspeed), 
     #xlim=c(-2.5, 5.5), 
     ylim=c(0, 25), 
     ylab="Poisson predicted fs'", xlab="Windspeed (m/s)")
points(exp(predict(m2, data.frame(HurComp=d$HurComp, locGrp=1, yr=2014),na.action=na.exclude))~d$HurComp, xlim=c(-2.5, 5.5), ylim=c(0, 25), col="blue")
points(exp(predict(m2, data.frame(HurComp=d$HurComp, locGrp=2, yr=2014),na.action=na.exclude))~d$HurComp, xlim=c(-2.5, 5.5), ylim=c(0, 25), col="red")
points(exp(predict(m2, data.frame(HurComp=d$HurComp, locGrp=3, yr=2014),na.action=na.exclude))~d$HurComp, xlim=c(-2.5, 5.5), ylim=c(0, 25), col="green")
points(exp(predict(m2, data.frame(HurComp=d$HurComp, locGrp=4, yr=2014),na.action=na.exclude))~d$HurComp, xlim=c(-2.5, 5.5), ylim=c(0, 25), col="yellow")
points(exp(predict(m2, data.frame(HurComp=d$HurComp, locGrp=5, yr=2014),na.action=na.exclude))~d$HurComp, xlim=c(-2.5, 5.5), ylim=c(0, 25), col="purple")

abline(v= ((33 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= ((63 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= ((82 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= ((95 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= ((112 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= ((136 - historical_mean_windspeed) / historical_sd_windspeed))
