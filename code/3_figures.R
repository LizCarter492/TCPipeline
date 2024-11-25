############################################################################################################################
#           TODO
#           *make this pretty (maybe an R Markdown?)
#           *much finishing of figures, see "TODOs" throughout
#           *description of each dataset loaded in .RData
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
library(emmeans)



############################################################################################################################
#
#     Figure 3: a) summary of mixed-effects poisson model parameters and fit, predicting frequency of 60 day post-TC pipeline failures above average annual failure rate given TC intensity,
#     with random intercepts on time and location, and a random slope on location, compared to null model; b) intercept shifts on location associated with prediction
#     skill of TC intensity as a parameter, with open circles representing null model intercept, closed circles representing 
#     TC intensity model intercept, and error bars representing 95% confidence intervals of intercept estimate; 3) as in b, but for 
#     time random intercept. 
############################################################################################################################
load("./model_outputs/model_outputs.RData")
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
locGrps = c("High Plains","Northeast", "MI Delta", "Southeast", "Midwest", "Texas Coast")


summary(ms1)
summary(ms2)
#a, model table
tab_model(ms1,ms2, show.r2=FALSE, show.aic=TRUE,  show.ci=FALSE, show.icc=FALSE, show.loglik=TRUE, show.est=TRUE, show.re.var=TRUE, auto.label = FALSE, 
          pred.labels=c("Intercept","TCI"),
          dv.labels=c("Null model", "TC Intensity"),
          p.style="stars")


#b, year caterpillar plot
x <- ranef(ms1,condVar=TRUE)[[1]]
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


x <- ranef(ms2,condVar=TRUE)[[1]]
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


svg(file="./figures/Figure3b.svg", width =9, height=4)
shapel = rep(c(1, 16), 46)
sizel = rep(c(4.5, 2.5),46)
p <- ggplot(test, aes(ID,y, group=grp)) #+ coord_flip()

p <- p + geom_hline(yintercept=0)
p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0.5, colour=rep(c("darkgray", "black"), 46))
p <- p + geom_point(aes(colour = ID), shape=shapel, colour=rep(topo.colors(46)[ord], each=2), size=sizel)
p <- p +  scale_fill_continuous(guide = guide_colourbar())
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "Year",
       y = "Intercept",
       color = "Legend") +
  scale_color_manual(values = terrain.colors(11), breaks= seq(1975,2021, by=5))

dev.off()

#c, location catepillar plot

x <- ranef(ms1,condVar=TRUE)[[2]]
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


x <- ranef(ms2,condVar=TRUE)[[2]][[1]]
pv   <- attr(ranef(ms2,condVar=TRUE)[[2]], "postVar")[1,1,]
cols <- 1:(length(pv))
se   <- sqrt(pv)

pDf2  <- data.frame(y=unlist(x)[ord],
                    ci=1.96*se[ord],
                    nQQ=rep(qnorm(ppoints(length(x))), 1),
                    ID=factor(rep(rownames(ranef(ms2,condVar=TRUE)[[2]]), 1)[ord], levels=rownames(ranef(ms2,condVar=TRUE)[[2]])[ord]),
                    ind=gl(1, length(x), labels=names(ranef(ms2,condVar=TRUE)[[2]])),
                    mod=rep("Poisson", length(x)))


test<-rbindlist(list(pDf, pDf2))[order(ID)]
test$grp <-paste(test$ID, test$mod, sep="_")
test$grp <- factor(test$grp, levels=test$grp)



shapel = rep(c(1, 16), 6)
sizel = rep(c(4.5, 2.5),6)

svg(file="./figures/Figure3c.svg", width =4.5, height=5)
p <- ggplot(test, aes(ID,y)) + coord_flip()
p <- p + theme(legend.position="none")
p <- p + geom_hline(yintercept=0)
p <- p + geom_point(size=sizel, shape=shapel, 
                    colour=rep(cmap[ord], each=2)) 
p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0.15, colour=rep(c("gray","black"), 6))#rep(terrain.colors(46)[rep(ord, each=2)], each=2))
p <- p + scale_x_discrete("location group", labels=locGrps[ord])
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))                      #  scale_fill_gradientn(
dev.off()

############################################################################################################################
#
#                                       Figure 4 a,b
#
############################################################################################################################
#Read in HURDAT2 raw data
h<-read.csv("./raw_data/NOAAHURDAT.csv")
h$nameyr<-paste(h$name,as.character(h$year), sep="")
h$minlat<-NA
for(i in 1:length(h$minlat)){
  h$minlat[i]<-min(h$lat[h$nameyr %in% h$nameyr[i]])
}
#Format dataframe of annual min/max values
time = seq(min(h$year), max(h$year))
cli = as.data.frame(time)
cli$ws<-rep(NA)
cli$z<-rep(NA)
cli$lat<-rep(NA)

for(i in 1:length(time)){
  cli$ws[i]<-max(h$wind[h$year==cli$time[i]], na.rm=T)
  cli$z[i]<-min(h$pressure[h$year==cli$time[i]], na.rm=T)
  cli$lat[i]<-min(h$minlat[h$year==cli$time[i]], na.rm=T)
}
#Create svg file to save image
svg(file = "./figures/Figure4ab.svg", width=5, height =6)
# set plot parameters
par(mfrow=c(3,1), mar=c(2,5,1,1))
# plot windspeed panel
# train linear regression of annual max windspeed by year
m1=lm(ws~time, data=cli)
# generate data for plotting trendline and prediction interval lines
l = as.data.frame(predict(m1, cli, interval = "prediction"))
l$time = time

# plot raw data
plot(cli$ws~cli$time, type="l", 
     ylim=c(min(l$lwr), max(l$upr)),
     ylab="annual max. windspeed \n(knots)",
     xlab="year")
#add trendline and prediction interval
points(l$fit~l$time, type="l", lwd=2)
points(l$lwr~l$time, type="l", lwd=1, lty=2)
points(l$upr~l$time, type="l", lwd=1, lty=2)
# add confidence interval
lc = as.data.frame(predict(m1, cli, interval = "confidence"))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=2, col="blue")
points(lc$upr~lc$time, type="l", lwd=1, lty=2, col="blue")
# print lm coefficient and confidence interval in legend
legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))

# plot pressure panel
# train linear regression of annual min pressure by year
m2=lm(z~time, data=cli)
l = as.data.frame(predict(m2, cli, interval = "prediction"))
l$time = time
plot(cli$z~cli$time, type="l", 
     ylim=c(min(l$lwr), max(l$upr)),
     ylab="annual min. pressure \n(kPa)",
     xlab="year")
points(l$fit~l$time, type="l", lwd=2)
points(l$lwr~l$time, type="l", lwd=1, lty=2)
points(l$upr~l$time, type="l", lwd=1, lty=2)
lc = as.data.frame(predict(m2, cli, interval = "confidence"))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=2, col="blue")
points(lc$upr~lc$time, type="l", lwd=1, lty=2, col="blue")
legend(x="bottomleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m2)$coefficients[2,1], digits=3)), "**", sep="")))

# plot latitude panel
# train linear regression of annual min latitude by year
m3=lm(lat~time, data=cli)
l = as.data.frame(predict(m3, cli, interval = "prediction"))
l$time = time
plot(cli$lat~cli$time, type="l", 
     ylim=c(min(l$lwr), max(cli$lat)),
     ylab="annual min. latitude \n(degrees)",
     xlab="year")
points(l$fit~l$time, type="l", lwd=2)
points(l$lwr~l$time, type="l", lwd=1, lty=2)
points(l$upr~l$time, type="l", lwd=1, lty=2)
lc = as.data.frame(predict(m3, cli, interval = "confidence"))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=2, col="blue")
points(lc$upr~lc$time, type="l", lwd=1, lty=2, col="blue")
legend(x="bottomleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m3)$coefficients[2,1], digits=3)), "***", sep="")))

dev.off()
hm1<-m1
hm2<-m2
hm3<-m3
save(h, cli, hm1,hm2,hm3, file="./model_outputs/hur_trends.RData")
#################################################################################
#
#                                        Fig 4c
#
################################################################################
rm(list=ls())
load("./model_outputs/hur_trends.RData")
load("./model_outputs/model_outputs.RData")
h<-read.csv("./raw_data/NOAAHURDAT.csv")
h_p<-read.csv("./analysis_data/analysis_ready_binomial_k6.csv")
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
locGrps = c("High Plains","Northeast", "MI Delta", "Southeast", "Midwest", "Texas Coast")

# Historical statistics 
u_ws <- u["wind"]
sd_ws <- sd["wind"]
u_z <- u["pressure"]
sd_z <- sd["pressure"]
u_lat<- u["minlat"]
sd_lat<- sd["minlat"]


time = c(1970, 2010, 2050)
pred = as.data.frame(time)
pred$ws<-pred$z<-pred$lat<-rep(NA)

pred$ws = predict(hm1, newdata=pred)
pred$z = predict(hm2, pred)
pred$lat = predict(hm3, pred)
pred

# Standardize predictions 
spred = pred
spred$ws = (pred$ws - u_ws)/sd_ws
spred$z = (pred$z - u_z)/sd_z
spred$lat = rep(0,3)
spred
names(spred)<-c("year", "minlat", "pressure", "wind")

# Transform predictions into pc space
PC_pred <- as.data.frame(predict(pca, newdata=spred))
PC_pred$year = spred$year
PC_pred


wind = (c(34, 63,82,95,112, 136)-u_ws)/sd_ws
pressure = (c(1020, 1000,980,965,945,920)-u_z)/sd_z
minlat = rep(0,6)
hur_cat=data.frame(wind, pressure, minlat)
PC_hurcat <-as.data.frame(predict(pca, newdata=hur_cat))

yr_s=2014
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
d <- d[order(d$HurComp),]
par(mfrow=c(1,1))
svg(file = "./figures/Figure4c.svg", width=6, height =6)

m1<-ms1
m2<-ms2
#plot data
HurComp<-seq(min(d$HurComp), 7, by = 0.05)
plot(exp(predict(ms2, data.frame(HurComp=d$HurComp[d$locGrp==0], locGrp=0, yr=yr_s),
                 na.action=na.exclude))~d$HurComp[d$locGrp==0], 
     type="l",
     lwd=4,
     col=cmap[1],
     xlim=c(-2.2, 6), 
     ylim=c(0, 35), 
     ylab="Poisson predicted fs'", xlab="TC intensity")

abline(v= PC_pred[1,1], lty =1, lwd=30, col="lightgoldenrod")#((136 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= PC_pred[2,1], lty =1, lwd=30, col="lightgoldenrod")#((136 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= PC_pred[3,1], lty =1, lwd=30, col="lightgoldenrod")#((136 - historical_mean_windspeed) / historical_sd_windspeed))

for(i in 0:5){
  p4p<-predict(ms2, data.frame(HurComp=HurComp, locGrp=i, yr=yr_s))
  se<-as.data.frame(emmeans(ms2, ~HurComp, at = list(x = HurComp, locGrp=i, yr=yr_s)))$SE
  
  polygon(x = c(HurComp, rev(HurComp)),
          y = c(exp(p4p) - 1.96*exp(se),
                rev(exp(p4p) +1.96*exp(se))),
          col =  adjustcolor(cmap[i+1], alpha.f = 0.1), border = NA)
  
  
}
points(exp(predict(ms2, data.frame(HurComp=d$HurComp[d$locGrp==1], locGrp=1, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==1],type="l", lwd=4,col=cmap[2])
points(exp(predict(ms2, data.frame(HurComp=d$HurComp[d$locGrp==2], locGrp=2, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==2],type="l", lwd=4,col=cmap[3])
points(exp(predict(ms2, data.frame(HurComp=d$HurComp[d$locGrp==3], locGrp=3, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==3],type="l", lwd=4,col=cmap[4])
points(exp(predict(ms2, data.frame(HurComp=d$HurComp[d$locGrp==4], locGrp=4, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==4],type="l", lwd=6,col=cmap[5])
points(exp(predict(ms2, data.frame(HurComp=d$HurComp[d$locGrp==5], locGrp=5, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==5],type="l", lwd=4,col=cmap[6])

points(exp(predict(ms2, data.frame(HurComp=HurComp, locGrp=0, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[1])
points(exp(predict(ms2, data.frame(HurComp=HurComp, locGrp=1, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[2])
points(exp(predict(ms2, data.frame(HurComp=HurComp, locGrp=2, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[3])
points(exp(predict(ms2, data.frame(HurComp=HurComp, locGrp=3, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[4])
points(exp(predict(ms2, data.frame(HurComp=HurComp, locGrp=4, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[5])
points(exp(predict(ms2, data.frame(HurComp=HurComp, locGrp=5, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[6])

x=data.frame(yr=c(1970,2010,2050))
for(i in 0:5){
  yr= seq(min(h_p$yr, na.rm=T), max(h_p$yr, na.rm=T))
  wind = rep(NA,length(yr))
  pressure = rep(NA,length(yr))
  for(j in 1:length(yr)){
    wind[j]<-max(h_p$wind[h_p$locGrp==i & h_p$yr==yr[j]], na.rm=T)
    pressure[j]<--1*max(h_p$pressure[h_p$locGrp==i & h_p$yr==yr[j]], na.rm=T)
  }
  wind[wind==-Inf]<-NA
  pressure[pressure==-Inf]<-NA
  df<-data.frame('wind'=wind, 'pressure'=pressure, minlat=rep(0, length(yr)))
  for(k in 1:(dim(df)[1])){
    df[k,]<-(df[k,]-u)/sd
  }
  df$minlat=rep(-3.5)
  df$HurComp <- predict(pca, df)[,1]
  
  m=lm(HurComp~yr, data=df)
  p = predict(m,newdata=x)
  print(i)
  print(exp(predict(ms2, data.frame(HurComp=p, locGrp=i, yr=yr_s),na.action=na.exclude)))
  points(exp(predict(ms2, data.frame(HurComp=p[1], locGrp=i, yr=yr_s),na.action=na.exclude))~p[1],type="p", pch = 25, cex=1.5,bg=cmap[i+1])
  points(exp(predict(ms2, data.frame(HurComp=p[2], locGrp=i, yr=yr_s),na.action=na.exclude))~p[2],type="p", pch = 23, cex=1.5,bg=cmap[i+1])
  points(exp(predict(ms2, data.frame(HurComp=p[3], locGrp=i, yr=yr_s),na.action=na.exclude))~p[3],type="p", pch = 24, cex=1.5,bg=cmap[i+1])

}

abline(v= PC_hurcat$Comp.1[1], lty =2, lwd=2, col="darkgray")#((33 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= PC_hurcat$Comp.1[2], lty =2, lwd=2, col="darkgray")#((63 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= PC_hurcat$Comp.1[3], lty =2, lwd=2, col="darkgray")#((82 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= PC_hurcat$Comp.1[4], lty =2, lwd=2, col="darkgray")#((95 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= PC_hurcat$Comp.1[5], lty =2, lwd=2, col="darkgray")#((112 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= PC_hurcat$Comp.1[6], lty =2, lwd=2, col="darkgray")#((136 - historical_mean_windspeed) / historical_sd_windspeed))



txt_y=rep(35, 7)
txt_x=c((min(HurComp)+PC_hurcat$Comp.1[1])/2,
        (PC_hurcat$Comp.1[1] + PC_hurcat$Comp.1[2])/2,
        (PC_hurcat$Comp.1[2] + PC_hurcat$Comp.1[3])/2,
        (PC_hurcat$Comp.1[3] + PC_hurcat$Comp.1[4])/2,
        (PC_hurcat$Comp.1[4] + PC_hurcat$Comp.1[5])/2,
        (PC_hurcat$Comp.1[5] + PC_hurcat$Comp.1[6])/2,
        ((PC_hurcat$Comp.1[6] + max(HurComp))/2-0.2))
labs=c(-1,0,1,2,3,4,5)
text( txt_x, txt_y,labs,
      cex = 1.5, pos = 1, col = "darkgray")

legend(x = "left", c(locGrps, "1970 ann.max", "2010 ann.max", "2050 ann.max"), col=c(cmap, rep("black", 3)), lty=c(rep(1,6), rep(NA,3)),lwd=c(rep(3,6),rep(NA,3)),pch=c(rep(NA,6),c( 25,23,24)))



dev.off()



############################################################################################################################
#
#                                       Figure S2
#
############################################################################################################################
#Figure S1a
load("./model_outputs/model_outputs.RData")
custom_breaks <- c(0,1, 6, 11, 16,21,26,31,36,41,46,51,56,61,66,71,76,81,86,91,96,101,106,111,116,121,126)
svg(file = "./figures/FigureS2.svg", width=8, height =5)
par(mfrow=c(1,2))
hist(d_fs$failTF_sum, breaks=custom_breaks, freq=TRUE, xlab="fs", main="", ylim=c(0,300))

#Figure S1b
hist(d$failTF_sum, breaks=custom_breaks, freq=TRUE, xlab="fs'", main="", ylim=c(0,300))
dev.off()
############################################################################################################################
#
#                                       Figure S3
#
############################################################################################################################
load("./model_outputs/model_outputs.RData")
corr_matrix <- cor(dn)
svg(file = "./figures/FigureS3.svg", width=3, height =3)
ggcorrplot(corr_matrix)
dev.off()

############################################################################################################################
#
#                                       Figure S4
#
############################################################################################################################
load("./model_outputs/model_outputs.RData")
#Use the "label" and "select.ind" arguments to select only larger storms or high impact storms to label
svg(file = "./figures/FigureS4.svg", width=8, height =5)
fviz_pca_var(pca, col.var = "cos2",
             repel = FALSE)
dev.off()

############################################################################################################################
#
#                                       Figure S5
#
############################################################################################################################
#Read in data
h<-read.csv("./analysis_data/analysis_ready_binomial_k6.csv")
h <- h[h$lat>0,]
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
locGrps = c("High Plains","Northeast", "MI Delta", "Southeast", "Midwest", "Texas Coast")


#Create dataframe for annual max/min values
#time = seq(min(h$year), max(h$year))
time = seq(min(h$yr, na.rm=T), max(h$yr, na.rm=T))

#Format annual max/min dataframes
ws = as.data.frame(time)
ws$g1<-ws$g2<-ws$g3<-ws$g4<-ws$g5<-ws$g6<-rep(NA)

for(i in 1:length(time)){
  ws$g1[i]<-max(h$wind[h$yr==ws$time[i] & h$locGrp==0], na.rm=T)
  ws$g2[i]<-max(h$wind[h$yr==ws$time[i] & h$locGrp==1], na.rm=T)
  ws$g3[i]<-max(h$wind[h$yr==ws$time[i] & h$locGrp==2], na.rm=T)
  ws$g4[i]<-max(h$wind[h$yr==ws$time[i] & h$locGrp==3], na.rm=T)
  ws$g5[i]<-max(h$wind[h$yr==ws$time[i] & h$locGrp==4], na.rm=T)
  ws$g6[i]<-max(h$wind[h$yr==ws$time[i] & h$locGrp==5], na.rm=T)
  
}
ws<-replace(ws, ws==-Inf, NA)
ws<-replace(ws, ws==Inf, NA)

z = as.data.frame(time)
z$g1<-z$g2<-z$g3<-z$g4<-z$g5<-z$g6<-rep(NA)
h$pressure <-h$pressure *-1
for(i in 1:length(time)){
  z$g1[i]<-min(h$pressure[h$yr==z$time[i] & h$locGrp==0], na.rm=T)
  z$g2[i]<-min(h$pressure[h$yr==z$time[i] & h$locGrp==1], na.rm=T)
  z$g3[i]<-min(h$pressure[h$yr==z$time[i] & h$locGrp==2], na.rm=T)
  z$g4[i]<-min(h$pressure[h$yr==z$time[i] & h$locGrp==3], na.rm=T)
  z$g5[i]<-min(h$pressure[h$yr==z$time[i] & h$locGrp==4], na.rm=T)
  z$g6[i]<-min(h$pressure[h$yr==z$time[i] & h$locGrp==5], na.rm=T)
  
}
z<-replace(z, z==-Inf, NA)
z<-replace(z, z==Inf, NA)

lat = as.data.frame(time)
lat$g1<-lat$g2<-lat$g3<-lat$g4<-lat$g5<-lat$g6<-rep(NA)

for(i in 1:length(time)){
  lat$g1[i]<-min(h$minlat[h$yr==lat$time[i] & h$locGrp==0], na.rm=T)
  lat$g2[i]<-min(h$minlat[h$yr==lat$time[i] & h$locGrp==1], na.rm=T)
  lat$g3[i]<-min(h$minlat[h$yr==lat$time[i] & h$locGrp==2], na.rm=T)
  lat$g4[i]<-min(h$minlat[h$yr==lat$time[i] & h$locGrp==3], na.rm=T)
  lat$g5[i]<-min(h$minlat[h$yr==lat$time[i] & h$locGrp==4], na.rm=T)
  lat$g6[i]<-min(h$minlat[h$yr==lat$time[i] & h$locGrp==5], na.rm=T)
  
}

lat<-replace(lat, lat==-Inf, NA)
lat<-replace(lat, lat==Inf, NA)

#Create pdf file to save image
svg(file = "./figures/FigureS5.svg", width=5, height =6)

# set plot parameters
par(mfrow=c(3,1), mar=c(2,5,1,10), xpd=TRUE)


####### Windspeed panel
leg_txt = c("g1", "g2", "g3", "g4", "g5", "g6")
coef = rep(NA, 6)
pval = rep(NA, 6)
i=1
col <- paste("g", as.character(i), sep="")
m=lm(get(col)~time, data=ws[complete.cases(ws),])
# generate data for plotting trendline and prediction interval lines
coef[leg_txt %in% col] = print(summary(m)$coefficients[2,1])
pval[leg_txt %in% col] = print(summary(m)$coefficients[2,4])
l = as.data.frame(predict(m, ws, interval = "prediction"))
l$time = time

# plot raw data
plot(ws[,col]~ws$time, type="l", col=cmap[i], 
     ylim=c(0,200),
     ylab="annual max. windspeed \n(knots)",
     xlab="year")
l = as.data.frame(predict(m, ws[,c(col, "time")], interval = "prediction"))
l$time = time

# plot raw data
points(ws[,col]~ws$time, type="l", col=cmap[i])
#add trendline and prediction interval
points(l$fit~l$time, type="l", lwd=1, col=cmap[i])
# add confidence interval
lc = as.data.frame(predict(m, ws[,c(col,"time")], interval = "confidence"))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
# print lm coefficient and confidence interval in legend
print(summary(m)$coefficients[2,c(1,4)])


for(i in 2:6){
  col <- paste("g", as.character(i), sep="")
  m=lm(get(col)~time, data=ws[complete.cases(ws),])
  coef[leg_txt %in% col] = print(summary(m)$coefficients[2,1])
  pval[leg_txt %in% col] = print(summary(m)$coefficients[2,4])
  
  # generate data for plotting trendline and prediction interval lines
  l = as.data.frame(predict(m, ws[,c(col, "time")], interval = "prediction"))
  l$time = time
  
  # plot raw data
  points(ws[,col]~ws$time, type="l", col=cmap[i])
  #add trendline and prediction interval
  points(l$fit~l$time, type="l", lwd=2, col=cmap[i])
  # add confidence interval
  lc = as.data.frame(predict(m, ws[,c(col,"time")], interval = "confidence"))
  lc$time = time
  points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  # print lm coefficient and confidence interval in legend
  print(summary(m)$coefficients[2,c(1,4)])
  
}

stats = rep(NA)
for(i in 1:6){
  if(pval[i] <= 0.01){
    stats[i] <- paste(as.character(round(coef[i], digits=3)), "***", sep="")
  } else if (pval[i]<=0.05){
    stats[i] <- paste(as.character(round(coef[i], digits=3)), "**", sep="")
  } else if (pval[i]<=0.1){
    stats[i] <- paste(as.character(round(coef[i], digits=3)), "*", sep="")
  } else{
    stats[i] <- as.character(round(coef[i], digits=3))
  }
}
legend(x = "topright",inset=c(-0.275,0), paste(expression(beta), "=", stats, sep=" "), col=cmap, lty=1)

####### Pressure panel
leg_txt = c("g1", "g2", "g3", "g4", "g5", "g6")
coef = rep(NA, 6)
pval = rep(NA, 6)
i=1
col <- paste("g", as.character(i), sep="")
m=lm(get(col)~time, data=z[complete.cases(z),])
# generate data for plotting trendline and prediction interval lines
coef[leg_txt %in% col] = print(summary(m)$coefficients[2,1])
pval[leg_txt %in% col] = print(summary(m)$coefficients[2,4])
l = as.data.frame(predict(m, z, interval = "prediction"))
l$time = time

# plot raw data
plot(z[,col]~z$time, type="l", col=cmap[i], 
     ylim=c(900,1100),
     ylab="annual min. pressure \n(kPa)",
     xlab="year")
l = as.data.frame(predict(m, z[,c(col, "time")], interval = "prediction"))
l$time = time

# plot raw data
points(z[,col]~z$time, type="l", col=cmap[i])
#add trendline and prediction interval
points(l$fit~l$time, type="l", lwd=1, col=cmap[i])
# add confidence interval
lc = as.data.frame(predict(m, z[,c(col,"time")], interval = "confidence"))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
# print lm coefficient and confidence interval in legend
print(summary(m)$coefficients[2,c(1,4)])


for(i in 2:6){
  col <- paste("g", as.character(i), sep="")
  m=lm(get(col)~time, data=z[complete.cases(z),])
  coef[leg_txt %in% col] = print(summary(m)$coefficients[2,1])
  pval[leg_txt %in% col] = print(summary(m)$coefficients[2,4])
  
  # generate data for plotting trendline and prediction interval lines
  l = as.data.frame(predict(m, z[,c(col, "time")], interval = "prediction"))
  l$time = time
  
  # plot raw data
  points(z[,col]~z$time, type="l", col=cmap[i])
  #add trendline and prediction interval
  points(l$fit~l$time, type="l", lwd=2, col=cmap[i])
  # add confidence interval
  lc = as.data.frame(predict(m, z[,c(col,"time")], interval = "confidence"))
  lc$time = time
  points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  # print lm coefficient and confidence interval in legend
  print(summary(m)$coefficients[2,c(1,4)])

}

stats = rep(NA)
for(i in 1:6){
  if(pval[i] <= 0.01){
    stats[i] <- paste(as.character(round(coef[i], digits=3)), "***", sep="")
  } else if (pval[i]<=0.05){
    stats[i] <- paste(as.character(round(coef[i], digits=3)), "**", sep="")
  } else if (pval[i]<=0.1){
    stats[i] <- paste(as.character(round(coef[i], digits=3)), "*", sep="")
  } else{
    stats[i] <- as.character(round(coef[i], digits=3))
  }
}
legend(x = "topright",inset=c(-0.275,0), paste(expression(beta), "=", stats, sep=" "), col=cmap, lty=1)

####### Latitude panel
leg_txt = c("g1", "g2", "g3", "g4", "g5", "g6")
coef = rep(NA, 6)
pval = rep(NA, 6)
i=1
col <- paste("g", as.character(i), sep="")
m=lm(get(col)~time, data=lat[complete.cases(lat),])
# generate data for plotting trendline and prediction interval lines
coef[leg_txt %in% col] = print(summary(m)$coefficients[2,1])
pval[leg_txt %in% col] = print(summary(m)$coefficients[2,4])
l = as.data.frame(predict(m, lat, interval = "prediction"))
l$time = time

# plot raw data
plot(lat[,col]~lat$time, type="l", col=cmap[i], 
     ylim=c(5,50),
     ylab="annual min. storm latitude \n(degrees)",
     xlab="year")
l = as.data.frame(predict(m, lat[,c(col, "time")], interval = "prediction"))
l$time = time

# plot raw data
points(lat[,col]~lat$time, type="l", col=cmap[i])
#add trendline and prediction interval
points(l$fit~l$time, type="l", lwd=1, col=cmap[i])
# add confidence interval
lc = as.data.frame(predict(m, lat[,c(col,"time")], interval = "confidence"))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
print(summary(m)$coefficients[2,c(1,4)])


for(i in 2:6){
  col <- paste("g", as.character(i), sep="")
  m=lm(get(col)~time, data=lat[complete.cases(lat),])
  coef[leg_txt %in% col] = print(summary(m)$coefficients[2,1])
  pval[leg_txt %in% col] = print(summary(m)$coefficients[2,4])
  
  # generate data for plotting trendline and prediction interval lines
  l = as.data.frame(predict(m, lat[,c(col, "time")], interval = "prediction"))
  l$time = time
  
  # plot raw data
  points(lat[,col]~lat$time, type="l", col=cmap[i])
  #add trendline and prediction interval
  points(l$fit~l$time, type="l", lwd=2, col=cmap[i])
  # add confidence interval
  lc = as.data.frame(predict(m, lat[,c(col,"time")], interval = "confidence"))
  lc$time = time
  points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  # print lm coefficient and confidence interval in legend
  print(summary(m)$coefficients[2,c(1,4)])

}

stats = rep(NA)
for(i in 1:6){
  if(pval[i] <= 0.01){
    stats[i] <- paste(as.character(round(coef[i], digits=3)), "***", sep="")
  } else if (pval[i]<=0.05){
    stats[i] <- paste(as.character(round(coef[i], digits=3)), "**", sep="")
  } else if (pval[i]<=0.1){
    stats[i] <- paste(as.character(round(coef[i], digits=3)), "*", sep="")
  } else{
    stats[i] <- as.character(round(coef[i], digits=3))
  }
}
legend(x = "topright",inset=c(-0.275,0), paste(expression(beta), "=", stats, sep=" "), col=cmap, lty=1)

dev.off()

#################################################################################################3###########################
#
#     Figure S6 & S7, see ipynb
#
#
#############################################################################################################################

plot(d_fs$failTF_sum~(as.numeric(as.character((d_fs$yr)))), xlab="year", ylab="fs", pch=10)
tl<-lm(d_fs$failTF_sum~(as.numeric(as.character((d_fs$yr)))))
abline(tl$coefficients, col="red",lwd=3)
############################################################################################################################
#
#     Figure S8: a) summary of mixed-effects binomial model parameters and fit, predicting likelihood of >=1 pipeline failures given TC intensity,
#     with random intercepts on time and location, compared to null model; b) intercept shifts on location associated with prediction
#     skill of TC intensity as a parameter, with open circles representing null model intercept, closed circles representing 
#     TC intensity model intercept, and error bars representing 95% confidence intervals of intercept estimate; 3) as in b, but for 
#     time random intercept. 
############################################################################################################################
load("./model_outputs/model_outputs.RData")
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
locGrps = c("High Plains","Northeast", "MI Delta", "Southeast", "Midwest", "Texas Coast")


summary(mb1)
summary(mb2)


#Figure 1
#a, model table
tab_model(mb1,mb2, show.r2=FALSE, show.aic=TRUE,  show.ci=FALSE, show.icc=FALSE, show.loglik=TRUE, show.est=TRUE, show.re.var=TRUE, auto.label = FALSE, 
          pred.labels=c("Intercept","TCI"),
          dv.labels=c("Null model", "TC Intensity"),
          p.style="stars")


#c, year caterpillar plot
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


svg(file="./figures/FigureS8b.svg", width =9, height=4)
shapel = rep(c(1, 16), 46)
sizel = rep(c(4.5, 2.5),46)
p <- ggplot(test, aes(ID,y, group=grp))# + coord_flip()

p <- p + geom_hline(yintercept=0)
p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0.5, colour=rep(c("darkgray", "black"), 46))
p <- p + geom_point(aes(colour = ID), shape=shapel, colour=rep(topo.colors(46)[ord], each=2), size=sizel)
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                       panel.grid.minor = element_blank(), 
                       axis.line = element_line(colour = "black")) +
  labs(x = "Year",
       y = "Intercept",
       color = "Legend") +
  scale_color_manual(values = terrain.colors(11), breaks= seq(1975,2021, by=5))
p <- p +  scale_fill_continuous(guide = guide_colourbar())
dev.off()

#b, location caterpillar plot

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
svg(file="./figures/FigureS8c.svg", width =4.5, height=5)
p <- ggplot(test, aes(ID,y)) + coord_flip()
p <- p + theme(legend.position="none")
p <- p + geom_hline(yintercept=0)
p <- p + geom_point(size=sizel, shape=shapel, 
                    colour=rep(cmap[ord], each=2)) 
p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0.15, colour=rep(c("gray","black"), 6))#rep(terrain.colors(46)[rep(ord, each=2)], each=2))
p <- p + scale_x_discrete("location group", labels=locGrps[ord])
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))                      #  scale_fill_gradientn(
dev.off()



############################################################################################################################
#
#     Figure S9: Marilyn add
#
################################################################################################################
load("./model_outputs/model_outputs.RData")

############################################################################################################################
#
#     Figure S10
#
################################################################################################################
h<-read.csv("./raw_data/NOAAHURDAT.csv")
h$nameyr<-paste(h$name,as.character(h$year), sep="")
h$minlat<-NA
h$maxlon<-NA
for(i in 1:length(h$minlat)){
  h$minlat[i]<-min(h$lat[h$nameyr %in% h$nameyr[i]])
  h$maxlon[i]<-max(h$long[h$nameyr %in% h$nameyr[i]])
}
#Format dataframe of annual min/max values
time = seq(min(h$year), max(h$year))
cli = as.data.frame(time)
cli$lat<-rep(NA)
cli$lon<-rep(NA)
for(i in 1:length(time)){
  cli$lat[i]<-min(h$minlat[h$year ==time[i]])
  cli$lon[i]<-max(h$maxlon[h$year ==time[i]])
}
#Create pdf file to save image

svg(file = "./figures/FigureS10.svg", width=5, height =6)
# set plot parameters
par(mfrow=c(2,1), mar=c(2,5,1,1))
# plot lat panel
# train linear regression of lat by year
m1=lm(lat~time, data=cli)
# generate data for plotting trendline and prediction interval lines
l = as.data.frame(predict(m1, cli, interval = "prediction"))
l$time = time

# plot raw data
plot(cli$lat~cli$time, type="l", 
     ylim=c(min(l$lwr), max(l$upr)),
     ylab="annual min. lattitude \n(degrees)",
     xlab="year")
#add trendline and prediction interval
points(l$fit~l$time, type="l", lwd=2)
points(l$lwr~l$time, type="l", lwd=1, lty=2)
points(l$upr~l$time, type="l", lwd=1, lty=2)
# add confidence interval
lc = as.data.frame(predict(m1, cli, interval = "confidence"))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=2, col="blue")
points(lc$upr~lc$time, type="l", lwd=1, lty=2, col="blue")
# print lm coefficient and confidence interval in legend
legend(x="bottomleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))

# plot pressure panel
# train linear regression of annual max longitude by year
m2=lm(lon~time, data=cli)
l = as.data.frame(predict(m2, cli, interval = "prediction"))
l$time = time
plot(cli$lon~cli$time, type="l", 
     ylim=c(min(l$lwr), max(l$upr)),
     ylab="annual max. longitude \n(degrees)",
     xlab="year")
points(l$fit~l$time, type="l", lwd=2)
points(l$lwr~l$time, type="l", lwd=1, lty=2)
points(l$upr~l$time, type="l", lwd=1, lty=2)
lc = as.data.frame(predict(m2, cli, interval = "confidence"))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=2, col="blue")
points(lc$upr~lc$time, type="l", lwd=1, lty=2, col="blue")
legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m2)$coefficients[2,1], digits=3)), "**", sep="")))

dev.off()


