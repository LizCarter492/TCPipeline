############################################################################################################################
#          Figures and table from Smith and Carter (2025) 
#############################################################################################################################
#Load required libraries####
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
library(ggplot2)
library(viridis)
library(gridExtra)
library(gtable)
library(webshot)

####################################################################################
#
#         Figure 1
#
##################################################################################
#Fig 1a, b####
#created in ArcPro

#Fig 1c####
d= read.csv( "./analysis_data/analysis_ready.csv")

#format data
d$name = as.factor(d$name)
d$pressure = d$pressure_min 
d$locGrp = as.factor(d$locGrp)
d$per = NA
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
#create data frame of year and locgrp combos
years <- 1975:2022
locGrps <- 0:5
fo <- expand.grid(yr = years, locGrp = locGrps)
fo$per = NA
for(i in 1:length(fo$per)){
  if(fo$yr[i] < 1980){
    fo$per[i] <-1
  } else if (fo$yr[i] < 1985){
    fo$per[i] <-2
  } else if (fo$yr[i] < 1990){
    fo$per[i] <-3
  } else if (fo$yr[i] < 1995){
    fo$per[i] <-4
  } else if (fo$yr[i] < 2000){
    fo$per[i] <-5
  } else if (fo$yr[i] < 2005){
    fo$per[i] <-6
  } else if (fo$yr[i] < 2010){
    fo$per[i] <-7
  } else if (fo$yr[i] < 2015){
    fo$per[i] <-8
  } else {
    fo$per[i] <-9
  }
}

fo$fo<-NA

#Calculate fo
for(y in 1:dim(fo)[1]){
      fo$fo[y] <-mean(d$failTF_sum[which(d$locGrp==fo$locGrp[y] & d$per==fo$per[y] & vars=="NoStorm")], na.rm=FALSE)/(0.5*6.083)
}

fo$fo[fo$locGrp==4 & fo$per==1]<-19.151734

# Define the color map and locGrp levels
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
locGrps = c("High Plains", "Northeast", "MI Delta", "Southeast", "Midwest", "Texas Coast")


#Plot fs and fo
svg(file="./figures/Figure1c.svg", width =5, height=9)
par(mfrow = c(6, 1),          # 6 rows, 1 column
    mar = c(1, 4, 0.5, 6), 
    oma = c(3, 0, 2, 2))

for(i in 0:5){  
  # Plot the data
  plot(d$failTF_sum[d$locGrp == i] ~ as.numeric(as.character(d$yr[d$locGrp == i])), 
       col = cmap[i + 1], 
       pch = 16, 
       ylim=c(0, 100),
       xlim=c(1975,2022),
       ylab = expression(f[s]), 
       xaxt = ifelse(i == 5, "s", "n"),
       xlab = ifelse(i == 5, "Year", ""),  # Add xlab only to the bottom plot
       axes = TRUE)            # Disable automatic axes for control
  
  # Draw custom axes (only bottom plot gets x-axis)
  if (i == 5) {
    axis(1, at = seq(1975, 2022, by = 5))   # Custom x-axis ticks (year)
  }
  axis(2)  # Add y-axis
  
  # Get valid years and corresponding f_o values
  valid_fo_values <- fo$fo[fo$locGrp == i]
  valid_years_fo <- fo$yr[fo$locGrp == i]
  
  # Ensure polygon does not extend beyond valid years
  valid_data_indices <- which(valid_years_fo >= 1975 & valid_years_fo <= 2022)
  
  
  # Find the first non-zero 'fo' value to adjust the upper left corner of the polygon
  first_non_zero_fo_index <- which(valid_fo_values > 0)[1]  # Find first non-zero value
  first_non_zero_fo_value <- valid_fo_values[first_non_zero_fo_index]
  
 
  
  # Plot the gray background below the fo$fo line
  polygon(c(valid_years_fo[valid_data_indices], rev(valid_years_fo[valid_data_indices])),
          c(rep(-1.5, length(valid_years_fo[valid_data_indices])), rev(valid_fo_values[valid_data_indices])),
          col = "gray80", border = NA)  # Fill area with gray color

#  Plot the line for f_o 
  points(fo$fo[fo$locGrp == i] ~ fo$yr[fo$locGrp == i], type = "l", col = cmap[i + 1], lwd=2)  # Line for f_o
  
# replot points  
  points(d$failTF_sum[d$locGrp == i] ~ as.numeric(as.character(d$yr[d$locGrp == i])), pch=16, col = cmap[i + 1])  # Line for f_o
}
dev.off()





# Full legend labels, including f_o (line) and f_s (point)
svg(file="./figures/Figure1c_legend.svg", width =3, height=3)

legend_labels <- c(locGrps, expression(f[o]), expression(f[s]))

# Define symbols: NA for entries without point, NA for entries without line
pch_vals <- c(rep(16, 6), NA, 16)   # f[o] has no point, f[s] has point
lty_vals <- c(rep(NA, 6), 1, NA)    # f[o] has line, f[s] no line
col_vals <- c(cmap, "black", "black")  # Use black for f[o] and f[s]


par(mfrow = c(1, 1),      # Single figure for the legend
    mar = c(0, 0, 0, 0),  # No margins around the legend plot
    oma = c(0, 0, 0, 0))  # No outer margin for legend

# Create a blank plot
plot.new()
# Add legend (adjust inset and position as needed)
legend("center", inset = c(-0.25, 0), 
       legend = legend_labels,
       col = col_vals,
       pch = pch_vals,
       lty = lty_vals,
       pt.cex = 1.5,
       title = "Location Groups",
       xpd = NA, bty = "n")
dev.off()



############################################################################################################################
#
#         Figure 2
#
############################################################################################################################
load("./model_outputs/model_outputs.RData")
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
locGrps = c("High Plains","Northeast", "MI Delta", "Southeast", "Midwest", "Texas Coast")


summary(ms1)
summary(ms2)
summary(ms3)
#Fig2a, model table####
tab_model(ms1,ms2, show.r2=FALSE, show.aic=TRUE,  show.ci=FALSE, show.icc=FALSE, show.loglik=TRUE, show.est=TRUE, show.re.var=TRUE, auto.label = FALSE, 
          pred.labels=c("Intercept","TCI"),
          dv.labels=c("Null model", "TCI"),
          p.style="stars",
          file = "./figures/Figure2a.html")

#Fig2b, year caterpillar plot####
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


svg(file="./figures/Figure2b.svg", width =9, height=4)
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

#Fig2c, location caterpillar plot####

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


x <- ranef(ms3,condVar=TRUE)[[2]][[1]]
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

svg(file="./figures/Figure2c.svg", width =4.5, height=5)
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
#         Calculation of numeric values in logistic regression discussion
#
###################################################################################################################
library(car)
load("./model_outputs/hur_trends.RData")
load("./model_outputs/model_outputs.RData")

# Historical statistics 
u_ws <- u["wind"]
sd_ws <- sd["wind"]
u_z <- u["pressure"]
sd_z <- sd["pressure"]


# Initialize  cutoffs for Saffir Simpson Hurricane scale
wind = (c(34, 63,82,95,112, 136)-u_ws)/sd_ws
pressure = (c(1020, 1000,980,965,945,920)-u_z)/sd_z
hur_cat=data.frame(wind, pressure)

PC_hurcat <-as.data.frame(predict(pca, newdata=hur_cat))

#Summary table interpretting logistic regression in terms of SS hurricane scale
BTC<-as.data.frame(PC_hurcat$Comp.1)
names(BTC)<-c("HurComp")

year_value <- 1979
locGrp_value <- 5

# Add columns to BTC
BTC$yr <- rep(year_value, nrow(BTC))
BTC$locGrp <- rep(locGrp_value, nrow(BTC))

Intercept<-fixef(mb2)[1]
BetaTCI<-fixef(mb2)[2]

BTC

BTC$P<-(exp(Intercept) + exp(BetaTCI*BTC$HurComp))/(1+(exp(Intercept) + exp(BetaTCI*BTC$HurComp)))


exp(predict(mb2, BTC))/(1+exp(predict(mb2, BTC)))

exp(predict(mb2, BTC))/(1+exp(predict(mb2, BTC)))
############################################################################################################################
#
#         Figure 3
#
############################################################################################################################
#Fig3a,b####
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

for(i in 1:length(time)){
  cli$ws[i]<-max(h$wind[h$year==cli$time[i]], na.rm=T)
  cli$z[i]<-min(h$pressure[h$year==cli$time[i]], na.rm=T)
}
#Create svg file to save image
svg(file = "./figures/Figure3ab.svg", width=5, height =6)
# set plot parameters
par(mfrow=c(2,1), mar=c(2,5,1,1))
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


dev.off()
hm1<-m1
hm2<-m2
hm3<-m3
save(h, cli, hm1,hm2,hm3, file="./model_outputs/hur_trends.RData")

#Fig 3c####
rm(list=ls())
load("./model_outputs/hur_trends.RData")
load("./model_outputs/model_outputs.RData")
h<-read.csv("./raw_data/NOAAHURDAT.csv")
h_p<-read.csv("./analysis_data/analysis_ready_binomial.csv")

# Set plot parameters
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
locGrps = c("High Plains","Northeast", "MI Delta", "Southeast", "Midwest", "Texas Coast")

# Historical statistics 
u_ws <- u["wind"]
sd_ws <- sd["wind"]
u_z <- u["pressure"]
sd_z <- sd["pressure"]

# Predictions for 1970, 2010, 2050
time = c(1970, 2010, 2050)
pred = as.data.frame(time)
pred$ws<-pred$z<-rep(NA)

pred$ws = predict(hm1, pred)
pred$z = predict(hm2, pred)
pred

# Standardize predictions 
spred = pred
spred$ws = (pred$ws - u_ws)/sd_ws
spred$z = (pred$z - u_z)/sd_z
names(spred)<-c("year", "pressure", "wind")
spred

# Transform predictions into pc space
PC_pred <- as.data.frame(predict(pca, newdata=spred))
PC_pred$year = spred$year
PC_pred

# Initialize  cutoffs for Saffir Simpson Hurricane scale
wind = (c(34, 63,82,95,112, 136)-u_ws)/sd_ws
pressure = (c(1020, 1000,980,965,945,920)-u_z)/sd_z
hur_cat=data.frame(wind, pressure)

PC_hurcat <-as.data.frame(predict(pca, newdata=hur_cat))
PC_hurcat

# Pick an intercept value
yr_s=2014

#
d <- d[order(d$HurComp),]
par(mfrow=c(1,1))
svg(file = "./figures/Figure3c.svg", width=6, height =6)


#plot data
HurComp<-seq(-2.5, 7, by = 0.05)

#Initialize plot with fs' projections for the historical range of max TCI for locGrp 0
plot(exp(predict(ms2, data.frame(HurComp=d$HurComp[d$locGrp==0], locGrp=0, yr=yr_s),
                 na.action=na.exclude))~d$HurComp[d$locGrp==0], 
     type="l",
     lwd=4,
     col=cmap[1],
     xlim=c(-2.2, 6), 
     ylim=c(0, 35), 
     ylab="Poisson predicted fs'", xlab="TC intensity")


#vertical bars representing 1970, 2010, and 2050 annual maximum TCI
abline(v= PC_pred[1,1], lty =1, lwd=30, col="lightgoldenrod")
abline(v= PC_pred[2,1], lty =1, lwd=30, col="lightgoldenrod")
abline(v= PC_pred[3,1], lty =1, lwd=30, col="lightgoldenrod")

# plot confidence interval of predicted fs'
for(i in 0:5){
  p4p<-predict(ms2, data.frame(HurComp=HurComp, locGrp=i, yr=yr_s))
  se<-as.data.frame(emmeans(ms2, ~HurComp, at = list(x = HurComp, locGrp=i, yr=yr_s)))$SE
  
  polygon(x = c(HurComp, rev(HurComp)),
          y = c(exp(p4p) - 1.96*exp(se),
                rev(exp(p4p) +1.96*exp(se))),
          col =  adjustcolor(cmap[i+1], alpha.f = 0.1), border = NA)
  
  
}

#add fs' projections for the historical range of max TCI for locGrp 1-5
points(exp(predict(ms2, data.frame(HurComp=d$HurComp[d$locGrp==1], locGrp=1, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==1],type="l", lwd=4,col=cmap[2])
points(exp(predict(ms2, data.frame(HurComp=d$HurComp[d$locGrp==2], locGrp=2, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==2],type="l", lwd=4,col=cmap[3])
points(exp(predict(ms2, data.frame(HurComp=d$HurComp[d$locGrp==3], locGrp=3, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==3],type="l", lwd=4,col=cmap[4])
points(exp(predict(ms2, data.frame(HurComp=d$HurComp[d$locGrp==4], locGrp=4, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==4],type="l", lwd=6,col=cmap[5])
points(exp(predict(ms2, data.frame(HurComp=d$HurComp[d$locGrp==5], locGrp=5, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==5],type="l", lwd=4,col=cmap[6])

#add fs' projections for the full global range of max TCI
points(exp(predict(ms2, data.frame(HurComp=HurComp, locGrp=0, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[1])
points(exp(predict(ms2, data.frame(HurComp=HurComp, locGrp=1, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[2])
points(exp(predict(ms2, data.frame(HurComp=HurComp, locGrp=2, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[3])
points(exp(predict(ms2, data.frame(HurComp=HurComp, locGrp=3, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[4])
points(exp(predict(ms2, data.frame(HurComp=HurComp, locGrp=4, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[5])
points(exp(predict(ms2, data.frame(HurComp=HurComp, locGrp=5, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[6])

#plot local fs' projections, based on observed linear trends in WS and z 
for(i in 0:5){
  yr= seq(min(h_p$yr, na.rm=T), max(h_p$yr, na.rm=T))
  h_p_loc<-h_p[h_p$locGrp==i,]
  wind <- h_p_loc %>%
    group_by(yr) %>%
    summarise(max_wind = max(wind, na.rm = TRUE))
  pressure <- h_p_loc %>%
    group_by(yr) %>%
    summarise(min_z = min(pressure, na.rm = TRUE))
  wind$max_wind[wind$max_wind==10] = NA
  pressure$min_z[pressure$min_z==1050] = NA
  
  wind$max_wind = (wind$max_wind - u_ws)/sd_ws
  pressure$min_z = (pressure$min_z-u_z)/sd_z
  
  hurI<-merge(wind, pressure, by="yr")
  names(hurI)<-c("yr", "wind", "pressure")
  
  mw<-lm(wind~yr, data=hurI)
  mp<-lm(pressure~yr, data=hurI)
  
  lpm<-data.frame(yr=c(1970, 2010, 2050), wind=c(NA, NA, NA), pressure=c(NA, NA, NA))
  lpm$wind<-predict(mw, newdata=lpm)
  lpm$pressure<-predict(mp, newdata=lpm)
  lpm$HurComp<-predict(pca, lpm)[,1]
  print(exp(predict(ms2, data.frame(HurComp=lpm$HurComp, locGrp=i, yr=yr_s),na.action=na.exclude)))
  points(exp(predict(ms2, data.frame(HurComp=lpm$HurComp[1], locGrp=i, yr=yr_s),na.action=na.exclude))~lpm$HurComp[1],type="p", pch = 25, cex=1.5,bg=cmap[i+1])
  points(exp(predict(ms2, data.frame(HurComp=lpm$HurComp[2], locGrp=i, yr=yr_s),na.action=na.exclude))~lpm$HurComp[2],type="p", pch = 23, cex=1.5,bg=cmap[i+1])
  points(exp(predict(ms2, data.frame(HurComp=lpm$HurComp[3], locGrp=i, yr=yr_s),na.action=na.exclude))~lpm$HurComp[3],type="p", pch = 24, cex=1.5,bg=cmap[i+1])
}  

#draw lines corresponding to cutoffs for Safir Simpson hurricane scale
abline(v= PC_hurcat$Comp.1[1], lty =2, lwd=2, col="darkgray")# Tropical Storm
abline(v= PC_hurcat$Comp.1[2], lty =2, lwd=2, col="darkgray")# Category 1
abline(v= PC_hurcat$Comp.1[3], lty =2, lwd=2, col="darkgray")# Category 2
abline(v= PC_hurcat$Comp.1[4], lty =2, lwd=2, col="darkgray")# Category 3
abline(v= PC_hurcat$Comp.1[5], lty =2, lwd=2, col="darkgray")# Category 4
abline(v= PC_hurcat$Comp.1[6], lty =2, lwd=2, col="darkgray")# Category 5



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
#Fig 3d####
# Define column and row names
locGrps = c("High Plains","Northeast", "MI Delta", "Southeast", "Midwest", "Texas Coast")
col_names <- c("1970", "2010", "2050")
row_names <- as.character(0:5)

# Create an empty data frame
T3d <- data.frame(matrix(NA, nrow = length(row_names), ncol = length(col_names)),
                 stringsAsFactors = FALSE)

# Assign column and row names
colnames(T3d) <- col_names
rownames(T3d) <- row_names

# View the empty data frame
T3d
for(i in row_names){
  T3d[i,]<-exp(predict(ms2, data.frame(HurComp=PC_pred[,1], locGrp=i, yr=yr_s),na.action=na.exclude))
}
T3d<-round(T3d[,])
rownames(T3d)<-locGrps

gt <- tableGrob(T3d, theme = ttheme_minimal(
  widths=list(3,1,1,1),
))

# Set up the SVG device to capture the plot
svg(file = "./figures/Figure3d.svg", width=6, height =6)

# Draw the table in the SVG file
grid.draw(gt)

# Close the SVG device (saves the file)
dev.off()
# Print the table to the plotting window
############################################################################################################################
#
#         Table S1
#
############################################################################################################################
tab_model(mb1,mb2, mb3, show.r2=FALSE, show.aic=TRUE,  show.ci=FALSE, show.icc=FALSE, show.loglik=TRUE, show.est=TRUE, show.re.var=TRUE, auto.label = FALSE, 
          pred.labels=c("Intercept","TCI", "WS"),
          dv.labels=c("Null model", "TC Intensity", "Max Windspeed"),
          p.style="stars",
          file = "./figures/TableS1.html")


############################################################################################################################
#
#         Table S2
#
############################################################################################################################
tab_model(ms1,ms2, ms3, show.r2=FALSE, show.aic=TRUE,  show.ci=FALSE, show.icc=FALSE, show.loglik=TRUE, show.est=TRUE, show.re.var=TRUE, auto.label = FALSE, 
          pred.labels=c("Intercept","TCI", "WS"),
          dv.labels=c("Null model", "TC Intensity", "Max Windspeed"),
          p.style="stars",
          file = "./figures/TableS2.html")


############################################################################################################################
#
#         Figure S1
#
############################################################################################################################
#See .ipynb
############################################################################################################################
#
#         Figure S2
#
############################################################################################################################
#See .ipynb
############################################################################################################################
#
#         Figure S3
#
############################################################################################################################
#See .ipynb

############################################################################################################################
#
#         Figure S4
#
############################################################################################################################
load("./model_outputs/model_outputs.RData")
corr_matrix <- cor(dn)
svg(file = "./figures/FigureS4.svg", width=3, height =3)
ggcorrplot(corr_matrix)
dev.off()

############################################################################################################################
#
#         Figure S5
#
############################################################################################################################
load("./model_outputs/model_outputs.RData")
#Use the "label" and "select.ind" arguments to select only larger storms or high impact storms to label
svg(file = "./figures/FigureS5.svg", width=8, height =5)
fviz_pca_var(pca, col.var = "cos2",
             repel = FALSE)
dev.off()

############################################################################################################################
#
#         Figure S6
#
############################################################################################################################
#Read in data
h<-read.csv("./analysis_data/analysis_ready_binomial.csv")
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

############################################################################################################################
#
#         Figure S7
#
############################################################################################################################
load("./model_outputs/model_outputs.RData")
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
locGrps = c("High Plains","Northeast", "MI Delta", "Southeast", "Midwest", "Texas Coast")


summary(mb1)
summary(mb2)


#FigS7a#####
tab_model(mb1,mb2, show.r2=FALSE, show.aic=TRUE,  show.ci=FALSE, show.icc=FALSE, show.loglik=TRUE, show.est=TRUE, show.re.var=TRUE, auto.label = FALSE, 
          pred.labels=c("Intercept","TCI"),
          dv.labels=c("Null model", "TC Intensity"),
          p.style="stars",
          file="./figures/FigureS7a.html")


#FigS7c, year caterpillar plot####
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


svg(file="./figures/FigureS7b.svg", width =9, height=4)
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

#FigS7b, location caterpillar plot####
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
svg(file="./figures/FigureS7c.svg", width =4.5, height=5)
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
#         Figure S8
#
################################################################################################################
rm(list=ls())
load("./model_outputs/hur_trends.RData")
load("./model_outputs/model_outputs.RData")
h<-read.csv("./raw_data/NOAAHURDAT.csv")
h_p<-read.csv("./analysis_data/analysis_ready_binomial.csv")
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
locGrps = c("High Plains","Northeast", "MI Delta", "Southeast", "Midwest", "Texas Coast")

# Historical statistics 
u_ws <- u["wind"]
sd_ws <- sd["wind"]


time = c(1970, 2010, 2050)
pred = as.data.frame(time)
pred$ws<-rep(NA)
pred$ws = predict(hm1, newdata=pred)
pred

# Standardize predictions 
spred = pred
spred$ws = (pred$ws - u_ws)/sd_ws
names(spred)<-c("year", "wind")
spred

sswind = (c(34, 63,82,95,112, 136)-u_ws)/sd_ws

yr_s=2014
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
d <- d[order(d$ws_std),]
par(mfrow=c(1,1))
svg(file = "./figures/FigureS8.svg", width=6, height =6)

#m1<-ms1
#m2<-ms3
#plot data
#HurComp<-seq(min(d$HurComp), 7, by = 0.05)
ws_std<-seq(min(d$ws_std), 4, by = 0.05)
plot(exp(predict(ms3, data.frame(ws_std=d$ws_std[d$locGrp==0], locGrp=0, yr=yr_s),
                 na.action=na.exclude))~d$ws_std[d$locGrp==0], 
     type="l",
     lwd=4,
     col=cmap[1],
     xlim=c(-1.7, 4), 
     ylim=c(0, 35), 
     ylab="Poisson predicted fs'", xlab=expression("maximum windspeed (" * sigma * ")"))

abline(v= spred[1,2], lty =1, lwd=30, col="lightgoldenrod")#((136 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= spred[2,2], lty =1, lwd=30, col="lightgoldenrod")#((136 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= spred[3,2], lty =1, lwd=30, col="lightgoldenrod")#((136 - historical_mean_windspeed) / historical_sd_windspeed))

for(i in 0:5){
  p4p<-predict(ms3, data.frame(ws_std=ws_std, locGrp=i, yr=yr_s))
  se<-as.data.frame(emmeans(ms3, ~ws_std, at = list(x = ws_std, locGrp=i, yr=yr_s)))$SE
  
  polygon(x = c(ws_std, rev(ws_std)),
          y = c(exp(p4p) - 1.96*exp(se),
                rev(exp(p4p) +1.96*exp(se))),
          col =  adjustcolor(cmap[i+1], alpha.f = 0.1), border = NA)
  
  
}
points(exp(predict(ms3, data.frame(ws_std=d$ws_std[d$locGrp==1], locGrp=1, yr=yr_s),na.action=na.exclude))~d$ws_std[d$locGrp==1],type="l", lwd=4,col=cmap[2])
points(exp(predict(ms3, data.frame(ws_std=d$ws_std[d$locGrp==2], locGrp=2, yr=yr_s),na.action=na.exclude))~d$ws_std[d$locGrp==2],type="l", lwd=4,col=cmap[3])
points(exp(predict(ms3, data.frame(ws_std=d$ws_std[d$locGrp==3], locGrp=3, yr=yr_s),na.action=na.exclude))~d$ws_std[d$locGrp==3],type="l", lwd=4,col=cmap[4])
points(exp(predict(ms3, data.frame(ws_std=d$ws_std[d$locGrp==4], locGrp=4, yr=yr_s),na.action=na.exclude))~d$ws_std[d$locGrp==4],type="l", lwd=6,col=cmap[5])
points(exp(predict(ms3, data.frame(ws_std=d$ws_std[d$locGrp==5], locGrp=5, yr=yr_s),na.action=na.exclude))~d$ws_std[d$locGrp==5],type="l", lwd=4,col=cmap[6])

points(exp(predict(ms3, data.frame(ws_std=ws_std, locGrp=0, yr=yr_s),na.action=na.exclude))~ws_std,type="l", lwd=2,lty=3,col=cmap[1])
points(exp(predict(ms3, data.frame(ws_std=ws_std, locGrp=1, yr=yr_s),na.action=na.exclude))~ws_std,type="l", lwd=2,lty=3,col=cmap[2])
points(exp(predict(ms3, data.frame(ws_std=ws_std, locGrp=2, yr=yr_s),na.action=na.exclude))~ws_std,type="l", lwd=2,lty=3,col=cmap[3])
points(exp(predict(ms3, data.frame(ws_std=ws_std, locGrp=3, yr=yr_s),na.action=na.exclude))~ws_std,type="l", lwd=2,lty=3,col=cmap[4])
points(exp(predict(ms3, data.frame(ws_std=ws_std, locGrp=4, yr=yr_s),na.action=na.exclude))~ws_std,type="l", lwd=2,lty=3,col=cmap[5])
points(exp(predict(ms3, data.frame(ws_std=ws_std, locGrp=5, yr=yr_s),na.action=na.exclude))~ws_std,type="l", lwd=2,lty=3,col=cmap[6])
h$yr<-h$year

abline(v= sswind[1], lty =2, lwd=2, col="darkgray")#((33 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= sswind[2], lty =2, lwd=2, col="darkgray")#((63 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= sswind[3], lty =2, lwd=2, col="darkgray")#((82 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= sswind[4], lty =2, lwd=2, col="darkgray")#((95 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= sswind[5], lty =2, lwd=2, col="darkgray")#((112 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= sswind[6], lty =2, lwd=2, col="darkgray")#((136 - historical_mean_windspeed) / historical_sd_windspeed))



txt_y=rep(35, 7)
txt_x=c((min(ws_std)+sswind[1])/2,
        (sswind[1] + sswind[2])/2,
        (sswind[2] + sswind[3])/2,
        (sswind[3] + sswind[4])/2,
        (sswind[4] + sswind[5])/2,
        (sswind[5] + sswind[6])/2,
        ((sswind[6] + max(ws_std))/2-0.2))
labs=c(-1,0,1,2,3,4,5)
text( txt_x, txt_y,labs,
      cex = 1.5, pos = 1, col = "darkgray")

legend(x = "left", c(locGrps), col=c(cmap), lty=c(rep(1,6)),lwd=c(rep(3,6)),pch=c(rep(NA,6)))



dev.off()


############################################################################################################################
#
#         Figure S9
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


