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
library(emmeans)


#Set working directory to model_outputs
#setwd("./model_outputs")

#load data generated in step 2 (./code/2_modelTrain.R)
load("../model_outputs/model_outputs.RData")

cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
locGrps = c("High Plains","Northeast", "MI Delta", "Southeast", "Midwest", "Texas Coast")


############################################################################################################################
#
#     Figure 2: a) summary of mixed-effects binomial model parameters and fit, predicting likelihood of >=1 pipeline failures given TC intensity,
#     with random intercepts on time and location, compared to null model; b) intercept shifts on location associated with prediction
#     skill of TC intensity as a parameter, with open circles representing null model intercept, closed circles representing 
#     TC intensity model intercept, and error bars representing 95% confidence intervals of intercept estimate; 3) as in b, but for 
#     time random intercept. 
############################################################################################################################
load("./model_outputs/model_outputs.RData")
#Figure 1, table component (TODO: format these model summaries into a table, add them as a panel to figure 2, you don't need to use R)
summary(mb1)
summary(mb2)


#Figure 1
#a, model table
tab_model(mb1,mb2, show.r2=FALSE, show.aic=TRUE,  show.ci=FALSE, show.icc=FALSE, show.loglik=TRUE, show.est=TRUE, show.re.var=TRUE, auto.label = FALSE, 
          pred.labels=c("Intercept","TCI"),
          dv.labels=c("Null model", "TC Intensity"),
          p.style="stars")


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


svg(file="../figures/Figure1b.svg", width =9, height=4)
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
svg(file="../figures/Figure1c.svg", width =4.5, height=5)
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
#     Figure 3: a) summary of mixed-effects poisson model parameters and fit, predicting frequency of 60 day post-TC pipeline failures above average annual failure rate given TC intensity,
#     with random intercepts on time and location, compared to null model; b) intercept shifts on location associated with prediction
#     skill of TC intensity as a parameter, with open circles representing null model intercept, closed circles representing 
#     TC intensity model intercept, and error bars representing 95% confidence intervals of intercept estimate; 3) as in b, but for 
#     time random intercept. 
############################################################################################################################
load("./model_outputs/model_outputs.RData")
summary(m1)
summary(m2)

#a, model table
tab_model(m1,m2, show.r2=FALSE, show.aic=TRUE,  show.ci=FALSE, show.icc=FALSE, show.loglik=TRUE, show.est=TRUE, show.re.var=TRUE, auto.label = FALSE, 
          pred.labels=c("Intercept","TCI"),
          dv.labels=c("Null model", "TC Intensity"),
          p.style="stars")


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


#png(file="../figures/Figure2b.png", width =900, height=400)
svg(file="../figures/Figure2b.svg", width =9, height=4)
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

#Figure 2, location catepillar plot

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

svg(file="../figures/Figure2c.svg", width =4.5, height=5)
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
#     Figure 3_2: a) summary of mixed-effects poisson model parameters and fit, predicting frequency of 60 day post-TC pipeline failures above average annual failure rate given TC intensity,
#     with random intercepts on time and location, and a random slope on location, compared to null model; b) intercept shifts on location associated with prediction
#     skill of TC intensity as a parameter, with open circles representing null model intercept, closed circles representing 
#     TC intensity model intercept, and error bars representing 95% confidence intervals of intercept estimate; 3) as in b, but for 
#     time random intercept. 
############################################################################################################################
load("./model_outputs/model_outputs.RData")
summary(ms1)
summary(ms2)
#a, model table
tab_model(ms1,ms2, show.r2=FALSE, show.aic=TRUE,  show.ci=FALSE, show.icc=FALSE, show.loglik=TRUE, show.est=TRUE, show.re.var=TRUE, auto.label = FALSE, 
          pred.labels=c("Intercept","TCI"),
          dv.labels=c("Null model", "TC Intensity"),
          p.style="stars")


#Figure 2a model table
tab_model(ms1, ms2, show.ICC=FALSE, show.r2=FALSE, show.aic=TRUE)
#Figure 2b year caterpillar plot
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


#png(file="../figures/Figure2b.png", width =900, height=400)
svg(file="../figures/Figure2_2b.svg", width =9, height=4)
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

#Figure 2, location catepillar plot

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

svg(file="../figures/Figure2_2c.svg", width =4.5, height=5)
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
#                                       Figure 4d
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
#                                       Figure 4 Option 1 (All HURDAT)
#
############################################################################################################################
#Read in HURDAT2 raw data#####
h<-read.csv("../raw_data/NOAAHURDAT.csv")
h$nameyr<-paste(h$name,as.character(h$year), sep="")
h$minlat<-NA
for(i in 1:length(h$minlat)){
  h$minlat[i]<-min(h$lat[h$nameyr %in% h$nameyr[i]])
}
#Format dataframe of annual min/max values####
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
#Create pdf file to save image####
#pdf(file = "../figures/Figure4_1.pdf", width=5, height =6)
svg(file = "../figures/Figure4_1.svg", width=5, height =6)
# set plot parameters####
par(mfrow=c(3,1), mar=c(2,5,1,1))
# plot windspeed panel####
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

# plot pressure panel####
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

# plot latitude panel ####
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
#####
dev.off()
hm1<-m1
hm2<-m2
hm3<-m3
save(h, cli, hm1,hm2,hm3, file="../model_outputs/hur_trends.RData")

############################################################################################################################
#
#                                       Figure S9
############################################################################################################################
#Read in HURDAT2 raw data#####
h<-read.csv("./raw_data/NOAAHURDAT.csv")
h$nameyr<-paste(h$name,as.character(h$year), sep="")
h$minlat<-NA
h$maxlon<-NA
for(i in 1:length(h$minlat)){
  h$minlat[i]<-min(h$lat[h$nameyr %in% h$nameyr[i]])
  h$maxlon[i]<-max(h$long[h$nameyr %in% h$nameyr[i]])
}
#Format dataframe of annual min/max values####
time = seq(min(h$year), max(h$year))
cli = as.data.frame(time)
cli$lat<-rep(NA)
cli$lon<-rep(NA)
for(i in 1:length(time)){
  cli$lat[i]<-min(h$minlat[h$year ==time[i]])
  cli$lon[i]<-max(h$maxlon[h$year ==time[i]])
}
#Create pdf file to save image####

svg(file = "./figures/FigureS9.svg", width=5, height =6)
# set plot parameters####
par(mfrow=c(2,1), mar=c(2,5,1,1))
# plot lat panel####
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

# plot pressure panel####
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
#####
dev.off()


############################################################################################################################
#
#                   Figure 4 Option 1 Panel 2
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


#Format dataframe of annual min/max values####

time = seq(min(as.numeric(as.character(d$yr))), max(as.numeric(as.character(d$yr))))
cli = as.data.frame(time)
cli$HC<-rep(NA)

for(i in 1:length(time)){
  cli$HC[i]<-max(d$HurComp[as.numeric(as.character(d$yr))==cli$time[i]], na.rm=T)
}
cli$HC[cli$HC==-Inf]<-NA

# set plot parameters####
par(mfrow=c(3,1), mar=c(2,5,1,1))
# plot windspeed panel####
# train linear regression of annual max windspeed by year
m4=lm(HC~time, data=cli)
# generate data for plotting trendline and prediction interval lines
l = as.data.frame(predict(m4, cli, interval = "prediction"))
l$time = time

# plot raw data
plot(cli$HC~cli$time, type="l", 
     ylim=c(min(l$lwr), max(l$upr)),
     ylab="annual max. TC intensity",
     xlab="year")
#add trendline and prediction interval
points(l$fit~l$time, type="l", lwd=2)
points(l$lwr~l$time, type="l", lwd=1, lty=2)
points(l$upr~l$time, type="l", lwd=1, lty=2)
# add confidence interval
lc = as.data.frame(predict(m4, cli, interval = "confidence"))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=2, col="blue")
points(lc$upr~lc$time, type="l", lwd=1, lty=2, col="blue")
# print lm coefficient and confidence interval in legend
legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m4)$coefficients[2,1], digits=3)), "*", sep="")))


############################################################################################################################
#
#                   Figure 4 Option 2
#
############################################################################################################################
#Read in data####
h<-read.csv("./analysis_data/analysis_ready_binomial_k6.csv")
h <- h[h$lat>0,]
h$year<-h$yr
h$pressure <- -1 * h$pressure
#Create dataframe for annual max/min values
time = seq(min(h$year), max(h$year))
#Format annual max/min values for plotting####
cli = as.data.frame(time)
cli$ws<-rep(NA)
cli$z<-rep(NA)
cli$lat<-rep(NA)

for(i in 1:length(time)){
  cli$ws[i]<-max(h$wind[h$year==cli$time[i]], na.rm=T)
  cli$z[i]<-min(h$pressure[h$year==cli$time[i]], na.rm=T)
  cli$lat[i]<-min(h$minlat[h$year==cli$time[i]], na.rm=T)

}
cli[cli==-Inf]<-NA
cli[cli==Inf]<-NA
#Create pdf file to save image###
#pdf(file = "../figures/Figure4_2.pdf", width=5, height =6)
svg(file = "../figures/Figure4_2.svg", width=5, height =6)

# set plot parameters####
par(mfrow=c(3,1), mar=c(2,5,1,1))
# train linear regression of annual max windspeed by year
m1=lm(ws~time, data=cli[complete.cases(cli),])
# generate data for plotting trendline and prediction interval lines
l = as.data.frame(predict(m1, cli, interval = "prediction"))
l$time = time

# plot winspeed panel ####
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
legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "", sep="")))

# plot pressure panel####
# train linear regression of annual min pressure by year
m2=lm(z~time, data=cli[complete.cases(cli),])
l = as.data.frame(predict(m2, cli, interval = "prediction"))
l$time = time
plot(cli$z~cli$time, type="l", 
     ylim=c(min(l$lwr), max(cli$z[complete.cases(cli)])),
     ylab="annual min. pressure \n(kPa)",
     xlab="year")
points(l$fit~l$time, type="l", lwd=2)
points(l$lwr~l$time, type="l", lwd=1, lty=2)
points(l$upr~l$time, type="l", lwd=1, lty=2)
lc = as.data.frame(predict(m2, cli, interval = "confidence"))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=2, col="blue")
points(lc$upr~lc$time, type="l", lwd=1, lty=2, col="blue")
legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m2)$coefficients[2,1], digits=3)), "", sep="")))

#plot latitude panel ####
# train linear regression of annual min latitude by year
m3=lm(lat~time, data=cli[complete.cases(cli),])
l = as.data.frame(predict(m3, cli, interval = "prediction"))
l$time = time
plot(cli$lat~cli$time, type="l", 
     ylim=c(min(l$lwr), max(l$upr)),
     ylab="annual min. latitude \n(degrees)",
     xlab="year")
points(l$fit~l$time, type="l", lwd=2)
points(l$lwr~l$time, type="l", lwd=1, lty=2)
points(l$upr~l$time, type="l", lwd=1, lty=2)
lc = as.data.frame(predict(m3, cli, interval = "confidence"))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=2, col="blue")
points(lc$upr~lc$time, type="l", lwd=1, lty=2, col="blue")
legend(x="bottomleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m3)$coefficients[2,1], digits=3)), "*", sep="")))
#####
dev.off()


############################################################################################################################
#
#                                       Figure 4 option 3
#
############################################################################################################################
#Read in data####

h<-read.csv("../analysis_data/analysis_ready_binomial_k6.csv")
h <- h[h$lat>0,]
#Create dataframe for annual max/min values
#time = seq(min(h$year), max(h$year))
time = seq(min(h$yr, na.rm=T), max(h$yr, na.rm=T))

#Format annual max/min dataframes####
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

#Create pdf file to save image#####
#pdf(file = "../figures/Figure4_3.pdf", width=5, height =6)
svg(file = "../figures/Figure4_3.svg", width=5, height =6)

# set plot parameters
par(mfrow=c(3,1), mar=c(2,5,1,10), xpd=TRUE)


####### Windspeed panel#####
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
#l$time = seq(min(h$year), max(h$year))
l$time = time

# plot raw data
plot(ws[,col]~ws$time, type="l", col=cmap[i], 
     ylim=c(0,200),
     ylab="annual max. windspeed \n(knots)",
     xlab="year")
l = as.data.frame(predict(m, ws[,c(col, "time")], interval = "prediction"))
#l$time = seq(min(h$year), max(h$year))
l$time = time

# plot raw data
points(ws[,col]~ws$time, type="l", col=cmap[i])
#add trendline and prediction interval
points(l$fit~l$time, type="l", lwd=1, col=cmap[i])
#points(l$lwr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
#points(l$upr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
# add confidence interval
lc = as.data.frame(predict(m, ws[,c(col,"time")], interval = "confidence"))
#lc$time = seq(min(h$year), max(h$year))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
# print lm coefficient and confidence interval in legend
#legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
print(summary(m)$coefficients[2,c(1,4)])


for(i in 2:6){
  col <- paste("g", as.character(i), sep="")
  m=lm(get(col)~time, data=ws[complete.cases(ws),])
  coef[leg_txt %in% col] = print(summary(m)$coefficients[2,1])
  pval[leg_txt %in% col] = print(summary(m)$coefficients[2,4])
  
  # generate data for plotting trendline and prediction interval lines
  l = as.data.frame(predict(m, ws[,c(col, "time")], interval = "prediction"))
  #l$time = seq(min(h$year), max(h$year))
  l$time = time
  
  # plot raw data
  points(ws[,col]~ws$time, type="l", col=cmap[i])
  #add trendline and prediction interval
  points(l$fit~l$time, type="l", lwd=2, col=cmap[i])
  #points(l$lwr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
  #points(l$upr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
  # add confidence interval
  lc = as.data.frame(predict(m, ws[,c(col,"time")], interval = "confidence"))
  #lc$time = seq(min(h$year), max(h$year))
  lc$time = time
  points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  # print lm coefficient and confidence interval in legend
  print(summary(m)$coefficients[2,c(1,4)])
  # legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
  
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

####### Pressure panel#####
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
#l$time = seq(min(h$year), max(h$year))
l$time = time

# plot raw data
plot(z[,col]~z$time, type="l", col=cmap[i], 
     ylim=c(900,1100),
     ylab="annual min. pressure \n(kPa)",
     xlab="year")
l = as.data.frame(predict(m, z[,c(col, "time")], interval = "prediction"))
#l$time = seq(min(h$year), max(h$year))
l$time = time

# plot raw data
points(z[,col]~z$time, type="l", col=cmap[i])
#add trendline and prediction interval
points(l$fit~l$time, type="l", lwd=1, col=cmap[i])
#points(l$lwr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
#points(l$upr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
# add confidence interval
lc = as.data.frame(predict(m, z[,c(col,"time")], interval = "confidence"))
#lc$time = seq(min(h$year), max(h$year))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
# print lm coefficient and confidence interval in legend
#legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
print(summary(m)$coefficients[2,c(1,4)])


for(i in 2:6){
  col <- paste("g", as.character(i), sep="")
  m=lm(get(col)~time, data=z[complete.cases(z),])
  coef[leg_txt %in% col] = print(summary(m)$coefficients[2,1])
  pval[leg_txt %in% col] = print(summary(m)$coefficients[2,4])
  
  # generate data for plotting trendline and prediction interval lines
  l = as.data.frame(predict(m, z[,c(col, "time")], interval = "prediction"))
  #l$time = seq(min(h$year), max(h$year))
  l$time = time
  
  # plot raw data
  points(z[,col]~z$time, type="l", col=cmap[i])
  #add trendline and prediction interval
  points(l$fit~l$time, type="l", lwd=2, col=cmap[i])
  #points(l$lwr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
  #points(l$upr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
  # add confidence interval
  lc = as.data.frame(predict(m, z[,c(col,"time")], interval = "confidence"))
  #lc$time = seq(min(h$year), max(h$year))
  lc$time = time
  points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  # print lm coefficient and confidence interval in legend
  print(summary(m)$coefficients[2,c(1,4)])
  # legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
  
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

####### Latitude panel#####
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
#l$time = seq(min(h$year), max(h$year))
l$time = time

# plot raw data
plot(lat[,col]~lat$time, type="l", col=cmap[i], 
     ylim=c(5,50),
     ylab="annual min. storm latitude \n(degrees)",
     xlab="year")
l = as.data.frame(predict(m, lat[,c(col, "time")], interval = "prediction"))
#l$time = seq(min(h$year), max(h$year))
l$time = time

# plot raw data
points(lat[,col]~lat$time, type="l", col=cmap[i])
#add trendline and prediction interval
points(l$fit~l$time, type="l", lwd=1, col=cmap[i])
#points(l$lwr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
#points(l$upr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
# add confidence interval
lc = as.data.frame(predict(m, lat[,c(col,"time")], interval = "confidence"))
#lc$time = seq(min(h$year), max(h$year))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
# print lm coefficient and confidence interval in legend
#legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
print(summary(m)$coefficients[2,c(1,4)])


for(i in 2:6){
  col <- paste("g", as.character(i), sep="")
  m=lm(get(col)~time, data=lat[complete.cases(lat),])
  coef[leg_txt %in% col] = print(summary(m)$coefficients[2,1])
  pval[leg_txt %in% col] = print(summary(m)$coefficients[2,4])
  
  # generate data for plotting trendline and prediction interval lines
  l = as.data.frame(predict(m, lat[,c(col, "time")], interval = "prediction"))
  #l$time = seq(min(h$year), max(h$year))
  l$time = time
  
  # plot raw data
  points(lat[,col]~lat$time, type="l", col=cmap[i])
  #add trendline and prediction interval
  points(l$fit~l$time, type="l", lwd=2, col=cmap[i])
  #points(l$lwr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
  #points(l$upr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
  # add confidence interval
  lc = as.data.frame(predict(m, lat[,c(col,"time")], interval = "confidence"))
  #lc$time = seq(min(h$year), max(h$year))
  lc$time = time
  points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  # print lm coefficient and confidence interval in legend
  print(summary(m)$coefficients[2,c(1,4)])
  # legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
  
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
#####
dev.off()
############################################################################################################################
#
#                                       Figure 4 option 4
#
############################################################################################################################
#Read in data####
#h<-read.csv("../raw_data/NOAAHURDAT.csv")

h<-read.csv("../analysis_data/analysis_ready_binomial_k6.csv")
h <- h[h$lat>0,]
#Create dataframe for annual max/min values
#time = seq(min(h$year), max(h$year))
time = seq(min(h$yr, na.rm=T), max(h$yr, na.rm=T)-1)

#Format annual max/min dataframes####
ws = as.data.frame(time)
ws$g1<-ws$g2<-ws$g3<-ws$g4<-ws$g5<-ws$g6<-rep(NA)

for(i in 1:length(time)){
  ws$g1[i]<-max(h$wind[h$yr==cli$time[i] & h$locGrp==0], na.rm=T)
  ws$g2[i]<-max(h$wind[h$yr==cli$time[i] & h$locGrp==1], na.rm=T)
  ws$g3[i]<-max(h$wind[h$yr==cli$time[i] & h$locGrp==2], na.rm=T)
  ws$g4[i]<-max(h$wind[h$yr==cli$time[i] & h$locGrp==3], na.rm=T)
  ws$g5[i]<-max(h$wind[h$yr==cli$time[i] & h$locGrp==4], na.rm=T)
  ws$g6[i]<-max(h$wind[h$yr==cli$time[i] & h$locGrp==5], na.rm=T)
  
}
ws<-replace(ws, ws==-Inf, NA)
ws<-replace(ws, ws==Inf, NA)

z = as.data.frame(time)
z$g1<-z$g2<-z$g3<-z$g4<-z$g5<-z$g6<-rep(NA)
h$pressure <-h$pressure *-1
for(i in 1:length(time)){
  z$g1[i]<-min(h$pressure[h$yr==cli$time[i] & h$locGrp==0], na.rm=T)
  z$g2[i]<-min(h$pressure[h$yr==cli$time[i] & h$locGrp==1], na.rm=T)
  z$g3[i]<-min(h$pressure[h$yr==cli$time[i] & h$locGrp==2], na.rm=T)
  z$g4[i]<-min(h$pressure[h$yr==cli$time[i] & h$locGrp==3], na.rm=T)
  z$g5[i]<-min(h$pressure[h$yr==cli$time[i] & h$locGrp==4], na.rm=T)
  z$g6[i]<-min(h$pressure[h$yr==cli$time[i] & h$locGrp==5], na.rm=T)
  
}
z<-replace(z, z==-Inf, NA)
z<-replace(z, z==Inf, NA)

lat = as.data.frame(time)
lat$g1<-lat$g2<-lat$g3<-lat$g4<-lat$g5<-lat$g6<-rep(NA)

for(i in 1:length(time)){
  lat$g1[i]<-min(h$minlat[h$yr==cli$time[i] & h$locGrp==0], na.rm=T)
  lat$g2[i]<-min(h$minlat[h$yr==cli$time[i] & h$locGrp==1], na.rm=T)
  lat$g3[i]<-min(h$minlat[h$yr==cli$time[i] & h$locGrp==2], na.rm=T)
  lat$g4[i]<-min(h$minlat[h$yr==cli$time[i] & h$locGrp==3], na.rm=T)
  lat$g5[i]<-min(h$minlat[h$yr==cli$time[i] & h$locGrp==4], na.rm=T)
  lat$g6[i]<-min(h$minlat[h$yr==cli$time[i] & h$locGrp==5], na.rm=T)
  
}

lat<-replace(lat, lat==-Inf, NA)
lat<-replace(lat, lat==Inf, NA)

#Create pdf file to save image#####
pdf(file = "../figures/Figure4_4.pdf", width=7, height =6)
# set plot parameters
par(mfrow=c(3,1), mar=c(2,5,1,10), xpd=TRUE)


####### Windspeed panel#####
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
#l$time = seq(min(h$year), max(h$year))
l$time = time

# plot raw data
plot(ws[,col]~ws$time, type="p", col=cmap[i], 
     ylim=c(0,200),
     ylab="annual max. windspeed \n(knots)",
     xlab="year")
l = as.data.frame(predict(m, ws[,c(col, "time")], interval = "prediction"))
#l$time = seq(min(h$year), max(h$year))
l$time = time

# plot raw data
points(ws[,col]~ws$time, type="p", col=cmap[i])
#add trendline and prediction interval
points(l$fit~l$time, type="l", lwd=2, col=cmap[i])
#points(l$lwr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
#points(l$upr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
# add confidence interval
lc = as.data.frame(predict(m, ws[,c(col,"time")], interval = "confidence"))
#lc$time = seq(min(h$year), max(h$year))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
# print lm coefficient and confidence interval in legend
#legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
print(summary(m)$coefficients[2,c(1,4)])


for(i in 2:6){
  col <- paste("g", as.character(i), sep="")
  m=lm(get(col)~time, data=ws[complete.cases(ws),])
  coef[leg_txt %in% col] = print(summary(m)$coefficients[2,1])
  pval[leg_txt %in% col] = print(summary(m)$coefficients[2,4])
  
  # generate data for plotting trendline and prediction interval lines
  l = as.data.frame(predict(m, ws[,c(col, "time")], interval = "prediction"))
  #l$time = seq(min(h$year), max(h$year))
  l$time = time
  
  # plot raw data
  points(ws[,col]~ws$time, type="p", col=cmap[i])
  #add trendline and prediction interval
  points(l$fit~l$time, type="l", lwd=2, col=cmap[i])
  #points(l$lwr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
  #points(l$upr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
  # add confidence interval
  lc = as.data.frame(predict(m, ws[,c(col,"time")], interval = "confidence"))
  #lc$time = seq(min(h$year), max(h$year))
  lc$time = time
  points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  # print lm coefficient and confidence interval in legend
  print(summary(m)$coefficients[2,c(1,4)])
  # legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
  
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

####### Pressure panel#####
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
#l$time = seq(min(h$year), max(h$year))
l$time = time

# plot raw data
plot(z[,col]~z$time, type="p", col=cmap[i], 
     ylim=c(900,1100),
     ylab="annual min. pressure \n(kPa)",
     xlab="year")
l = as.data.frame(predict(m, z[,c(col, "time")], interval = "prediction"))
#l$time = seq(min(h$year), max(h$year))
l$time = time

# plot raw data
points(z[,col]~z$time, type="p", col=cmap[i])
#add trendline and prediction interval
points(l$fit~l$time, type="l", lwd=2, col=cmap[i])
#points(l$lwr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
#points(l$upr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
# add confidence interval
lc = as.data.frame(predict(m, z[,c(col,"time")], interval = "confidence"))
#lc$time = seq(min(h$year), max(h$year))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
# print lm coefficient and confidence interval in legend
#legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
print(summary(m)$coefficients[2,c(1,4)])


for(i in 2:6){
  col <- paste("g", as.character(i), sep="")
  m=lm(get(col)~time, data=z[complete.cases(z),])
  coef[leg_txt %in% col] = print(summary(m)$coefficients[2,1])
  pval[leg_txt %in% col] = print(summary(m)$coefficients[2,4])
  
  # generate data for plotting trendline and prediction interval lines
  l = as.data.frame(predict(m, z[,c(col, "time")], interval = "prediction"))
  #l$time = seq(min(h$year), max(h$year))
  l$time = time
  
  # plot raw data
  points(z[,col]~z$time, type="p", col=cmap[i])
  #add trendline and prediction interval
  points(l$fit~l$time, type="l", lwd=2, col=cmap[i])
  #points(l$lwr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
  #points(l$upr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
  # add confidence interval
  lc = as.data.frame(predict(m, z[,c(col,"time")], interval = "confidence"))
  #lc$time = seq(min(h$year), max(h$year))
  lc$time = time
  points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  # print lm coefficient and confidence interval in legend
  print(summary(m)$coefficients[2,c(1,4)])
  # legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
  
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

####### Latitude panel#####
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
#l$time = seq(min(h$year), max(h$year))
l$time = time

# plot raw data
plot(lat[,col]~lat$time, type="p", col=cmap[i], 
     ylim=c(5,50),
     ylab="annual min. storm latitude \n(degrees)",
     xlab="year")
l = as.data.frame(predict(m, lat[,c(col, "time")], interval = "prediction"))
#l$time = seq(min(h$year), max(h$year))
l$time = time

# plot raw data
points(lat[,col]~lat$time, type="p", col=cmap[i])
#add trendline and prediction interval
points(l$fit~l$time, type="l", lwd=2, col=cmap[i])
#points(l$lwr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
#points(l$upr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
# add confidence interval
lc = as.data.frame(predict(m, lat[,c(col,"time")], interval = "confidence"))
#lc$time = seq(min(h$year), max(h$year))
lc$time = time
points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
# print lm coefficient and confidence interval in legend
#legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
print(summary(m)$coefficients[2,c(1,4)])


for(i in 2:6){
  col <- paste("g", as.character(i), sep="")
  m=lm(get(col)~time, data=lat[complete.cases(lat),])
  coef[leg_txt %in% col] = print(summary(m)$coefficients[2,1])
  pval[leg_txt %in% col] = print(summary(m)$coefficients[2,4])
  
  # generate data for plotting trendline and prediction interval lines
  l = as.data.frame(predict(m, lat[,c(col, "time")], interval = "prediction"))
  #l$time = seq(min(h$year), max(h$year))
  l$time = time
  
  # plot raw data
  points(lat[,col]~lat$time, type="p", col=cmap[i])
  #add trendline and prediction interval
  points(l$fit~l$time, type="l", lwd=2, col=cmap[i])
  #points(l$lwr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
  #points(l$upr~l$time, type="l", lwd=1, lty=2, col=cmap[i])
  # add confidence interval
  lc = as.data.frame(predict(m, lat[,c(col,"time")], interval = "confidence"))
  #lc$time = seq(min(h$year), max(h$year))
  lc$time = time
  points(lc$lwr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  points(lc$upr~lc$time, type="l", lwd=1, lty=3, col=cmap[i])
  # print lm coefficient and confidence interval in legend
  print(summary(m)$coefficients[2,c(1,4)])
  # legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
  
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
#####
dev.off()


############################################################################################################################
#
#                                       Figure 4 option 5
#
############################################################################################################################
#Read in data####
#h<-read.csv("../raw_data/NOAAHURDAT.csv")

h<-read.csv("../analysis_data/analysis_ready_binomial_k6.csv")
h <- h[h$lat>0,]
h$pressure = h$pressure*-1
#Create dataframe for annual max/min values
time = seq(min(h$yr), max(h$yr))

#Format annual max/min dataframes####
#Create pdf file to save image#####
pdf(file = "../figures/Figure4_5.pdf", width=7, height =6)
# set plot parameters
par(mfrow=c(3,1), mar=c(2,5,1,10), xpd=TRUE)


####### Windspeed panel#####
coef = rep(NA, 6)
pval = rep(NA, 6)
yr= seq(min(h$yr, na.rm=T), max(h$yr, na.rm=T))
i=0
m=lm(wind~yr, data=h[h$locGrp==i,])
# generate data for plotting trendline and prediction interval lines
coef[i+1] = print(summary(m)$coefficients[2,1])
pval[i+1] = print(summary(m)$coefficients[2,4])
d = as.data.frame(yr)
d$wind = rep(NA)
l = as.data.frame(predict(m, d, interval = "prediction"))

# plot raw data
plot(h[h$locGrp==i,"wind"]~h[h$locGrp==i,"yr"], type="p", col=cmap[i+1], 
     ylim=c(-10,160),
     xlim=c(1974,2022),
     ylab="averate TC windspeed \n(knots)",
     xlab="year")
#add trendline and prediction interval
points(l$fit~d$yr, type="l", lwd=2, col=cmap[i+1])
# add confidence interval
points(l$lwr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
points(l$upr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
# print lm coefficient and confidence interval in legend
#legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
print(summary(m)$coefficients[2,c(1,4)])


for(i in 1:5){
  m=lm(wind~yr, data=h[h$locGrp==i,])
  # generate data for plotting trendline and prediction interval lines
  coef[i+1] = print(summary(m)$coefficients[2,1])
  pval[i+1] = print(summary(m)$coefficients[2,4])
  d = as.data.frame(yr)
  d$wind = rep(NA)
  l = as.data.frame(predict(m, d, interval = "prediction"))
  
  # plot raw data
  points(h[h$locGrp==i,"wind"]~h[h$locGrp==i,"yr"], type="p", col=cmap[i+1])
  #add trendline and prediction interval
  points(l$fit~d$yr, type="l", lwd=2, col=cmap[i+1])
  # add confidence interval
  points(l$lwr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
  points(l$upr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
  # print lm coefficient and confidence interval in legend
  #legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
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
legend(x = "topright",inset=c(-0.21,0), paste(expression(beta), "=", stats, sep=" "), col=cmap, lty=1)

####### Pressure panel#####
coef = rep(NA, 6)
pval = rep(NA, 6)
yr= seq(min(h$yr, na.rm=T), max(h$yr, na.rm=T))
i=0
m=lm(pressure~yr, data=h[h$locGrp==i,])
# generate data for plotting trendline and prediction interval lines
coef[i+1] = print(summary(m)$coefficients[2,1])
pval[i+1] = print(summary(m)$coefficients[2,4])
d = as.data.frame(yr)
d$pressure = rep(NA)
l = as.data.frame(predict(m, d, interval = "prediction"))

# plot raw data
plot(h[h$locGrp==i,"pressure"]~h[h$locGrp==i,"yr"], type="p", col=cmap[i+1], 
     ylim=c(910,1050),
     xlim=c(1974,2022),
     ylab="average TC pressure \n(kPa)",
     xlab="year")
#add trendline and prediction interval
points(l$fit~d$yr, type="l", lwd=2, col=cmap[i+1])
# add confidence interval
points(l$lwr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
points(l$upr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
# print lm coefficient and confidence interval in legend
#legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
print(summary(m)$coefficients[2,c(1,4)])


for(i in 1:5){
  m=lm(pressure~yr, data=h[h$locGrp==i,])
  # generate data for plotting trendline and prediction interval lines
  coef[i+1] = print(summary(m)$coefficients[2,1])
  pval[i+1] = print(summary(m)$coefficients[2,4])
  d = as.data.frame(yr)
  d$pressure = rep(NA)
  l = as.data.frame(predict(m, d, interval = "prediction"))
  
  # plot raw data
  points(h[h$locGrp==i,"pressure"]~h[h$locGrp==i,"yr"], type="p", col=cmap[i+1])
  #add trendline and prediction interval
  points(l$fit~d$yr, type="l", lwd=2, col=cmap[i+1])
  # add confidence interval
  points(l$lwr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
  points(l$upr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
  # print lm coefficient and confidence interval in legend
  #legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
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
legend(x = "topright",inset=c(-0.21,0), paste(expression(beta), "=", stats, sep=" "), col=cmap, lty=1)
####### Latitude panel#####
coef = rep(NA, 6)
pval = rep(NA, 6)
yr= seq(min(h$yr, na.rm=T), max(h$yr, na.rm=T))
i=0
m=lm(minlat~yr, data=h[h$locGrp==i,])
# generate data for plotting trendline and prediction interval lines
coef[i+1] = print(summary(m)$coefficients[2,1])
pval[i+1] = print(summary(m)$coefficients[2,4])
d = as.data.frame(yr)
d$minlat = rep(NA)
l = as.data.frame(predict(m, d, interval = "prediction"))

# plot raw data
plot(h[h$locGrp==i,"minlat"]~h[h$locGrp==i,"yr"], type="p", col=cmap[i+1], 
     ylim=c(0,45),
     xlim=c(1974,2022),
     ylab="average TC starting latitude \n(degrees)",
     xlab="year")
#add trendline and prediction interval
points(l$fit~d$yr, type="l", lwd=2, col=cmap[i+1])
# add confidence interval
points(l$lwr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
points(l$upr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
# print lm coefficient and confidence interval in legend
#legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
print(summary(m)$coefficients[2,c(1,4)])


for(i in 1:5){
  m=lm(minlat~yr, data=h[h$locGrp==i,])
  # generate data for plotting trendline and prediction interval lines
  coef[i+1] = print(summary(m)$coefficients[2,1])
  pval[i+1] = print(summary(m)$coefficients[2,4])
  d = as.data.frame(yr)
  d$minlat = rep(NA)
  l = as.data.frame(predict(m, d, interval = "prediction"))
  
  # plot raw data
  points(h[h$locGrp==i,"minlat"]~h[h$locGrp==i,"yr"], type="p", col=cmap[i+1])
  #add trendline and prediction interval
  points(l$fit~d$yr, type="l", lwd=2, col=cmap[i+1])
  # add confidence interval
  points(l$lwr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
  points(l$upr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
  # print lm coefficient and confidence interval in legend
  #legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
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
legend(x = "topright",inset=c(-0.21,0), paste(expression(beta), "=", stats, sep=" "), col=cmap, lty=1)

#####
dev.off()

############################################################################################################################
#
#                                       Figure 4 option 6
#
############################################################################################################################

#Read in data####
#h<-read.csv("../raw_data/NOAAHURDAT.csv")

h<-read.csv("../analysis_data/analysis_ready_binomial_k6.csv")
h <- h[h$lat>0,]
h$pressure = h$pressure*-1
#Create dataframe for annual max/min values
time = seq(min(h$yr), max(h$yr))

#Create pdf file to save image#####
pdf(file = "../figures/Figure4_6.pdf", width=7, height =6)
# set plot parameters
par(mfrow=c(3,1), mar=c(2,5,1,10), xpd=TRUE)


####### Windspeed panel#####
coef = rep(NA, 6)
pval = rep(NA, 6)
yr= seq(min(h$yr, na.rm=T), max(h$yr, na.rm=T))
i=0
m=lm(wind~yr, data=h[h$locGrp==i,])
# generate data for plotting trendline and prediction interval lines
coef[i+1] = print(summary(m)$coefficients[2,1])
pval[i+1] = print(summary(m)$coefficients[2,4])
d = as.data.frame(yr)
d$wind = rep(NA)
l = as.data.frame(predict(m, d, interval = "prediction"))

# plot raw data
plot(h[h$locGrp==i,"wind"]~h[h$locGrp==i,"yr"], type="p", col="white", 
     ylim=c(-10,160),
     xlim=c(1974,2022),
     ylab="averate TC windspeed \n(knots)",
     xlab="year")
#add trendline and prediction interval
points(l$fit~d$yr, type="l", lwd=2, col=cmap[i+1])
# add confidence interval
points(l$lwr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
points(l$upr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
# print lm coefficient and confidence interval in legend
#legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
print(summary(m)$coefficients[2,c(1,4)])


for(i in 1:5){
  m=lm(wind~yr, data=h[h$locGrp==i,])
  # generate data for plotting trendline and prediction interval lines
  coef[i+1] = print(summary(m)$coefficients[2,1])
  pval[i+1] = print(summary(m)$coefficients[2,4])
  d = as.data.frame(yr)
  d$wind = rep(NA)
  l = as.data.frame(predict(m, d, interval = "prediction"))
  
  # plot raw data
  #points(h[h$locGrp==i,"wind"]~h[h$locGrp==i,"yr"], type="p", col="white")
  #add trendline and prediction interval
  points(l$fit~d$yr, type="l", lwd=2, col=cmap[i+1])
  # add confidence interval
  points(l$lwr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
  points(l$upr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
  # print lm coefficient and confidence interval in legend
  #legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
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
legend(x = "topright",inset=c(-0.21,0), paste(expression(beta), "=", stats, sep=" "), col=cmap, lty=1)

####### Pressure panel#####
coef = rep(NA, 6)
pval = rep(NA, 6)
yr= seq(min(h$yr, na.rm=T), max(h$yr, na.rm=T))
i=0
m=lm(pressure~yr, data=h[h$locGrp==i,])
# generate data for plotting trendline and prediction interval lines
coef[i+1] = print(summary(m)$coefficients[2,1])
pval[i+1] = print(summary(m)$coefficients[2,4])
d = as.data.frame(yr)
d$pressure = rep(NA)
l = as.data.frame(predict(m, d, interval = "prediction"))

# plot raw data
plot(h[h$locGrp==i,"pressure"]~h[h$locGrp==i,"yr"], type="p", col="white", 
     ylim=c(910,1050),
     xlim=c(1974,2022),
     ylab="average TC pressure \n(kPa)",
     xlab="year")
#add trendline and prediction interval
points(l$fit~d$yr, type="l", lwd=2, col=cmap[i+1])
# add confidence interval
points(l$lwr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
points(l$upr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
# print lm coefficient and confidence interval in legend
print(summary(m)$coefficients[2,c(1,4)])


for(i in 1:5){
  m=lm(pressure~yr, data=h[h$locGrp==i,])
  # generate data for plotting trendline and prediction interval lines
  coef[i+1] = print(summary(m)$coefficients[2,1])
  pval[i+1] = print(summary(m)$coefficients[2,4])
  d = as.data.frame(yr)
  d$pressure = rep(NA)
  l = as.data.frame(predict(m, d, interval = "prediction"))
  
  # plot raw data
  #points(h[h$locGrp==i,"pressure"]~h[h$locGrp==i,"yr"], type="p", col="white")
  #add trendline and prediction interval
  points(l$fit~d$yr, type="l", lwd=2, col=cmap[i+1])
  # add confidence interval
  points(l$lwr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
  points(l$upr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
  # print lm coefficient and confidence interval in legend
  #legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
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
legend(x = "topright",inset=c(-0.21,0), paste(expression(beta), "=", stats, sep=" "), col=cmap, lty=1)
####### Latitude panel#####
coef = rep(NA, 6)
pval = rep(NA, 6)
int = rep(NA, 6)
yr= seq(min(h$yr, na.rm=T), max(h$yr, na.rm=T))
i=0
m=lm(minlat~yr, data=h[h$locGrp==i,])
# generate data for plotting trendline and prediction interval lines
coef[i+1] = print(summary(m)$coefficients[2,1])
pval[i+1] = print(summary(m)$coefficients[2,4])
d = as.data.frame(yr)
d$minlat = rep(NA)
l = as.data.frame(predict(m, d, interval = "prediction"))

# plot raw data
plot(h[h$locGrp==i,"minlat"]~h[h$locGrp==i,"yr"], type="p", col="white", 
     ylim=c(0,45),
     xlim=c(1974,2022),
     ylab="average TC starting latitude \n(degrees)",
     xlab="year")
#add trendline and prediction interval
points(l$fit~d$yr, type="l", lwd=2, col=cmap[i+1])
# add confidence interval
points(l$lwr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
points(l$upr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
# print lm coefficient and confidence interval in legend
#legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
print(summary(m)$coefficients[2,c(1,4)])


for(i in 1:5){
  m=lm(minlat~yr, data=h[h$locGrp==i,])
  # generate data for plotting trendline and prediction interval lines
  coef[i+1] = print(summary(m)$coefficients[2,1])
  pval[i+1] = print(summary(m)$coefficients[2,4])
  d = as.data.frame(yr)
  d$minlat = rep(NA)
  l = as.data.frame(predict(m, d, interval = "prediction"))
  
  # plot raw data
  #points(h[h$locGrp==i,"minlat"]~h[h$locGrp==i,"yr"], type="p", col="white")
  #add trendline and prediction interval
  points(l$fit~d$yr, type="l", lwd=2, col=cmap[i+1])
  # add confidence interval
  points(l$lwr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
  points(l$upr~d$yr, type="l", lwd=1, lty=2, col=cmap[i+1])
  # print lm coefficient and confidence interval in legend
  #legend(x="topleft", legend=c(expression(paste(beta, " =")), paste(as.character(round(summary(m1)$coefficients[2,1], digits=3)), "*", sep="")))
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
legend(x = "topright",inset=c(-0.21,0), paste(expression(beta), "=", stats, sep=" "), col=cmap, lty=1)

#####
dev.off()

############################################################################################################################
#
#                     Figure 5 
#
############################################################################################################################
rm(list=ls())
load("./model_outputs/hur_trends.RData")
load("./model_outputs/model_outputs.RData")
h<-read.csv("./raw_data/NOAAHURDAT.csv")
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
locGrps = c("High Plains","Northeast", "MI Delta", "Southeast", "Midwest", "Texas Coast")

# Historical statistics 
u_ws <- u["wind"]
sd_ws <- sd["wind"]
u_z <- u["pressure"]
sd_z <- sd["pressure"]
u_lat<- u["minlat"]
sd_lat<- sd["minlat"]

#
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
spred$lat = rep(0,3)#(pred$lat - u_lat)/sd_lat
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
# 
# hur_cat$minlat_min = rep(0)
# # hur_cat$pressure_min = rep(0)
# for(i in 1:length(wind_max)){
# #   hur_cat$minlat_min[i] = mean(h$lat[h$wind > wind_max[i]-4 & wind_max[i]+4])
# # hur_cat$pressure_min[i] = mean(h$pressure[h$wind > wind_max[i]-4 & wind_max[i]+4])
# }
# hur_cat$wind_max = (hur_cat$wind_max - u_ws)/sd_ws
# #hur_cat$pressure_min = (hur_cat$pressure_min - u_z)/40#sd_z
# hur_cat$pressure_min = (hur_cat$pressure_min - u_z)/sd_z
# hur_cat$minlat_min = (hur_cat$minlat_min - u_lat)/sd_lat
# 
# #hur_cat$pressure_min = 1
# hur_cat$minlat_min = -4
# 
# 
# PC_hurcat <-as.data.frame(predict(pca, newdata=hur_cat))

yr_s=2014
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
d <- d[order(d$HurComp),]
par(mfrow=c(1,1))
#pdf(file = "../figures/Figure5.pdf", width=6, height =6)
svg(file = "./figures/Figure5_2.svg", width=6, height =6)

m1<-ms1
m2<-ms2
#plot data####
HurComp<-seq(min(d$HurComp), 7, by = 0.05)
plot(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==0], locGrp=0, yr=yr_s),
                 na.action=na.exclude))~d$HurComp[d$locGrp==0], 
     type="l",
     lwd=4,
     col=cmap[1],
     xlim=c(-2.2, 7), 
     ylim=c(0, 30), 
     ylab="Poisson predicted fs'", xlab="TC intensity")

for(i in 0:5){
p4p<-predict(m2, data.frame(HurComp=HurComp, locGrp=i, yr=yr_s))
se<-as.data.frame(emmeans(m2, ~HurComp, at = list(x = HurComp, locGrp=i, yr=yr_s)))$SE

polygon(x = c(HurComp, rev(HurComp)),
        y = c(exp(p4p) - 1.96*exp(se),
              rev(exp(p4p) +1.96*exp(se))),
        col =  adjustcolor(cmap[i+1], alpha.f = 0.1), border = NA)


}
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==1], locGrp=1, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==1],type="l", lwd=4,col=cmap[2])
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==2], locGrp=2, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==2],type="l", lwd=4,col=cmap[3])
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==3], locGrp=3, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==3],type="l", lwd=4,col=cmap[4])
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==4], locGrp=4, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==4],type="l", lwd=6,col=cmap[5])
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==5], locGrp=5, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==5],type="l", lwd=4,col=cmap[6])

points(exp(predict(m2, data.frame(HurComp=HurComp, locGrp=0, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[1])
points(exp(predict(m2, data.frame(HurComp=HurComp, locGrp=1, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[2])
points(exp(predict(m2, data.frame(HurComp=HurComp, locGrp=2, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[3])
points(exp(predict(m2, data.frame(HurComp=HurComp, locGrp=3, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[4])
points(exp(predict(m2, data.frame(HurComp=HurComp, locGrp=4, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[5])
points(exp(predict(m2, data.frame(HurComp=HurComp, locGrp=5, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[6])

points(exp(predict(m2, data.frame(HurComp=PC_pred[1,1], locGrp=0, yr=yr_s),na.action=na.exclude))~PC_pred[1,1],type="p", pch = 25, cex=1.5,bg=cmap[1])
points(exp(predict(m2, data.frame(HurComp=PC_pred[1,1], locGrp=1, yr=yr_s),na.action=na.exclude))~PC_pred[1,1],type="p", pch = 25, cex=1.5,bg=cmap[2])
points(exp(predict(m2, data.frame(HurComp=PC_pred[1,1], locGrp=2, yr=yr_s),na.action=na.exclude))~PC_pred[1,1],type="p", pch = 25, cex=1.5,bg=cmap[3])
points(exp(predict(m2, data.frame(HurComp=PC_pred[1,1], locGrp=3, yr=yr_s),na.action=na.exclude))~PC_pred[1,1],type="p", pch = 25, cex=1.5,bg=cmap[4])
points(exp(predict(m2, data.frame(HurComp=PC_pred[1,1], locGrp=4, yr=yr_s),na.action=na.exclude))~PC_pred[1,1],type="p", pch = 25, cex=1.5,bg=cmap[5])
points(exp(predict(m2, data.frame(HurComp=PC_pred[1,1], locGrp=5, yr=yr_s),na.action=na.exclude))~PC_pred[1,1],type="p", pch = 25, cex=1.5,bg=cmap[6])

points(exp(predict(m2, data.frame(HurComp=PC_pred[2,1], locGrp=0, yr=yr_s),na.action=na.exclude))~PC_pred[2,1],type="p", pch = 23, cex=1.5,bg=cmap[1])
points(exp(predict(m2, data.frame(HurComp=PC_pred[2,1], locGrp=1, yr=yr_s),na.action=na.exclude))~PC_pred[2,1],type="p", pch = 23, cex=1.5,bg=cmap[2])
points(exp(predict(m2, data.frame(HurComp=PC_pred[2,1], locGrp=2, yr=yr_s),na.action=na.exclude))~PC_pred[2,1],type="p", pch = 23, cex=1.5,bg=cmap[3])
points(exp(predict(m2, data.frame(HurComp=PC_pred[2,1], locGrp=3, yr=yr_s),na.action=na.exclude))~PC_pred[2,1],type="p", pch = 23, cex=1.5,bg=cmap[4])
points(exp(predict(m2, data.frame(HurComp=PC_pred[2,1], locGrp=4, yr=yr_s),na.action=na.exclude))~PC_pred[2,1],type="p", pch = 23, cex=1.5,bg=cmap[5])
points(exp(predict(m2, data.frame(HurComp=PC_pred[2,1], locGrp=5, yr=yr_s),na.action=na.exclude))~PC_pred[2,1],type="p", pch = 23, cex=1.5,bg=cmap[6])

points(exp(predict(m2, data.frame(HurComp=PC_pred[3,1], locGrp=0, yr=yr_s),na.action=na.exclude))~PC_pred[3,1],type="p", pch = 24, cex=1.5,bg=cmap[1])
points(exp(predict(m2, data.frame(HurComp=PC_pred[3,1], locGrp=1, yr=yr_s),na.action=na.exclude))~PC_pred[3,1],type="p", pch = 24, cex=1.5,bg=cmap[2])
points(exp(predict(m2, data.frame(HurComp=PC_pred[3,1], locGrp=2, yr=yr_s),na.action=na.exclude))~PC_pred[3,1],type="p", pch = 24, cex=1.5,bg=cmap[3])
points(exp(predict(m2, data.frame(HurComp=PC_pred[3,1], locGrp=3, yr=yr_s),na.action=na.exclude))~PC_pred[3,1],type="p", pch = 24, cex=1.5,bg=cmap[4])
points(exp(predict(m2, data.frame(HurComp=PC_pred[3,1], locGrp=4, yr=yr_s),na.action=na.exclude))~PC_pred[3,1],type="p", pch = 24, cex=1.5,bg=cmap[5])
points(exp(predict(m2, data.frame(HurComp=PC_pred[3,1], locGrp=5, yr=yr_s),na.action=na.exclude))~PC_pred[3,1],type="p", pch = 24, cex=1.5,bg=cmap[6])

abline(v= PC_hurcat$Comp.1[1], lty =2, lwd=2, col="darkgray")#((33 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= PC_hurcat$Comp.1[2], lty =2, lwd=2, col="darkgray")#((63 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= PC_hurcat$Comp.1[3], lty =2, lwd=2, col="darkgray")#((82 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= PC_hurcat$Comp.1[4], lty =2, lwd=2, col="darkgray")#((95 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= PC_hurcat$Comp.1[5], lty =2, lwd=2, col="darkgray")#((112 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= PC_hurcat$Comp.1[6], lty =2, lwd=2, col="darkgray")#((136 - historical_mean_windspeed) / historical_sd_windspeed))

txt_y=rep(50, 7)
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

print(exp(predict(m2, data.frame(HurComp=PC_pred[,1], locGrp=0, yr=yr_s),na.action=na.exclude)))
print(exp(predict(m2, data.frame(HurComp=PC_pred[,1], locGrp=1, yr=yr_s),na.action=na.exclude)))
print(exp(predict(m2, data.frame(HurComp=PC_pred[,1], locGrp=2, yr=yr_s),na.action=na.exclude)))
print(exp(predict(m2, data.frame(HurComp=PC_pred[,1], locGrp=3, yr=yr_s),na.action=na.exclude)))
print(exp(predict(m2, data.frame(HurComp=PC_pred[,1], locGrp=4, yr=yr_s),na.action=na.exclude)))
print(exp(predict(m2, data.frame(HurComp=PC_pred[,1], locGrp=5, yr=yr_s),na.action=na.exclude)))

#####
dev.off()

#################################################################################
#
#       Fig 5_2
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

#
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
#pdf(file = "../figures/Figure5.pdf", width=6, height =6)
svg(file = "./figures/Figure5_2.svg", width=6, height =6)

m1<-ms1
m2<-ms2
#plot data####
HurComp<-seq(min(d$HurComp), 7, by = 0.05)
plot(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==0], locGrp=0, yr=yr_s),
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
  p4p<-predict(m2, data.frame(HurComp=HurComp, locGrp=i, yr=yr_s))
  se<-as.data.frame(emmeans(m2, ~HurComp, at = list(x = HurComp, locGrp=i, yr=yr_s)))$SE
  
  polygon(x = c(HurComp, rev(HurComp)),
          y = c(exp(p4p) - 1.96*exp(se),
                rev(exp(p4p) +1.96*exp(se))),
          col =  adjustcolor(cmap[i+1], alpha.f = 0.1), border = NA)
  
  
}
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==1], locGrp=1, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==1],type="l", lwd=4,col=cmap[2])
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==2], locGrp=2, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==2],type="l", lwd=4,col=cmap[3])
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==3], locGrp=3, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==3],type="l", lwd=4,col=cmap[4])
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==4], locGrp=4, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==4],type="l", lwd=6,col=cmap[5])
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==5], locGrp=5, yr=yr_s),na.action=na.exclude))~d$HurComp[d$locGrp==5],type="l", lwd=4,col=cmap[6])

points(exp(predict(m2, data.frame(HurComp=HurComp, locGrp=0, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[1])
points(exp(predict(m2, data.frame(HurComp=HurComp, locGrp=1, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[2])
points(exp(predict(m2, data.frame(HurComp=HurComp, locGrp=2, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[3])
points(exp(predict(m2, data.frame(HurComp=HurComp, locGrp=3, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[4])
points(exp(predict(m2, data.frame(HurComp=HurComp, locGrp=4, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[5])
points(exp(predict(m2, data.frame(HurComp=HurComp, locGrp=5, yr=yr_s),na.action=na.exclude))~HurComp,type="l", lwd=2,lty=3,col=cmap[6])

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
  print(exp(predict(m2, data.frame(HurComp=p, locGrp=i, yr=yr_s),na.action=na.exclude)))
  points(exp(predict(m2, data.frame(HurComp=p[1], locGrp=i, yr=yr_s),na.action=na.exclude))~p[1],type="p", pch = 25, cex=1.5,bg=cmap[i+1])
  points(exp(predict(m2, data.frame(HurComp=p[2], locGrp=i, yr=yr_s),na.action=na.exclude))~p[2],type="p", pch = 23, cex=1.5,bg=cmap[i+1])
  points(exp(predict(m2, data.frame(HurComp=p[3], locGrp=i, yr=yr_s),na.action=na.exclude))~p[3],type="p", pch = 24, cex=1.5,bg=cmap[i+1])

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


#####
dev.off()



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

###################################################################################################################
#
#                 HurInt by Regions
#
##################################################################################################################
cmap = c("#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#bcbd22", "#17becf")
locGrps = c("High Plains","Northeast", "MI Delta", "Southeast", "Midwest", "Texas Coast")
h<-read.csv("./analysis_data/analysis_ready_binomial_k6.csv")
h <- h[h$lat>0,]
h$pressure = h$pressure*-1
#Create dataframe for annual max/min values
time = seq(min(h$yr), max(h$yr))

plot(h$lat~h$lon, col="white")
points(h$lat[h$locGrp==0]~h$lon[h$locGrp==0], col=cmap[1])
points(h$lat[h$locGrp==1]~h$lon[h$locGrp==1], col=cmap[2])
points(h$lat[h$locGrp==2]~h$lon[h$locGrp==2], col=cmap[3])
points(h$lat[h$locGrp==3]~h$lon[h$locGrp==3], col=cmap[4])
points(h$lat[h$locGrp==4]~h$lon[h$locGrp==4], col=cmap[5])
points(h$lat[h$locGrp==5]~h$lon[h$locGrp==5], col=cmap[6])

#####################################################################################################################################################
#
#                   End script (clean up below)
#
###################################################################################################################




# Ensure locGrp in newdata matches levels in original model
levels(newdata$locGrp) <- levels(d$locGrp)  

# # Make predictions using the binomial model
# binom_predictions <- predict(mb2, newdata = newdata,  type = "response", allow.new.levels = TRUE)
# 
# # Make predictions using the Poisson model
# poisson_predictions <- exp(predict(m2, newdata = newdata,  type = "response", allow.new.levels = TRUE))
# 
# # Print the predictions
# cat("Binomial model predictions (probability of failure):\n")
# print(binom_predictions)
# 
# cat("\nPoisson model predictions (expected number of failures):\n")
# print(poisson_predictions)
# 
# # Plotting the predictions
# plot(exp(predict(m2, newdata = newdata, allow.new.levels = TRUE)), type = "p", col = "blue", pch = 19,
#      xlab = "Location Group", ylab = "Predicted Failures",
#      main = "Predicted Pipeline Failures by Location Group (2050)")
# 

#explore predictions
plot(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==0], locGrp=0, yr=2014),
                 na.action=na.exclude))~d$wind_max[d$locGrp==0], 
     xlim=c(10, 150), 
     ylim=c(0, 25), 
     ylab="Poisson predicted fs'", xlab="Windspeed (m/s)")
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==1], locGrp=1, yr=2014),na.action=na.exclude))~d$wind_max[d$locGrp==1], xlim=c(-2.5, 5.5), ylim=c(0, 25), col="blue")
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==2], locGrp=2, yr=2014),na.action=na.exclude))~d$wind_max[d$locGrp==2], xlim=c(-2.5, 5.5), ylim=c(0, 25), col="red")
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==3], locGrp=3, yr=2014),na.action=na.exclude))~d$wind_max[d$locGrp==3], xlim=c(-2.5, 5.5), ylim=c(0, 25), col="green")
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==4], locGrp=4, yr=2014),na.action=na.exclude))~d$wind_max[d$locGrp==4], xlim=c(-2.5, 5.5), ylim=c(0, 25), col="yellow")
points(exp(predict(m2, data.frame(HurComp=d$HurComp[d$locGrp==5], locGrp=5, yr=2014),na.action=na.exclude))~d$wind_max[d$locGrp==5], xlim=c(-2.5, 5.5), ylim=c(0, 25), col="purple")

abline(v= 33)#((33 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= 63)#((63 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= 82)#((82 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= 95)#((95 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= 112)#((112 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= 136)#((136 - historical_mean_windspeed) / historical_sd_windspeed))


plot(d$failTF_sum[d$locGrp==0]~d$wind_max[d$locGrp==0], 
     col=cmap[1],
     xlim=c(10, 150), 
     ylim=c(0, 25), 
     ylab="Poisson predicted fs'", xlab="Windspeed (m/s)")
points(d$failTF_sum[d$locGrp==1]~d$wind_max[d$locGrp==1], col=cmap[2])
points(d$failTF_sum[d$locGrp==2]~d$wind_max[d$locGrp==2], col=cmap[3])
points(d$failTF_sum[d$locGrp==3]~d$wind_max[d$locGrp==3], col=cmap[4])
points(d$failTF_sum[d$locGrp==4]~d$wind_max[d$locGrp==4], col=cmap[5])
points(d$failTF_sum[d$locGrp==5]~d$wind_max[d$locGrp==5], col=cmap[6])

abline(v= 33)#((33 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= 63)#((63 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= 82)#((82 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= 95)#((95 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= 112)#((112 - historical_mean_windspeed) / historical_sd_windspeed))
abline(v= 136)#((136 - historical_mean_windspeed) / historical_sd_windspeed))
