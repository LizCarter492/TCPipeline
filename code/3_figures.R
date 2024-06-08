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
#TBD