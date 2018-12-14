library(ggplot2)
library(cowplot)
library(utils)
library(geosphere)
library(sp)
library(rgeos)
library(downloader)
library(maptools)
library(mapmisc)
library(sgeostat)
library(psych)

# read in data
dstate=read.csv("data/dstate.csv")
bbs=read.csv("data/BBS_Annual_Indices_Estimates_2015_7-29-2016.csv") # read in BBS data
colnames(bbs)<-c("sp","state","year","ind","cred")
head(bbs)

g=subset(bbs,sp=="s05460"); head(g)
levels(bbs$state)
gg=subset(g,state %in% c("CON", "NH ", "MIS", "NY ","PA ","NJ ","DEL","MD ","WV ",
                         "VA ","NC ","KY ","TEN", "OHI", "IND", "ILL","MIC", 
                         "ALA", "GA ", "SC ", "WIS")) # 


head(gg); tail(gg)
gg$r = 100*c(NA,diff(gg$ind))/gg$ind
gg=droplevels(gg)
head(gg,10)
gg1=subset(gg,year<1991); head(gg1)
gg2=subset(gg,year>1990); head(gg2)


# checking data graphically and average index values using tapply
ggplot(subset(gg,state %in% c("ILL","IND"))) + geom_line(aes(x=year,y=ind,color=state),size=2)
ggplot(subset(gg,year!=1966 & state %in% c("ILL","IND"))) + geom_line(aes(x=year,y=r,color=state),size=2)
tapply(gg$ind,gg$state,mean, na.rm=T) 
# MA, RI, VT, ME have no data data; exclude CT, NH, and MS later due to small or 0 pops
# better way to make it interchangeable with any species is to merge with state name lookup table
# need to make that

# reshaping data to wide format for correlations
g.r.matrix=reshape(subset(gg,year!=1966)[c("state","year", "r")],idvar="year",timevar = "state", direction = "wide")
colnames(g.r.matrix) <- c("year","AL","CT","DE","GA","IL","IN","KY","MD",
                    "MI","MS","NC","NH","NJ","NY","OH","PA","SC","TN","VA","WI","WV")
head(g.r.matrix); tail(g.r.matrix); ncol(g.r.matrix)
g.r.matrix$CT[41]<-NaN

g.test=corr.test(g.r.matrix[c("AL","CT","DE","GA","IL","IN","KY","MD",
                              "MI","MS","NC","NH","NJ","NY","OH","PA","SC","TN","VA","WI","WV")], adjust = "none")
g.cor = data.frame(reshape2::melt(g.test$r), reshape2::melt(g.test$p)[3]); head(g.cor); nrow(g.cor)
colnames(g.cor) = c("state1","state2","cor","pval"); head(g.cor)
g.cor$pairid = paste(as.character(g.cor$state1),"-",as.character(g.cor$state2),sep=""); head(g.cor); max(g.cor$cor);nrow(g.cor)
g.cor=g.cor[g.cor$cor!=1,]; head(g.cor); nrow(g.cor)
g.cor[230:235,]


#### combining GRSP correlation and interstate distances
head(g.cor); nrow(g.cor)
head(dstate); nrow(dstate)

gcordist = merge(g.cor,dstate,by='pairid',all.x=T); head(gcordist); nrow(gcordist)
colnames(gcordist)[c(2,3)] <- c("state1","state2"); head(gcordist)
gcordist = gcordist[,c("pairid","state1","state2","cor","pval","dist")]; head(gcordist)
gcordist$sig = ifelse(gcordist$pval<0.05,"Y","N"); head(gcordist)

#remove three states with very low populations
gcordist= subset(gcordist,!(state1 %in% c("CT","MS","NH")) & !(state2 %in% c("CT","MS","NH")))

### SCATTERPLOT OF SYNCHRONY VS. DISTANCE
X11(13,9)
ggplot(subset(gcordist,!(state1 %in% c("CT","MS","NH")) & !(state2 %in% c("CT","MS","NH")))) + 
  geom_point(aes(x=dist,y=cor,color=sig),size=4) + 
  geom_smooth(aes(x=dist,y=cor),method="lm",color="black") 
#ggsave("GRSP_cor_1966-2015.png")

# with significant pairs labeled
ggplot(subset(gcordist,!(state1 %in% c("CT","MS","NH")) & !(state2 %in% c("CT","MS","NH")))) + 
  geom_point(aes(x=dist,y=cor,color=sig),size=4) + 
  geom_smooth(aes(x=dist,y=cor),method="lm",color="black") +
  geom_label(data=subset(gcordist,sig=="Y"),aes(x=dist,y=cor,label=pairid))

###
### NOTE PAIRS ARE IN DUPLICATE AT THE MOMENT.  FIGURE OUT HOW TO REMOVE LATER!!!!
### % SIGNIFICANT IN BINS WILL NOT CHANGE, BUT SAMPLE SIZE IS 1/2 OF WHAT IT APPEARS

###### plotting % signficant by distance bins 
# (to improve: put distance bin into dstate.csv to make easier to switch spp)

# create 200m bins
binlabels1=c("100-200","201-400","401-600","601-800","801-1000",
                   "1001-1200","1201-1400","1401-1600","1601-1800","1801-2000")
binlabels2=c("100-200","201-300","301-400","401-500","501-600","601-700",
             "701-800","801-900","901-1000",">1000")
binlabels3=c("100-300","301-500","501-700","701-900","901-1100","1101-1300",
             "1301-1500",">1500")
gcordist$dcat1=cut(gcordist$dist,c(0,200,400,600,800,1000,1200,1400,1600,1800,2000),
    labels=binlabels1)
gcordist$dcat2=cut(gcordist$dist,c(0,200,300,400,500,600,700,800,900,1000,2000),
                   labels=binlabels2)
gcordist$dcat3=cut(gcordist$dist,c(0,300,500,700,900,1100,1300,1500,2000),
                   labels=binlabels3)
gcordist$sig1=ifelse(gcordist$sig=="Y",1,0)
head(gcordist)

gbin1 = cbind.data.frame(pctsig=100*tapply(gcordist$sig1,gcordist$dcat1,mean),
      n=tapply(gcordist$sig1,gcordist$dcat1,length))
gbin1$bin=factor(rownames(gbin1),levels=binlabels1); gbin1

gbin2 = cbind.data.frame(pctsig=100*tapply(gcordist$sig1,gcordist$dcat2,mean),
      n=tapply(gcordist$sig1,gcordist$dcat2,length))
gbin2$bin=factor(rownames(gbin2),levels=binlabels2); gbin2

gbin3 = cbind.data.frame(pctsig=100*tapply(gcordist$sig1,gcordist$dcat3,mean),
      n=tapply(gcordist$sig1,gcordist$dcat3,length))
gbin3$bin=factor(rownames(gbin3),levels=binlabels3); gbin3

x11(13,8.5)
ggplot(gbin3) + geom_col(aes(y=pctsig,x=bin)) + labs(x="distance (km)",y="% signifcant correlations")
#ggsave("GRSP_bincor_1966-2015.png")


## NEXT: need to create n = 1000 null (random) distributions for each bin like Martin et al. 
## to test spatial extent of synchrony



# nex: correlate with distance and repeat for historic vs. recent periods
# note: try correlating average cor with yield? avg growth rate?
# do it separately for early and late period
#write.csv(data.frame(tapply(ggall$cor,ggall$state2,mean)),"grsp.mean.cor.1967-2015.csv") # changed name of gg above to get early and late

# relating historic vs. recent period spatial synchrony
late=read.csv("grsp.mean.cor.1991-2015.csv"); colnames(late)=c("state","late.cor"); late
early=read.csv("grsp.mean.cor.1967-1990.csv"); colnames(early)=c("state","early.cor"); early
late.early=cbind.data.frame(late,early)[c(1,2,4)]
ggplot(late.early) + geom_point(aes(x=early.cor, y=late.cor,color=state),size=4) + geom_abline(slope=1,intercept=0)

### MAKE A GRAPH SHOWING A SERIES OF LINES CONNECTING HISTORIC/RECENT BY STATE
late.early.melt=reshape2::melt(late.early)
ggplot(late.early.melt) + geom_line(aes(x=state,y=value)) # not working