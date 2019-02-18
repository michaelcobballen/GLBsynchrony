####
#### WRITE CODE TO EVALUATE SPATIAL SYNCHRONY WHILE WAITING FOR CORRECT DATA (BBSBAYES OR CAR MODEL)
####

library(dplyr)
library(ggplot2)
library(psych)

g = read.csv("data/dgrid.csv")

# make a data frame of all possible grid/year combinations
grid.list = unique(grid.year$grid)
year.list = 1966:2017
grid.year.list = expand.grid(grid = grid.list,year = year.list) %>%
  arrange(grid,year) %>%
  mutate(grid.year = paste(grid,year,sep="."))

# attach the bird data to the complete list and format for correlation analysis
grid.year.full = grid.year.list %>%
  left_join(grid.year, by = "grid.year") %>%
  select(grid=grid.x, year=year.x, grid.year, mean.count) %>%
  mutate(log.count = log(mean.count+1), logcountlag1 = lag(log.count)) %>% # log+1 transform count ala Hanski and create lag for AR and differencing analyses
  mutate(logcountlag1 = as.numeric(ifelse(year==1966,"NA",logcountlag1))) %>% # zero out 1966 which should have no lagged value
  filter(is.na(log.count) == F & is.na(logcountlag1) == F) %>% # temporarily remove NAs to get residuals and first diffs of logged counts
  mutate(change = log.count - logcountlag1) %>% # first diff of logged counts, ~ % change
  group_by(grid) %>%
  mutate(AR1.R = resid(lm(log.count~logcountlag1))) %>% # calculate residuals from first order autoregression (Hanski-style)
  ungroup() %>%
  right_join(grid.year.list, by = "grid.year") %>% # merge back with the full grid/year list to allow matching in correlation analyses
  select(grid = grid.y, year = year.y, grid.year, mean.count, logcountlag1, change, AR1.R)

# could filter bird data grid cells here, e.g., 
# filter(count>0) %>%  # if using % difference method (growth rate)
# could then exclude grid cells with fewer than 10 year-pairs? if use % difference model
# otherwise try 1st order autoregressive model like Hanski etc., Michel?
# exclude cells with > XX years below a certain pop size? 10 years of zeros?

#gridbirdsynch = function(yearstart, yearend) {
yearstart = 1966; yearend = 2017  # for testing the function

grid.bird.mat = 
  grid.year.full %>%
  dplyr::select(grid,year,AR1.R) %>%
  spread(grid,AR1.R) %>%
  filter(year > yearstart-1 & year < yearend+1)

grid.bird.mat = grid.bird.mat[,2:ncol(grid.bird.mat)]


# create matrix of Pearson correlations among states
grid.bird.test = 
  grid.bird.mat %>%
  corr.test(adjust="none")

grid.bird.cor = 
  data.frame(reshape2::melt(grid.bird.test$r), 
             reshape2::melt(grid.bird.test$p)[3],
             reshape2::melt(grid.bird.test$n)[3]) %>%
  select(grid1=Var1, grid2=Var2, cor=value, pval=value.1, n=value.2) %>%
  mutate(pairid = paste(grid1,".",grid2,sep="")) %>%
  filter(n>19, grid1!=grid2) %>% # excludes grid pairs with less than 7 years in common (Hanski's criteria)
  filter(grid1!=grid2) # excludes self-correlations
  #filter(is.na(pval)==F) # excludes NA correlations resulting from too many zeros - should be removed earlier via some form of rule (ie. < XX 0's)
  #distinct(pval, .keep_all = TRUE) # should cut number in half, but no...
  # still need to figure out why some cors are NA; seems to be because of zeros...

grid.bird.cordist = grid.bird.cor %>%
  left_join(dgrid,by='pairid') %>%
  select(pairid, grid1=grid1.x, grid2=grid2.x, cor,pval,n,dist,dcat) %>%
  mutate(sig = ifelse(pval < 0.05,"Y","N"), sig1 = as.numeric(ifelse(pval < 0.05,"1","0"))) %>%
  distinct(cor, .keep_all = TRUE)
# note: distinct' above excluded also the ones with NAs; need to understand why they were NAs eventually

#return(grid.bird.cordist)
#}
 
X11(13,9)
ggplot(grid.bird.cordist) + 
#  geom_point(aes(x=dist,y=cor,color=sig),size=0.1,alpha=0.5) +
 geom_smooth(aes(x=dist,y=cor),method="loess",color="black", se=F) 
#ggsave("figures/grid_sync_scatter/GRSP_cor_grid_1967-2017_24n+_lineonly.jpg")


######
###### plotting % signficant correlations by distance bins 
######

binlabels = c("<100","100-300","301-500","501-700","701-900",
              "901-1100","1101-1300","1301-1500",">1500")
corbin = cbind.data.frame(pctsig=100*tapply(grid.bird.cordist$sig1,grid.bird.cordist$dcat,mean,na.rm=T),
                          n=tapply(grid.bird.cordist$sig1,grid.bird.cordist$dcat,length))
corbin$bin=factor(rownames(corbin),levels=binlabels)

#x11(13,8.5)
ggplot(corbin) + 
  geom_col(aes(y=pctsig,x=bin)) + 
  labs(x="distance (km)",y="% of correlations with P < 0.05") #+
#  ylim(c(0,35))
#ggsave("figures/bin_cor/zBOBO_bincor_all_1967-2015.png")
#ggsave("figures/bin_cor/BOBO_bincor_all_1967-1990.png")
#ggsave("figures/bin_cor/BOBO_bincor_all_1991-2015.png")


## NEXT: need to create n = 1000 null (random) distributions for each bin like Martin et al. 
## to test spatial extent of synchrony


