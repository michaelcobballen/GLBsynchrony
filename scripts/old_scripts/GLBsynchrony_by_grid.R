####
#### WRITE CODE TO EVALUATE SPATIAL SYNCHRONY WHILE WAITING FOR CORRECT DATA (BBSBAYES OR CAR MODEL)
####

library(dplyr)
library(ggplot2)
library(psych)
library(tidyr)

g = read.csv("data/dgrid.92.17.csv")

# make a data frame of all possible grid/year combinations
grid.list.92.17 = unique(grid.year.92.17$grid)
year.list.92.17 = 1992:2017
grid.year.list.92.17 = expand.grid(grid = grid.list.92.17,year = year.list.92.17) %>%
  arrange(grid,year) %>%
  mutate(grid.year = paste(grid,year,sep="."))

# attach the bird data to the complete list and format for correlation analysis
grid.year.full.92.17 = grid.year.list.92.17 %>%
  left_join(grid.year.92.17, by = "grid.year") %>%
  select(grid=grid.x, year=year.x, grid.year, mean.count) %>%
  mutate(log.count = log(mean.count+1), logcountlag1 = lag(log.count)) %>% # log+1 transform count ala Hanski and create lag for AR and differencing analyses
  mutate(logcountlag1 = as.numeric(ifelse(year==1992,"NA",logcountlag1))) %>% # assign "NA" to 1992 which should have no lagged value
  filter(is.na(log.count) == F & is.na(logcountlag1) == F) %>% # temporarily remove NAs to get residuals and first diffs of logged counts
  mutate(change = log.count - logcountlag1) %>% # first diff of logged counts, ~ % change
  group_by(grid) %>%
  mutate(AR1.R = resid(lm(log.count~logcountlag1))) %>% # calculate residuals from first order autoregression (Hanski-style)
  ungroup() %>%
  right_join(grid.year.list.92.17, by = "grid.year") %>% # merge back with the full grid/year list to allow matching in correlation analyses
  select(grid = grid.y, year = year.y, grid.year, mean.count, log.count,logcountlag1, change, AR1.R)

### make a correlation matrix of bird pop data
grid.bird.mat.92.17 = 
  grid.year.full.92.17 %>%
  dplyr::select(grid,year,AR1.R) %>%
  spread(grid,AR1.R) %>%
  filter(year > 1992 & year < 2018)

grid.bird.mat.92.17 = grid.bird.mat.92.17[,2:ncol(grid.bird.mat.92.17)]

# create matrix of Pearson correlations among states
grid.bird.test.92.17 = 
  grid.bird.mat.92.17 %>%
  corr.test(adjust="none")

grid.bird.cor.92.17 = 
  data.frame(reshape2::melt(grid.bird.test.92.17$r), 
             reshape2::melt(grid.bird.test.92.17$p)[3],
             reshape2::melt(grid.bird.test.92.17$n)[3]) %>%
  select(gridA=Var1, gridB=Var2, cor=value, pval=value.1, n=value.2) %>%
  mutate(pairid = paste(gridA,".",gridB,sep="")) %>%
  filter(n>19) %>% # excludes grid pairs with less than XX years in common (Hanski's criteria was 7)
  filter(gridA!=gridB) %>% # excludes self-correlations
  filter(is.na(pval)==F) %>% # excludes NA correlations resulting from too many zeros if these weren't removed earlier - they were! So nothing changes.
  distinct(pval, .keep_all = TRUE) # should cut number in half, and yes it does!

grid.bird.cordist.92.17 = grid.bird.cor.92.17 %>%
  left_join(g,by='pairid') %>%
  select(pairid, gridA=gridA.x, gridB=gridB.x, cor,pval,n,dist) %>%
  mutate(sig = ifelse(pval < 0.05,"Y","N"), sig1 = as.numeric(ifelse(pval < 0.05,"1","0"))) %>%
  distinct(cor, .keep_all = TRUE) # not sure why I have this here, but it doesn't seem to change anything

#######################################
###### plotting correlations by distance - continuous 
#######################################

X11(13,9)
ggplot(grid.bird.cordist.92.17) + 
#  geom_point(aes(x=dist,y=cor,color=sig),size=0.1,alpha=0.5) +
  geom_smooth(aes(x=dist,y=cor),method="auto",color="black", se=T) +
  xlim(c(0,2000))
#ggsave("figures/grid_sync_scatter/GRSP_cor_grid_92.17_lines_AR_5zeros_n20.jpg")
    # suggestion to use method="auto" as "loess" didn't have enough memory
    # https://groups.google.com/forum/#!topic/ggplot2/enavD18MmkY


#######################################
###### plotting % signficant correlations by distance bins 
#######################################

binlabels = c("<100","101-200","201-300","301-400","401-500","501-600",
              "601-700","701-800","801-900","901-1000","1001-1100","1101-1200",
              "1201-1300","1301-1400","1401-1500","1501-1600","1601-1700",
              "1701-1800","1801-1900","1901-2000",">2000")
binlabels2 = c("<200","201-400","401-600","601-800","801-1000","1001-1200",
              "1201-1400","1401-1600","1601-1800","1801-2000",">2000")

grid.bird.cordist.92.17 = grid.bird.cordist.92.17 %>%
    mutate(dcat=cut(grid.bird.cordist.92.17$dist,c(0,100,200,300,400,500,600,700,800,900,
                                        1000,1100,1200,1300,1400,1500,1600,1700,
                                        1800,1900,2000,100000),
                    labels=binlabels)) %>%
    mutate(dcat2 = cut(grid.bird.cordist.92.17$dist,c(0,200,400,600,800,1000,
                                        1200,1400,1600,1800,2000,100000),
                     labels=binlabels2))
    


corbin1 = grid.bird.cordist.92.17 %>%
  group_by(dcat) %>%
  summarise(pctsig = mean(sig1)*100, n = length(sig1)) %>%
  ungroup()
  
corbin2 = grid.bird.cordist.92.17 %>%
  group_by(dcat2) %>%
  summarise(pctsig = mean(sig1)*100, n = length(sig1)) %>%
  ungroup()


x11(13,8.5)
ggplot(corbin1) + 
  geom_col(aes(y=pctsig,x=dcat)) + 
  labs(x="distance (km)",y="% of correlations with P < 0.05")
#ggsave("figures/grid_bin_cor/GRSPgrid_100bin_92.17_AR_5zeros_n20.png")

x11(13,8.5)
ggplot(corbin2) + 
  geom_col(aes(y=pctsig,x=dcat2)) + 
  labs(x="distance (km)",y="% of correlations with P < 0.05") 
#ggsave("figures/grid_bin_cor/GRSPgrid_200bin_92.17_AR_5zeros_n20.png")

## NEXT: need to create n = 1000 null (random) distributions for each bin like Martin et al. 
## to test spatial extent of synchrony


