library(ggplot2)
library(dplyr)
library(tidyr)
library(psych)

h = read.csv("data/AllHay_Yield_survey_allstates.csv")
h = subset(h,Period=="YEAR")
head(h); nrow(h) # 5241

detr.yld = 
  h %>% 
  filter(Year>1965 & Year<2017) %>%
  group_by(State) %>%
  mutate(ryld = resid(lm(Value~Year)))

ggplot(detr.yld) + 
  geom_line(aes(x=Year,y=Value)) + 
  geom_line(aes(x=Year,y=ryld)) +
  facet_wrap(~State)

yld.mat = 
  detr.yld %>%
  select(State,Year,ryld) %>%
  spread(State,ryld)

colnames(yld.mat) = c("year","AL","AK","AZ","AR","CA","CO","CT","DE","FL","GA","ID","IL",
                     "IN","IA","KS","KY","LA","ME","MD","MA","MI","MN","MS","MO","MT",
                     "NE","NV","NH","NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI",
                     "SC","SD","TN","TX","UT","VT","VA","WA","WV","WI","WY")
yld.mat=cormat[,c("AL","AK","AZ","AR","CA","CO","CT","DE","FL","GA","ID","IL",
                 "IN","IA","KS","KY","LA","ME","MD","MA","MI","MN","MS","MO","MT",
                 "NE","NV","NH","NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI",
                 "SC","SD","TN","TX","UT","VT","VA","WA","WV","WI","WY")]

yld.test = 
  cormat %>%
  corr.test(adjust="none")
  
yld.cor = 
  data.frame(reshape2::melt(yld.test$r), 
             reshape2::melt(yld.test$p)[3])
colnames(yld.cor) = c("state1","state2","cor","pval"); head(yld.cor)
yld.cor$pairid = paste(as.character(yld.cor$state1),"-",as.character(yld.cor$state2),sep=""); head(yld.cor); nrow(yld.cor)
yld.cor=yld.cor[yld.cor$cor!=1,]; head(yld.cor); nrow(yld.cor); max(yld.cor$cor)
yld.cor[230:235,] 

# subsetting the region of interest (east of Mississippi River)
yld.cor = subset(yld.cor,state1 %in% c("NH","NY","PA","NJ","ME","CT","RI","MA","VT",
                                  "DE","MD","WV","VA","NC","KY","TN", "OH","AL",
                                  "MS","SC","VA","NC", "MI","GA","IN","IL","WI") &
                   state2 %in% c("NH","NY","PA","NJ","ME","CT","RI","MA","VT",
                                 "DE","MD","WV","VA","NC","KY","TN", "OH","AL",
                                 "MS","SC","VA","NC", "MI","GA","IN","IL","WI") ); nrow(yld.cor)

