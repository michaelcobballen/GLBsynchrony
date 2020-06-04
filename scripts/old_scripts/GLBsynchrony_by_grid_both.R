### TO PLOT 1X1 GRID AND 2X2 GRID SYNCHRONY ON SAME GRAPH
### NEED TO RUN THE FOLLOWING SCRIPTS FIRST:
    ### bbs_to_grid4.R, GLBsynchrony_by_grid4.R, bbs_to_grid.R, GLBsynchrony_by_grid.R 
library(cowplot)
library(ggplot2)

X11(13,9)
ggplot(grid4.bird.cordist.92.17) + 
  geom_smooth(aes(x=dist,y=cor),method="auto",color="orange", se=T) +
  geom_smooth(data=grid.bird.cordist.92.17, 
              aes(x=dist,y=cor),method="auto",color="darkgreen", se=T) +
  xlim(c(0,2000)) +
  labs(title="GRSP spatial synchrony - 1 vs. 4 deg. grid cells", x = "Distance (km)", y = "Pearson's correlation") +
  theme_bw() +
  theme(text = element_text(size=20)) 
  
#ggsave("figures/GRSP_cor_grid1&4_92.17_lines_AR_5zeros_n20.jpg")




#######################################
###### plotting % signficant correlations by distance bins 
#######################################

binlabels2 = c("<200","201-400","401-600","601-800","801-1000","1001-1200",
               "1201-1400","1401-1600","1601-1800","1801-2000",">2000")

### GRID 4 INFO
grid4.bird.cordist.92.17 = grid4.bird.cordist.92.17 %>%
  mutate(dcat2 = cut(grid4.bird.cordist.92.17$dist,c(0,200,400,600,800,1000,
                                                     1200,1400,1600,1800,2000,100000),
                     labels=binlabels2))

corbin2.4 = grid4.bird.cordist.92.17 %>%
  group_by(dcat2) %>%
  summarise(pctsig = mean(sig1)*100, n = length(sig1)) %>%
  ungroup()

### GRID 1 INFO
grid.bird.cordist.92.17 = grid.bird.cordist.92.17 %>%
  mutate(dcat2 = cut(grid.bird.cordist.92.17$dist,c(0,200,400,600,800,1000,
                                                     1200,1400,1600,1800,2000,100000),
                     labels=binlabels2))

corbin2 = grid.bird.cordist.92.17 %>%
  group_by(dcat2) %>%
  summarise(pctsig = mean(sig1)*100, n = length(sig1)) %>%
  ungroup()

### COMBINDED PLOTS

x11(13,8.5)
ggplot(corbin2.4) + 
  geom_col(aes(y=pctsig,x=dcat2), fill = "orange", alpha=.5) + 
  geom_col(data = corbin2, aes(y=pctsig,x=dcat2), fill="darkgreen", alpha = .5) +
  labs(x="Distance (km)",y="% of correlations with P < 0.05", title="GRSP spatial synchrony - 1 vs. 4 deg. grid cells") +
  theme_bw() +
  theme(text = element_text(size=20), axis.text.x = element_text(size=12)) 
#ggsave("figures/GRSPgrid1&4_200bin_92.17_AR_5zeros_n20.png")

