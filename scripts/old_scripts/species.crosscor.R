#######################################
###### cross-correlation among species
#######################################

grsp.grid.year.92.17 = bbsync("5460",1992,2017,"AR","grid.year")
bobo.grid.year.92.17 = bbsync("4940",1992,2017,"AR","grid.year")
eame.grid.year.92.17 = bbsync("5010",1992,2017,"AR","grid.year")
savs.grid.year.92.17 = bbsync("5420",1992,2017,"AR","grid.year")

grsp.bobo.92.17 = grsp.grid.year.92.17 %>%
  left_join(bobo.grid.year.92.17, by = "grid4.year") %>%
  filter(is.na(AR.x)==F,is.na(AR.y)==F) %>%
  group_by(grid4.x) %>%
  summarize(ccor = cor(AR.x, AR.y)) %>%
  mutate(grid4_lat = as.numeric(substr(grid4.x,1,2))) %>%
  mutate(grid4_lon = -1*as.numeric(substr(grid4.x,4,10))) %>%
  dplyr::select(grid4 = grid4.x, ccor,grid4_lat,grid4_lon) %>%
  mutate(ccorcat = cut(ccor, breaks = c(-1,0,0.1,.2,.3,.4,.5,.6)))

grsp.eame.92.17 = grsp.grid.year.92.17 %>%
  left_join(eame.grid.year.92.17, by = "grid4.year") %>%
  filter(is.na(AR.x)==F,is.na(AR.y)==F) %>%
  group_by(grid4.x) %>%
  summarize(ccor = cor(AR.x, AR.y)) %>%
  mutate(grid4_lat = as.numeric(substr(grid4.x,1,2))) %>%
  mutate(grid4_lon = -1*as.numeric(substr(grid4.x,4,10))) %>%
  dplyr::select(grid4 = grid4.x, ccor,grid4_lat,grid4_lon) %>%
  mutate(ccorcat = cut(ccor, breaks = c(-1,0,0.1,.2,.3,.4,.5,1)))

grsp.savs.92.17 = grsp.grid.year.92.17 %>%
  left_join(savs.grid.year.92.17, by = "grid4.year") %>%
  filter(is.na(AR.x)==F,is.na(AR.y)==F) %>%
  group_by(grid4.x) %>%
  summarize(ccor = cor(AR.x, AR.y)) %>%
  mutate(grid4_lat = as.numeric(substr(grid4.x,1,2))) %>%
  mutate(grid4_lon = -1*as.numeric(substr(grid4.x,4,10))) %>%
  dplyr::select(grid4 = grid4.x, ccor,grid4_lat,grid4_lon) %>%
  mutate(ccorcat = cut(ccor, breaks = c(-1,0,0.1,.2,.3,.4,.5,1)))

bobo.savs.92.17 = bobo.grid.year.92.17 %>%
  left_join(savs.grid.year.92.17, by = "grid4.year") %>%
  filter(is.na(AR.x)==F,is.na(AR.y)==F) %>%
  group_by(grid4.x) %>%
  summarize(ccor = cor(AR.x, AR.y)) %>%
  mutate(grid4_lat = as.numeric(substr(grid4.x,1,2))) %>%
  mutate(grid4_lon = -1*as.numeric(substr(grid4.x,4,10))) %>%
  dplyr::select(grid4 = grid4.x, ccor,grid4_lat,grid4_lon) %>%
  mutate(ccorcat = cut(ccor, breaks = c(-1,0,0.1,.2,.3,.4,.5,1)))

bobo.eame.92.17 = bobo.grid.year.92.17 %>%
  left_join(eame.grid.year.92.17, by = "grid4.year") %>%
  filter(is.na(AR.x)==F,is.na(AR.y)==F) %>%
  group_by(grid4.x) %>%
  summarize(ccor = cor(AR.x, AR.y)) %>%
  mutate(grid4_lat = as.numeric(substr(grid4.x,1,2))) %>%
  mutate(grid4_lon = -1*as.numeric(substr(grid4.x,4,10))) %>%
  dplyr::select(grid4 = grid4.x, ccor,grid4_lat,grid4_lon) %>%
  mutate(ccorcat = cut(ccor, breaks = c(-1,0,0.1,.2,.3,.4,.5,1)))

eame.savs.92.17 = eame.grid.year.92.17 %>%
  left_join(savs.grid.year.92.17, by = "grid4.year") %>%
  filter(is.na(AR.x)==F,is.na(AR.y)==F) %>%
  group_by(grid4.x) %>%
  summarize(ccor = cor(AR.x, AR.y)) %>%
  mutate(grid4_lat = as.numeric(substr(grid4.x,1,2))) %>%
  mutate(grid4_lon = -1*as.numeric(substr(grid4.x,4,10))) %>%
  dplyr::select(grid4 = grid4.x, ccor,grid4_lat,grid4_lon) %>%
  mutate(ccorcat = cut(ccor, breaks = c(-1,0,0.1,.2,.3,.4,.5,1)))

#######################################
###### mapping cross-correlation among species
#######################################

#grsp v. bobo
x11(10,10)
na + # na is defined in "mapping mean synchrony" above
  geom_point(aes(x = grid4_lon, y = grid4_lat, color = ccorcat), size = 2, alpha = .7, data = grsp.bobo.92.17) + 
  scale_colour_manual(values = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  #  facet_wrap(~sp, ncol=2) +
  labs(colour="Cross\ncorrelation", x="", y="", title="GRSP v. BOBO, 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))

# grsp v. eame
x11(10,10)
na + # na is defined in "mapping mean synchrony" above
  geom_point(aes(x = grid4_lon, y = grid4_lat, color = ccorcat), size = 2, alpha = .7, data = grsp.eame.92.17) + 
  scale_colour_manual(values = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  #  facet_wrap(~sp, ncol=2) +
  labs(colour="Cross\ncorrelation", x="", y="", title="GRSP v. EAME, 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))

# grsp v. savs
x11(10,10)
na + # na is defined in "mapping mean synchrony" above
  geom_point(aes(x = grid4_lon, y = grid4_lat, color = ccorcat), size = 2, alpha = .7, data = grsp.savs.92.17) + 
  scale_colour_manual(values = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  #  facet_wrap(~sp, ncol=2) +
  labs(colour="Cross\ncorrelation", x="", y="", title="GRSP v. SAVS, 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))

# bobo v. savs
x11(10,10)
na + # na is defined in "mapping mean synchrony" above
  geom_point(aes(x = grid4_lon, y = grid4_lat, color = ccorcat), size = 2, alpha = .7, data = bobo.savs.92.17) + 
  scale_colour_manual(values = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  #  facet_wrap(~sp, ncol=2) +
  labs(colour="Cross\ncorrelation", x="", y="", title="BOBO v. SAVS, 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))

# bobo v. eame
x11(10,10)
na + # na is defined in "mapping mean synchrony" above
  geom_point(aes(x = grid4_lon, y = grid4_lat, color = ccorcat), size = 2, alpha = .7, data = bobo.eame.92.17) + 
  scale_colour_manual(values = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  #  facet_wrap(~sp, ncol=2) +
  labs(colour="Cross\ncorrelation", x="", y="", title="BOBO v. EAME, 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))

# eame v. savs
x11(10,10)
na + # na is defined in "mapping mean synchrony" above
  geom_point(aes(x = grid4_lon, y = grid4_lat, color = ccorcat), size = 2, alpha = .7, data = eame.savs.92.17) + 
  scale_colour_manual(values = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  #  facet_wrap(~sp, ncol=2) +
  labs(colour="Cross\ncorrelation", x="", y="", title="EAME v. SAVS, 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))