###############################################################################################
###
### This is an attempt to get 1-degree-grid-cell index values via the bbsBayes package
### It is buggy at the moment and producing non-sensical results compared with USGS published indices
### Would be easiest if it can be made to work
###
###############################################################################################

# To install v1.0.0-beta.2 from Github:
#install.packages("devtools")
#library(devtools)
#devtools::install_github("BrandonEdwards/bbsBayes", ref = "v1.0.0-beta.2", dependencies = TRUE)

# To install the development version from GitHub:
#install.packages("devtools")
#library(devtools)
#devtools::install_github("BrandonEdwards/bbsBayes")

library(bbsBayes)
library(ggplot2)
library(geofacet)

bbs_data <- fetch_bbs_data()

bbs_sample <- fetch_sample_data()

strat_data <- stratify(bbs_data, stratify_by = "state")

jags_data <- prepare_jags_data(strat_data, 
                               species_to_run = "Grasshopper Sparrow", 
                               model = "firstdiff")
mod <- run_model(jags_data = jags_data)

#mod <- run_model(jags_data = jags_data,
#                 n_burnin = 1000,
#                 n_iter=1000,
#                 n_adapt = 500)

strat_indices <- generate_strata_indices(mod)
strat_trend <- generate_strata_trends(indices = strat_indices)

generate_map(strat_trend, stratify_by = "state")

#write.csv(strat_indices, "data/bbsBayes_state_indices_slope.csv")
#write.csv(strat_indices, "data/bbsBayes_state_indices_firstdiff_v1.csv")


grsp_state_slope = read.csv("data/bbsBayes_state_indices_slope.csv")
grsp_state_firstdiff = read.csv("data/bbsBayes_state_indices_firstdiff.csv")
grsp_state_firstdiff_v1 = read.csv("data/bbsBayes_state_indices_firstdiff_v1.csv")

#X11(26,18)
ggplot(grsp_state_firstdiff_v1) + 
  geom_line(aes(x=Year,y=Index))  +
  facet_geo(~Stratum, scales = "free_y") +
  theme_bw() +
  xlab("Year") + ylab("BBS Index") +
  theme(axis.text=element_text(size=3.25))
# ggsave("figures/bbs_indices/GRSP_bbsBayes_state_slope.png")
# ggsave("figures/bbs_indices/GRSP_bbsBayes_state_firstdiff.png")
# ggsave("figures/bbs_indices/GRSP_bbsBayes_state_firstdiff_v1.png")

