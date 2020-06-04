{  
  ## Annual index model for the NABBS, with spatial conditional autoregressive (CAR) structure on stratum-level annual indices
  ## based on the BBS trend model described in Sauer and Link 2011, and including the modification described in Smith et al. 2014
  ## this model differs from the BBS trend models in the following ways:
  ## no beta term - no estimated continuous long-term trend
  ## year effects are fit separately for each year, with a spatial CAR structure among neighbouring strata
  ## CAR year effects are constrained to sum to zero, so a continent wide mean yeareffect term is also fit
  ## stratum-specific intercepts are not estimated, they are supplied as data (effectively an off-set) - log(average observed count across all routes, observers, and years in each stratum)
  ## the supplying the intercepts as offsets ensures that the spatial year effects are only modeling annual departures from the stratum-level average abundance, and not the variation in average abundance among strata.
  
  ## data to supply:
  #count - (vector, length ncounts) total number of birds observed on route in a given year 
  #strat - (vector, length ncounts) analytical stratum
  #obser - (vector, length ncounts) unique combination of route and observer
  #firstyr - (vector, length ncounts) indicator variable (0,1) identifying the first year of a given observer on a given route
  #year - (vector, length ncounts) year (1,number of years in analysis)
  
  #ncounts - (vector, length 1) number of route*years included in the analysis
  #nstrata - (vector, length 1) number of analytical strata included in the analysis
  #ymin - (vector, length 1) earliest year included in the analysis (almost always = 1)
  #ymax - (vector, length 1) latest year included in the analysis
  #nonzeroweight - (vector, length nstrata) proportion of routes in each stratum that are included in the analysis
  #nobservers - (vector, length nstrata) number of unique Observer-Route combinations in each stratum
  #strata - (vector, length nstrata) log mean count (across routes and strata) for each stratum
# data specific to the car.normal() function (see the WinBUGS map documentation for more details)
  #sumNumNeigh - (vector, length 1) total number of neighbour relationships
  #adj - (vector, length sumNumNeigh) neighbouring strata lists for each stratum
  #num - (vector, length nstrata) number of neighbour relationships for each stratum
  
  
  
  
  
  
  #### main model statement, counts, overdispersion effects, and fit statistics  ######
  
  
  for( k in 1 : ncounts ) {
    log(lambda[k]) <-  obs[strat[k],obser[k]] + strata[strat[k]] + eta*firstyr[k] + noise[k] + yeareffect[year[k],strat[k]]
    noise[k] ~ dnorm(0.0, taunoise)
    count[k] ~ dpois(lambda[k])
    fcount[k] ~ dpois(lambda[k])
    err[k] <- pow(count[k]-lambda[k],2)/lambda[k]
    ferr[k] <- pow(fcount[k]-lambda[k],2)/lambda[k]
  }
  
  gof <- sum(err[1:ncounts])
  fgof <- sum(ferr[1:ncounts])
  diffgof <- gof-fgof
  posdiff <- step(diffgof)
  
  taunoise ~ dgamma(0.001,0.001)
  sdnoise <- 1 / pow(taunoise, 0.5)
  
  #### end main model statement, counts, overdispersion effects, and fit statistics  ######
  
  
  
  #### observer effect priors  ######
  ############ alternative, informative priors that can help avoid numerical overlow in BUGS program
  #mulogtauobs ~ dnorm(0,4)
  #taulogtauobs ~ dgamma(2,0.2)
  ############ 
  
  mulogtauobs ~ dnorm(0.0,1.0E-6) 
  taulogtauobs ~ dgamma(0.001,0.001) 
  
  eta ~ dnorm( 0.0,1.0E-6) # first-year observer effect
  
  #### end observer effect priors  ######
  
  
  
  #### stratum-level effects  ######
  for( s in 1 : nstrata ) {
    
    #### stratum-specific observer effects #### 
    for( i in 1 : nobservers[s] ) {
      obs[s,i] ~ dnorm( 0.0,tauobs[s])
    }
    
    log(tauobs[s]) <- logtauobs[s]
    logtauobs[s] ~ dnorm(mulogtauobs,taulogtauobs)
    sdobs[s] <- 1 / pow(tauobs[s], 0.5)
    #### end stratum-specific observer effects #### 
    
    
    #### Sum of year effect and stratum*year CAR year effects ####   
    for( y in ymin : ymax) {
      yeareffect[y,s] <- yearcar[y,s] + betacar[y]
    } 
    
    #### end Sum of year effect and stratum*year CAR year effects ####   
  }    
  
  #### end stratum-level effects  ######
  
  
  #### CAR year effects  ######
  for( y in ymin : ymax) {
    yearcar[y,1:nstrata] ~ car.normal(adj[], weights[], num[], tauyear[y]) 
    tauyear[y] ~ dgamma(0.5,0.0005)
    betacar[y] ~ dflat()
  }
  for(k in 1:sumNumNeigh) {weights[k] <- 1}
  #### end CAR year effects  ######
  
  
  #### summary statistics - annual indices = n[s,t]  ######
  for( s in 1 : nstrata ) {
    for( t in ymin : ymax) {
      n[s,t] <- nonzeroweight[s]*exp(yeareffect[t,s]+ strata[s] + 0.5*sdnoise*sdnoise+ 0.5*sdobs[s]*sdobs[s])
    }
  }
  
  #### end summary statistics ######
}


