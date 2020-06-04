incclust = c("42_102", "40_102", "38_102", "36_102", "36_104", "38_100", "36_100", 
             "34_100", "32_100", "34_98", "32_98", "30_98")

grsp.ts = bbsync("5460", 1966, 2017, "AR", "grid.year") %>%
  dplyr::filter(grid4 %in% incclust) %>%
  mutate(sp = "GRSP")

casp.ts = bbsync("5780", 1966, 2017, "AR", "grid.year") %>%
  filter(grid4 %in% incclust) %>%
  mutate(sp = "CASP")

hola.ts = bbsync("4740", 1966, 2017, "AR", "grid.year") %>%
  filter(grid4 %in% incclust) %>%
  mutate(sp = "HOLA")

weme.ts = bbsync("5011", 1966, 2017, "AR", "grid.year") %>%
  filter(grid4 %in% incclust) %>%
  mutate(sp = "WEME")

dick.ts = bbsync("6040", 1966, 2017, "AR", "grid.year") %>%
  filter(grid4 %in% incclust) %>%
  mutate(sp = "DICK")

larb.ts = bbsync("6050", 1966, 2017, "AR", "grid.year") %>%
  filter(grid4 %in% incclust) %>%
  mutate(sp = "LARB")

# bind together time series data for all 6 species for cross-synchrony calculations
cross.data = full_join(grsp.ts, casp.ts, by = "grid4.year") %>%
  dplyr::select(grid4.year, GRSP.ar = AR.x, GRSP.meancount = mean.count.x, CASP.ar = AR.y, CASP.meancount = mean.count.y) %>%
  full_join(hola.ts, by = "grid4.year") %>%
  dplyr::select(grid4.year, GRSP.ar, GRSP.meancount, CASP.ar, CASP.meancount, HOLA.ar = AR, HOLA.meancount = mean.count) %>%
  full_join(weme.ts, by = "grid4.year") %>%
  dplyr::select(grid4.year, GRSP.ar, GRSP.meancount, CASP.ar, CASP.meancount, HOLA.ar, HOLA.meancount, WEME.ar = AR, WEME.meancount = mean.count) %>%
  full_join(dick.ts, by = "grid4.year") %>%
  dplyr::select(grid4.year, GRSP.ar, GRSP.meancount, CASP.ar, CASP.meancount, HOLA.ar, HOLA.meancount, WEME.ar, WEME.meancount, DICK.ar = AR, DICK.meancount = mean.count) %>%
  full_join(larb.ts, by = "grid4.year") %>%
  dplyr::select(grid4.year, GRSP.ar, GRSP.meancount, CASP.ar, CASP.meancount, HOLA.ar, HOLA.meancount, WEME.ar, WEME.meancount, DICK.ar, DICK.meancount, LARB.ar = AR, LARB.meancount = mean.count) %>%
  mutate(grid4 = ifelse(nchar(grid4.year)==10,substr(grid4.year,1,5),substr(grid4.year,1,6))) %>%
  mutate(year = ifelse(nchar(grid4.year)==10,substr(grid4.year,7,10),substr(grid4.year,8,11)))