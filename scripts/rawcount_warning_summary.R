# make a list of files containing the raw count & warningsummary data sets
rawcount.files = sprintf(
  "data//raw_counts//%s.rawcount.warnings.csv",
  paste0(
    species$spcode,
    ".",
    substr(species$startyr, 3, 4),
    ".",
    substr(species$endyr, 3, 4)
  )
)

# read in all raw count data into list form
rawcount.list = lapply(rawcount.files, function(x)
  read.csv(x))
names(rawcount.list) <-   paste0(species$spcode,
                                 ".",
                                 substr(species$startyr, 3, 4),
                                 ".",
                                 substr(species$endyr, 3, 4))

# create data frames of mean raw count summaries and Jags warnings
rawsum = do.call("rbind", rawcount.list) %>%
  mutate(id = substr(row.names(do.call("rbind", rawcount.list)), 1, 10))

warnings = rawsum %>%
  dplyr::select(id, warning1, warning2) %>%
  distinct() %>%
  arrange(id)

rawsum = rawsum %>%
  dplyr::select(-warning1, -warning2)

rm(rawcount.list)
rm(rawcount.files)