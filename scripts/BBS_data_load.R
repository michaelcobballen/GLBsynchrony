path <- "data//states//"
files <- list.files(path=path, pattern="*.csv")
for(file in files)
{
  perpos <- which(strsplit(file, "")[[1]]==".")
  assign(
    gsub(" ","",substr(file, 1, perpos-1)), 
    read.csv(paste(path,file,sep=""),header=T,skip=0))
}

# combine data files for birds (by state), and read in route and weather data
# BBS fields: ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/NABBS_DataFieldDefinitions_1996-2017.xml
b = rbind.data.frame(Alabama,Alaska,Alberta,Arizona,Arkansa,BritCol,Califor,Colorad,Connect,Delawar,
                     Florida,Georgia,Idaho,Illinoi,Indiana,Iowa,Kansas,Kentuck, Louisia,Maine,
                     Manitob,Marylan,Massach,Michiga,Minneso,Mississ,Missour,Montana,NBrunsw,
                     NCaroli,NDakota,Nebrask,Nevada,Newfoun,NHampsh,NJersey,NMexico,NovaSco,Nunavut,
                     NWTerri,NYork,Ohio,Oklahom,Ontario,Oregon,PEI,Pennsyl,Quebec,RhodeIs,
                     Saskatc,SCaroli,SDakota,Tenness,Texas,Utah,Vermont,Virgini,W_Virgi,Washing,
                     Wiscons,Wyoming,Yukon) %>% 
  mutate(rid = paste(StateNum,".",Route,sep=""))
r = read.csv("data/routes.csv") %>% 
  mutate(statenum=as.character(statenum),rid = paste(statenum,".",Route,sep=""))
w = read.csv("data/weather.csv") %>% 
  mutate(rid = paste(StateNum,".",Route,sep=""))