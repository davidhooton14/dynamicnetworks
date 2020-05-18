library(readr)
library(dplyr)
library(igraph)
library(maps)
library(countrycode)
library(matrixcalc)
library(abind)
library(amen)
library(reshape2)
library(ggplot2)
library(splines)
library(quantmod)





# DOT trade data
dot2 = read.csv("data/imf/DOT2.csv",stringsAsFactors=F)
ccodes = read.csv("data/imf/countrycodes.csv",stringsAsFactors=F)

dot2 = dot2 %>% subset(Country.Code %in% ccodes$IMF.Code & Counterpart.Country.Code %in% ccodes$IMF.Code)
dot2[,8:ncol(dot2)] = apply(dot2[,3:ncol(dot2)],2,as.numeric)
# will introduce NAs (deliberately) and give warnings that can be ignored

# inflation adjustment data 
getSymbols("CPIAUCSL", src='FRED') #Consumer Price Index for All Urban Consumers: All Items
avg.cpi = apply.yearly(CPIAUCSL, mean)
cf = as.numeric(avg.cpi/as.numeric(avg.cpi['2020']))
names(cf) = 1947:2020

dotcountries = dot2[,c(2,6)] %>% left_join(ccodes,by=c("Country.Code"="IMF.Code")) %>%
  left_join(ccodes,by=c("Counterpart.Country.Code"="IMF.Code"))
colnames(dotcountries) = c("fromIMF","toIMF","fromISO","fromName","toISO","toName")
dot2 = cbind(dotcountries,dot2[,8:ncol(dot2)])
colnames(dot2)[7:ncol(dot2)] = as.character(1948:2018)
# adjust for inflation
dot2[,7:ncol(dot2)] = dot2[,7:ncol(dot2)]/cf[as.character(1948:2018)]





# ATOP alliance data
atopdyads = read.csv("data/atop_v4.01_data__csv_/ATOP V4.01 Data (csv)/atop4_01dy.csv")
cowcodes = read.csv("data/atop_v4.01_data__csv_/ATOP V4.01 Data (csv)/COW country codes.csv")

atopdyads$fromcode = atopdyads$dyad %/% 1000
atopdyads$tocode = atopdyads$dyad %% 1000
atopdyads = atopdyads %>% 
  left_join(cowcodes,by=c("fromcode"="CCode")) %>%
  rename(fromAbb = StateAbb, fromName = StateNme) %>%
  left_join(cowcodes,by=c("tocode"="CCode")) %>%
  rename(toAbb = StateAbb, toName = StateNme)





# COLDAT colony data
COLDAT_dyads <- read.csv("data/colony/COLDAT_dyads.csv", stringsAsFactors=FALSE)
Cdyads = COLDAT_dyads %>% subset(col==1)
write.csv(Cdyads,"data/colony/Cdyads.csv")





# get capital city positions (to be used as country locations)
caps = maps::world.cities %>% subset(capital==1) %>% arrange(country.etc)
caps = dotcountries[match(unique(dotcountries$fromISO),dotcountries$fromISO),c("fromName","fromISO")] %>%
  left_join(caps,by=c("fromName"="country.etc")) %>% subset(!(fromISO %in% c("KOS","SSD")))
# one of many times where country or capital names were slightly incorrect across datasets and had to be matched manually
manual = data.frame(matrix(c("BHS","Nassau","BRN","Bandar Seri Begawan","HKG","Hong Kong","COD","Kinshasa","COG","Brazzaville",
                             "CIV","Yamoussoukro","CPV","Praia","FJI","Suva","GMB","Banjul","KOR","Soul","KGZ","Biskek","LAO","Vientiane","MKD","Skopje",
                             "MNE","Podgorica","STP","Sao Tome","SRB","Belgrade","SVK","Bratislava","KNA","Basseterre","LCA","Castries","VCT","Kingstown",
                             "TLS","Dili","GBR","London","USA","Washington"),ncol=2,byrow=T),stringsAsFactors=F)
manual = manual[order(manual[,1]),]
manual = left_join(manual,maps::world.cities,by=c("X2"="name"))[-c(2,9,10,25,26),]
caps[caps$fromISO %in% manual[,1],] = left_join(caps[caps$fromISO %in% manual[,1],1:2],manual[,-c(3)],by=c("fromISO"="X1"))[,]
caps[35,c(5,6)] = c(22.32,114.17)





# creating trade network graphs from the data - an inefficient process

# full graphs
dotgraphs = list() # graphs by year
dotyears = list() # data by year
dotgls = list() # map layouts (countries and capital positions) by year

for(y in 1953:2018){
  doty = dot2[,c(1:6,y-1941)]
  doty = doty[!is.na(doty[,7]),] %>%
    subset(!(fromISO %in% c("TWN","KOS","SSD")) & !(toISO %in% c("TWN","KOS","SSD")))
  doty = do.call(rbind,lapply(split(doty,doty$fromISO), function(df) df[order(-df[,7]),]))
  doty$ranks = unlist(lapply(split(doty,doty$fromISO),function(df) order(-df[,7])))
  dotyears[[as.character(y)]] = doty
  
  
  dotgy = graph_from_data_frame(doty[,c(3,5,7,8)],directed=T) %>% simplify(edge.attr.comb="sum")
  E(dotgy)$width = ((log(edge.attributes(dotgy)[[1]])-10)/7)^2
  dotgraphs[[as.character(y)]] = dotgy
  
  dotgly = as.matrix(caps[match(V(dotgy)$name,caps$fromISO),c("long","lat")])
  dotgls[[as.character(y)]] = norm_coords(dotgly,xmin=-2,xmax=2,ymin=-2/3,ymax=1)
}

# data/graphs where only each node's top K trade links are included
dotgraphsK = list()
dotyearsK = list()
dotglsK = list()
K = 5

for(y in 1953:2018){
  doty = dot2[,c(1:6,y-1941)]
  doty = doty[!is.na(doty[,7]),] %>%
    subset(!(fromISO %in% c("TWN","KOS","SSD")) & !(toISO %in% c("TWN","KOS","SSD")))
  # select top K trade partners for each country
  doty = do.call(rbind,lapply(split(doty,doty$fromISO), function(df) head(df[order(-df[,7]),],K)))
  doty$ranks = unlist(lapply(split(doty,doty$fromISO),function(df) order(-df[,7])))
  dotyearsK[[as.character(y)]] = doty
  
  dotgy = graph_from_data_frame(doty[,c(3,5,7,8)],directed=F) %>% simplify(edge.attr.comb="mean")
  E(dotgy)$width = ((log(edge.attributes(dotgy)[[1]])-10)/7)^2
  dotgraphsK[[as.character(y)]] = dotgy
  
  dotgly = as.matrix(caps[match(V(dotgy)$name,caps$fromISO),c("long","lat")])
  dotglsK[[as.character(y)]] = norm_coords(dotgly,xmin=-2,xmax=2,ymin=-2/3,ymax=1)
}


# small graphs: top M countries by total cumulative exports
M = 40
dottotals = dot2[,-c(1,2,4:6)] %>% group_by(fromISO) %>% summarise_all(sum,na.rm=T)
dottotals$all = rowSums(dottotals[,-1])
largestMexporters = dottotals$fromISO[order(-dottotals$all)][1:M]

dotgssmall = list()
dotyssmall = list()
dotglssmall = list()
K = 5

for(y in 1953:2018){
  doty = dot2[dot2$fromISO %in% largestMexporters & dot2$toISO %in% largestMexporters,c(1:6,y-1941)]
  doty = doty[!is.na(doty[,7]),] %>%
    subset(!(fromISO %in% c("TWN","KOS","SSD")) & !(toISO %in% c("TWN","KOS","SSD")))
  # select top K trade partners for each country
  # doty = do.call(rbind,lapply(split(doty,doty$fromISO), function(df) head(df[order(-df[,7]),],K)))
  doty$ranks = unlist(lapply(split(doty,doty$fromISO),function(df) order(-df[,7])))/39
  dotyssmall[[as.character(y)]] = doty
  
  dotgy = graph_from_data_frame(doty[,c(3,5,7,8)],directed=F) %>% simplify(edge.attr.comb="mean")
  E(dotgy)$width = ((log(edge.attributes(dotgy)[[1]])-10)/7)^2
  dotgssmall[[as.character(y)]] = dotgy
  
  dotgly = as.matrix(caps[match(V(dotgy)$name,caps$fromISO),c("long","lat")])
  dotglssmall[[as.character(y)]] = norm_coords(dotgly,xmin=-2,xmax=2,ymin=-2/3,ymax=1)
}





# more data processing for the alliance data

atopcodes = read.csv("data/atop_v4.01_data__csv_/atopcodes.csv",stringsAsFactors=F)

caps2 = atopcodes %>% left_join(caps[,c(2,3,5,6)],by="fromISO")
manual2 = data.frame(matrix(c("CUB","Havanna","CZS","Prague","GDR","Berlin","GFR","Bonn",
                              "LIE","Vaduz","MNC","Monte Carlo","NGR","Abuja","PRK","Pyongyang","SOM","Mogadishu",
                              "TAW","Taipei","YUG","Podgorica"),ncol=2,byrow=T),stringsAsFactors=F)
manual2 = manual2[order(manual2[,1]),]
manual2 = left_join(manual2,maps::world.cities,by=c("X2"="name"))[-c(3,5),]
caps2[caps2$fromISO %in% manual2[,1],4:6] = left_join(caps2[caps2$fromISO %in% manual2[,1],c(1,3)],manual2,by=c("fromISO"="X1"))[,c(3,6:7)]
caps2 = caps2[apply(caps2,1,function(v) !any(is.na(v))),]

atopdyads = atopdyads[(atopdyads$fromAbb %in% caps2$fromAbb)
                      & (atopdyads$toAbb %in% caps2$fromAbb),]

# creating alliance network graphs from the data

atopgraphs = list()
atopyears = list()
atopgls = list()

for(y in 1900:2016){
  atopy = atopdyads[atopdyads$year==y,c(2,27:32)] %>% 
    left_join(atopcodes[,2:3],by="fromAbb") %>%
    left_join(atopcodes[,2:3],by=c("toAbb"="fromAbb")) %>%
    rename(fromISO=fromISO.x,toISO=fromISO.y) # can ignore warnings
  atopyears[[as.character(y)]] = atopy
  atopgy = graph_from_data_frame(atopy[,8:9]) %>% 
    as.undirected(mode="collapse") %>% simplify()
  E(atopgy)$width = rep(2,length(E(atopgy)))
  E(atopgy)$ranks = rep(1,length(E(atopgy)))
  atopgraphs[[as.character(y)]] = atopgy
  atopgl = sweep(as.matrix(caps2[match(V(atopgy)$name,caps2$fromISO),c("long","lat")]),
                 2,c(150,75),FUN='/')
  atopgls[[as.character(y)]] = atopgl
}






