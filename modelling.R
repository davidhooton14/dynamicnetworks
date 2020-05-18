# data processing to create the trade value/eigenvector centrality rank plots

# total trade per year
totalexportsy = dot2[,-c(1,2,4,5,6)] %>% group_by(fromISO) %>% summarise_at(vars(-1),sum,na.rm=T) 
totalexports = totalexportsy %>% transmute(ISO=fromISO,totalexp=rowSums(.[-1],na.rm=T))
totalimportsy = dot2[,-c(1,2,3,4,6)] %>% group_by(toISO) %>% summarise_at(vars(-1),sum,na.rm=T) %>% transmute(ISO=toISO,totalimp=rowSums(.[-1],na.rm=T))
totalimports = totalimportsy %>% transmute(ISO=toISO,totalexp=rowSums(.[-1],na.rm=T))
totaltradey = cbind(totalexportsy[,1],totalexportsy[,-c(1:6)]+totalimportsy[totalimportsy$toISO!="TWN",-c(1:6)])

# collate by decade
ttdecranked = totaltradey %>% transmute(ISO=fromISO,'1950s' = -rowSums(.[2:8]),'1960s' = -rowSums(.[9:18]),'1970s' = -rowSums(.[19:28]),
                                        '1980s' = -rowSums(.[29:38]),'1990s' = -rowSums(.[39:48]),'2000s' = -rowSums(.[49:58]),'2010s' = -rowSums(.[59:67])) %>%
  mutate_at(vars(-1),rank)
rownames(ttdecranked) = ttdecranked$ISO
ttdecranked$'1940s' = ttdecranked$'1950s'; ttdecranked$'2020s' = ttdecranked$'2010s'
ttdecranked = ttdecranked[,c(1,9,2:8,10)]

totaltradeyordered = totaltradey %>% mutate_at(vars(-1),order)
totaltrade = full_join(totalexports,totalimports,by="ISO") %>% mutate(total = log(rowSums(cbind(totalexp,totalimp),na.rm=T))/10 -1.5 ) %>% as.data.frame
rownames(totaltrade) = totaltrade$ISO

ttdecranked_plot = t(ttdecranked[apply(ttdecranked[,-1],1,function(r) any(r<11)),-1]) %>% as.data.frame
ttdecranked_plot$Decade = seq(1945,2025,10)
ttdecranked_plot = ttdecranked_plot %>% melt(id.vars="Decade") %>% rename(Country=variable,Rank=value)
ttdecranked_plot$Rank[ttdecranked_plot$Rank > 10] = rep(11.5,sum(ttdecranked_plot$Rank > 10))

colstoplot1 = c("#86bf86","#4c8732","#d63e52","#db0400","#c75306","#c9c920","#5113a1","#39a7e3","#0dbf99","#1ab317","#e83397","#fcba03","#247319","#e3e645","#e00d0d")
traderankplot = ggplot(ttdecranked_plot,aes(x=Decade,y=Rank,group=Country)) + geom_line(size=1.2,aes(color=Country)) + 
  scale_color_manual(values=colstoplot1,guide=FALSE) + 
  coord_cartesian(xlim=c(1955,2015),ylim=c(0,10)) + 
  scale_x_continuous(breaks=seq(1955,2015,10),labels=colnames(ttdecranked)[3:9],minor_breaks=NULL) +
  scale_y_reverse(name="Country ranked by total trade",breaks=1:10,minor_breaks=NULL,labels=ttdecranked$ISO[order(ttdecranked$`1950s`)][1:10],
                  sec.axis=sec_axis(~.,breaks=1:10,labels=ttdecranked$ISO[order(ttdecranked$`2010s`)][1:10])) + 
  theme_bw()
ggsave("plots/traderankplot.pdf",traderankplot,units="in",width=7,height=2.9)


# eigenvector centralities by year 
ecs = lapply(dotgraphs,function(g) sort(eigen_centrality(g,weights=log(edge.attributes(g)[[1]]+1))$'vector',decreasing=T)) # (weighted)
ecnames = V(simplify(do.call(union,dotgraphs)))$name
ecdf = matrix(nrow=length(ecnames),ncol=1+length(ecs),dimnames=list(ecnames,c("ISO",names(ecs)))) %>% as.data.frame
for(i in 1:length(ecs)){
  ecdf[names(ecs[[i]]),i+1] = as.numeric(ecs[[i]])
}
ecdf[,1] = rownames(ecdf)
# collate by decade
ecdfranked = ecdf %>% transmute(ISO=ISO,'1950s' = -rowSums(.[2:8]),'1960s' = -rowSums(.[9:18]),'1970s' = -rowSums(.[19:28]),
                                '1980s' = -rowSums(.[29:38]),'1990s' = -rowSums(.[39:48]),'2000s' = -rowSums(.[49:58]),'2010s' = -rowSums(.[59:67])) %>%
  mutate_at(vars(-1),rank)
ecdfranked$'1940s' = ecdfranked$'1950s'; ecdfranked$'2020s' = ecdfranked$'2010s'
ecdfranked = ecdfranked[order(ecdfranked$ISO),c(1,9,2:8,10)]
rownames(ecdfranked) = ecdfranked$ISO

ecdfranked_plot = t(ecdfranked[apply(ecdfranked[,-1],1,function(r) any(r<11)),-1]) %>% as.data.frame
ecdfranked_plot$Decade = seq(1945,2025,10)
ecdfranked_plot = ecdfranked_plot %>% melt(id.vars="Decade") %>% rename(Country=variable,Rank=value)
ecdfranked_plot$Rank[ecdfranked_plot$Rank > 10] = rep(11.5,sum(ecdfranked_plot$Rank > 10))

colstoplot2 = c("#d63e52","#8a8a8a","#c75306","#c9c920","#53edcc","#ede213","#5113a1","#39a7e3","#1ab317","#e83397","#fcba03","#e3e645","#e00d0d")
colstoplot3 = c("#917200","#d63e52","#8a8a8a","#c75306","#c9c920","#ede213","#5113a1","#39a7e3","#247319","#1ab317","#e83397","#fcba03","#e3e645","#e00d0d")
ecrankplot1 = ggplot(ecdfranked_plot,aes(x=Decade,y=Rank,group=Country)) + geom_line(size=1.2,aes(color=Country)) + 
  scale_color_manual(values=colstoplot2,guide=FALSE) + 
  coord_cartesian(xlim=c(1955,2015),ylim=c(0,10)) + 
  scale_x_continuous(breaks=seq(1955,2015,10),labels=colnames(ecdfranked)[3:9],minor_breaks=NULL) +
  scale_y_reverse(name="Country ranked by eigenvector \ncentrality (unweighted)",breaks=1:10,minor_breaks=NULL,labels=ecdfranked$ISO[order(ecdfranked$`1950s`)][1:10],
                  sec.axis=sec_axis(~.,breaks=1:10,labels=ecdfranked$ISO[order(ecdfranked$`2010s`)][1:10])) + 
  theme_bw()
ggsave("plots/ecrankplot.pdf",ecrankplot1,units="in",width=7,height=2.9)




# work out intercountry distances
dotdistancematrix = matrix(nrow=nrow(caps),ncol=nrow(caps),
                           dimnames=list(caps$fromISO,caps$fromISO))
for(i in 1:nrow(caps)){
  dotdistancematrix[i,] = sapply(1:nrow(caps),function(j) 
    distGeo(caps[i,c("long","lat")],caps[j,c("long","lat")])/1000)
  # great circle distance in km
}





# AME models
Yadjs = list()
Yadjbigs = list()
Xds = list()
dotames = list()
simpdotames = list()
for(y in seq(1960,2010,2)){
  # create matrices of predictor variables
  # distances
  Xd1_ = sqrt(dotdistancematrix)/10 # gets coefficient on same scale
  rs2 = intersect(row.names(Xd1_),unique(atopyears[[as.character(y)]]$fromISO))
  Xd2_ = matrix(rep(0,prod(dim(Xd1_))),nrow=nrow(Xd1_),ncol=ncol(Xd1_),
                dimnames = list(row.names(Xd1_),row.names(Xd1_)))
  Xd2_[rs2,rs2] = as.matrix(as_adj(atopgraphs[[as.character(y)]]))[rs2, rs2]
  rs3 = intersect(row.names(Xd1_),row.names(colonyadj))
  Xd3_ = matrix(rep(0,prod(dim(Xd1_))),nrow=nrow(Xd1_),ncol=ncol(Xd1_),
                dimnames = list(row.names(Xd1_),row.names(Xd1_)))
  Xd3_[rs3,rs3] = colonyadj[rs3,rs3]
  
  Xd_ = abind(Xd1_,Xd2_,Xd3_,along=3)
  dimnames(Xd_)[[3]] = c("sqrtdist","alliance","colony")
  Xds[[as.character(y)]] = Xd_
  
  Yadjbig = matrix(rep(0,prod(dim(Xd1_))),nrow=nrow(Xd1_),ncol=ncol(Xd1_),
                dimnames = list(row.names(Xd1_),row.names(Xd1_)))
  Yadj = as.matrix(as_adj(dotgraphs[[as.character(y)]],attr=as.character(y)))
  Yadj = log(1+Yadj/1e3) # trade now in $1000s and log transformed
  Yadjbig[row.names(Yadj),row.names(Yadj)] = Yadj 
  Yadjbigs[[as.character(y)]] = Yadjbig
  Yadjs[[as.character(y)]] = Yadj

  dotamey = ame(Y=Yadj>0,Xdyad=Xd_[row.names(Yadj),row.names(Yadj),],Xr=rowMeans(Yadj),Xc=colMeans(Yadj),
                model="bin",dcor=T,burn=150,nscan=2000,odens=10,R=2)
  dotames[[as.character(y)]] = dotamey
  
  simpdotamey = ame(Y=Yadj>0,Xdyad=Xd_[row.names(Yadj),row.names(Yadj),1],Xr=rowMeans(Yadj),Xc=colMeans(Yadj),
                    model="bin",dcor=T,burn=150,nscan=2000,odens=10,R=2)
  simpdotames[[as.character(y)]] = simpdotamey
}
summary(dotames[[1]])
summary(dotamey)

sapply(1:maxT,function(i) sapply(abs(colMeans(dotames[[i]]$BETA)/apply(dotames[[i]]$BETA,2,sd)),dnorm))

plotcoefs = t(dotamecoefs[4:6,]/dotamecoefs[4:6,2]) %>% as.data.frame %>%
  cbind(year=seq(1960,2010,2)) %>%
  melt(id.vars="year") %>%
  rename(Coefficient=variable)

dotcoefplot = ggplot(plotcoefs,aes(x=year,group=Coefficient,y=value)) + 
#  geom_hline(aes(yintercept=0),col=rgb(0.2,0.2,0.2,0.7)) +
  geom_line(aes(color=Coefficient),size=1.4,show.legend=T) + 
  theme_bw() + labs(x="Year",y="Coefficient (relative to 1962)") + 
  scale_color_manual(values=c("#39A780",rgb(55, 110, 191,maxColorValue=255),
                              rgb(235,128,27,maxColorValue=255)),labels=c("Distance","Alliance","Colony"))
ggsave("plots/dotcoefplot.pdf",dotcoefplot,width=6,height=3,units="in")

# attempt at a time-dependent model
# diagnostics are terrible
Yadjt = abind(Yadjbigs,along=3)
maxT = dim(Yadjt)[3]
Xdt = abind(Xds,along=4)
Ylag = array(Yadjt[,,1:(maxT-1)],dim=c(dim(Yadjt)[1:2],1,maxT-1))
tYlag = array(apply(Ylag[,,1,],3,t),dim=dim(Ylag))
Yrmeans = array(apply(Yadjt,3,rowMeans),dim=c(ncol(Yadjt),1,maxT))
Ycmeans = array(apply(Yadjt,3,colMeans),dim=c(ncol(Yadjt),1,maxT))
amet = ame_rep(Yadjt[,,2:T],Xd=abind(Xdt[,,,2:6],"Ylag"=Ylag,"tYlag"=tYlag,along=3),Xr=Yrmeans,Xc=Ycmeans,  
               model="nrm",burn=100,nscan=1000,odens=10,R=maxT)
summary(amet)




