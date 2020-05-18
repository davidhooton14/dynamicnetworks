# please run setup.R to load all packages and data into the environment
# the following functions have not been exhaustively tested, use wisely

# layout phase of the algorithm
layout_dynamically = function(graphlist,updatespy=0,niter=1,wattr=NULL,initl=NULL,method="fr",wtmultiplier=1){
  # parameters: graphlist (required)
  # updatespy: number of force-directed updates for each input graph
  # niter: number of iterations to run each force-directed algorithm for on each update - gets scaled depending on the method
  # wattr: string for the edge attribute name that should be used as edge weights - if left as NULL will result in no weights
  # initl: matrix of initial node positions - will be calculated by 'method' if not supplied
  # method: method to use for layout updates - must be one of 'fr' or 'kk'
  # wtmultiplier: multiplicative scaling factor for the edge weights
  
  retls = list()
  
  fullgraph = simplify(do.call(union,graphlist))
  K = length(graphlist)
  allnames = V(fullgraph)$name
  n = length(allnames)
  
  # matrix to store all layouts
  full_layout = matrix(nrow=n,ncol=2*(updatespy+1)*(K-1)+2)
  row.names(full_layout) = allnames
  
  tmpg = graphlist[[1]]
  if(!is.null(wattr)){
    wts = edge_attr(delete.vertices(tmpg,degree(tmpg)==0),wattr)
  } else { wts = NULL }
  isolated = V(tmpg)$name[which(degree(tmpg)==0)]
  
  # create an initial layout if one hasn't been provided
  if(is.null(initl)){
    if(method=="fr"){
      tmpl = layout_with_fr(delete.vertices(tmpg,degree(tmpg)==0),weights=wts)
    } else if(method=="kk"){
      tmpl = layout_with_kk(delete.vertices(tmpg,degree(tmpg)==0),weights=wts)
    } else {print("invalid method");return}
  } else { tmpl = initl }
  
  # add initial layout to full layout matrix
  rownames(tmpl) = setdiff(V(tmpg)$name,isolated)
  retls[[names(graphlist)[1]]] = tmpl
  if(length(isolated)>0){
    isomat = matrix(rep(0,2*length(isolated)),ncol=2)
    row.names(isomat) = isolated
    tmpl = rbind(tmpl,isomat)
  }
  tmpl = tmpl[order(row.names(tmpl)),]
  
  full_layout[row.names(tmpl),1:2] = tmpl
  colnames(full_layout) = c(paste0(names(graphlist)[1],"_x"),paste0(names(graphlist)[1],"_y"),
                            as.vector(sapply(names(graphlist)[-1],function(y) 
                              sapply(0:updatespy,function(z) c(paste0(y,"_x_",z+1),paste0(y,"_y_",z+1)) ) ) ))
  
  for(k in 2:K){
    print(k)
    tmpg2 = graphlist[[k]]
    # calculate new vertices whose positions need initialising
    newvs = setdiff(names(V(tmpg2)),names(V(tmpg)))
    oldvs = setdiff(names(V(tmpg)),names(V(tmpg2)))
    connected = names(V(tmpg))
    isolated = setdiff(allnames,names(which(degree(tmpg2)!=0)))
    
    tmpl2 = full_layout[,c(2*(updatespy+1)*(k-2)+1,2*(updatespy+1)*(k-2)+2)]
    
    # add new vertices according to their neighbours positions in the previous
    # layout, or according to their neighbours positions in the full graph if
    # there are none, or in neither apply then place in the corner
    newvnbrs = lapply(newvs,function(v) 
      intersect(neighbors(tmpg2,v)$name,connected))
    newvpos = lapply(newvnbrs,function(js) colMeans(matrix(tmpl[as.character(js),],ncol=2)))
    # may warn about NAs being introduced by coercion
    
    switch = c(-1,-1) # decides the corner, switches each time used
    for(i in seq_along(newvs)){
      j = newvs[i]
      if(any(is.na(newvpos[[i]]))){
        newvnbrs[[i]] = intersect(neighbors(fullgraph,j)$name,connected)
        newvpos[[i]] = colMeans(tmpl[newvnbrs[[i]],,drop=F])
      }
      if(any(is.na(newvpos[[i]]))){
        lims = apply(abs(tmpl),2,function(v) max(v,na.rm=T))
        newvpos[[i]] = switch*lims
        switch = switch %*% matrix(c(0,1,-1,0),nrow=2)
      }
    }
    
    # put the new positions in the matrix
    if(length(newvpos)>0){
      newvposmat = matrix(unlist(newvpos),ncol=2,byrow=T)
      tmpl2[newvs,] = newvposmat
    }
    
    for(i in 0:updatespy){
      posmat = tmpl2[V(tmpg2)$name,]
      
      tmpg2l = delete.vertices(tmpg2,degree(tmpg2)==0)
      if(!is.null(wattr)){
        wts = edge_attr(tmpg2l,wattr)
      } else { wts = NULL }
      
      # call a force directed algorithm for a layout update
      if(method=="fr"){
        retlk = layout_with_fr(tmpg2l,niter=niter,start.temp=vcount(tmpg2l)^0.2/4,
                               coords=posmat,weights=wts*wtmultiplier)
      } else if(method=="kk"){
        retlk = layout_with_kk(tmpg2l,maxiter=niter*sqrt(vcount(tmpg2l)),kkconst=sqrt(vcount(tmpg2l))*0.3,
                               coords=posmat,weights=wtmultiplier/(wts+0.1))
      } else {print("invalid method");break}
      
      # center the layout by subtracting the median node position
      retlk = sweep(retlk,2,apply(retlk,2,median,na.rm=T))
      rownames(retlk) = names(V(tmpg2l))
      retls[[paste0(names(graphlist)[k],"_",i+1)]] = retlk
      full_layout[V(tmpg2)$name,c(2*(updatespy+1)*(k-2)+3+2*i,2*(updatespy+1)*(k-2)+4+2*i)] = retlk
      tmpl2[V(tmpg2)$name,] = retlk
    }
    
    tmpl = tmpl2
    tmpg = tmpg2
  }
  
  retls[['full']] = full_layout
  return(retls)
}



# smoothing phase of the algorithm
# designed to work with output from the previous function
# layout list must contain a "full" layout
# names of graphlist and layoutlist must be integers (in string format) within a range given by startyear and endyear
splinelayoutsmoothing = function(graphlist,layoutlist,mindf=8,maxdf=24,startyear=1953,endyear=2018,graphsperyear=3){ 
  # parameters: graphlist, layoutlist (required)
  # mindf, maxdf: min and max degrees of freedom for the spline regression 
  # startyear, endyear: start and end years of the data in the graphs 
  #                       - workaround: name the lists 1:N, set startyear=1, endyear=N
  
  xsplines = list(); ysplines = list(); allcnbools = list()
  for(cn in names(V(simplify(do.call(igraph::union,graphlist))))){
    # a nodes dof will depend on how many graphs it appears in to try and prevent overfitting when there isnt data
    cnbools = sapply(graphlist,function(g) cn %in% names(V(g)))
    allcnbools[[cn]] = cnbools
    cnsplinedf = max(mindf,maxdf*sum(cnbools)/length(cnbools))
    
    # fit the splines
    k = graphsperyear
    cnxspline = lm(layoutlist[['full']][cn,c(T,F)] ~ bs(seq(startyear+(k-1)/k,endyear+(k-1)/k,1/k),df=cnsplinedf))
    cnyspline = lm(layoutlist[['full']][cn,c(F,T)] ~ bs(seq(startyear+(k-1)/k,endyear+(k-1)/k,1/k),df=cnsplinedf))
    xsplines[[cn]] = cnxspline; ysplines[[cn]] = cnyspline
  }
  
  # get fitted values from the splines and create new layouts
  newxvals = sapply(xsplines,function(s) predict(s,newdata=list()))
  newyvals = sapply(ysplines,function(s) predict(s,newdata=list()))
  newvals = cbind(newxvals[i,],newyvals[i,])
  retls = lapply(1:nrow(newxvals),function(i) 
    t(rbind(newxvals[i,rownames(layoutlist[[i]])],newyvals[i,rownames(layoutlist[[i]])])) )
  
  return(retls)
}




# interpolation phase of the algorithm
interpolateframes = function(l1,l2,n){
  # creates a list of layouts with node positions linearly interpolated between l1 and l2
  # considers a layout on the union of vertices in l1 and l2
  
  l1tmp = l1
  l1 = rbind(l1,l2[setdiff(row.names(l2),row.names(l1)),,drop=F])
  l2 = rbind(l2,l1tmp[setdiff(row.names(l1tmp),row.names(l2)),,drop=F])
  vec1 = as.vector(l1); vec2 = as.vector(l2[row.names(l1),])
  return(lapply(0:(n-1),function(x) matrix(vec1 + (x/n)*(vec2-vec1),ncol=2,
                                           dimnames=list(rownames(l1)))    ))
}




# plots all the animation frames (or slices) to be converted to video
makegraphframes = function(graphlist,layoutlist,datalist=NULL,path_="plots/slicespng/slice",
                           no_iframes=6,vsf=1,vaf=1,vtf=1,res_=80,edge_alpha=1,fullxlims=NA,fullylims=NA,names_to_highlight=NULL){
  # parameters: graphlist, layoutlist, datalist (unused), path
  # no_iframes: number of interpolation frames between each layout
  # vsf: multiplicative scaling factor for vertex size, vaf: same for opacity, vtf: same for label text (stacks with vsf)
  # res: png resolution, edge_alpha: scaling factor for edge opacity
  # fullxlims, fullylims: override the automatically set frame boundaries
  # names_to_highlight: vector of vertex names which will be highlighted in a different colour
  
  # calculate frame boundaries if none were manually supplied
  # will throw up warnings which can be ignored 
  if(is.na(fullxlims)){
    fullxlims = c(quantile(sapply(layoutlist,function(l) quantile(l[,1],0.01,na.rm=T)),0.01),
                  quantile(sapply(layoutlist,function(l) quantile(l[,1],0.99,na.rm=T)),0.99))*1.05 }
  if(is.na(fullylims)){
    fullylims = c(quantile(sapply(layoutlist,function(l) quantile(l[,2],0.01,na.rm=T)),0.01),
                  quantile(sapply(layoutlist,function(l) quantile(l[,2],0.99,na.rm=T)),0.99))*1.05 }
  
  frameno = 1 # frame counter
  for(i in 1:(length(graphlist))){
    print(paste("i:",sprintf("%02d",i),"  year:",names(graphlist)[i]))
    # this function takes a long time to render all the outputs, some printing is needed for reassurance that it's working
    g = graphlist[[i]]
    # opacity scales exponentially with edge attribute rank
    # where there are no ranks in the graph I set a rank parameter that is 1 for every edge - workaround
    edgealphas = 200*(2/3)^(edge_attr(g)$ranks-1)
    
    for(j in 0:2){
      k = i*3+j-4 + 2*(i==1)
      if(k==length(layoutlist)-1){print("end");break}
      l = layoutlist[[k]]
      l2 = layoutlist[[k+1]]
      
      # compute and print out of bounds vertices
      badvs = l[,1] < fullxlims[1] | l[,1] > fullxlims[2] |
        l[,2] < fullylims[1] | l[,2] > fullylims[2]
      if(sum(badvs)>0){print(paste0("oob vertices for k = ",k,": ",paste(names(V(g))[badvs],collapse=" ")))}
      
      # create interpolated layouts
      interfs = interpolateframes(l,l2,no_iframes)
      
      if(j < 2 && i !=1){ # layout steps where vertices cannot be added/removed
        for(f in interfs){
          # if applicable, get numeric vertex ids of vertices to be highlighted
          hl_vs = as.numeric(sapply(names(V(g)),function(x) x %in% names_to_highlight))+1
          
          # render
          png(paste0(path_,sprintf("%04d",frameno),".png"),width=12.5,height=7.5,units="in",res=res_)
          plot(g,edge.color=rgb(235,128,27,edge_alpha*edgealphas,maxColorValue=255),
               vertex.color=rgb(c(55,224)[hl_vs],c(110,79)[hl_vs],c(191,117)[hl_vs],rep(250*vaf,vcount(g)),maxColorValue=255),
               vertex.frame.color=gray(rep(0.2,vcount(g))),vertex.label.color=gray(rep(0.1,vcount(g))),
               layout=f,rescale=F,edge.width=E(g)$width,xlim=fullxlims*1.02,ylim=fullylims*1.02,
               vertex.label.family="sans",vertex.size=25*vsf+5*hl_vs,vertex.label.cex=0.5*vtf*vsf+0.05*hl_vs,asp=0)
          text(fullxlims[1]*0.92,fullylims[2]*0.92,names(graphlist)[i],cex=1)
          dev.off()
          print(frameno)
          frameno=frameno+1
        }
      } else if((i < length(graphlist) && j == 2) || i == 1){ 
        # layout steps where vertices could be added/removed as the underlying graph is changing
        
        g2 = graphlist[[i+1]]
        gu = igraph::union(g,g2) # we require both graphs in order to transition between them
        # get old and new vertex/edge opacities
        eas = 200*(2/3)^(edge_attr(gu)$ranks_1-1); eas[is.na(eas)] = 0
        eas2 = 200*(2/3)^(edge_attr(gu)$ranks_2-1); eas2[is.na(eas2)] = 0
        ews = edge_attr(gu)$width_1; ews[is.na(ews)] = 0
        ews2 = edge_attr(gu)$width_2; ews2[is.na(ews2)] = 0
        oldvs = setdiff(names(V(g)),names(V(g2))); oldvbs = names(V(gu)) %in% oldvs
        newvs = setdiff(names(V(g2)),names(V(g))); newvbs = names(V(gu)) %in% newvs
        oldes = incident(gu,oldvs,"all"); oldebs = E(gu) %in% oldes
        newes = incident(gu,newvs,"all"); newebs = E(gu) %in% newes
        
        for(x in seq_along(interfs)){
          f = interfs[[x]]
          f = f[names(V(gu)),]
          # interpolate between old and new vertex/edge opacities
          newedgealphas = eas2*x/length(interfs) + eas*(1-x/length(interfs))
          newedgewidths = ews2*x/length(interfs) + ews*(1-x/length(interfs))
          vertexalphas = (!newvbs & !oldvbs) + (oldvbs*(1-x/length(interfs))) + (newvbs*x/length(interfs))
          hl_vs = as.numeric(sapply(names(V(gu)),function(x) x %in% names_to_highlight))+1
          
          # render
          png(paste0(path_,sprintf("%04d",frameno),".png"),width=12.5,height=7.5,units="in",res=res_)
          plot(gu,edge.color=rgb(235,128,27,edge_alpha*newedgealphas,maxColorValue=255),
               vertex.color=rgb(c(55,224)[hl_vs],c(110,79)[hl_vs],c(191,117)[hl_vs],250*vaf*vertexalphas,maxColorValue=255),
               vertex.frame.color=gray(rep(0.2,length(vertexalphas)),alpha=vertexalphas),
               vertex.label.color=gray(rep(0.1,length(vertexalphas)),alpha=vertexalphas),layout=f,
               rescale=F,edge.width=newedgewidths,xlim=fullxlims*1.02,ylim=fullylims*1.02,
               vertex.label.family="sans",vertex.size=25*vsf+5*hl_vs,vertex.label.cex=0.5*vtf*vsf+0.05*hl_vs,asp=0)
          text(fullxlims[1]*0.92,fullylims[2]*0.92,names(graphlist)[i],cex=1)
          dev.off()
          print(frameno)
          frameno=frameno+1
        }
      }
      
      if(i == 1){break}
    }
  }
  # function end
}

# example run
dotdynls = layout_dynamically(dotgraphs,niter=2,updatespy=2,wattr="width",initl=dotgls[[1]]*2)
dotdynlssmooth = splinelayoutsmoothing(dotgraphs,dotdynls,mindf=6,maxdf=16)
makegraphframes(dotgraphs,dotdynlssmooth,datalist=dotyears,path_="plots/slicespng/slice",res_=80,vsf=0.9,vaf=0.9,edge_alpha=0.75,fullxlims=c(-4.5,4.5),fullylims=c(-3.5,3.5))


# run on command line in the folder the images are in to generate video (needs ffmpeg installed)
# ffmpeg -y -framerate 24 -i slice%04d.png -c:v libx264 -pix_fmt yuv420p -an -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" slicevideo24.mp4




