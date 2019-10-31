StanPosterior.attribute<-function(stan.model,all.only=T){
  nc=length(StanPosterior.class(stan.model))
  class.name<-unique(apply(combn(rep(c(0,1),log2(nc)),log2(nc)),2,function(x){paste(x,collapse = "")}))
  posteriorAttr<-get_posterior_mean(stan.model,pars = c("posteriorPC"))
  rowNames<-row.names(posteriorAttr)
  newRowNames<-rowNames
  for(i in 1: length(newRowNames)){
    id<-str_match_all(rowNames[i],"[0-9]+")[[1]][,1]
    newRowNames[i]<-paste(paste("p",id[1],sep=''),paste("c",class.name[as.numeric(as.character(id[2]))],sep=''),sep='')
  } 
  np=nrow(posteriorAttr)/nc
  if(all.only){
    posteriorAttr<-unlist(lapply(1:np,function(x){which.max(matrix(posteriorAttr[,ncol(posteriorAttr)],np,nc,byrow = T)[x,])}))
    
    posteriorAttr<-class.name[posteriorAttr]
  }
  if(!all.only){
    row.names(posteriorAttr)<-newRowNames
    
  }
  posteriorAttr
}