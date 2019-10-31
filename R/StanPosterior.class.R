StanPosterior.class<-function(stan.model,all.only=T){
  res<-get_posterior_mean(stan.model,pars = c("Vc"))
  nc<-nrow(res)
  class.name<-unique(apply(combn(rep(c(0,1),log2(nc)),log2(nc)),2,function(x){paste(x,collapse = "")}))
  rownames(res)<-class.name
  if(all.only){
    res<-res[,ncol(res)]
    
  }
  res
}
