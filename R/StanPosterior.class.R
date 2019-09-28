StanPosterior.class<-function(stan.model,all.only=T){
  res<-get_posterior_mean(stan.model,pars = c("Vc"))
  if(all.only){
    res<-res[,ncol(res)]
  }
  res
}