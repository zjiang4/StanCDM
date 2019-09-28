StanPosterior.attribute<-function(stan.model,all.only=T){
  nc=length(StanPosterior.class(stan.model))
  posteriorAttr<-get_posterior_mean(stan.model,pars = c("posteriorPC"))
  np=nrow(posteriorAttr)/nc
  if(all.only){
    posteriorAttr<-unlist(lapply(1:np,function(x){which.max(matrix(posteriorAttr[,ncol(posteriorAttr)],np,nc,byrow = T)[x,])}))
  }
  posteriorAttr
}