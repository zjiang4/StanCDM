
StanLCDM_gC.run<-function(Qmatrix,response.matrix,SubjectID=NA,TimestampID=NA,quad.structure=F,variance.equal=F,
                          script.path=NA,save.path=getwd(),save.name="LCDM_gC_uninf",iter=1000,warmup = 0,
                          chain.num=3,init.list='random',control.list=NA){
  
  
  
  data.list<-Generate.datalist(Qmatrix,respMatrix,SubjectID=SubjectID,TimestampID=TimestampID)
  time.vector<-data.list$TimeVec
  if(is.na(control.list)){control.list<-list(adapt_delta=0.82)}
  if(is.na(script.path)==T){
    options(warn=-1)
    StanLCDM_gC.script(Qmatrix,time.vector=time.vector,quad.structure=quad.structure,variance.equal=variance.equal, save.path=save.path,save.name="LCDM_gC_uninf")
    script.path<-paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
    options(warn=0)
    compiled_model<-stan_model(script.path)
  }else{
    compiled_model<-stan_model(script.path)
  }
  
  estimated_model<-tryCatch(sampling(compiled_model,
                                     data = data.list,
                                     iter = iter,
                                     init = init.list,
                                     warmup = warmup,
                                     chains=chain.num,
                                     control=control.list),
                            error=function(e){"The estimation process is terminated with errors"})
  estimated_model
  
}