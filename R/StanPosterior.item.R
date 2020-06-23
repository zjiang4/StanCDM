StanPosterior.item<-function(stan.model,lambda.parm=TRUE){
  x<-stan.model
  summaryStan<-summary(x)$'summary'
  temp<-!grepl("PImat|log_lik|contributionsI|posteriorIC|Vc|posteriorPC|lp",rownames(summaryStan))  
  summaryStan<-summaryStan[temp,]
  if(lambda.parm==TRUE){
    summaryStan[grepl("l",rownames(summaryStan)),]
  }else{
    summaryStan[!grepl("l",rownames(summaryStan)),]
  }
}

