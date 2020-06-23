#targetPath<-"\\\\libfs1\\users$\\home\\zjiang17\\My Documents\\GitHub\\StanDCM\\R\\";setwd(targetPath)
#AllFiles<-list.files(targetPath)
#AllFiles<-AllFiles[AllFiles!="StanDCM.R"]

#for(i in 1:length(AllFiles)){
 #  line <- readLines(AllFiles[i])
#  write(line, "StanDCM.txt",append = T)
#  write("\n#SEPERATION#", "StanDCM.txt",append = T)
#}

#SEPERATION#
#' @title Generate data list
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param response.matrix save the .stan file to somewhere; the default path is getwd()
#' @param GroupID name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\data.RData") ;Qmatrix<-cbind(Qmatrix,rep(1,9));Qmatrix[1,1]<-0



Generate.datalist<-function(Qmatrix,response.matrix,GroupID=NA,SubjectID=NA,TimestampID=NA,quad.structure=F){
  if(quad.structure==T){Nz=3}else{Nz=2}
  if(!is.na(TimestampID)[1]&!is.na(SubjectID)[1]){
    if(sum(is.na(response.matrix))!=0){print('Stop!The response dataset contains missing value(s)')}else{
      Generate.dataList<-list(Y=response.matrix,
                              Nr=nrow(response.matrix),
                              Na=ncol(Qmatrix),
                              Np = length(unique(SubjectID)),
                              Ni = ncol(response.matrix),
                              Nz = Nz,
                              No = length(unique(TimestampID)),
                              SubjectID = SubjectID,
                              OccasionID = rank(unique(TimestampID))[match(TimestampID,unique(TimestampID))],
                              TimeVec = unique(TimestampID)[order(unique(TimestampID))]
      )
    }
  }else{
    if(sum(is.na(response.matrix))!=0){print('Stop!The response dataset contains missing value(s)')}else{
      if( length((unique(unlist(unique(c(response.matrix))))))>2 ){
        Generate.dataList<-list(Y=response.matrix, Na=ncol(Qmatrix),
                                Np = nrow(response.matrix),
                                Ni = ncol(response.matrix),Nc=2^(ncol(Qmatrix)),
                                Ns = length((unique(unlist(unique(c(response.matrix)))))))
      }else{
        Generate.dataList<-list(Y=response.matrix, Na=ncol(Qmatrix),
                                Np = nrow(response.matrix),
                                Ni = ncol(response.matrix),Nc=2^(ncol(Qmatrix)))
      }
    }
    if(!is.na(GroupID)[1]){Generate.dataList$GroupID=GroupID}
    
    
    
    
  }
  if(quad.structure==T){Generate.dataList$TimeQuadVec=
    (unique(TimestampID)[order(unique(TimestampID))])^2}
  Generate.dataList
}

#SEPERATION#
#' @title Install and load packages.
#'
#' @description
#' \code{Install.package} allows other functions to install and load packages.
#'
#' @param needed_packages the packages to be installed and loaded.
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @examples
#' \dontrun{
#' Install.package("ggplot2")
#' }
#'

Install.package = function(needed_packages){
  for (i in 1:length(needed_packages)){
    haspackage = require(needed_packages[i], character.only = TRUE)
    if (haspackage==FALSE){
      install.packages(needed_packages[i])
      require(needed_packages[i], character.only = TRUE)
    }
  }
}


StanPosterior.item<-function(stan.model,lambda.parm=TRUE){
  x<-stan.model
  summaryStan<-summary(x)$'summary'
  temp<-!grepl("PImat|log_lik|contributionsI|posteriorIC|Vc|posteriorPC|lp",rownames(summaryStan))  
  summaryStan<-summaryStan[temp,]
  if(lambda.parm=TRUE){
    summaryStan[grepl("l",rownames(summaryStan)),]
  }else{
    summaryStan[!grepl("l",rownames(summaryStan)),]
  }
}



#SEPERATION#
#' @title Generate Stan code and Run the estimation for LCDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param savepath save the .stan file to somewhere; the default path is getwd()
#' @param savename name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\data.RData") ;Qmatrix<-cbind(Qmatrix,rep(1,9));Qmatrix[1,1]<-0


Parm.name<-function(Qmatrix){

  #Load packages
  Install.package("plyr")
  Install.package('stringr')


  nc<-ncol(Qmatrix)
  nr<-nrow(Qmatrix)
  temp.table.col<-unique(apply(combn(rep(c(0,1),nc),nc),2,function(x){paste(x,collapse = "")}))
  temp.table.col<-temp.table.col[order(temp.table.col)]
  temp.table<-matrix(0,nr,length(temp.table.col))
  colnames(temp.table)<-temp.table.col
  rownames(temp.table)<-paste('item',c(1:nr),sep='')
  temp.table<-as.data.frame(temp.table)
  for (i in 1:nr){
    temp.table[i,]<-paste('l',i,'_0',sep='')
  }
  intercept<-temp.table[,1]

  #Generate attribute combinations
  comb.generator<-function(x.vector){
    if(length(x.vector)>1){
      temp.attr<-x.vector
      temp.attr.sav<-NULL
      for(i in 1:length(temp.attr)){
        temp.1<-combn(temp.attr,i)
        temp.2<-apply(temp.1,2,function(x){paste(x,collapse = "")})
        temp.attr.sav<-c(temp.attr.sav,temp.2)
      }
    }
    if(length(x.vector)==1){temp.attr.sav<-x.vector}
    temp.attr.sav
  }
  #vectors needed for combination.generator
  Item.load.id<-list()
  for ( i in 1:nr){
    Item.load.id[[i]]<-grep('1',Qmatrix[i,])}

  Attr.load.id<-list()
  attr.load.id<-matrix(0,length(temp.table.col),nc)
  for ( i in 1:length(temp.table.col)){
    attr.load.id[i,]<-unlist(strsplit(temp.table.col[i],split=''))
    Attr.load.id[[i]]<-grep('1',attr.load.id[i,])
  }

  #Generate Combination for both Item.load and Attr.load
  Item.Comb<-list()
  for ( i in 1:nr){
    Item.Comb[[i]]<-comb.generator(Item.load.id[[i]])
  }
  Attr.Comb<-list()
  for ( i in 2:length(temp.table.col)){
    Attr.Comb[[1]]<-0
    Attr.Comb[[i]]<-comb.generator(Attr.load.id[[i]])
  }
  constraints.list<-list()
  nway.inter.list<-list()
  for(i in 1:nr){
    for(a in 2:length(temp.table.col)){
      ifzero<-as.numeric(paste(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],collapse=''))
      if((!is.na(ifzero))){
        temp.table[i,a]<-paste(c(temp.table[i,a],
                                 paste("S","l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='',collapse='')
        ),collapse='')
        if(a==length(temp.table.col)){
          nway.inter.list[[i]]<-nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])])
          constraints.list[[i]]<-paste("l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='')
        }
      }
    }
  }

  #Create Lambda Table
  Lamda.Table<-temp.table
  for(i in 1:nr){
    for(a in 1:length(Lamda.Table)){
      t.ref<-unique(as.character(Lamda.Table[i,]))
      pos<-c(1:length(t.ref))[Lamda.Table[i,a]==t.ref]
      temp.table[i,a]<-paste("t",i,"_",pos,sep='')}}

  #Generate LCDM specification
  out<-list()
  out[[1]]<-Lamda.Table
  out[[2]]<-temp.table
  out[[3]]<-constraints.list
  out[[4]]<-nway.inter.list
  out[[5]]<-intercept
  OUTPUT<-out
  nclass<-ncol(OUTPUT[[1]]);Nc<-nclass

  #Produce kernel expressions across items and attributes
  Kernel.exp<-OUTPUT[[1]]
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])}
    }
  }

  Classp.exp1<-colnames(OUTPUT[[1]])
  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                        nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))
  if(length(name.inter)!=0){
    for (inter in 1: length(name.inter)){
      temp.nw<-numway.inter[inter]
      temp.nm<-name.inter[inter]
      temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
      temp.sel<-paste(unlist(strsplit(temp.nm,split = '_'))[1],"_",(1:(temp.nw-1)),sep='')
      first.sel<-unlist(OUTPUT[[3]])[grep(paste((temp.sel),collapse="|"),unlist(OUTPUT[[3]]))]
      second.sel<-sub(".*_.", "", first.sel)
      for (sel in 1:length(temp.subnm)){
        SEL<-second.sel[sel]
        Constrain.List1<-rbind(
          paste(temp.nm,">-(0", paste("+",first.sel[grep(SEL,second.sel)],
                                      sep='',collapse=''),")",sep=''),Constrain.List1)
      }
    }
    Constrain.List1<-as.character(Constrain.List1)
  }else{
    Constrain.List1<-NULL
  }



  itemParmName<-c(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==2],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==3],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==4],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==5],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==6],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==7],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==8],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==9],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==10],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==11],OUTPUT[[5]])
  numMainEffect<-length(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1])
  Constrain.List<-paste('  real<lower=0>',itemParmName[1:numMainEffect],';\n ')
  Unconstrain.List<-paste('  real',itemParmName[-(1:numMainEffect)],';\n ')
  Reparm<-as.data.frame(matrix(0,nr,nclass))

  trueParmName<-itemParmName
  out.list<-list()
  out.list$parm.name<-trueParmName
  out.list$intercept.name<-OUTPUT[[5]]
  out.list$main.name<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1]
  out.list$interaction.name<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])!=1]
  out.list$class.expression<-Classp.exp1
  out.list
}




#SEPERATION#
#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanCRUM.run<-function(Qmatrix,response.matrix,script.path=NA,save.path=getwd(),save.name="CRUM_uninf",iter=1000,warmup = 0,
                       chain.num=3,init.list='random',control.list=NA){
  rstan.detect<-tryCatch(library("rstan"),error=function(e){"rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details."})
  if(length(rstan.detect)==1){
    #break
    stop()
  }
  Cdm.init<-F
  if(init.list=='cdm'){
    Cdm.init<-T
    Install.package(c("CDM","stringr"))
    trueParmName<-Parm.name(Qmatrix=Qmatrix)$parm.name
    Classp.exp1<-Parm.name(Qmatrix=Qmatrix)$class.expression
    mod1<-gdina( data =respMatrix, q.matrix = Qmatrix , maxit=700,link = "logit",progress=F)
    CDMresult<-as.data.frame(coef(mod1))
    library(stringr)
    CDM.parm.name<-paste(paste(paste('l',CDMresult[,3],sep=''),'_',sep=''),str_count(CDMresult$partype.attr,"Attr"),sep='')
    CDM.parm.name<-paste(CDM.parm.name,
                         unlist(lapply(strsplit(unlist(lapply(strsplit(CDMresult$partype.attr, 'Attr', fixed=FALSE),function(x){paste(x,collapse="")})),'-'),function(x){paste(x,collapse="")})),
                         sep='')
    CDM.parm.est<-CDMresult$est
    parm.ini<-round(CDM.parm.est[match(trueParmName,CDM.parm.name)],4)
    CDM.prop.est<-mod1$attribute.patt
    prop.ini<-CDM.prop.est[match(Classp.exp1,rownames(CDM.prop.est)),1]
    inilist1<-paste('list(',paste(noquote(paste(noquote(unlist(list(
      paste(paste('Vc=c(',paste((prop.ini),collapse=','),')',collapse=','))))
    ))),collapse=',') ,')',collapse='')

    inilist1<-eval(parse(n =2000000 ,text=inilist1))
    for( i in 2:chain.num){
      temp.text<-paste('inilist',i,"<-inilist1",sep='')
      eval(parse(text=(temp.text)))
    }
    temp.text<-paste('init.list<-list(',paste(paste('inilist',1:chain.num,sep=''),collapse = ","),')',sep='')
    eval(parse(text=(temp.text)))
  }
  data.list<-Generate.datalist(Qmatrix,response.matrix)

  if(is.na(control.list)){control.list<-list(adapt_delta=0.82)}
  if(is.na(script.path)==T){
    options(warn=-1)
    StanCRUM.script(Qmatrix,save.path=save.path,save.name=save.name)
    script.path<-paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
    options(warn=0)
    compiled_model<-stan_model(script.path)
  }else{
    compiled_model<-stan_model(script.path)
  }
  if(Cdm.init==T){
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})
  }else{
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})

  }

  estimated_model
}


#SEPERATION#
#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanCRUM.script<-function(Qmatrix,save.path=getwd(),save.name="CRUM_uninf"){
  #Load packages
  Install.package("plyr")
  Install.package('stringr')

  nc<-ncol(Qmatrix)
  nr<-nrow(Qmatrix)
  temp.table.col<-unique(apply(combn(rep(c(0,1),nc),nc),2,function(x){paste(x,collapse = "")}))
  temp.table.col<-temp.table.col[order(temp.table.col)]
  temp.table<-matrix(0,nr,length(temp.table.col))
  colnames(temp.table)<-temp.table.col
  rownames(temp.table)<-paste('item',c(1:nr),sep='')
  temp.table<-as.data.frame(temp.table)
  for (i in 1:nr){
    temp.table[i,]<-paste('l',i,'_0',sep='')
  }
  intercept<-temp.table[,1]

  #Generate attribute combinations
  comb.generator<-function(x.vector){
    if(length(x.vector)>1){
      temp.attr<-x.vector
      temp.attr.sav<-NULL
      for(i in 1:length(temp.attr)){
        temp.1<-combn(temp.attr,i)
        temp.2<-apply(temp.1,2,function(x){paste(x,collapse = "")})
        temp.attr.sav<-c(temp.attr.sav,temp.2)
      }
    }
    if(length(x.vector)==1){temp.attr.sav<-x.vector}
    temp.attr.sav
  }
  #vectors needed for combination.generator
  Item.load.id<-list()
  for ( i in 1:nr){
    Item.load.id[[i]]<-grep('1',Qmatrix[i,])}

  Attr.load.id<-list()
  attr.load.id<-matrix(0,length(temp.table.col),nc)
  for ( i in 1:length(temp.table.col)){
    attr.load.id[i,]<-unlist(strsplit(temp.table.col[i],split=''))
    Attr.load.id[[i]]<-grep('1',attr.load.id[i,])
  }

  #Generate Combination for both Item.load and Attr.load
  Item.Comb<-list()
  for ( i in 1:nr){
    Item.Comb[[i]]<-comb.generator(Item.load.id[[i]])
  }
  Attr.Comb<-list()
  for ( i in 2:length(temp.table.col)){
    Attr.Comb[[1]]<-0
    Attr.Comb[[i]]<-comb.generator(Attr.load.id[[i]])
  }
  constraints.list<-list()
  nway.inter.list<-list()
  for(i in 1:nr){
    for(a in 2:length(temp.table.col)){
      ifzero<-as.numeric(paste(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],collapse=''))
      if((!is.na(ifzero))){
        temp.table[i,a]<-paste(c(temp.table[i,a],
                                 paste("S","l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='',collapse='')
        ),collapse='')
        if(a==length(temp.table.col)){
          nway.inter.list[[i]]<-nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])])
          constraints.list[[i]]<-paste("l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='')
        }
      }
    }
  }

  #Create Lambda Table
  Lamda.Table<-temp.table
  for(i in 1:nr){
    for(a in 1:length(Lamda.Table)){
      t.ref<-unique(as.character(Lamda.Table[i,]))
      pos<-c(1:length(t.ref))[Lamda.Table[i,a]==t.ref]
      temp.table[i,a]<-paste("t",i,"_",pos,sep='')}}

  #Generate LCDM specification
  out<-list()
  out[[1]]<-Lamda.Table
  out[[2]]<-temp.table
  out[[3]]<-constraints.list
  out[[4]]<-nway.inter.list
  out[[5]]<-intercept
  OUTPUT<-out
  nclass<-ncol(OUTPUT[[1]]);Nc<-nclass

  #Produce kernel expressions across items and attributes
  Kernel.exp<-OUTPUT[[1]]
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])}
    }
  }


  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                        nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))

  if(length(name.inter)!=0){
    for (inter in 1: length(name.inter)){
      temp.nw<-numway.inter[inter]
      temp.nm<-name.inter[inter]
      temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
      temp.sel<-paste(unlist(strsplit(temp.nm,split = '_'))[1],"_",(1:(temp.nw-1)),sep='')
      first.sel<-unlist(OUTPUT[[3]])[grep(paste((temp.sel),collapse="|"),unlist(OUTPUT[[3]]))]
      second.sel<-sub(".*_.", "", first.sel)
      for (sel in 1:length(temp.subnm)){
        SEL<-second.sel[sel]
        Constrain.List1<-rbind(
          paste(temp.nm,">-(0", paste("+",first.sel[grep(SEL,second.sel)],
                                      sep='',collapse=''),")",sep=''),Constrain.List1)
      }
    }
    Constrain.List1<-as.character(Constrain.List1)
  }else{
    Constrain.List1<-NULL
  }

  itemParmName<-c(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==2],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==3],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==4],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==5],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==6],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==7],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==8],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==9],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==10],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==11],OUTPUT[[5]])
  numMainEffect<-length(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1])
  Constrain.List<-paste('  real<lower=0>',itemParmName[1:numMainEffect],';\n ')
  Unconstrain.List<-paste('  real',intercept,';\n ')
  Reparm<-as.data.frame(matrix(0,nr,nclass))


  #Produce Stan code for PImat parameter
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      Reparm[loopi,loopc]<-paste('  PImat[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp[loopi,loopc]),');\n',sep='')
    }
  }

  ####################053119update####################
  for(loopi in 1:nr){
    for(loopc in 1:Nc){
      Reparm[loopi,loopc]<- str_replace_all(Reparm[loopi,loopc],paste((name.inter),collapse="|"),"0")
    }
  }
  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',c(itemParmName[1:numMainEffect],intercept),'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))
  ####################053119update END  ##############

  #Likelihood Stan code
  Likelihood<-'
  \n
  //Likelihood
  for (iterp in 1:Np){
    for (iterc in 1:Nc){
      for (iteri in 1:Ni){
        if (Y[iterp,iteri] == 1)
          contributionsI[iteri]=bernoulli_lpmf(1|PImat[iteri,iterc]);
        else
          contributionsI[iteri]=bernoulli_lpmf(0|PImat[iteri,iterc]);
      }
      contributionsC[iterc]=log(Vc[iterc])+sum(contributionsI);
    }
  target+=log_sum_exp(contributionsC);
  }
  '


  #Data Specification
  data.spec<-'
data{
  int Np;
  int Ni;
  int Nc;
  matrix[Np, Ni] Y;
}
  '

#Parameter Specification
parm.spec<-paste(c('
parameters{
  simplex[Nc] Vc;\n ',paste0(Constrain.List),paste0(Unconstrain.List),
                   '}\n'),collapse='')

#Reparameter Specification
transparm.spec<-paste(c('
  transformed parameters{
  matrix[Ni, Nc] PImat;\n',
                        paste0(unlist(Reparm)),'}\n'),collapse='')

#Model Specification
model.spec<-paste(c('\nmodel {\n',paste(c(Modelcontainer,Parmprior,Likelihood),sep=''),'\n}',sep=''))
model.spec<-model.spec[!startsWith(str_remove_all(model.spec," "),"~")]

#Generated Quantities Specification
generatedQuantities.spec<-'
  \n
generated quantities {

 vector[Ni] log_lik[Np];
 vector[Ni] contributionsI;
 matrix[Ni,Nc] contributionsIC;
 
 matrix[Ni,Nc] posteriorIC;
 matrix[Np,Nc] posteriorPC;



 //Posterior
 for (iterp in 1:Np){
   for (iteri in 1:Ni){
     for (iterc in 1:Nc){
       if (Y[iterp,iteri] == 1)
          contributionsI[iteri]=bernoulli_lpmf(1|PImat[iteri,iterc]);
       else
           contributionsI[iteri]=bernoulli_lpmf(0|PImat[iteri,iterc]);
       contributionsIC[iteri,iterc]=log(Vc[iterc])+contributionsI[iteri];
       posteriorIC[iteri,iterc]=contributionsI[iteri];
      }
      log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
    }
   for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
  }
}
'
if (.Platform$OS.type == "unix") {
  filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
}else{
  filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
}

sink(file= filename,append=FALSE)
cat(
  paste(c('   ',
          data.spec,parm.spec,transparm.spec,model.spec,generatedQuantities.spec)
  ))
sink(NULL)

}


#SEPERATION#
#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanCRUM_mG.run<-function(Qmatrix,
                          response.matrix,
                          GroupID,
                          fixeditem.vector=NA,
                          class.equal=T,
                          script.path=NA,save.path=getwd(),save.name="CRUM_uninf_multiG",
                          iter=1000,warmup = 0,
                          chain.num=3,init.list='random',control.list=NA){
  group.num<-length(unique(GroupID))
  rstan.detect<-tryCatch(library("rstan"),error=function(e){"rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details."})
  if(length(rstan.detect)==1){
    stop()
  }
  Cdm.init<-F
  if(init.list=='cdm'){
    Cdm.init<-T
    Install.package(c("CDM","stringr"))
    trueParmName<-Parm.name(Qmatrix=Qmatrix)$parm.name
    Classp.exp1<-Parm.name(Qmatrix=Qmatrix)$class.expression
    mod1<-gdina( data =response.matrix, q.matrix = Qmatrix , maxit=700,link = "logit",progress=F)
    CDMresult<-as.data.frame(coef(mod1))
    library(stringr)
    CDM.parm.name<-paste(paste(paste('l',CDMresult[,3],sep=''),'_',sep=''),str_count(CDMresult$partype.attr,"Attr"),sep='')
    CDM.parm.name<-paste(CDM.parm.name,
                         unlist(lapply(strsplit(unlist(lapply(strsplit(CDMresult$partype.attr, 'Attr', fixed=FALSE),function(x){paste(x,collapse="")})),'-'),function(x){paste(x,collapse="")})),
                         sep='')
    CDM.parm.est<-CDMresult$est
    parm.ini<-round(CDM.parm.est[match(trueParmName,CDM.parm.name)],4)
    CDM.prop.est<-mod1$attribute.patt
    prop.ini<-CDM.prop.est[match(Classp.exp1,rownames(CDM.prop.est)),1]
    inilist1<-paste('list(',paste(noquote(paste(noquote(unlist(list(
      paste(paste('Vc=c(',paste((prop.ini),collapse=','),')',collapse=','))))
    ))),collapse=',')  ,')',collapse='')

    IniList1<-NULL
    if(!class.equal){
      temp.inilist1<-eval(parse(n =2000000 ,text=inilist1))
      eval(parse(n =2000000 ,text=paste(paste('IniList1$Vc_g',1:group.num,sep=''),"<-temp.inilist1$Vc",sep='') ))
      inilist1<-IniList1
    }else{
      inilist1<-eval(parse(n =2000000 ,text=inilist1))}


    for( i in 2:chain.num){
      temp.text<-paste('inilist',i,"<-inilist1",sep='')
      eval(parse(text=(temp.text)))
    }
    temp.text<-paste('init.list<-list(',paste(paste('inilist',1:chain.num,sep=''),collapse = ","),')',sep='')
    eval(parse(text=(temp.text)))
  }
  data.list<-Generate.datalist(Qmatrix,response.matrix,GroupID)

  if(is.na(control.list)){control.list<-list(adapt_delta=0.82)}
  if(is.na(script.path)==T){
    options(warn=-1)
    #Need to update script
    StanCRUM_mG.script(Qmatrix=Qmatrix,
                       group.num=group.num,
                       fixeditem.vector=fixeditem.vector,
                       class.equal=class.equal,
                       save.path=save.path,save.name=save.name)
    script.path<-paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
    options(warn=0)
    compiled_model<-stan_model(script.path)}
  else{
    compiled_model<-stan_model(script.path)
  }
  if(Cdm.init==T){
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})
  }else{
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})

  }

  estimated_model
}

#SEPERATION#
#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanCRUM_mG.script<-function(Qmatrix,
                             group.num,
                             fixeditem.vector=NA,
                             class.equal=T,
                             save.path=getwd(),save.name="CRUM_uninf_multiG"){

  #Load packages
  Install.package("plyr")
  Install.package('stringr')

  nc<-ncol(Qmatrix)
  nr<-nrow(Qmatrix)
  temp.table.col<-unique(apply(combn(rep(c(0,1),nc),nc),2,function(x){paste(x,collapse = "")}))
  temp.table.col<-temp.table.col[order(temp.table.col)]
  temp.table<-matrix(0,nr,length(temp.table.col))
  colnames(temp.table)<-temp.table.col
  rownames(temp.table)<-paste('item',c(1:nr),sep='')
  temp.table<-as.data.frame(temp.table)
  for (i in 1:nr){
    temp.table[i,]<-paste('l',i,'_0',sep='')
  }
  intercept<-temp.table[,1]

  #Generate attribute combinations
  comb.generator<-function(x.vector){
    if(length(x.vector)>1){
      temp.attr<-x.vector
      temp.attr.sav<-NULL
      for(i in 1:length(temp.attr)){
        temp.1<-combn(temp.attr,i)
        temp.2<-apply(temp.1,2,function(x){paste(x,collapse = "")})
        temp.attr.sav<-c(temp.attr.sav,temp.2)
      }
    }
    if(length(x.vector)==1){temp.attr.sav<-x.vector}
    temp.attr.sav
  }
  #vectors needed for combination.generator
  Item.load.id<-list()
  for ( i in 1:nr){
    Item.load.id[[i]]<-grep('1',Qmatrix[i,])}

  Attr.load.id<-list()
  attr.load.id<-matrix(0,length(temp.table.col),nc)
  for ( i in 1:length(temp.table.col)){
    attr.load.id[i,]<-unlist(strsplit(temp.table.col[i],split=''))
    Attr.load.id[[i]]<-grep('1',attr.load.id[i,])
  }

  #Generate Combination for both Item.load and Attr.load
  Item.Comb<-list()
  for ( i in 1:nr){
    Item.Comb[[i]]<-comb.generator(Item.load.id[[i]])
  }
  Attr.Comb<-list()
  for ( i in 2:length(temp.table.col)){
    Attr.Comb[[1]]<-0
    Attr.Comb[[i]]<-comb.generator(Attr.load.id[[i]])
  }
  constraints.list<-list()
  nway.inter.list<-list()
  for(i in 1:nr){
    for(a in 2:length(temp.table.col)){
      ifzero<-as.numeric(paste(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],collapse=''))
      if((!is.na(ifzero))){
        temp.table[i,a]<-paste(c(temp.table[i,a],
                                 paste("S","l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='',collapse='')
        ),collapse='')
        if(a==length(temp.table.col)){
          nway.inter.list[[i]]<-nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])])
          constraints.list[[i]]<-paste("l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='')
        }
      }
    }
  }

  #Create Lambda Table
  Lamda.Table<-temp.table
  for(i in 1:nr){
    for(a in 1:length(Lamda.Table)){
      t.ref<-unique(as.character(Lamda.Table[i,]))
      pos<-c(1:length(t.ref))[Lamda.Table[i,a]==t.ref]
      temp.table[i,a]<-paste("t",i,"_",pos,sep='')}}

  #Generate CRUM specification
  out<-list()
  out[[1]]<-Lamda.Table
  out[[2]]<-temp.table
  out[[3]]<-constraints.list
  out[[4]]<-nway.inter.list
  out[[5]]<-intercept
  OUTPUT<-out
  nclass<-ncol(OUTPUT[[1]]);Nc<-nclass

  #Produce kernel expressions across items and attributes
  Kernel.exp<-OUTPUT[[1]]
  Kernel.exp.detect<-OUTPUT[[1]] #052719updates
  Kernel.exp.CRUM<-OUTPUT[[1]] #052719updates
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])
      Kernel.exp.detect[i,j]<-NA} #052719updates
    }
  }
  for (i in 1:nrow(OUTPUT[[1]])){ #052719updates
    theClosestEffect<-which(is.na(Kernel.exp.detect[i,]))[1] #052719updates
    useToReplaceLonger<-Kernel.exp[i,theClosestEffect] #052719updates
    Kernel.exp.CRUM[i,is.na(Kernel.exp.detect[i,])]<-useToReplaceLonger #052719updates
  } #052719updates

  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                        nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))

  if(length(name.inter)!=0){
    for (inter in 1: length(name.inter)){
      temp.nw<-numway.inter[inter]
      temp.nm<-name.inter[inter]
      temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
      temp.sel<-paste(unlist(strsplit(temp.nm,split = '_'))[1],"_",(1:(temp.nw-1)),sep='')
      first.sel<-unlist(OUTPUT[[3]])[grep(paste((temp.sel),collapse="|"),unlist(OUTPUT[[3]]))]
      second.sel<-sub(".*_.", "", first.sel)
      for (sel in 1:length(temp.subnm)){
        SEL<-second.sel[sel]
        Constrain.List1<-rbind(
          paste(temp.nm,">-(0", paste("+",first.sel[grep(SEL,second.sel)],
                                      sep='',collapse=''),")",sep=''),Constrain.List1)
      }
    }
    Constrain.List1<-as.character(Constrain.List1)
  }else{
    Constrain.List1<-NULL
  }

  itemParmName<-c(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==2],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==3],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==4],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==5],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==6],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==7],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==8],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==9],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==10],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==11],OUTPUT[[5]])
  numMainEffect<-length(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1])
  Constrain.List<-paste('  real<lower=0>',itemParmName[1:numMainEffect],';\n ')
  Unconstrain.List<-paste('  real',itemParmName[-(1:numMainEffect)],';\n ')
  Reparm<-as.data.frame(matrix(0,nr,nclass))

  #############################################################
  #######060319update: Multiple Group##########################
  #############################################################
  if(!is.na(fixeditem.vector)[1]){
    fixedItem.vector<-c(1:nr)[-fixeditem.vector]
  }else{fixedItem.vector<-c(1:nr)}

  #############################################################
  #######060719update: Multiple Group##########################
  #############################################################
  Kernel.exp.CRUM<-Kernel.exp
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      Reparm[loopi,loopc]<-paste('  PImat[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp[loopi,loopc]),');\n',sep='')
    }
  }
  for(loopi in 1:nr){
    for(loopc in 1:Nc){
      Reparm[loopi,loopc]<- str_replace_all(Reparm[loopi,loopc],paste((name.inter),collapse="|"),"0")
    }
  }
  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',c(itemParmName[1:numMainEffect],intercept),'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))
  #############################################################
  #######060719update: Multiple Group End######################
  #############################################################

  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))
  #############################################################


  Kernel.exp.CRUM.groupName<-paste("Kernel.exp.CRUM_g",c(1:group.num),sep='')
  for(i in 1:group.num){
    tempfill.Kernel.exp.CRUM<-Kernel.exp.CRUM
    temp.Kernel.exp.CRUM<-Kernel.exp.CRUM[fixedItem.vector,]
    for(j in 1:nrow(temp.Kernel.exp.CRUM)){
      for(z in 1:ncol(temp.Kernel.exp.CRUM)){
        temp.Kernel.exp.CRUM[j,z]<-paste(temp.Kernel.exp.CRUM[j,z],'_g',i,sep='')
      }
    }
    for(j in 1:nrow(temp.Kernel.exp.CRUM)){
      for(z in 1:ncol(temp.Kernel.exp.CRUM)){
        temp.Kernel.exp.CRUM[j,z]<-str_replace_all(temp.Kernel.exp.CRUM[j,z],"\\+",paste("_g",i,"+",sep=''))
      }
    }
    tempfill.Kernel.exp.CRUM[fixedItem.vector,]<-temp.Kernel.exp.CRUM
    assign(Kernel.exp.CRUM.groupName[i],tempfill.Kernel.exp.CRUM)
  }
  Kernel.exp.CRUM.list<-list()
  for(i in 1:group.num){Kernel.exp.CRUM.list[[i]]<-eval(parse(text=paste("Kernel.exp.CRUM_g",i,sep='')))}
  #############################################################
  ##########060719update: Multiple Group End###################
  #############################################################

  #############################################################
  ##########060719update: Multiple Group#######################
  #############################################################

  PImat.groupName<-paste("PImat_g",c(1:group.num),sep='')
  Reparm.multigroup<-array(0,dim = c(nr,nclass,group.num))
  #Produce Stan code for PImat parameter
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      for (loopg in 1:group.num){
        Reparm.multigroup[loopi,loopc,loopg]<-paste('  PImat_g',loopg,'[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp.CRUM.list[[loopg]][loopi,loopc]),');\n',sep='')

      }
    }
  }

  #############################################################
  ##########060719update: Multiple Group ######################
  #############################################################

  intercept.multigroup<-array(0,dim = c(nr,1,group.num))
  mainEff.multigroup<-array(0,dim = c(numMainEffect,1,group.num))
  #Group Invariant Parameter Name
  fixedParmName<-NULL
  if(!is.na(fixeditem.vector)[1]){
    for(i in fixeditem.vector){
      fixedParmName<-c(fixedParmName,out[[3]][[i]])
    }
    fixedParmName<-c(fixedParmName,out[[5]][fixeditem.vector])
  }
  #Group Variant Parameter Name
  freeParmName<-c(out[[5]],itemParmName)[!c(out[[5]],itemParmName)%in%fixedParmName]

  for(i in 1:group.num){
    tempfill.intercept<-out[[5]]
    tempfill.mainEff<-itemParmName[1:numMainEffect]


    temp.intercept<-tempfill.intercept[fixedItem.vector]
    temp.mainEff<-tempfill.mainEff[!tempfill.mainEff%in%fixedParmName]

    for(j in 1:length(temp.intercept)){
      temp.intercept[j]<-paste(temp.intercept[j],"_g",i,sep='')
    }

    for(j in 1:length(temp.mainEff)){
      temp.mainEff[j]<-paste(temp.mainEff[j],"_g",i,sep='')
    }


    tempfill.intercept[fixedItem.vector]<-temp.intercept
    tempfill.mainEff[!tempfill.mainEff%in%fixedParmName]<-temp.mainEff

    intercept.multigroup[,,i]<-tempfill.intercept
    mainEff.multigroup[,,i]<-tempfill.mainEff
  }

  Constrain.List<-paste('  real<lower=0>',unique(mainEff.multigroup),';\n ')
  Unconstrain.List<-paste('  real',unique(c(intercept.multigroup)),';\n ')
  #############################################################
  ##########060319update: Multiple Group End###################
  #############################################################

  #############################################################
  #####060719update:Change from update.Parmprior&fix     ######
  #############################################################
  update.Parmprior.multiGroup<-paste(paste('   ',unique(c(intercept.multigroup)),'~normal(0,5)',';\n',sep=''),
                                     paste('   ',unique(c(mainEff.multigroup)),'~normal(0,5)',';\n',sep=''))

  update.Parmprior.multiGroup<-unique(update.Parmprior.multiGroup)
  update.Parmprior.multiGroup<-c("   //Prior\n ",paste('   ',fixedParmName,'~normal(0,5)',';\n',sep=''),update.Parmprior.multiGroup )
  if(class.equal){
    update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='') )
  }else{
    for(i in 1:group.num){
      update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,
                                     paste('   Vc_g',i,'~dirichlet(rep_vector(2.0, Nc));\n',sep='') )
    }
  }
  ##therefore we can use: fix.Parmprior,update.Parmprior
  #############################################################
  #############################################################

  #############################################################
  #####060419update:Likelihood Add PImat_g#####################
  #############################################################
  PImat.likelihood1<-NULL
  PImat.likelihood0<-NULL
  Vc.likelihood<-NULL
  for(loopg in 1:group.num){
    temp.PImat.likelihood1<-paste(paste('          if (GroupID[iterp]==',loopg,')',sep=''),
                                  paste('            contributionsI[iteri]=bernoulli_lpmf(1|PImat_g',loopg,'[iteri,iterc]);',sep=''),
                                  sep='\n')
    PImat.likelihood1<-paste(PImat.likelihood1,temp.PImat.likelihood1,sep='\n')
    temp.PImat.likelihood0<-paste(paste('          if (GroupID[iterp]==',loopg,')',sep=''),
                                  paste('            contributionsI[iteri]=bernoulli_lpmf(0|PImat_g',loopg,'[iteri,iterc]);',sep=''),
                                  sep='\n')
    PImat.likelihood0<-paste(PImat.likelihood0,temp.PImat.likelihood0,sep='\n')
  }
  if(!class.equal){
    for(loopg in 1:group.num){
      temp.Vc.likelihood<-paste(paste('       if (GroupID[iterp]==',loopg,')',sep=''),
                                paste('         contributionsC[iterc]=log(Vc_g',loopg,'[iterc])+sum(contributionsI);',sep=''),
                                sep='\n')
      Vc.likelihood<-paste(Vc.likelihood,temp.Vc.likelihood,sep='\n')
    }
  }

  #Likelihood Stan code
  if(class.equal){
    Likelihood<-paste('
  \n
  //Likelihood
  for (iterp in 1:Np){
    for (iterc in 1:Nc){
      for (iteri in 1:Ni){
        if (Y[iterp,iteri] == 1)'
                      ,PImat.likelihood1,'\n',
                      '        else'
                      ,PImat.likelihood0,
                      '}
      contributionsC[iterc]=log(Vc[iterc])+sum(contributionsI);
    }
  target+=log_sum_exp(contributionsC);
  }
  ',sep='')}else{
    Likelihood<-paste('
  \n
  //Likelihood
  for (iterp in 1:Np){
    for (iterc in 1:Nc){
      for (iteri in 1:Ni){
        if (Y[iterp,iteri] == 1)'
                      ,PImat.likelihood1,'\n',
                      '        else'
                      ,PImat.likelihood0,
                      '}\n',
                      Vc.likelihood
                      ,'


  }
  target+=log_sum_exp(contributionsC);
 }
                      ',sep='')

  }
  #############################################################
  #####060419update:Likelihood Add PImat_g  end################
  #############################################################

  #Data Specification
  data.spec<-'
  data{
  int Np;
  int Ni;
  int Nc;
  matrix[Np, Ni] Y;
  vector[Np] GroupID;
  }
  '


  #############################################################
  #####060419update: Stan script##############################
  #############################################################

  #Parameter Specification
  if(class.equal){parm.spec<-paste(c('
  parameters{
  simplex[Nc] Vc;\n ',paste0(Constrain.List),paste0(Unconstrain.List),
                                     '}\n'),collapse='')}else{
                                       parm.spec<-paste(c('
  parameters{\n ',
                                                          paste(paste('   simplex[Nc] Vc_g',1:group.num, ";",sep=''),"\n"),
                                                          paste0(Constrain.List),paste0(Unconstrain.List),
                                                          '}\n'),collapse='')


                                     }

  #Reparameter Specification
  transparm.spec<-paste(c('
  transformed parameters{

                          ',
                          paste('  matrix[Ni, Nc] PImat_g',1:group.num,';\n',sep=''),
                          paste0(c(Reparm.multigroup)),'}\n'),collapse='')

  #Model Specification update052619
  model.spec<-paste(c('\nmodel {\n',paste(c(Modelcontainer,update.Parmprior.multiGroup,Likelihood),sep=''),'\n}',sep=''))
  model.spec<-model.spec[!startsWith(str_remove_all(model.spec," "),"~")]
  #Generated Quantities Specification
  IC.generatedquantities<-NULL
  if(!class.equal){
    for(loopg in 1:group.num){
      temp.IC.generatedquantities<-paste(paste('       if (GroupID[iterp]==',loopg,')',sep=''),
                                         paste('         contributionsIC[iteri,iterc]=log(Vc_g',loopg,'[iterc])+contributionsI[iteri];',sep=''),
                                         sep='\n')
      IC.generatedquantities<-paste(IC.generatedquantities,temp.IC.generatedquantities,sep='\n')
    }
  }


  if(class.equal){
    generatedQuantities.spec<-paste('
  \n
  generated quantities {
  vector[Ni] log_lik[Np];
  vector[Ni] contributionsI;
  matrix[Ni,Nc] contributionsIC;

  matrix[Ni,Nc] posteriorIC;
  matrix[Np,Nc] posteriorPC;


  //Posterior
  for (iterp in 1:Np){
    for (iteri in 1:Ni){
      for (iterc in 1:Nc){
        if (Y[iterp,iteri] == 1)'
                                    ,PImat.likelihood1,'\n',
                                    '        else'
                                    ,PImat.likelihood0,
                                    '
          contributionsIC[iteri,iterc]=log(Vc[iterc])+contributionsI[iteri];
          posteriorIC[iteri,iterc]=contributionsI[iteri];
        }
      log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
    }
    for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
  }
  }
  ',sep='')}else{
    generatedQuantities.spec <- paste( '
    \n
    generated quantities {
    vector[Ni] log_lik[Np];
    vector[Ni] contributionsI;
    matrix[Ni,Nc] contributionsIC;

    matrix[Ni,Nc] posteriorIC;
    matrix[Np,Nc] posteriorPC;


    //Posterior
    for (iterp in 1:Np){
      for (iteri in 1:Ni){
        for (iterc in 1:Nc){
          if (Y[iterp,iteri] == 1)',
                                       PImat.likelihood1,'\n',
                                       '        else',
                                       PImat.likelihood0,'\n',
                                       "   ",IC.generatedquantities,'\n',
                             '          posteriorIC[iteri,iterc]=contributionsI[iteri];','\n
      }\n',
                                       '     log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
     }
     for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
   }
   }
      ',
                                       sep = ''
    )
  }




  if (.Platform$OS.type == "unix") {
    filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
  }else{
    filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
  }

  sink(file=filename, append=FALSE)
  cat(
    paste(c('   ',
            data.spec,parm.spec,transparm.spec,model.spec,generatedQuantities.spec)
    ))
  sink(NULL)

}


#SEPERATION#
#' @title A function to generate leave-one-out cross-validation for Stan Model
#'
#' @description
#' The StanLCDM.loofit Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param savepath save the .stan file to somewhere; the default path is getwd()
#' @param savename name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export

StanDCM.fit <- function(x, pars ="log_lik", cores = 2, save_psis =TRUE) {
  loo1 <- loo(x = x, pars = pars, save_psis = save_psis, cores = cores)
  loo1
}

#SEPERATION#
#' @title A function to calculate pecentage of extreme p-values for sumscore distribution
#'
#' @description
#' The StanDCM.ppmc Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param stan.model A rStan object
#' @param response.matrix the response matrix used by RStan Object
#' @param n.sim number of simulations for Posterior Predictive Model Checking
#' @param n.burnin number of burn-ins
#' @param plot.option logical. whether to provide a plot for ppmc
#' @return p-values tables
#'
#' @author {Jihong Zhang, University of Iowa, \email{jihong-zhang@uiowa.edu}}
#'
#' @export
#' @examples
#' load("data.RData")
#' Qmatrix<-cbind(Qmatrix,rep(1,9)); Qmatrix[1,1]<-0
#' dim(respMatrix)
#' misspecifiedQmatrix <- Qmatrix
#' misspecifiedQmatrix[1:6,] <- 1-Qmatrix[1:6,]
#' misspecifiedQmatrix[1,3] = 0
#' mod2 <- StanDINA.run(misspecifiedQmatrix,response.matrix = respMatrix,iter=100,init.list='cdm', chain.num = 3, warmup = 20)
#' StanDCM.ppmc(stan.model = mod2, response.matrix = respMatrix, n.sim = 1000, n.burnin = 1, plot.option = FALSE)
#' end - start


StanDCM.ppmc <- function(stan.model, response.matrix, n.sim = 6000, n.burnin = 20, plot.option=FALSE) {
  Install.package(c("Rlab", "MCMCpack", "tidyr", "dplyr", "pbapply"))
  if(plot.option == TRUE) Install.package("ggplot2")

  mod1 <- stan.model
  posterior.matrix <- as.matrix(mod1)
  Ni <- ncol(response.matrix)
  Np <- nrow(response.matrix)
  PImat.name <- grep(names(mod1), pattern = "PImat", value = T)
  Nc <-  length(PImat.name) / Ni
  Vc.name <- grep(names(mod1), pattern = "Vc", value = T)
  n.sim = n.sim
  n.iter = mod1@stan_args[[1]]$iter

  if (n.burnin >= n.iter) {
    stop('Number of iteration less than number of burn-ins')
  }

  Vc.posterior.matrix <- posterior.matrix[(n.burnin +1): n.iter, Vc.name]
  PImat.posterior.matrix <- posterior.matrix[(n.burnin +1): n.iter, PImat.name]
  # simulate response matrix based on PImat and Vc
  pseudo.sumscore.dist.matrix = NULL
  if (is.na(n.sim)) {
    n.sim = Np
  }
  time = 0
  # create progress bar
  # pb <- txtProgressBar(min = 0, max = n.sim, style = 3)
  pseudo.sumscore.extract <- function() {
    time <<- time + 1
    iter <- sample(1:(n.iter-n.burnin), 1, replace = TRUE)
    # select Predicted Probability Matrix
    PImat.select <- matrix(PImat.posterior.matrix[iter,], Ni, Nc)

    # draw from dirichlet distribution and generate class for each person
    Vc.draw <- rdirichlet(1, rep(1, Nc))
    pseudo.class.vector<-apply(rmultinom(Np, 1, Vc.draw), 2, function(x){which(x==1)})
    pseudo.response.matrix <- t(PImat.select[, pseudo.class.vector])
    pseudo.response.matrix <- t(apply(pseudo.response.matrix, 1, function(x) rbinom(Ni, 1, x)))
    pseudo.sumscore.vector <- apply(pseudo.response.matrix, 1, sum)
    pseudo.sumscore.dist.vector <- data.frame(table(pseudo.sumscore.vector))
    pseudo.sumscore.dist.vector$time = time
    pseudo.sumscore.dist.vector
  }
  # (2) Set up the style for progree bar
  pboptions(type = "txt", style = 3, char = "=")
  #（1）extract the pseudo sumscore for n.sim times
  pseudo.sumscore.dist.matrix = pbreplicate(n.sim, pseudo.sumscore.extract(), simplify = FALSE)
  pseudo.sumscore.dist.matrix.long <- do.call(rbind,pseudo.sumscore.dist.matrix)
  pseudo.sumscore.dist.matrix.wide <- spread(key = pseudo.sumscore.vector, value = Freq , pseudo.sumscore.dist.matrix.long)
  pseudo.sumscore.dist.df <- pseudo.sumscore.dist.matrix.wide[,-1]
  #(2) Compute the observed raw score distribution
  rep.quantile <- apply(pseudo.sumscore.dist.df, 2, quantile,na.rm = TRUE, probs = c(0.05, 0.5, 0.95))
  rep.quantile <- data.frame(t(rep.quantile))
  colnames(rep.quantile) <- c("Q5", "Q50", "Q95")
  rep.quantile$score <- as.numeric(rownames(rep.quantile))
  response.sumscore.obs <- apply(response.matrix, 1, sum)
  observed <- data.frame(table(response.sumscore.obs))
  colnames(observed) <- c("score", "obs")
  observed$score <- as.numeric(as.character(observed$score))
  # (3) combine them together
  df <- left_join(observed, rep.quantile, by = "score")
  df$'WithinPosterior(95%)' <- ifelse(df$obs > df$Q95 | df$obs < df$Q5, FALSE, TRUE)


  if(plot.option == TRUE){
    p <- ggplot(df, aes(x = score + 1, y = obs)) +
          theme_bw() +
          labs(x = "Raw score", y = "Number of examinees", title = "The observed and replicated raw score distributions")

    # (2) Add the entire replicated raw score distributions using violin plot
    # and jittered data points
    p <- p + geom_violin(data = pseudo.sumscore.dist.matrix.long, aes(x = pseudo.sumscore.vector, y = Freq))
    p <- p + geom_jitter(data = pseudo.sumscore.dist.matrix.long, aes(x = pseudo.sumscore.vector, y = Freq, color = pseudo.sumscore.vector),
                         position = position_jitter(width = 0.12), alpha = 0.15) + theme(legend.position = "none")

    # (3) Overlay 5%, 50% (median), and 95% quantiles of each replicated raw
    # score distribution
    p <- p +
      geom_point(aes(y = Q50), size = 4, shape = 2) +
      geom_line(aes(y = Q50),size = 1, linetype = "dashed")
    p <- p + geom_line(aes(y = Q5), size = 1, linetype = "dotted")
    p <- p + geom_line(aes(y = Q95), size = 1, linetype = "dotted")

    # (4) Overlay the observed raw score distribution
    p <- p + geom_point(size = 4) + geom_line(size = 1)
    print(p)
  }

  df
}








#SEPERATION#
#' @title DCM calibration under the DINA model via Stan
#'
#'
#' @description
#' \code{StanDINA} uses Stan program to calibrate the deterministic inputs, noisy and gate model for dichotomous responses, and its
#' extension
#'
#' In addition, users are allowed to specify design matrix and link function for each item, and distinct models may be used
#' in a single test for different items. The attributions can be either dichomous or polytomous.

#' @usage
#' StanDINA.run<-function(Qmatrix,response.matrix,script.path=NA,save.path=getwd(),save.name="DINA_uninf",iter=1000,warmup = floor(iter/2),
#' chain.num=3,init.list='random',control.list=NA)
#'
#' @param Qmatrix A required matrix
#' @param response.matrix save the .stan file to somewhere; the default path is getwd()
#' @param script.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name the name of saved stan file
#' @param iter number of iteration of MCMC estimation. defalts to 1000.
#' @param chain.num number of MCMC chain.num.
#' @param init.list the initial values. 'random' or 'CDM'
#' @param control.list the controlled parameters
#'
#' @return StanDINA returens an object of class StanDINA. Methods for StanDINA objects include
#' \code{\link{extract}} for extract for extracting various components, \code{\link{coef}} for
#' extracting strctural parameters.
#'
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu} \cr Jihong Zhang, University of Iowa, \email{jihong-zhang@uiowa.edu}}
#'
#' @export
#' @examples
#' \dontrun{
#' #----------- DINA model-----------#
#' mod1<-StanDINA.run(Qmatrix, respMatrix, iter=20, init.list='cdm', chain.num = 3)
#' summary(mod1)
#' }


StanDINA.run<-function(Qmatrix,response.matrix,
                       script.path=NA,save.path=getwd(),save.name="DINA_uninf",
                       iter=1000,warmup = 0,
                       chain.num=3,init.list='random',control.list=NA){
  rstan.detect<-tryCatch(library("rstan"),error=function(e){"rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details."})
  if(length(rstan.detect)==1){
    stop()
  }
  Cdm.init<-F
  if(init.list=='cdm'){
    Cdm.init<-T
    Install.package(c("CDM","stringr"))
    trueParmName<-Parm.name(Qmatrix=Qmatrix)$parm.name
    Classp.exp1<-Parm.name(Qmatrix=Qmatrix)$class.expression
    mod1<-gdina( data =respMatrix, q.matrix = Qmatrix , maxit=700,link = "logit",progress=F)
    CDMresult<-as.data.frame(coef(mod1))
    library(stringr)
    CDM.parm.name<-paste(paste(paste('l',CDMresult[,3],sep=''),'_',sep=''),str_count(CDMresult$partype.attr,"Attr"),sep='')
    CDM.parm.name<-paste(CDM.parm.name,
                         unlist(lapply(strsplit(unlist(lapply(strsplit(CDMresult$partype.attr, 'Attr', fixed=FALSE),function(x){paste(x,collapse="")})),'-'),function(x){paste(x,collapse="")})),
                         sep='')
    CDM.parm.est<-CDMresult$est
    parm.ini<-round(CDM.parm.est[match(trueParmName,CDM.parm.name)],4)
    CDM.prop.est<-mod1$attribute.patt
    prop.ini<-CDM.prop.est[match(Classp.exp1,rownames(CDM.prop.est)),1]
    inilist1<-paste('list(',paste(noquote(paste(noquote(unlist(list(
      paste(paste('Vc=c(',paste((prop.ini),collapse=','),')',collapse=','))))
    ))),collapse=',')  ,')',collapse='')

    inilist1<-eval(parse(n =2000000 ,text=inilist1))
    for( i in 2:chain.num){
      temp.text<-paste('inilist',i,"<-inilist1",sep='')
      eval(parse(text=(temp.text)))
    }
    temp.text<-paste('init.list<-list(',paste(paste('inilist',1:chain.num,sep=''),collapse = ","),')',sep='')
    eval(parse(text=(temp.text)))
  }
  data.list<-Generate.datalist(Qmatrix,response.matrix)

  if(is.na(control.list)){control.list<-list(adapt_delta=0.82)}
  if(is.na(script.path)==T){
    options(warn=-1)
    StanDINA.script(Qmatrix,save.path=save.path,save.name=save.name)
    script.path<-paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
    options(warn=0)
    compiled_model<-stan_model(script.path)
  }else{
    compiled_model<-stan_model(script.path)
  }
  if(Cdm.init==T){
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})
  }else{
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})

  }

  estimated_model
}

#SEPERATION#
#' @title A function to generate leave-one-out cross-validation for Stan Model
#'
#' @description
#' The StanLCDM.loofit Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages

StanDINA.script<-function(Qmatrix,save.path=getwd(),save.name="DINA_uninf"){
  #Load packages
  Install.package("plyr")
  Install.package('stringr')

  nc<-ncol(Qmatrix)
  nr<-nrow(Qmatrix)
  temp.table.col<-unique(apply(combn(rep(c(0,1),nc),nc),2,function(x){paste(x,collapse = "")}))
  temp.table.col<-temp.table.col[order(temp.table.col)]
  temp.table<-matrix(0,nr,length(temp.table.col))
  colnames(temp.table)<-temp.table.col
  rownames(temp.table)<-paste('item',c(1:nr),sep='')
  temp.table<-as.data.frame(temp.table)
  for (i in 1:nr){
    temp.table[i,]<-paste('l',i,'_0',sep='')
  }
  intercept<-temp.table[,1]

  #Generate attribute combinations
  comb.generator<-function(x.vector){
    if(length(x.vector)>1){
      temp.attr<-x.vector
      temp.attr.sav<-NULL
      for(i in 1:length(temp.attr)){
        temp.1<-combn(temp.attr,i)
        temp.2<-apply(temp.1,2,function(x){paste(x,collapse = "")})
        temp.attr.sav<-c(temp.attr.sav,temp.2)
      }
    }
    if(length(x.vector)==1){temp.attr.sav<-x.vector}
    temp.attr.sav
  }
  #vectors needed for combination.generator
  Item.load.id<-list()
  for ( i in 1:nr){
    Item.load.id[[i]]<-grep('1',Qmatrix[i,])}

  Attr.load.id<-list()
  attr.load.id<-matrix(0,length(temp.table.col),nc)
  for ( i in 1:length(temp.table.col)){
    attr.load.id[i,]<-unlist(strsplit(temp.table.col[i],split=''))
    Attr.load.id[[i]]<-grep('1',attr.load.id[i,])
  }

  #Generate Combination for both Item.load and Attr.load
  Item.Comb<-list()
  for ( i in 1:nr){
    Item.Comb[[i]]<-comb.generator(Item.load.id[[i]])
  }
  Attr.Comb<-list()
  for ( i in 2:length(temp.table.col)){
    Attr.Comb[[1]]<-0
    Attr.Comb[[i]]<-comb.generator(Attr.load.id[[i]])
  }
  constraints.list<-list()
  nway.inter.list<-list()
  for(i in 1:nr){
    for(a in 2:length(temp.table.col)){
      ifzero<-as.numeric(paste(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],collapse=''))
      if((!is.na(ifzero))){
        temp.table[i,a]<-paste(c(temp.table[i,a],
                                 paste("S","l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='',collapse='')
        ),collapse='')
        if(a==length(temp.table.col)){
          nway.inter.list[[i]]<-nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])])
          constraints.list[[i]]<-paste("l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='')
        }
      }
    }
  }

  #Create Lambda Table
  Lamda.Table<-temp.table
  for(i in 1:nr){
    for(a in 1:length(Lamda.Table)){
      t.ref<-unique(as.character(Lamda.Table[i,]))
      pos<-c(1:length(t.ref))[Lamda.Table[i,a]==t.ref]
      temp.table[i,a]<-paste("t",i,"_",pos,sep='')}}

  #Generate LCDM specification
  out<-list()
  out[[1]]<-Lamda.Table
  out[[2]]<-temp.table
  out[[3]]<-constraints.list
  out[[4]]<-nway.inter.list
  out[[5]]<-intercept
  OUTPUT<-out
  nclass<-ncol(OUTPUT[[1]]);Nc<-nclass

  #Produce kernel expressions across items and attributes
  Kernel.exp<-OUTPUT[[1]]
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])}
    }
  }


  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                        nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))

  if(length(name.inter)!=0){
    for (inter in 1: length(name.inter)){
      temp.nw<-numway.inter[inter]
      temp.nm<-name.inter[inter]
      temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
      temp.sel<-paste(unlist(strsplit(temp.nm,split = '_'))[1],"_",(1:(temp.nw-1)),sep='')
      first.sel<-unlist(OUTPUT[[3]])[grep(paste((temp.sel),collapse="|"),unlist(OUTPUT[[3]]))]
      second.sel<-sub(".*_.", "", first.sel)
      for (sel in 1:length(temp.subnm)){
        SEL<-second.sel[sel]
        Constrain.List1<-rbind(
          paste(temp.nm,">-(0", paste("+",first.sel[grep(SEL,second.sel)],
                                      sep='',collapse=''),")",sep=''),Constrain.List1)
      }
    }
    Constrain.List1<-as.character(Constrain.List1)
  }else{
    Constrain.List1<-NULL
  }

  itemParmName<-c(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==2],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==3],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==4],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==5],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==6],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==7],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==8],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==9],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==10],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==11],OUTPUT[[5]])
  numMainEffect<-length(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1])
  Constrain.List<-paste('  real<lower=0>',itemParmName[1:numMainEffect],';\n ')
  Unconstrain.List<-paste('  real',itemParmName[-(1:numMainEffect)],';\n ')
  Reparm<-as.data.frame(matrix(0,nr,nclass))

  #############################################################
  ###########052619update:The highest-interactionn############
  hi.interaction<-rep(1,nr)
  zero.list<-constraints.list
  fixparm.vec<-NULL
  for(i in 1:nr){
    hi.interaction[i]<-constraints.list[[i]][length(constraints.list[[i]])]
    if(length(zero.list[[i]])==1){zero.list[[i]]=NA}else{
      zero.list[[i]]<-constraints.list[[i]][
        1:(length(constraints.list[[i]])-1)
        ]
      fixparm.vec<-c(fixparm.vec,zero.list[[i]])
    }
  } #intercept, hi.interaction,zero.list/fixparm.vec are what we need

  Constrain.List<-paste('  real<lower=0>',hi.interaction,';\n ')
  Unconstrain.List<-paste('  real',intercept,';\n ')
  #############################################################
  #############################################################

  #Produce Stan code for PImat parameter
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      Reparm[loopi,loopc]<-paste('  PImat[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp[loopi,loopc]),');\n',sep='')
    }
  }

  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))
  #############################################################
  ###########052619update:The highest-interactionn############
  update.Parmprior<-Parmprior
  fix.Parmprior<-NULL
  for(i in 1:length(Parmprior)){
    if(grepl(paste(fixparm.vec, collapse = "|"), Parmprior[i])){
      update.Parmprior[i]<-""
    }
  }
  fix.Parmprior<-c(paste('  real',fixparm.vec,';\n '),
                   paste(' ',fixparm.vec,"=0",';\n ')
  )

  #########052719update:create g and s parameters
  gParm<-rep(0,nr)
  sParm<-rep(0,nr)
  for(loopi in 1:nr){
    gParm[loopi]<-paste('  gParm[',loopi,']=inv_logit(',paste(Kernel.exp[loopi,1]),');\n',sep='')
    sParm[loopi]<-paste('  sParm[',loopi,']=1-inv_logit(',paste(Kernel.exp[loopi,nclass]),');\n',sep='')
  }
  ##therefore we can use: fix.Parmprior,update.Parmprior
  #############################################################
  #############################################################

  #Likelihood Stan code
  Likelihood<-'
  \n
  //Likelihood
  for (iterp in 1:Np){
    for (iterc in 1:Nc){
      for (iteri in 1:Ni){
        if (Y[iterp,iteri] == 1)
          contributionsI[iteri]=bernoulli_lpmf(1|PImat[iteri,iterc]);
        else
          contributionsI[iteri]=bernoulli_lpmf(0|PImat[iteri,iterc]);
      }
      contributionsC[iterc]=log(Vc[iterc])+sum(contributionsI);
    }
  target+=log_sum_exp(contributionsC);
  }
  '


  #Data Specification
  data.spec<-'
data{
  int Np;
  int Ni;
  int Nc;
  matrix[Np, Ni] Y;
}
  '
#Parameter Specification
parm.spec<-paste(c('
parameters{
  simplex[Nc] Vc;\n ',paste0(Constrain.List),paste0(Unconstrain.List),
                   '}\n'),collapse='')

#Reparameter Specification
transparm.spec<-paste(c('
 transformed parameters{
 matrix[Ni, Nc] PImat;
 vector[Ni] gParm;
 vector[Ni] sParm;\n',
                        fix.Parmprior,
                        gParm, #052719update
                        sParm, #052719update
                        paste0(unlist(Reparm)),'}\n'),collapse='')

#Model Specification update052619
model.spec<-paste(c('\nmodel {\n',paste(c(Modelcontainer,update.Parmprior,Likelihood),sep=''),'\n}',sep=''))
model.spec<-model.spec[!startsWith(str_remove_all(model.spec," "),"~")]
#Generated Quantities Specification
generatedQuantities.spec<-'
  \n
generated quantities {

 vector[Ni] log_lik[Np];
 vector[Ni] contributionsI;
 matrix[Ni,Nc] contributionsIC;
 
 matrix[Ni,Nc] posteriorIC;
 matrix[Np,Nc] posteriorPC;



 //Posterior
 for (iterp in 1:Np){
   for (iteri in 1:Ni){
     for (iterc in 1:Nc){
       if (Y[iterp,iteri] == 1)
          contributionsI[iteri]=bernoulli_lpmf(1|PImat[iteri,iterc]);
       else
           contributionsI[iteri]=bernoulli_lpmf(0|PImat[iteri,iterc]);
       contributionsIC[iteri,iterc]=log(Vc[iterc])+contributionsI[iteri];
       posteriorIC[iteri,iterc]=contributionsI[iteri];
      }
      log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
    }
   for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
  }
}
'

if (.Platform$OS.type == "unix") {
  filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
}else{
  filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
}

sink(file=filename,append=FALSE)
cat(
  paste(c('   ',
          data.spec,parm.spec,transparm.spec,model.spec,generatedQuantities.spec)
  ))
sink(NULL)

}


#SEPERATION#
#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanDINA_mG.run<-function(Qmatrix,
                          response.matrix,
                          GroupID,
                          fixeditem.vector=NA,
                          class.equal=T,
                          script.path=NA,save.path=getwd(),save.name="DINA_uninf_multiG",
                          iter=1000,warmup = 0,
                          chain.num=3,init.list='random',control.list=NA){
  group.num<-length(unique(GroupID))
  rstan.detect<-tryCatch(library("rstan"),error=function(e){"rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details."})
  if(length(rstan.detect)==1){
    stop()
  }
  Cdm.init<-F
  if(init.list=='cdm'){
    Cdm.init<-T
    Install.package(c("CDM","stringr"))
    trueParmName<-Parm.name(Qmatrix=Qmatrix)$parm.name
    Classp.exp1<-Parm.name(Qmatrix=Qmatrix)$class.expression
    mod1<-gdina( data =response.matrix, q.matrix = Qmatrix , maxit=700,link = "logit",progress=F)
    CDMresult<-as.data.frame(coef(mod1))
    library(stringr)
    CDM.parm.name<-paste(paste(paste('l',CDMresult[,3],sep=''),'_',sep=''),str_count(CDMresult$partype.attr,"Attr"),sep='')
    CDM.parm.name<-paste(CDM.parm.name,
                         unlist(lapply(strsplit(unlist(lapply(strsplit(CDMresult$partype.attr, 'Attr', fixed=FALSE),function(x){paste(x,collapse="")})),'-'),function(x){paste(x,collapse="")})),
                         sep='')
    CDM.parm.est<-CDMresult$est
    parm.ini<-round(CDM.parm.est[match(trueParmName,CDM.parm.name)],4)
    CDM.prop.est<-mod1$attribute.patt
    prop.ini<-CDM.prop.est[match(Classp.exp1,rownames(CDM.prop.est)),1]
    inilist1<-paste('list(',paste(noquote(paste(noquote(unlist(list(
      paste(paste('Vc=c(',paste((prop.ini),collapse=','),')',collapse=','))))
    ))),collapse=',')  ,')',collapse='')

    IniList1<-NULL
    if(!class.equal){
      temp.inilist1<-eval(parse(n =2000000 ,text=inilist1))
      eval(parse(n =2000000 ,text=paste(paste('IniList1$Vc_g',1:group.num,sep=''),"<-temp.inilist1$Vc",sep='') ))
      inilist1<-IniList1
    }else{
      inilist1<-eval(parse(n =2000000 ,text=inilist1))}


    for( i in 2:chain.num){
      temp.text<-paste('inilist',i,"<-inilist1",sep='')
      eval(parse(text=(temp.text)))
    }
    temp.text<-paste('init.list<-list(',paste(paste('inilist',1:chain.num,sep=''),collapse = ","),')',sep='')
    eval(parse(text=(temp.text)))
  }
  data.list<-Generate.datalist(Qmatrix,response.matrix,GroupID)

  if(is.na(control.list)){control.list<-list(adapt_delta=0.82)}
  if(is.na(script.path)==T){
    options(warn=-1)
    #Need to update script
    StanDINA_mG.script(Qmatrix=Qmatrix,
                       group.num=group.num,
                       fixeditem.vector=fixeditem.vector,
                       class.equal=class.equal,
                       save.path=save.path,save.name=save.name)
    script.path<-paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
    options(warn=0)
    compiled_model<-stan_model(script.path)}
  else{
    compiled_model<-stan_model(script.path)
  }
  if(Cdm.init==T){
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})
  }else{
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})

  }

  estimated_model
}


#SEPERATION#
#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanDINA_mG.script<-function(Qmatrix,
                             group.num,
                             fixeditem.vector=NA,
                             class.equal=T,
                             save.path=getwd(),save.name="DINA_uninf_multiG"){

  #Load packages
  Install.package("plyr")
  Install.package('stringr')

  nc<-ncol(Qmatrix)
  nr<-nrow(Qmatrix)
  temp.table.col<-unique(apply(combn(rep(c(0,1),nc),nc),2,function(x){paste(x,collapse = "")}))
  temp.table.col<-temp.table.col[order(temp.table.col)]
  temp.table<-matrix(0,nr,length(temp.table.col))
  colnames(temp.table)<-temp.table.col
  rownames(temp.table)<-paste('item',c(1:nr),sep='')
  temp.table<-as.data.frame(temp.table)
  for (i in 1:nr){
    temp.table[i,]<-paste('l',i,'_0',sep='')
  }
  intercept<-temp.table[,1]

  #Generate attribute combinations
  comb.generator<-function(x.vector){
    if(length(x.vector)>1){
      temp.attr<-x.vector
      temp.attr.sav<-NULL
      for(i in 1:length(temp.attr)){
        temp.1<-combn(temp.attr,i)
        temp.2<-apply(temp.1,2,function(x){paste(x,collapse = "")})
        temp.attr.sav<-c(temp.attr.sav,temp.2)
      }
    }
    if(length(x.vector)==1){temp.attr.sav<-x.vector}
    temp.attr.sav
  }
  #vectors needed for combination.generator
  Item.load.id<-list()
  for ( i in 1:nr){
    Item.load.id[[i]]<-grep('1',Qmatrix[i,])}

  Attr.load.id<-list()
  attr.load.id<-matrix(0,length(temp.table.col),nc)
  for ( i in 1:length(temp.table.col)){
    attr.load.id[i,]<-unlist(strsplit(temp.table.col[i],split=''))
    Attr.load.id[[i]]<-grep('1',attr.load.id[i,])
  }

  #Generate Combination for both Item.load and Attr.load
  Item.Comb<-list()
  for ( i in 1:nr){
    Item.Comb[[i]]<-comb.generator(Item.load.id[[i]])
  }
  Attr.Comb<-list()
  for ( i in 2:length(temp.table.col)){
    Attr.Comb[[1]]<-0
    Attr.Comb[[i]]<-comb.generator(Attr.load.id[[i]])
  }
  constraints.list<-list()
  nway.inter.list<-list()
  for(i in 1:nr){
    for(a in 2:length(temp.table.col)){
      ifzero<-as.numeric(paste(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],collapse=''))
      if((!is.na(ifzero))){
        temp.table[i,a]<-paste(c(temp.table[i,a],
                                 paste("S","l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='',collapse='')
        ),collapse='')
        if(a==length(temp.table.col)){
          nway.inter.list[[i]]<-nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])])
          constraints.list[[i]]<-paste("l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='')
        }
      }
    }
  }

  #Create Lambda Table
  Lamda.Table<-temp.table
  for(i in 1:nr){
    for(a in 1:length(Lamda.Table)){
      t.ref<-unique(as.character(Lamda.Table[i,]))
      pos<-c(1:length(t.ref))[Lamda.Table[i,a]==t.ref]
      temp.table[i,a]<-paste("t",i,"_",pos,sep='')}}

  #Generate LCDM specification
  out<-list()
  out[[1]]<-Lamda.Table
  out[[2]]<-temp.table
  out[[3]]<-constraints.list
  out[[4]]<-nway.inter.list
  out[[5]]<-intercept
  OUTPUT<-out
  nclass<-ncol(OUTPUT[[1]]);Nc<-nclass

  #Produce kernel expressions across items and attributes
  Kernel.exp<-OUTPUT[[1]]
  Kernel.exp.detect<-OUTPUT[[1]] #052719updates
  Kernel.exp.dina<-OUTPUT[[1]] #052719updates
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])
      Kernel.exp.detect[i,j]<-NA} #052719updates
    }
  }
  for (i in 1:nrow(OUTPUT[[1]])){ #052719updates
    theClosestEffect<-which(is.na(Kernel.exp.detect[i,]))[1] #052719updates
    useToReplaceLonger<-Kernel.exp[i,theClosestEffect] #052719updates
    Kernel.exp.dina[i,is.na(Kernel.exp.detect[i,])]<-useToReplaceLonger #052719updates
  } #052719updates

  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                        nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))

  if(length(name.inter)!=0){
    for (inter in 1: length(name.inter)){
      temp.nw<-numway.inter[inter]
      temp.nm<-name.inter[inter]
      temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
      temp.sel<-paste(unlist(strsplit(temp.nm,split = '_'))[1],"_",(1:(temp.nw-1)),sep='')
      first.sel<-unlist(OUTPUT[[3]])[grep(paste((temp.sel),collapse="|"),unlist(OUTPUT[[3]]))]
      second.sel<-sub(".*_.", "", first.sel)
      for (sel in 1:length(temp.subnm)){
        SEL<-second.sel[sel]
        Constrain.List1<-rbind(
          paste(temp.nm,">-(0", paste("+",first.sel[grep(SEL,second.sel)],
                                      sep='',collapse=''),")",sep=''),Constrain.List1)
      }
    }
    Constrain.List1<-as.character(Constrain.List1)
  }else{
    Constrain.List1<-NULL
  }

  itemParmName<-c(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==2],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==3],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==4],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==5],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==6],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==7],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==8],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==9],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==10],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==11],OUTPUT[[5]])
  numMainEffect<-length(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1])
  Constrain.List<-paste('  real<lower=0>',itemParmName[1:numMainEffect],';\n ')
  Unconstrain.List<-paste('  real',itemParmName[-(1:numMainEffect)],';\n ')
  Reparm<-as.data.frame(matrix(0,nr,nclass))


  #############################################################
  ###########052619update:The highest-interactionn############
  hi.interaction<-rep(1,nr)
  zero.list<-constraints.list
  fixparm.vec<-NULL
  for(i in 1:nr){
    hi.interaction[i]<-constraints.list[[i]][length(constraints.list[[i]])]
    if(length(zero.list[[i]])==1){zero.list[[i]]=NA}else{
      zero.list[[i]]<-constraints.list[[i]][
        1:(length(constraints.list[[i]])-1)
        ]
      fixparm.vec<-c(fixparm.vec,zero.list[[i]])
    }
  } #intercept, hi.interaction,zero.list/fixparm.vec are what we need

  Constrain.List<-paste('  real<lower=0>',hi.interaction,';\n ')
  Unconstrain.List<-paste('  real',intercept,';\n ')

  #############################################################
  #######060319update: Multiple Group##########################
  #############################################################
  if(!is.na(fixeditem.vector)[1]){
    fixedItem.vector<-c(1:nr)[-fixeditem.vector]
  }else{fixedItem.vector<-c(1:nr)}

  #############################################################
  #######060719update: Multiple Group##########################
  #############################################################
  Kernel.exp.dina<-Kernel.exp
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      Reparm[loopi,loopc]<-paste('  PImat[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp[loopi,loopc]),');\n',sep='')
    }
  }
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      Kernel.exp.dina[loopi,loopc]<-str_replace_all(Kernel.exp.dina[loopi,loopc],paste(fixparm.vec, collapse = "|"),"0")
      Kernel.exp.dina[loopi,loopc]<-str_replace_all(Kernel.exp.dina[loopi,loopc],"\\+0","")
      Reparm[loopi,loopc]<-str_replace_all(Reparm[loopi,loopc],paste(fixparm.vec, collapse = "|"),"0")
      Reparm[loopi,loopc]<-str_replace_all(Reparm[loopi,loopc],"\\+0","")
    }
  }
  #############################################################
  #######060719update: Multiple Group End######################
  #############################################################

  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))  #############################################################


  Kernel.exp.dina.groupName<-paste("Kernel.exp.dina_g",c(1:group.num),sep='')
  for(i in 1:group.num){
    tempfill.Kernel.exp.dina<-Kernel.exp.dina
    temp.Kernel.exp.dina<-Kernel.exp.dina[fixedItem.vector,]
    for(j in 1:nrow(temp.Kernel.exp.dina)){
      for(z in 1:ncol(temp.Kernel.exp.dina)){
        temp.Kernel.exp.dina[j,z]<-paste(temp.Kernel.exp.dina[j,z],'_g',i,sep='')
      }
    }
    for(j in 1:nrow(temp.Kernel.exp.dina)){
      for(z in 1:ncol(temp.Kernel.exp.dina)){
        temp.Kernel.exp.dina[j,z]<-str_replace_all(temp.Kernel.exp.dina[j,z],"\\+",paste("_g",i,"+",sep=''))
      }
    }
    tempfill.Kernel.exp.dina[fixedItem.vector,]<-temp.Kernel.exp.dina
    assign(Kernel.exp.dina.groupName[i],tempfill.Kernel.exp.dina)
  }
  Kernel.exp.dina.list<-list()
  for(i in 1:group.num){Kernel.exp.dina.list[[i]]<-eval(parse(text=paste("Kernel.exp.dina_g",i,sep='')))}
  #############################################################
  ##########060719update: Multiple Group End###################
  #############################################################

  #############################################################
  ##########060719update: Multiple Group#######################
  #############################################################

  PImat.groupName<-paste("PImat_g",c(1:group.num),sep='')
  Reparm.multigroup<-array(0,dim = c(nr,nclass,group.num))
  #Produce Stan code for PImat parameter
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      for (loopg in 1:group.num){
        Reparm.multigroup[loopi,loopc,loopg]<-paste('  PImat_g',loopg,'[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp.dina.list[[loopg]][loopi,loopc]),');\n',sep='')

      }
    }
  }
  #############################################################
  ##########060319update: Multiple Group End###################
  #############################################################
  #########052719update:create g and s parameters
  #########052719update:create g and s parameters
  gParm<-rep(0,nr)
  sParm<-rep(0,nr)
  for(loopi in 1:nr){
    gParm[loopi]<-paste('  gParm[',loopi,']=inv_logit(',paste(Kernel.exp.dina[loopi,1]),');\n',sep='')
    sParm[loopi]<-paste('  sParm[',loopi,']=1-inv_logit(',paste(Kernel.exp.dina[loopi,nclass]),');\n',sep='')
  }
  #############################################################
  ##########060319update: Multiple Group ######################
  #############################################################
  gParm.multigroup<-array(0,dim = c(nr,1,group.num))
  sParm.multigroup<-array(0,dim = c(nr,1,group.num))
  for(i in 1:group.num){
    tempfill.gParm<-gParm
    tempfill.sParm<-sParm

    tempfill.gParm<-str_replace_all(tempfill.gParm,"gParm",paste("gParm_g",i,sep=''))
    tempfill.sParm<-str_replace_all(tempfill.sParm,"sParm",paste("sParm_g",i,sep=''))

    temp.gParm<-gParm[fixedItem.vector]
    temp.sParm<-sParm[fixedItem.vector]

    for(j in 1:length(temp.gParm)){
      temp.gParm[j]<-str_replace_all(temp.gParm[j],"\\)",paste("_g",i,")",sep=''))
      temp.sParm[j]<-str_replace_all(temp.sParm[j],"\\)",paste("_g",i,")",sep=''))
    }
    for(j in 1:length(temp.gParm)){
      temp.gParm[j]<-str_replace_all(temp.gParm[j],"\\+",paste("_g",i,")",sep=''))
      temp.sParm[j]<-str_replace_all(temp.sParm[j],"\\+",paste("_g",i,"+",sep=''))
    }

    temp.gParm<-str_replace_all(temp.gParm,"gParm",paste("gParm_g",i,sep=''))
    temp.sParm<-str_replace_all(temp.sParm,"sParm",paste("sParm_g",i,sep=''))

    tempfill.gParm[fixedItem.vector]<-temp.gParm
    tempfill.sParm[fixedItem.vector]<-temp.sParm

    gParm.multigroup[,,i]<-tempfill.gParm
    sParm.multigroup[,,i]<-tempfill.sParm

  }

  #############################################################
  ##########060319update: Multiple Group ######################
  #############################################################
  keep.oneMainEffect.multigroup<-NULL
  intercept.multigroup<-array(0,dim = c(nr,1,group.num))
  mainEff.multigroup<-array(0,dim = c(numMainEffect,1,group.num))
  interaction.multigroup<-array(0,dim = c(length(name.inter),1,group.num))
  #Group Invariant Parameter Name
  fixedParmName<-NULL
  if(!is.na(fixeditem.vector)[1]){
    for(i in fixeditem.vector){
      fixedParmName<-c(fixedParmName,out[[3]][[i]])
    }
    fixedParmName<-c(fixedParmName,out[[5]][fixeditem.vector])
  }
  #Group Variant Parameter Name
  freeParmName<-c(out[[5]],itemParmName)[!c(out[[5]],itemParmName)%in%fixedParmName]

  for(i in 1:group.num){
    tempfill.intercept<-out[[5]]
    tempfill.mainEff<-itemParmName[1:numMainEffect]
    tempfill.interaction<-name.inter
    tempfill.keep.oneMainEffect<-hi.interaction

    temp.intercept<-tempfill.intercept[fixedItem.vector]
    temp.mainEff<-tempfill.mainEff[!tempfill.mainEff%in%fixedParmName]
    temp.interaction<-tempfill.interaction[!tempfill.interaction%in%fixedParmName]
    temp.keep.oneMainEffect<-tempfill.keep.oneMainEffect[!tempfill.keep.oneMainEffect%in%fixedParmName]

    for(j in 1:length(temp.intercept)){
      temp.intercept[j]<-paste(temp.intercept[j],"_g",i,sep='')
    }

    for(j in 1:length(temp.mainEff)){
      temp.mainEff[j]<-paste(temp.mainEff[j],"_g",i,sep='')
    }

    for(j in 1:length(temp.interaction)){
      temp.interaction[j]<-paste(temp.interaction[j],"_g",i,sep='')
    }

    for(j in 1:length(temp.keep.oneMainEffect)){
      temp.keep.oneMainEffect[j]<-paste(temp.keep.oneMainEffect[j],"_g",i,sep='')
    }

    tempfill.intercept[fixedItem.vector]<-temp.intercept
    tempfill.mainEff[!tempfill.mainEff%in%fixedParmName]<-temp.mainEff
    tempfill.interaction[!tempfill.interaction%in%fixedParmName]<-temp.interaction
    tempfill.keep.oneMainEffect[!tempfill.keep.oneMainEffect%in%fixedParmName]<-temp.keep.oneMainEffect

    intercept.multigroup[,,i]<-tempfill.intercept
    mainEff.multigroup[,,i]<-tempfill.mainEff
    interaction.multigroup[,,i]<-tempfill.interaction
    keep.oneMainEffect.multigroup<-c(keep.oneMainEffect.multigroup,tempfill.keep.oneMainEffect)
  }

  Constrain.List<-paste('  real<lower=0>',unique(keep.oneMainEffect.multigroup),';\n ')
  Unconstrain.List<-paste('  real',unique(c(intercept.multigroup)),';\n ')
  #############################################################
  ##########060319update: Multiple Group End###################
  #############################################################

 
  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))
  update.Parmprior<-Parmprior
  fix.Parmprior<-NULL
  for(i in 1:length(Parmprior)){
    if(grepl(paste(fixparm.vec, collapse = "|"), Parmprior[i])){
      update.Parmprior[i]<-""
    }
  }
  update.Parmprior<-update.Parmprior[update.Parmprior!='']


  fix.Parmprior<-c(paste('  real',fixparm.vec,';\n '),
                   paste(' ',fixparm.vec,"=0",';\n ')
  )
  #############################################################
  #####060319update:Change from update.Parmprior&fix     ######
  #############################################################
  update.Parmprior.multiGroup<-NULL
  for(i in 1:length(update.Parmprior)){
    for (j in 1:length(freeParmName)){
      for (z in 1:group.num){
        if(sum(grepl(freeParmName[j],update.Parmprior[i]))>=1){
          temp.update.Parmprior<-str_replace_all(update.Parmprior[i],freeParmName[j],paste(freeParmName[j],"_g",z,sep=''))
          update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,
                                         temp.update.Parmprior)
        }
      }
    }
  }
  update.Parmprior.multiGroup<-unique(update.Parmprior.multiGroup)
  update.Parmprior.multiGroup<-c("   //Prior\n",paste('   ',fixedParmName,'~normal(0,5)',';\n',sep=''),update.Parmprior.multiGroup )
  if(class.equal){
    update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='') )
  }else{
    for(i in 1:group.num){
      update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,
                                     paste('   Vc_g',i,'~dirichlet(rep_vector(2.0, Nc));\n',sep='') )
    }
  }

  ##therefore we can use: fix.Parmprior,update.Parmprior
  #############################################################
  #############################################################

  #############################################################
  #####060419update:Likelihood Add PImat_g#####################
  #############################################################
  PImat.likelihood1<-NULL
  PImat.likelihood0<-NULL
  Vc.likelihood<-NULL
  for(loopg in 1:group.num){
    temp.PImat.likelihood1<-paste(paste('          if (GroupID[iterp]==',loopg,')',sep=''),
                                  paste('            contributionsI[iteri]=bernoulli_lpmf(1|PImat_g',loopg,'[iteri,iterc]);',sep=''),
                                  sep='\n')
    PImat.likelihood1<-paste(PImat.likelihood1,temp.PImat.likelihood1,sep='\n')
    temp.PImat.likelihood0<-paste(paste('          if (GroupID[iterp]==',loopg,')',sep=''),
                                  paste('            contributionsI[iteri]=bernoulli_lpmf(0|PImat_g',loopg,'[iteri,iterc]);',sep=''),
                                  sep='\n')
    PImat.likelihood0<-paste(PImat.likelihood0,temp.PImat.likelihood0,sep='\n')
  }
  if(!class.equal){
    for(loopg in 1:group.num){
      temp.Vc.likelihood<-paste(paste('       if (GroupID[iterp]==',loopg,')',sep=''),
                                paste('         contributionsC[iterc]=log(Vc_g',loopg,'[iterc])+sum(contributionsI);',sep=''),
                                sep='\n')
      Vc.likelihood<-paste(Vc.likelihood,temp.Vc.likelihood,sep='\n')
    }
  }

  #Likelihood Stan code
  if(class.equal){
    Likelihood<-paste('
  \n
  //Likelihood
  for (iterp in 1:Np){
    for (iterc in 1:Nc){
      for (iteri in 1:Ni){
        if (Y[iterp,iteri] == 1)'
                      ,PImat.likelihood1,'\n',
                      '        else'
                      ,PImat.likelihood0,
                      '}
      contributionsC[iterc]=log(Vc[iterc])+sum(contributionsI);
    }
  target+=log_sum_exp(contributionsC);
  }
  ',sep='')}else{
    Likelihood<-paste('
  \n
  //Likelihood
  for (iterp in 1:Np){
    for (iterc in 1:Nc){
      for (iteri in 1:Ni){
        if (Y[iterp,iteri] == 1)'
                      ,PImat.likelihood1,'\n',
                      '        else'
                      ,PImat.likelihood0,
                      '}\n',
                      Vc.likelihood
                      ,'


  }
  target+=log_sum_exp(contributionsC);
 }
                      ',sep='')

  }
  #############################################################
  #####060419update:Likelihood Add PImat_g  end################
  #############################################################

  #Data Specification
  data.spec<-'
  data{
  int Np;
  int Ni;
  int Nc;
  matrix[Np, Ni] Y;
  vector[Np] GroupID;
  }
  '


  #############################################################
  #####060419update: Stan script##############################
  #############################################################

  #Parameter Specification
  if(class.equal){parm.spec<-paste(c('
  parameters{
  simplex[Nc] Vc;\n ',paste0(Constrain.List),paste0(Unconstrain.List),
                                     '}\n'),collapse='')}else{
                                       parm.spec<-paste(c('
  parameters{\n ',
                                                          paste(paste('   simplex[Nc] Vc_g',1:group.num, ";",sep=''),"\n"),
                                                          paste0(Constrain.List),paste0(Unconstrain.List),
                                                          '}\n'),collapse='')


                                     }

  #Reparameter Specification
  transparm.spec<-paste(c('
  transformed parameters{

                          ',
                          paste('  matrix[Ni, Nc] PImat_g',1:group.num,';\n',sep=''),
                          paste('  vector[Ni] gParm_g',1:group.num,';\n',sep=''),
                          paste('  vector[Ni] sParm_g',1:group.num,';\n',sep=''),

                          c(gParm.multigroup), #060419update
                          c(sParm.multigroup), #060419update
                          paste0(c(Reparm.multigroup)),'}\n'),collapse='')

  #Model Specification update052619
  model.spec<-paste(c('\nmodel {\n',paste(c(Modelcontainer,update.Parmprior.multiGroup,Likelihood),sep=''),'\n}',sep=''))
  model.spec<-model.spec[!startsWith(str_remove_all(model.spec," "),"~")]
  #Generated Quantities Specification
  IC.generatedquantities<-NULL
  if(!class.equal){
    for(loopg in 1:group.num){
      temp.IC.generatedquantities<-paste(paste('       if (GroupID[iterp]==',loopg,')',sep=''),
                                         paste('         contributionsIC[iteri,iterc]=log(Vc_g',loopg,'[iterc])+contributionsI[iteri];',sep=''),
                                         sep='\n')
      IC.generatedquantities<-paste(IC.generatedquantities,temp.IC.generatedquantities,sep='\n')
    }
  }


   if(class.equal){
    generatedQuantities.spec<-paste('
  \n
  generated quantities {
  vector[Ni] log_lik[Np];
  vector[Ni] contributionsI;
  matrix[Ni,Nc] contributionsIC;

  matrix[Ni,Nc] posteriorIC;
  matrix[Np,Nc] posteriorPC;


  //Posterior
  for (iterp in 1:Np){
    for (iteri in 1:Ni){
      for (iterc in 1:Nc){
        if (Y[iterp,iteri] == 1)'
                                    ,PImat.likelihood1,'\n',
                                    '        else'
                                    ,PImat.likelihood0,
                                    '
          contributionsIC[iteri,iterc]=log(Vc[iterc])+contributionsI[iteri];
          posteriorIC[iteri,iterc]=contributionsI[iteri];
        }
      log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
    }
    for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
  }
  }
  ',sep='')}else{
    generatedQuantities.spec <- paste( '
    \n
    generated quantities {
    vector[Ni] log_lik[Np];
    vector[Ni] contributionsI;
    matrix[Ni,Nc] contributionsIC;

    matrix[Ni,Nc] posteriorIC;
    matrix[Np,Nc] posteriorPC;


    //Posterior
    for (iterp in 1:Np){
      for (iteri in 1:Ni){
        for (iterc in 1:Nc){
          if (Y[iterp,iteri] == 1)',
                                       PImat.likelihood1,'\n',
                                       '        else',
                                       PImat.likelihood0,'\n',
                                       "   ",IC.generatedquantities,'\n',
                             '          posteriorIC[iteri,iterc]=contributionsI[iteri];','\n
      }\n',
                                       '     log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
     }
     for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
   }
   }
      ',
                                       sep = ''
    )
  }



  if (.Platform$OS.type == "unix") {
    filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
  }else{
    filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
  }

  sink(file=filename, append=FALSE)
  cat(
    paste(c('   ',
            data.spec,parm.spec,transparm.spec,model.spec,generatedQuantities.spec)
    ))
  sink(NULL)

}


#SEPERATION#
#' @title DCM calibration under the DINA model via Stan
#'
#'
#' @description
#' \code{StanDINA} uses Stan program to calibrate the deterministic inputs, noisy and gate model for dichotomous responses, and its
#' extension
#'
#' In addition, users are allowed to specify design matrix and link function for each item, and distinct models may be used
#' in a single test for different items. The attributions can be either dichomous or polytomous.

#' @usage
#' StanDINA.run<-function(Qmatrix,response.matrix,script.path=NA,save.path=getwd(),save.name="DINA_uninf",iter=1000,warmup = floor(iter/2),
#' chain.num=3,init.list='random',control.list=NA)
#'
#' @param Qmatrix A required matrix
#' @param response.matrix save the .stan file to somewhere; the default path is getwd()
#' @param script.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name the name of saved stan file
#' @param iter number of iteration of MCMC estimation. defalts to 1000.
#' @param chain.num number of MCMC chain.num.
#' @param init.list the initial values. 'random' or 'CDM'
#' @param control.list the controlled parameters
#'
#' @return StanDINA returens an object of class StanDINA. Methods for StanDINA objects include
#' \code{\link{extract}} for extract for extracting various components, \code{\link{coef}} for
#' extracting strctural parameters.
#'
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu} \cr Jihong Zhang, University of Iowa, \email{jihong-zhang@uiowa.edu}}
#'
#' @export
#' @examples
#' \dontrun{
#' #----------- DINO model-----------#
#' mod1<-StanDINO.run(Qmatrix, respMatrix, iter=20, init.list='cdm', chain.num = 3)
#' summary(mod1)
#' }

StanDINO.run<-function(Qmatrix,response.matrix,script.path=NA,save.path=getwd(),save.name="DINO_uninf",iter=1000,warmup = 0,
                       chain.num=3,init.list='random',control.list=NA){
  rstan.detect<-tryCatch(library("rstan"),error=function(e){"rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details."})
  if(length(rstan.detect)==1){
    stop()
  }
  Cdm.init<-F
  if(init.list=='cdm'){
    Cdm.init<-T
    Install.package(c("CDM","stringr"))
    trueParmName<-Parm.name(Qmatrix=Qmatrix)$parm.name
    Classp.exp1<-Parm.name(Qmatrix=Qmatrix)$class.expression
    mod1<-gdina( data =respMatrix, q.matrix = Qmatrix , maxit=700,link = "logit",progress=F)
    CDMresult<-as.data.frame(coef(mod1))
    library(stringr)
    CDM.parm.name<-paste(paste(paste('l',CDMresult[,3],sep=''),'_',sep=''),str_count(CDMresult$partype.attr,"Attr"),sep='')
    CDM.parm.name<-paste(CDM.parm.name,
                         unlist(lapply(strsplit(unlist(lapply(strsplit(CDMresult$partype.attr, 'Attr', fixed=FALSE),function(x){paste(x,collapse="")})),'-'),function(x){paste(x,collapse="")})),
                         sep='')
    CDM.parm.est<-CDMresult$est
    parm.ini<-round(CDM.parm.est[match(trueParmName,CDM.parm.name)],4)
    CDM.prop.est<-mod1$attribute.patt
    prop.ini<-CDM.prop.est[match(Classp.exp1,rownames(CDM.prop.est)),1]
    inilist1<-paste('list(',paste(noquote(paste(noquote(unlist(list(
      paste(paste('Vc=c(',paste((prop.ini),collapse=','),')',collapse=','))))
    ))),collapse=',') ,')',collapse='')

    inilist1<-eval(parse(n =2000000 ,text=inilist1))
    for( i in 2:chain.num){
      temp.text<-paste('inilist',i,"<-inilist1",sep='')
      eval(parse(text=(temp.text)))
    }
    temp.text<-paste('init.list<-list(',paste(paste('inilist',1:chain.num,sep=''),collapse = ","),')',sep='')
    eval(parse(text=(temp.text)))
  }
  data.list<-Generate.datalist(Qmatrix,response.matrix)

  if(is.na(control.list)){control.list<-list(adapt_delta=0.82)}
  if(is.na(script.path)==T){
    options(warn=-1)
    StanDINO.script(Qmatrix,save.path=save.path,save.name=save.name)
    script.path<-paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
    options(warn=0)
    compiled_model<-stan_model(script.path)
  }else{
    compiled_model<-stan_model(script.path)
  }
  if(Cdm.init==T){
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})
  }else{
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})

  }

  estimated_model
}

#SEPERATION#
#' @title Generate Stan code and Run the estimation for DINO model
#'
#'
#' @description
#' The \code{StanDINO.run} Function allows to automate Stan code geneartion for DINO model with binary resposnes
#'
#' @usage
#' StanDINA.run<-function(Qmatrix,response.matrix,script.path=NA,save.path=getwd(),
#' save.name="DINA_uninf",iter=1000,warmup = floor(iter/2),
#' chain.num=3,init.list='random',control.list=NA)
#'
#' @param Qmatrix A required matrix
#' @param response.matrix save the .stan file to somewhere; the default path is getwd()
#' @param script.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @param iter name the .stan
#' @param chain.num name the .stan
#' @param init.list name the .stan
#' @param control.list name the .stan
#'
#' @return StanDINA returens an object of class StanDINA. Methods for StanDINA objects include
#' \code{\link{extract}} for extract for extracting various components, \code{\link{coef}} for
#' extracting strctural parameters.
#'
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu} \cr Jihong Zhang, University of Iowa, \email{jihong-zhang@uiowa.edu}}
#'
#' @export
#' @examples
#' \dontrun{
#' #----------- DINO model-----------#
#' mod1<-StanDINO.run(Qmatrix, respMatrix, iter=20, init.list='cdm', chain.num = 3)
#' summary(mod1)
#' }

 
StanDINO.script<-function(Qmatrix,save.path=getwd(),save.name="DINO_uninf"){

  #Load packages
  Install.package("plyr")
  Install.package('stringr')

  nc<-ncol(Qmatrix)
  nr<-nrow(Qmatrix)
  temp.table.col<-unique(apply(combn(rep(c(0,1),nc),nc),2,function(x){paste(x,collapse = "")}))
  temp.table.col<-temp.table.col[order(temp.table.col)]
  temp.table<-matrix(0,nr,length(temp.table.col))
  colnames(temp.table)<-temp.table.col
  rownames(temp.table)<-paste('item',c(1:nr),sep='')
  temp.table<-as.data.frame(temp.table)
  for (i in 1:nr){
    temp.table[i,]<-paste('l',i,'_0',sep='')
  }
  intercept<-temp.table[,1]

  #Generate attribute combinations
  comb.generator<-function(x.vector){
    if(length(x.vector)>1){
      temp.attr<-x.vector
      temp.attr.sav<-NULL
      for(i in 1:length(temp.attr)){
        temp.1<-combn(temp.attr,i)
        temp.2<-apply(temp.1,2,function(x){paste(x,collapse = "")})
        temp.attr.sav<-c(temp.attr.sav,temp.2)
      }
    }
    if(length(x.vector)==1){temp.attr.sav<-x.vector}
    temp.attr.sav
  }
  #vectors needed for combination.generator
  Item.load.id<-list()
  for ( i in 1:nr){
    Item.load.id[[i]]<-grep('1',Qmatrix[i,])}

  Attr.load.id<-list()
  attr.load.id<-matrix(0,length(temp.table.col),nc)
  for ( i in 1:length(temp.table.col)){
    attr.load.id[i,]<-unlist(strsplit(temp.table.col[i],split=''))
    Attr.load.id[[i]]<-grep('1',attr.load.id[i,])
  }

  #Generate Combination for both Item.load and Attr.load
  Item.Comb<-list()
  for ( i in 1:nr){
    Item.Comb[[i]]<-comb.generator(Item.load.id[[i]])
  }
  Attr.Comb<-list()
  for ( i in 2:length(temp.table.col)){
    Attr.Comb[[1]]<-0
    Attr.Comb[[i]]<-comb.generator(Attr.load.id[[i]])
  }
  constraints.list<-list()
  nway.inter.list<-list()
  for(i in 1:nr){
    for(a in 2:length(temp.table.col)){
      ifzero<-as.numeric(paste(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],collapse=''))
      if((!is.na(ifzero))){
        temp.table[i,a]<-paste(c(temp.table[i,a],
                                 paste("S","l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='',collapse='')
        ),collapse='')
        if(a==length(temp.table.col)){
          nway.inter.list[[i]]<-nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])])
          constraints.list[[i]]<-paste("l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='')
        }
      }
    }
  }

  #Create Lambda Table
  Lamda.Table<-temp.table
  for(i in 1:nr){
    for(a in 1:length(Lamda.Table)){
      t.ref<-unique(as.character(Lamda.Table[i,]))
      pos<-c(1:length(t.ref))[Lamda.Table[i,a]==t.ref]
      temp.table[i,a]<-paste("t",i,"_",pos,sep='')}}

  #Generate LCDM specification
  out<-list()
  out[[1]]<-Lamda.Table
  out[[2]]<-temp.table
  out[[3]]<-constraints.list
  out[[4]]<-nway.inter.list
  out[[5]]<-intercept
  OUTPUT<-out
  nclass<-ncol(OUTPUT[[1]]);Nc<-nclass

  #Produce kernel expressions across items and attributes
  Kernel.exp<-OUTPUT[[1]]
  Kernel.exp.detect<-OUTPUT[[1]] #052719updates
  Kernel.exp.dino<-OUTPUT[[1]] #052719updates
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])
      Kernel.exp.detect[i,j]<-NA} #052719updates
    }
  }
  for (i in 1:nrow(OUTPUT[[1]])){ #052719updates
    theClosestEffect<-which(is.na(Kernel.exp.detect[i,]))[1] #052719updates
    useToReplaceLonger<-Kernel.exp[i,theClosestEffect] #052719updates
    Kernel.exp.dino[i,is.na(Kernel.exp.detect[i,])]<-useToReplaceLonger #052719updates
  } #052719updates

  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                        nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))

  if(length(name.inter)!=0){
    for (inter in 1: length(name.inter)){
      temp.nw<-numway.inter[inter]
      temp.nm<-name.inter[inter]
      temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
      temp.sel<-paste(unlist(strsplit(temp.nm,split = '_'))[1],"_",(1:(temp.nw-1)),sep='')
      first.sel<-unlist(OUTPUT[[3]])[grep(paste((temp.sel),collapse="|"),unlist(OUTPUT[[3]]))]
      second.sel<-sub(".*_.", "", first.sel)
      for (sel in 1:length(temp.subnm)){
        SEL<-second.sel[sel]
        Constrain.List1<-rbind(
          paste(temp.nm,">-(0", paste("+",first.sel[grep(SEL,second.sel)],
                                      sep='',collapse=''),")",sep=''),Constrain.List1)
      }
    }
    Constrain.List1<-as.character(Constrain.List1)
  }else{
    Constrain.List1<-NULL
  }

  itemParmName<-c(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==2],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==3],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==4],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==5],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==6],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==7],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==8],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==9],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==10],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==11],OUTPUT[[5]])
  numMainEffect<-length(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1])
  Constrain.List<-paste('  real<lower=0>',itemParmName[1:numMainEffect],';\n ')
  Unconstrain.List<-paste('  real',itemParmName[-(1:numMainEffect)],';\n ')
  Reparm<-as.data.frame(matrix(0,nr,nclass))

  #############################################################
  ###########052719update: Only one main effect is needed #####
  keep.oneMainEffect<-rep(1,nr)   #052719update:
  zero.list<-constraints.list
  fixparm.vec<-NULL #052719update: later they will be fixed to zero

  for(i in 1:nr){
    keep.oneMainEffect[i]<-OUTPUT[[3]][[i]][max(which(OUTPUT[[4]][[i]]==1))]#052719update:
    if(length(zero.list[[i]])==1){zero.list[[i]]=NA}else{#052719update:
      zero.list[[i]]<-constraints.list[[i]][-1]#052719update:
      fixparm.vec<-c(fixparm.vec,zero.list[[i]])#052719update:
    }
  } #intercept, hi.interaction,zero.list/fixparm.vec are what we need
  fixparm.vec<-unlist(constraints.list)[!unlist(constraints.list)%in%keep.oneMainEffect]
  Constrain.List<-paste('  real<lower=0>',keep.oneMainEffect,';\n ')#052719update:
  Unconstrain.List<-paste('  real',intercept,';\n ')
  #############################################################
  #############################################################

  #Produce Stan code for PImat parameter
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      Reparm[loopi,loopc]<-paste('  PImat[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp.dino[loopi,loopc]),');\n',sep='')
    }
  }

  #########052719update:create g and s parameters
  gParm<-rep(0,nr)
  sParm<-rep(0,nr)
  for(loopi in 1:nr){
    gParm[loopi]<-paste('  gParm[',loopi,']=inv_logit(',paste(Kernel.exp.dino[loopi,1]),');\n',sep='')
    sParm[loopi]<-paste('  sParm[',loopi,']=1-inv_logit(',paste(Kernel.exp.dino[loopi,nclass]),');\n',sep='')
  }


  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))
  #############################################################
  ###########052719update: Only one main effect is needed #####
  update.Parmprior<-Parmprior
  fix.Parmprior<-NULL
  for(i in 1:length(Parmprior)){
    if(grepl(paste(fixparm.vec, collapse = "|"), Parmprior[i])){
      update.Parmprior[i]<-""
    }
  }
  update.Parmprior<-update.Parmprior[update.Parmprior!='']
  fix.Parmprior<-c(paste('  real',fixparm.vec,';\n '),
                   paste(' ',fixparm.vec,"=0",';\n ')
  )
  ##therefore we can use: fix.Parmprior,update.Parmprior
  #############################################################
  #############################################################

  #Likelihood Stan code
  Likelihood<-'
  \n
  //Likelihood
  for (iterp in 1:Np){
    for (iterc in 1:Nc){
      for (iteri in 1:Ni){
        if (Y[iterp,iteri] == 1)
          contributionsI[iteri]=bernoulli_lpmf(1|PImat[iteri,iterc]);
        else
          contributionsI[iteri]=bernoulli_lpmf(0|PImat[iteri,iterc]);
      }
      contributionsC[iterc]=log(Vc[iterc])+sum(contributionsI);
    }
  target+=log_sum_exp(contributionsC);
  }
  '


  #Data Specification
  data.spec<-'
  data{
  int Np;
  int Ni;
  int Nc;
  matrix[Np, Ni] Y;
  }
  '
  #Parameter Specification
  parm.spec<-paste(c('
  parameters{
  simplex[Nc] Vc;\n ',paste0(Constrain.List),paste0(Unconstrain.List),
                     '}\n'),collapse='')

  #Reparameter Specification
  transparm.spec<-paste(c('
  transformed parameters{
  matrix[Ni, Nc] PImat;
  vector[Ni] gParm;
  vector[Ni] sParm;\n',
                          gParm, #052719update
                          sParm, #052719update
                          paste0(unlist(Reparm)),'}\n'),collapse='')

  #Model Specification update052619
  model.spec<-paste(c('\nmodel {\n',paste(c(Modelcontainer,update.Parmprior,Likelihood),sep=''),'\n}',sep=''))
  model.spec<-model.spec[!startsWith(str_remove_all(model.spec," "),"~")]
  #Generated Quantities Specification
 generatedQuantities.spec<-'
  \n
generated quantities {

 vector[Ni] log_lik[Np];
 vector[Ni] contributionsI;
 matrix[Ni,Nc] contributionsIC;
 
 matrix[Ni,Nc] posteriorIC;
 matrix[Np,Nc] posteriorPC;



 //Posterior
 for (iterp in 1:Np){
   for (iteri in 1:Ni){
     for (iterc in 1:Nc){
       if (Y[iterp,iteri] == 1)
          contributionsI[iteri]=bernoulli_lpmf(1|PImat[iteri,iterc]);
       else
           contributionsI[iteri]=bernoulli_lpmf(0|PImat[iteri,iterc]);
       contributionsIC[iteri,iterc]=log(Vc[iterc])+contributionsI[iteri];
       posteriorIC[iteri,iterc]=contributionsI[iteri];
      }
      log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
    }
   for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
  }
}
'


  if (.Platform$OS.type == "unix") {
    filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
  }else{
    filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
  }

  sink(file=filename, append=FALSE)
  cat(
    paste(c('   ',
            data.spec,parm.spec,transparm.spec,model.spec,generatedQuantities.spec)
    ))
  sink(NULL)

}

#SEPERATION#
#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanDINO_mG.run<-function(Qmatrix,
                          response.matrix,
                          GroupID,
                          fixeditem.vector=NA,
                          class.equal=T,
                          script.path=NA,save.path=getwd(),save.name="DINO_uninf_multiG",
                          iter=1000,warmup = 0,
                          chain.num=3,init.list='random',control.list=NA){
  group.num<-length(unique(GroupID))
  rstan.detect<-tryCatch(library("rstan"),error=function(e){"rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details."})
  if(length(rstan.detect)==1){
    stop()
  }
  Cdm.init<-F
  if(init.list=='cdm'){
    Cdm.init<-T
    Install.package(c("CDM","stringr"))
    trueParmName<-Parm.name(Qmatrix=Qmatrix)$parm.name
    Classp.exp1<-Parm.name(Qmatrix=Qmatrix)$class.expression
    mod1<-gdina( data =respMatrix, q.matrix = Qmatrix , maxit=700,link = "logit",progress=F)
    CDMresult<-as.data.frame(coef(mod1))
    library(stringr)
    CDM.parm.name<-paste(paste(paste('l',CDMresult[,3],sep=''),'_',sep=''),str_count(CDMresult$partype.attr,"Attr"),sep='')
    CDM.parm.name<-paste(CDM.parm.name,
                         unlist(lapply(strsplit(unlist(lapply(strsplit(CDMresult$partype.attr, 'Attr', fixed=FALSE),function(x){paste(x,collapse="")})),'-'),function(x){paste(x,collapse="")})),
                         sep='')
    CDM.parm.est<-CDMresult$est
    parm.ini<-round(CDM.parm.est[match(trueParmName,CDM.parm.name)],4)
    CDM.prop.est<-mod1$attribute.patt
    prop.ini<-CDM.prop.est[match(Classp.exp1,rownames(CDM.prop.est)),1]
    inilist1<-paste('list(',paste(noquote(paste(noquote(unlist(list(
      paste(paste('Vc=c(',paste((prop.ini),collapse=','),')',collapse=','))))
    ))),collapse=',')  ,')',collapse='')

    IniList1<-NULL
    if(!class.equal){
      temp.inilist1<-eval(parse(n =2000000 ,text=inilist1))
      eval(parse(n =2000000 ,text=paste(paste('IniList1$Vc_g',1:group.num,sep=''),"<-temp.inilist1$Vc",sep='') ))
      inilist1<-IniList1
    }else{
      inilist1<-eval(parse(n =2000000 ,text=inilist1))}


    for( i in 2:chain.num){
      temp.text<-paste('inilist',i,"<-inilist1",sep='')
      eval(parse(text=(temp.text)))
    }
    temp.text<-paste('init.list<-list(',paste(paste('inilist',1:chain.num,sep=''),collapse = ","),')',sep='')
    eval(parse(text=(temp.text)))
  }
  data.list<-Generate.datalist(Qmatrix,response.matrix,GroupID)

  if(is.na(control.list)){control.list<-list(adapt_delta=0.82)}
  if(is.na(script.path)==T){
    options(warn=-1)
    #Need to update script
    StanDINO_mG.script(Qmatrix=Qmatrix,
                       group.num=group.num,
                       fixeditem.vector=fixeditem.vector,
                       class.equal=class.equal,
                       save.path=save.path,save.name=save.name)
    script.path<-paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
    options(warn=0)
    compiled_model<-stan_model(script.path)}
  else{
    compiled_model<-stan_model(script.path)
  }
  if(Cdm.init==T){
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})
  }else{
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})

  }

  estimated_model
}


#SEPERATION#
#' @title Generate Stan code and Run the estimation for StanLCDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanDINO_mG.script<-function(Qmatrix,
                             group.num,
                             fixeditem.vector=NA,
                             class.equal=T,
                             save.path=getwd(),save.name="DINO_uninf_multiG"){

  #Load packages
  Install.package("plyr")
  Install.package('stringr')

  nc<-ncol(Qmatrix)
  nr<-nrow(Qmatrix)
  temp.table.col<-unique(apply(combn(rep(c(0,1),nc),nc),2,function(x){paste(x,collapse = "")}))
  temp.table.col<-temp.table.col[order(temp.table.col)]
  temp.table<-matrix(0,nr,length(temp.table.col))
  colnames(temp.table)<-temp.table.col
  rownames(temp.table)<-paste('item',c(1:nr),sep='')
  temp.table<-as.data.frame(temp.table)
  for (i in 1:nr){
    temp.table[i,]<-paste('l',i,'_0',sep='')
  }
  intercept<-temp.table[,1]

  #Generate attribute combinations
  comb.generator<-function(x.vector){
    if(length(x.vector)>1){
      temp.attr<-x.vector
      temp.attr.sav<-NULL
      for(i in 1:length(temp.attr)){
        temp.1<-combn(temp.attr,i)
        temp.2<-apply(temp.1,2,function(x){paste(x,collapse = "")})
        temp.attr.sav<-c(temp.attr.sav,temp.2)
      }
    }
    if(length(x.vector)==1){temp.attr.sav<-x.vector}
    temp.attr.sav
  }
  #vectors needed for combination.generator
  Item.load.id<-list()
  for ( i in 1:nr){
    Item.load.id[[i]]<-grep('1',Qmatrix[i,])}

  Attr.load.id<-list()
  attr.load.id<-matrix(0,length(temp.table.col),nc)
  for ( i in 1:length(temp.table.col)){
    attr.load.id[i,]<-unlist(strsplit(temp.table.col[i],split=''))
    Attr.load.id[[i]]<-grep('1',attr.load.id[i,])
  }

  #Generate Combination for both Item.load and Attr.load
  Item.Comb<-list()
  for ( i in 1:nr){
    Item.Comb[[i]]<-comb.generator(Item.load.id[[i]])
  }
  Attr.Comb<-list()
  for ( i in 2:length(temp.table.col)){
    Attr.Comb[[1]]<-0
    Attr.Comb[[i]]<-comb.generator(Attr.load.id[[i]])
  }
  constraints.list<-list()
  nway.inter.list<-list()
  for(i in 1:nr){
    for(a in 2:length(temp.table.col)){
      ifzero<-as.numeric(paste(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],collapse=''))
      if((!is.na(ifzero))){
        temp.table[i,a]<-paste(c(temp.table[i,a],
                                 paste("S","l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='',collapse='')
        ),collapse='')
        if(a==length(temp.table.col)){
          nway.inter.list[[i]]<-nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])])
          constraints.list[[i]]<-paste("l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='')
        }
      }
    }
  }

  #Create Lambda Table
  Lamda.Table<-temp.table
  for(i in 1:nr){
    for(a in 1:length(Lamda.Table)){
      t.ref<-unique(as.character(Lamda.Table[i,]))
      pos<-c(1:length(t.ref))[Lamda.Table[i,a]==t.ref]
      temp.table[i,a]<-paste("t",i,"_",pos,sep='')}}

  #Generate LCDM specification
  out<-list()
  out[[1]]<-Lamda.Table
  out[[2]]<-temp.table
  out[[3]]<-constraints.list
  out[[4]]<-nway.inter.list
  out[[5]]<-intercept
  OUTPUT<-out
  nclass<-ncol(OUTPUT[[1]]);Nc<-nclass

  #Produce kernel expressions across items and attributes
  Kernel.exp<-OUTPUT[[1]]
  Kernel.exp.detect<-OUTPUT[[1]] #052719updates
  Kernel.exp.dino<-OUTPUT[[1]] #052719updates
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])
      Kernel.exp.detect[i,j]<-NA} #052719updates
    }
  }
  for (i in 1:nrow(OUTPUT[[1]])){ #052719updates
    theClosestEffect<-which(is.na(Kernel.exp.detect[i,]))[1] #052719updates
    useToReplaceLonger<-Kernel.exp[i,theClosestEffect] #052719updates
    Kernel.exp.dino[i,is.na(Kernel.exp.detect[i,])]<-useToReplaceLonger #052719updates
  } #052719updates

  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                        nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))

  for (inter in 1: length(name.inter)){
    temp.nw<-numway.inter[inter]
    temp.nm<-name.inter[inter]
    temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
    temp.sel<-paste(unlist(strsplit(temp.nm,split = '_'))[1],"_",(1:(temp.nw-1)),sep='')
    first.sel<-unlist(OUTPUT[[3]])[grep(paste((temp.sel),collapse="|"),unlist(OUTPUT[[3]]))]
    second.sel<-sub(".*_.", "", first.sel)
    for (sel in 1:length(temp.subnm)){
      SEL<-second.sel[sel]
      Constrain.List1<-rbind(
        paste(temp.nm,">-(0", paste("+",first.sel[grep(SEL,second.sel)],
                                    sep='',collapse=''),")",sep=''),Constrain.List1)
    }
  }
  Constrain.List1<-as.character(Constrain.List1)

  itemParmName<-c(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==2],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==3],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==4],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==5],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==6],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==7],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==8],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==9],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==10],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==11],OUTPUT[[5]])
  numMainEffect<-length(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1])
  Constrain.List<-paste('  real<lower=0>',itemParmName[1:numMainEffect],';\n ')
  Unconstrain.List<-paste('  real',itemParmName[-(1:numMainEffect)],';\n ')
  Reparm<-as.data.frame(matrix(0,nr,nclass))

  #############################################################
  ###########052719update: Only one main effect is needed #####
  keep.oneMainEffect<-rep(1,nr)   #052719update:
  zero.list<-constraints.list
  fixparm.vec<-NULL #052719update: later they will be fixed to zero

  for(i in 1:nr){
    keep.oneMainEffect[i]<-OUTPUT[[3]][[i]][max(which(OUTPUT[[4]][[i]]==1))]#052719update:
    if(length(zero.list[[i]])==1){zero.list[[i]]=NA}else{#052719update:
      zero.list[[i]]<-constraints.list[[i]][-1]#052719update:
      fixparm.vec<-c(fixparm.vec,zero.list[[i]])#052719update:
    }
  } #intercept, hi.interaction,zero.list/fixparm.vec are what we need
  fixparm.vec<-unlist(constraints.list)[!unlist(constraints.list)%in%keep.oneMainEffect]
  Constrain.List<-paste('  real<lower=0>',keep.oneMainEffect,';\n ')#052719update:
  Unconstrain.List<-paste('  real',intercept,';\n ')


  #############################################################
  #######060319update: Multiple Group##########################
  #############################################################
  if(!is.na(fixeditem.vector)[1]){
    fixedItem.vector<-c(1:nr)[-fixeditem.vector]
  }else{fixedItem.vector<-c(1:nr)}
  #############################################################
  #######060319update: Multiple Group##########################
  #############################################################

  Kernel.exp.dino.groupName<-paste("Kernel.exp.dino_g",c(1:group.num),sep='')
  for(i in 1:group.num){
    tempfill.Kernel.exp.dino<-Kernel.exp.dino
    temp.Kernel.exp.dino<-Kernel.exp.dino[fixedItem.vector,]
    for(j in 1:nrow(temp.Kernel.exp.dino)){
      for(z in 1:ncol(temp.Kernel.exp.dino)){
        temp.Kernel.exp.dino[j,z]<-paste(temp.Kernel.exp.dino[j,z],'_g',i,sep='')
      }
    }
    for(j in 1:nrow(temp.Kernel.exp.dino)){
      for(z in 1:ncol(temp.Kernel.exp.dino)){
        temp.Kernel.exp.dino[j,z]<-str_replace_all(temp.Kernel.exp.dino[j,z],"\\+",paste("_g",i,"+",sep=''))
      }
    }
    tempfill.Kernel.exp.dino[fixedItem.vector,]<-temp.Kernel.exp.dino
    assign(Kernel.exp.dino.groupName[i],tempfill.Kernel.exp.dino)
  }
  Kernel.exp.dino.list<-list()
  for(i in 1:group.num){Kernel.exp.dino.list[[i]]<-eval(parse(text=paste("Kernel.exp.dino_g",i,sep='')))}
  #############################################################
  ##########060319update: Multiple Group End###################
  #############################################################

  #############################################################
  ##########060319update: Multiple Group#######################
  #############################################################

  PImat.groupName<-paste("PImat_g",c(1:group.num),sep='')
  Reparm.multigroup<-array(0,dim = c(nr,nclass,group.num))
  #Produce Stan code for PImat parameter
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      Reparm[loopi,loopc]<-paste('  PImat[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp.dino[loopi,loopc]),');\n',sep='')
    }
  }
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      for (loopg in 1:group.num){
        Reparm.multigroup[loopi,loopc,loopg]<-paste('  PImat_g',loopg,'[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp.dino.list[[loopg]][loopi,loopc]),');\n',sep='')

      }
    }
  }
  #############################################################
  ##########060319update: Multiple Group End###################
  #############################################################
  #########052719update:create g and s parameters
  gParm<-rep(0,nr)
  sParm<-rep(0,nr)
  for(loopi in 1:nr){
    gParm[loopi]<-paste('  gParm[',loopi,']=inv_logit(',paste(Kernel.exp.dino[loopi,1]),');\n',sep='')
    sParm[loopi]<-paste('  sParm[',loopi,']=1-inv_logit(',paste(Kernel.exp.dino[loopi,nclass]),');\n',sep='')
  }
  #############################################################
  ##########060319update: Multiple Group ######################
  #############################################################
  gParm.multigroup<-array(0,dim = c(nr,1,group.num))
  sParm.multigroup<-array(0,dim = c(nr,1,group.num))
  for(i in 1:group.num){
    tempfill.gParm<-gParm
    tempfill.sParm<-sParm

    tempfill.gParm<-str_replace_all(tempfill.gParm,"gParm",paste("gParm_g",i,sep=''))
    tempfill.sParm<-str_replace_all(tempfill.sParm,"sParm",paste("sParm_g",i,sep=''))

    temp.gParm<-gParm[fixedItem.vector]
    temp.sParm<-sParm[fixedItem.vector]

    for(j in 1:length(temp.gParm)){
      temp.gParm[j]<-str_replace_all(temp.gParm[j],"\\)",paste("_g",i,")",sep=''))
      temp.sParm[j]<-str_replace_all(temp.sParm[j],"\\)",paste("_g",i,")",sep=''))
    }
    for(j in 1:length(temp.gParm)){
      temp.gParm[j]<-str_replace_all(temp.gParm[j],"\\+",paste("_g",i,")",sep=''))
      temp.sParm[j]<-str_replace_all(temp.sParm[j],"\\+",paste("_g",i,"+",sep=''))
    }

    temp.gParm<-str_replace_all(temp.gParm,"gParm",paste("gParm_g",i,sep=''))
    temp.sParm<-str_replace_all(temp.sParm,"sParm",paste("sParm_g",i,sep=''))

    tempfill.gParm[fixedItem.vector]<-temp.gParm
    tempfill.sParm[fixedItem.vector]<-temp.sParm

    gParm.multigroup[,,i]<-tempfill.gParm
    sParm.multigroup[,,i]<-tempfill.sParm

  }

  #############################################################
  ##########060319update: Multiple Group ######################
  #############################################################
  keep.oneMainEffect.multigroup<-NULL
  intercept.multigroup<-array(0,dim = c(nr,1,group.num))
  mainEff.multigroup<-array(0,dim = c(numMainEffect,1,group.num))
  interaction.multigroup<-array(0,dim = c(length(name.inter),1,group.num))
  #Group Invariant Parameter Name
  fixedParmName<-NULL
  if(!is.na(fixeditem.vector)[1]){
    for(i in fixeditem.vector){
      fixedParmName<-c(fixedParmName,out[[3]][[i]])
    }
    fixedParmName<-c(fixedParmName,out[[5]][fixeditem.vector])
  }
  #Group Variant Parameter Name
  freeParmName<-c(out[[5]],itemParmName)[!c(out[[5]],itemParmName)%in%fixedParmName]

  for(i in 1:group.num){
    tempfill.intercept<-out[[5]]
    tempfill.mainEff<-itemParmName[1:numMainEffect]
    tempfill.interaction<-name.inter
    tempfill.keep.oneMainEffect<-keep.oneMainEffect

    temp.intercept<-tempfill.intercept[fixedItem.vector]
    temp.mainEff<-tempfill.mainEff[!tempfill.mainEff%in%fixedParmName]
    temp.interaction<-tempfill.interaction[!tempfill.interaction%in%fixedParmName]
    temp.keep.oneMainEffect<-tempfill.keep.oneMainEffect[!tempfill.keep.oneMainEffect%in%fixedParmName]

    for(j in 1:length(temp.intercept)){
      temp.intercept[j]<-paste(temp.intercept[j],"_g",i,sep='')
    }

    for(j in 1:length(temp.mainEff)){
      temp.mainEff[j]<-paste(temp.mainEff[j],"_g",i,sep='')
    }

    for(j in 1:length(temp.interaction)){
      temp.interaction[j]<-paste(temp.interaction[j],"_g",i,sep='')
    }

    for(j in 1:length(temp.keep.oneMainEffect)){
      temp.keep.oneMainEffect[j]<-paste(temp.keep.oneMainEffect[j],"_g",i,sep='')
    }

    tempfill.intercept[fixedItem.vector]<-temp.intercept
    tempfill.mainEff[!tempfill.mainEff%in%fixedParmName]<-temp.mainEff
    tempfill.interaction[!tempfill.interaction%in%fixedParmName]<-temp.interaction
    tempfill.keep.oneMainEffect[!tempfill.keep.oneMainEffect%in%fixedParmName]<-temp.keep.oneMainEffect

    intercept.multigroup[,,i]<-tempfill.intercept
    mainEff.multigroup[,,i]<-tempfill.mainEff
    interaction.multigroup[,,i]<-tempfill.interaction
    keep.oneMainEffect.multigroup<-c(keep.oneMainEffect.multigroup,tempfill.keep.oneMainEffect)
  }

  Constrain.List<-paste('  real<lower=0>',unique(keep.oneMainEffect.multigroup),';\n ')
  Unconstrain.List<-paste('  real',unique(c(intercept.multigroup)),';\n ')
  #############################################################
  ##########060319update: Multiple Group End###################
  #############################################################


  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))
  update.Parmprior<-Parmprior
  fix.Parmprior<-NULL
  for(i in 1:length(Parmprior)){
    if(grepl(paste(fixparm.vec, collapse = "|"), Parmprior[i])){
      update.Parmprior[i]<-""
    }
  }
  update.Parmprior<-update.Parmprior[update.Parmprior!='']


  fix.Parmprior<-c(paste('  real',fixparm.vec,';\n '),
                   paste(' ',fixparm.vec,"=0",';\n ')
  )
  #############################################################
  #####060319update:Change from update.Parmprior&fix     ######
  #############################################################
  update.Parmprior.multiGroup<-NULL
  for(i in 1:length(update.Parmprior)){
    for (j in 1:length(freeParmName)){
      for (z in 1:group.num){
        if(sum(grepl(freeParmName[j],update.Parmprior[i]))>=1){
          temp.update.Parmprior<-str_replace_all(update.Parmprior[i],freeParmName[j],paste(freeParmName[j],"_g",z,sep=''))
          update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,
                                         temp.update.Parmprior)
        }
      }
    }
  }
  update.Parmprior.multiGroup<-unique(update.Parmprior.multiGroup)
  update.Parmprior.multiGroup<-c("   //Prior\n",paste('   ',fixedParmName,'~normal(0,5)',';\n',sep=''),update.Parmprior.multiGroup )
  if(class.equal){
    update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='') )
  }else{
    for(i in 1:group.num){
      update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,
                                     paste('   Vc_g',i,'~dirichlet(rep_vector(2.0, Nc));\n',sep='') )
    }
  }

  ##therefore we can use: fix.Parmprior,update.Parmprior
  #############################################################
  #############################################################

  #############################################################
  #####060419update:Likelihood Add PImat_g#####################
  #############################################################
  PImat.likelihood1<-NULL
  PImat.likelihood0<-NULL
  Vc.likelihood<-NULL
  for(loopg in 1:group.num){
    temp.PImat.likelihood1<-paste(paste('          if (GroupID[iterp]==',loopg,')',sep=''),
                                  paste('            contributionsI[iteri]=bernoulli_lpmf(1|PImat_g',loopg,'[iteri,iterc]);',sep=''),
                                  sep='\n')
    PImat.likelihood1<-paste(PImat.likelihood1,temp.PImat.likelihood1,sep='\n')
    temp.PImat.likelihood0<-paste(paste('          if (GroupID[iterp]==',loopg,')',sep=''),
                                  paste('            contributionsI[iteri]=bernoulli_lpmf(0|PImat_g',loopg,'[iteri,iterc]);',sep=''),
                                  sep='\n')
    PImat.likelihood0<-paste(PImat.likelihood0,temp.PImat.likelihood0,sep='\n')
  }
  if(!class.equal){
    for(loopg in 1:group.num){
      temp.Vc.likelihood<-paste(paste('       if (GroupID[iterp]==',loopg,')',sep=''),
                                paste('         contributionsC[iterc]=log(Vc_g',loopg,'[iterc])+sum(contributionsI);',sep=''),
                                sep='\n')
      Vc.likelihood<-paste(Vc.likelihood,temp.Vc.likelihood,sep='\n')
    }
  }

  #Likelihood Stan code
  if(class.equal){
    Likelihood<-paste('
  \n
  //Likelihood
  for (iterp in 1:Np){
    for (iterc in 1:Nc){
      for (iteri in 1:Ni){
        if (Y[iterp,iteri] == 1)'
                      ,PImat.likelihood1,'\n',
                      '        else'
                      ,PImat.likelihood0,
                      '}
      contributionsC[iterc]=log(Vc[iterc])+sum(contributionsI);
    }
  target+=log_sum_exp(contributionsC);
  }
  ',sep='')}else{
    Likelihood<-paste('
  \n
  //Likelihood
  for (iterp in 1:Np){
    for (iterc in 1:Nc){
      for (iteri in 1:Ni){
        if (Y[iterp,iteri] == 1)'
                      ,PImat.likelihood1,'\n',
                      '        else'
                      ,PImat.likelihood0,
                      '}\n',
                      Vc.likelihood
                      ,'


  }
  target+=log_sum_exp(contributionsC);
 }
                      ',sep='')

  }
  #############################################################
  #####060419update:Likelihood Add PImat_g  end################
  #############################################################

  #Data Specification
  data.spec<-'
  data{
  int Np;
  int Ni;
  int Nc;
  matrix[Np, Ni] Y;
  vector[Np] GroupID;
  }
  '


  #############################################################
  #####060419update: Stan script##############################
  #############################################################

  #Parameter Specification
  if(class.equal){parm.spec<-paste(c('
  parameters{
  simplex[Nc] Vc;\n ',paste0(Constrain.List),paste0(Unconstrain.List),
                                     '}\n'),collapse='')}else{
                                       parm.spec<-paste(c('
  parameters{\n ',
                                                          paste(paste('   simplex[Nc] Vc_g',1:group.num, ";",sep=''),"\n"),
                                                          paste0(Constrain.List),paste0(Unconstrain.List),
                                                          '}\n'),collapse='')


                                     }

  #Reparameter Specification
  transparm.spec<-paste(c('
  transformed parameters{

                          ',
                          paste('  matrix[Ni, Nc] PImat_g',1:group.num,';\n',sep=''),
                          paste('  vector[Ni] gParm_g',1:group.num,';\n',sep=''),
                          paste('  vector[Ni] sParm_g',1:group.num,';\n',sep=''),

                          c(gParm.multigroup), #060419update
                          c(sParm.multigroup), #060419update
                          paste0(c(Reparm.multigroup)),'}\n'),collapse='')

  #Model Specification update052619
  model.spec<-paste(c('\nmodel {\n',paste(c(Modelcontainer,update.Parmprior.multiGroup,Likelihood),sep=''),'\n}',sep=''))

  #Generated Quantities Specification
  IC.generatedquantities<-NULL
  if(!class.equal){
    for(loopg in 1:group.num){
      temp.IC.generatedquantities<-paste(paste('       if (GroupID[iterp]==',loopg,')',sep=''),
                                         paste('         contributionsIC[iteri,iterc]=log(Vc_g',loopg,'[iterc])+contributionsI[iteri];',sep=''),
                                         sep='\n')
      IC.generatedquantities<-paste(IC.generatedquantities,temp.IC.generatedquantities,sep='\n')
    }
  }


   if(class.equal){
    generatedQuantities.spec<-paste('
  \n
  generated quantities {
  vector[Ni] log_lik[Np];
  vector[Ni] contributionsI;
  matrix[Ni,Nc] contributionsIC;

  matrix[Ni,Nc] posteriorIC;
  matrix[Np,Nc] posteriorPC;


  //Posterior
  for (iterp in 1:Np){
    for (iteri in 1:Ni){
      for (iterc in 1:Nc){
        if (Y[iterp,iteri] == 1)'
                                    ,PImat.likelihood1,'\n',
                                    '        else'
                                    ,PImat.likelihood0,
                                    '
          contributionsIC[iteri,iterc]=log(Vc[iterc])+contributionsI[iteri];
          posteriorIC[iteri,iterc]=contributionsI[iteri];
        }
      log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
    }
    for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
  }
  }
  ',sep='')}else{
    generatedQuantities.spec <- paste( '
    \n
    generated quantities {
    vector[Ni] log_lik[Np];
    vector[Ni] contributionsI;
    matrix[Ni,Nc] contributionsIC;

    matrix[Ni,Nc] posteriorIC;
    matrix[Np,Nc] posteriorPC;


    //Posterior
    for (iterp in 1:Np){
      for (iteri in 1:Ni){
        for (iterc in 1:Nc){
          if (Y[iterp,iteri] == 1)',
                                       PImat.likelihood1,'\n',
                                       '        else',
                                       PImat.likelihood0,'\n',
                                       "   ",IC.generatedquantities,'\n',
                             '          posteriorIC[iteri,iterc]=contributionsI[iteri];','\n
      }\n',
                                       '     log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
     }
     for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
   }
   }
      ',
                                       sep = ''
    )
  }




  if (.Platform$OS.type == "unix") {
    filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
  }else{
    filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
  }

  sink(file=filename, append=FALSE)
  cat(
    paste(c('   ',
            data.spec,parm.spec,transparm.spec,model.spec,generatedQuantities.spec)
    ))
  sink(NULL)

}

#SEPERATION#
#' @title Generate Stan code and Run the estimation for LCDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\data.RData")


StanLCDM.run<-function(Qmatrix,response.matrix,script.path=NA,save.path=getwd(),save.name="LCDM_uninf",iter=1000,warmup = 0,
                       chain.num=3, init.list='random',control.list=NA){
  rstan.detect<-tryCatch(library("rstan"),error=function(e){"rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details."})
  if(length(rstan.detect)==1){
    stop()
  }
  Cdm.init = F
  if(init.list=='cdm'){
    Cdm.init = T
    Install.package(c("CDM","stringr"))
    trueParmName<-Parm.name(Qmatrix=Qmatrix)$parm.name
    Classp.exp1<-Parm.name(Qmatrix=Qmatrix)$class.expression
    mod1<-gdina( data =response.matrix, q.matrix = Qmatrix , maxit=700,link = "logit",progress=F)
    CDMresult<-as.data.frame(coef(mod1))
    library(stringr)
    CDM.parm.name<-paste(paste(paste('l',CDMresult[,3],sep=''),'_',sep=''),str_count(CDMresult$partype.attr,"Attr"),sep='')
    CDM.parm.name<-paste(CDM.parm.name,
                         unlist(lapply(strsplit(unlist(lapply(strsplit(CDMresult$partype.attr, 'Attr', fixed=FALSE),function(x){paste(x,collapse="")})),'-'),function(x){paste(x,collapse="")})),
                         sep='')
    CDM.parm.est<-CDMresult$est
    # remove nagetive values in initial values
    CDM.parm.est[CDM.parm.est < 0]  <-  0.01
    parm.ini<-round(CDM.parm.est[match(trueParmName,CDM.parm.name)],4)
    parm.ini[abs(parm.ini)>10] <- 10*sign(parm.ini[abs(parm.ini)>10])
    parm.ini[parm.ini==0] <- 0.05

    CDM.prop.est<-mod1$attribute.patt
    prop.ini<-CDM.prop.est[match(Classp.exp1,rownames(CDM.prop.est)),1]
    inilist1<-paste('list(',paste(noquote(paste(noquote(unlist(list(paste(trueParmName,'=',round(parm.ini,2)),
                                                                    paste(paste('Vc=c(',paste((prop.ini),collapse=','),')',collapse=','))))
    ))),collapse=',')  ,')',collapse='')

    inilist1<-eval(parse(n =2000000 ,text=inilist1))
    for( i in 2:chain.num){
      temp.text<-paste('inilist',i,"<-inilist1",sep='')
      eval(parse(text=(temp.text)))
    }
    temp.text<-paste('init.list<-list(',paste(paste('inilist',1:chain.num,sep=''),collapse = ","),')',sep='')
    eval(parse(text=(temp.text)))
  }
  data.list<-Generate.datalist(Qmatrix,response.matrix)

  if(is.na(control.list)){control.list<-list(adapt_delta=0.82)}
  if(is.na(script.path)==T){
    options(warn=-1)
    StanLCDM.script(Qmatrix,save.path=save.path,save.name=save.name)
    if (.Platform$OS.type == "unix") {
      filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
    }else{
      filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
    }
    script.path<-filename
    options(warn=0)
    compiled_model<-stan_model(script.path)
  }else{
    compiled_model<-stan_model(script.path)
  }
  if(Cdm.init == T){
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains= chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})
  }else{
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains= chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})

  }

  estimated_model
}



#SEPERATION#
#' @title Generate Stan code and Run the estimation for LCDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanLCDM.script<-function(Qmatrix,save.path=getwd(),save.name="LCDM_uninf"){
  #Load packages
  Install.package("plyr")
  Install.package('stringr')

  nc<-ncol(Qmatrix)
  nr<-nrow(Qmatrix)
  temp.table.col<-unique(apply(combn(rep(c(0,1),nc),nc),2,function(x){paste(x,collapse = "")}))
  temp.table.col<-temp.table.col[order(temp.table.col)]
  temp.table<-matrix(0,nr,length(temp.table.col))
  colnames(temp.table)<-temp.table.col
  rownames(temp.table)<-paste('item',c(1:nr),sep='')
  temp.table<-as.data.frame(temp.table)
  for (i in 1:nr){
    temp.table[i,]<-paste('l',i,'_0',sep='')
  }
  intercept<-temp.table[,1]

  #Generate attribute combinations
  comb.generator<-function(x.vector){
    if(length(x.vector)>1){
      temp.attr<-x.vector
      temp.attr.sav<-NULL
      for(i in 1:length(temp.attr)){
        temp.1<-combn(temp.attr,i)
        temp.2<-apply(temp.1,2,function(x){paste(x,collapse = "")})
        temp.attr.sav<-c(temp.attr.sav,temp.2)
      }
    }
    if(length(x.vector)==1){temp.attr.sav<-x.vector}
    temp.attr.sav
  }
  #vectors needed for combination.generator
  Item.load.id<-list()
  for ( i in 1:nr){
    Item.load.id[[i]]<-grep('1',Qmatrix[i,])}

  Attr.load.id<-list()
  attr.load.id<-matrix(0,length(temp.table.col),nc)
  for ( i in 1:length(temp.table.col)){
    attr.load.id[i,]<-unlist(strsplit(temp.table.col[i],split=''))
    Attr.load.id[[i]]<-grep('1',attr.load.id[i,])
  }

  #Generate Combination for both Item.load and Attr.load
  Item.Comb<-list()
  for ( i in 1:nr){
    Item.Comb[[i]]<-comb.generator(Item.load.id[[i]])
  }
  Attr.Comb<-list()
  for ( i in 2:length(temp.table.col)){
    Attr.Comb[[1]]<-0
    Attr.Comb[[i]]<-comb.generator(Attr.load.id[[i]])
  }
  constraints.list<-list()
  nway.inter.list<-list()
  for(i in 1:nr){
    for(a in 2:length(temp.table.col)){
      ifzero<-as.numeric(paste(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],collapse=''))
      if((!is.na(ifzero))){
        temp.table[i,a]<-paste(c(temp.table[i,a],
                                 paste("S","l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='',collapse='')
        ),collapse='')
        if(a==length(temp.table.col)){
          nway.inter.list[[i]]<-nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])])
          constraints.list[[i]]<-paste("l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='')
        }
      }
    }
  }

  #Create Lambda Table
  Lamda.Table<-temp.table
  for(i in 1:nr){
    for(a in 1:length(Lamda.Table)){
      t.ref<-unique(as.character(Lamda.Table[i,]))
      pos<-c(1:length(t.ref))[Lamda.Table[i,a]==t.ref]
      temp.table[i,a]<-paste("t",i,"_",pos,sep='')}}

  #Generate LCDM specification
  out<-list()
  out[[1]]<-Lamda.Table
  out[[2]]<-temp.table
  out[[3]]<-constraints.list
  out[[4]]<-nway.inter.list
  out[[5]]<-intercept
  OUTPUT<-out
  nclass<-ncol(OUTPUT[[1]]);Nc<-nclass

  #Produce kernel expressions across items and attributes
  Kernel.exp<-OUTPUT[[1]]
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])}
    }
  }


  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                      nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))

  if(length(name.inter)!=0){
    for (inter in 1: length(name.inter)){
      temp.nw<-numway.inter[inter]
      temp.nm<-name.inter[inter]
      temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
      temp.sel<-paste(unlist(strsplit(temp.nm,split = '_'))[1],"_",(1:(temp.nw-1)),sep='')
      first.sel<-unlist(OUTPUT[[3]])[grep(paste((temp.sel),collapse="|"),unlist(OUTPUT[[3]]))]
      second.sel<-sub(".*_.", "", first.sel)
      for (sel in 1:length(temp.subnm)){
        SEL<-second.sel[sel]
        Constrain.List1<-rbind(
          paste(temp.nm,">-(0", paste("+",first.sel[grep(SEL,second.sel)],
                                      sep='',collapse=''),")",sep=''),Constrain.List1)
      }
    }
    Constrain.List1<-as.character(Constrain.List1)
  }else{
    Constrain.List1<-NULL
  }

  itemParmName<-c(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==2],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==3],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==4],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==5],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==6],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==7],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==8],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==9],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==10],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==11],OUTPUT[[5]])
  numMainEffect<-length(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1])
  Constrain.List<-paste('  real<lower=0>',itemParmName[1:numMainEffect],';\n ')
  Unconstrain.List<-paste('  real',itemParmName[-(1:numMainEffect)],';\n ')
  Reparm<-as.data.frame(matrix(0,nr,nclass))


  #Produce Stan code for PImat parameter
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      Reparm[loopi,loopc]<-paste('  PImat[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp[loopi,loopc]),');\n',sep='')
    }
  }

  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))
  #Likelihood Stan code
  Likelihood<-'
  \n
  //Likelihood
  for (iterp in 1:Np){
    for (iterc in 1:Nc){
      for (iteri in 1:Ni){
        if (Y[iterp,iteri] == 1)
          contributionsI[iteri]=bernoulli_lpmf(1|PImat[iteri,iterc]);
        else
          contributionsI[iteri]=bernoulli_lpmf(0|PImat[iteri,iterc]);
      }
      contributionsC[iterc]=log(Vc[iterc])+sum(contributionsI);
    }
  target+=log_sum_exp(contributionsC);
  }
  '


  #Data Specification
  data.spec<-'
data{
  int Np;
  int Ni;
  int Nc;
  matrix[Np, Ni] Y;
}
  '

  #Parameter Specification
  parm.spec<-paste(c('
parameters{
  simplex[Nc] Vc;\n ',paste0(Constrain.List),paste0(Unconstrain.List),
'}\n'),collapse='')

  #Reparameter Specification
  transparm.spec<-paste(c('
  transformed parameters{
  matrix[Ni, Nc] PImat;\n',
  paste0(unlist(Reparm)),'}\n'),collapse='')

  #Model Specification
  model.spec<-paste(c('\nmodel {\n',paste(c(Modelcontainer,Parmprior,Likelihood),sep=''),'\n}',sep=''))
  model.spec<-model.spec[!startsWith(str_remove_all(model.spec," "),"~")]

  #Generated Quantities Specification

generatedQuantities.spec<-'
  \n
generated quantities {

 vector[Ni] log_lik[Np];
 vector[Ni] contributionsI;
 matrix[Ni,Nc] contributionsIC;
 
 matrix[Ni,Nc] posteriorIC;
 matrix[Np,Nc] posteriorPC;



 //Posterior
 for (iterp in 1:Np){
   for (iteri in 1:Ni){
     for (iterc in 1:Nc){
       if (Y[iterp,iteri] == 1)
          contributionsI[iteri]=bernoulli_lpmf(1|PImat[iteri,iterc]);
       else
           contributionsI[iteri]=bernoulli_lpmf(0|PImat[iteri,iterc]);
       contributionsIC[iteri,iterc]=log(Vc[iterc])+contributionsI[iteri];
       posteriorIC[iteri,iterc]=contributionsI[iteri];
      }
      log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
    }
   for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
  }
}
'
  if (.Platform$OS.type == "unix") {
    filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
  }else{
    filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
  }

  sink(file= filename,append=FALSE)
    cat(
    paste(c('   ',
    data.spec,parm.spec,transparm.spec,model.spec,generatedQuantities.spec)
  ))
  sink(NULL)

}


#SEPERATION#

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

#SEPERATION#

StanLCDM_gC.script<-function(Qmatrix,time.vector,quad.structure=F,variance.equal=F, save.path=getwd(),save.name="LCDM_gC_uninf"){
  #Load packages
  Install.package("plyr")
  Install.package('stringr')
  nc<-ncol(Qmatrix);Na<-nc
  nr<-nrow(Qmatrix);
  
  
  #R commend
  attrSamplingStan<-NULL
  for(j in 1:nc){
    attrSamplingStan<-c(attrSamplingStan,
                        paste("A",j,"[iterp,itero] ~ bernoulli_logit(Beta0[",j,"]+Beta1[",j,"]*Theta[iterp,itero])",sep=''))
  }
  attrSamplingStan<-paste(attrSamplingStan      ,"; \n")
  
  
  attrPosteriorSamplingStan<-NULL
  for(j in 1:nc){
    attrPosteriorSamplingStan<-c(attrPosteriorSamplingStan,
                                 paste("A",j,"[iterp,itero] <- bernoulli_rng(inv_logit(Beta0[",j,"]+Beta1[",j,"]*Theta[iterp,itero]))",sep=''))
  }
  attrPosteriorSamplingStan<-paste(attrPosteriorSamplingStan      ,"; \n")
  
  
  defineA<-paste(paste(paste('matrix[Np,No] A',c(1:Na),sep=''),";",sep=''))
  defineAlpha<-paste(paste(paste('matrix[Np,No] Alpha',c(1:Na),sep=''),";",sep=''))
  
  
  
  defineAThreshold<-paste(paste("real <lower=0, upper=1> AThres",1:Na,sep=""),"[No];",sep='')
  AThresholdPrior<-paste(paste(" AThres",1:Na,sep=''),"~beta(2,2);",sep="")
  
  ####0617update: New parameters for GC#########
  occasion.num<-length(time.vector)
  group.num<-occasion.num
  if(quad.structure){
    thridOrder.num<-3
    timeQuad.vector<-(time.vector)^2
  }else{(thridOrder.num<-2)
    timeQuad.vector<-NA
  }
  No<-occasion.num
  
  temp.table.col<-unique(apply(combn(rep(c(0,1),nc),nc),2,function(x){paste(x,collapse = "")}))
  temp.table.col<-temp.table.col[order(temp.table.col)]
  temp.table<-matrix(0,nr,length(temp.table.col))
  colnames(temp.table)<-temp.table.col
  rownames(temp.table)<-paste('item',c(1:nr),sep='')
  temp.table<-as.data.frame(temp.table)
  for (i in 1:nr){
    temp.table[i,]<-paste('l',i,'_0',sep='')
  }
  intercept<-temp.table[,1]
  
  #Generate attribute combinations
  comb.generator<-function(x.vector){
    if(length(x.vector)>1){
      temp.attr<-x.vector
      temp.attr.sav<-NULL
      for(i in 1:length(temp.attr)){
        temp.1<-combn(temp.attr,i)
        temp.2<-apply(temp.1,2,function(x){paste(x,collapse = "")})
        temp.attr.sav<-c(temp.attr.sav,temp.2)
      }
    }
    if(length(x.vector)==1){temp.attr.sav<-x.vector}
    temp.attr.sav
  }
  #vectors needed for combination.generator
  Item.load.id<-list()
  for ( i in 1:nr){
    Item.load.id[[i]]<-grep('1',Qmatrix[i,])}
  
  Attr.load.id<-list()
  attr.load.id<-matrix(0,length(temp.table.col),nc)
  for ( i in 1:length(temp.table.col)){
    attr.load.id[i,]<-unlist(strsplit(temp.table.col[i],split=''))
    Attr.load.id[[i]]<-grep('1',attr.load.id[i,])
  }
  
  #Generate Combination for both Item.load and Attr.load
  Item.Comb<-list()
  for ( i in 1:nr){
    Item.Comb[[i]]<-comb.generator(Item.load.id[[i]])
  }
  Attr.Comb<-list()
  for ( i in 2:length(temp.table.col)){
    Attr.Comb[[1]]<-0
    Attr.Comb[[i]]<-comb.generator(Attr.load.id[[i]])
  }
  constraints.list<-list()
  nway.inter.list<-list()
  for(i in 1:nr){
    for(a in 2:length(temp.table.col)){
      ifzero<-as.numeric(paste(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],collapse=''))
      if((!is.na(ifzero))){
        temp.table[i,a]<-paste(c(temp.table[i,a],
                                 paste("S","l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='',collapse='')
        ),collapse='')
        if(a==length(temp.table.col)){
          nway.inter.list[[i]]<-nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])])
          constraints.list[[i]]<-paste("l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='')
        }
      }
    }
  }
  
  #Create Lambda Table
  Lamda.Table<-temp.table
  for(i in 1:nr){
    for(a in 1:length(Lamda.Table)){
      t.ref<-unique(as.character(Lamda.Table[i,]))
      pos<-c(1:length(t.ref))[Lamda.Table[i,a]==t.ref]
      temp.table[i,a]<-paste("t",i,"_",pos,sep='')}}
  
  #Generate LCDM specification
  out<-list()
  out[[1]]<-Lamda.Table
  out[[2]]<-temp.table
  out[[3]]<-constraints.list
  out[[4]]<-nway.inter.list
  out[[5]]<-intercept
  OUTPUT<-out
  nclass<-ncol(OUTPUT[[1]]);Nc<-nclass
   #Produce kernel expressions across items and attributes
  Kernel.exp<-OUTPUT[[1]]
  Kernel.exp.detect<-OUTPUT[[1]] #052719updates
  Kernel.exp.LCDM<-OUTPUT[[1]] #052719updates
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])
      Kernel.exp.detect[i,j]<-NA} #052719updates
    }
  }
 
  Kernel.exp.LCDM<-Kernel.exp
  
  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                        nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))
  
  
  if(length(name.inter)!=0){
    for (inter in 1: length(name.inter)){
      temp.nw<-numway.inter[inter]
      temp.nm<-name.inter[inter]
      temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
      temp.sel<-paste(unlist(strsplit(temp.nm,split = '_'))[1],"_",(1:(temp.nw-1)),sep='')
      first.sel<-unlist(OUTPUT[[3]])[grep(paste((temp.sel),collapse="|"),unlist(OUTPUT[[3]]))]
      second.sel<-sub(".*_.", "", first.sel)
      for (sel in 1:length(temp.subnm)){
        SEL<-second.sel[sel]
        Constrain.List1<-rbind(
          paste(temp.nm,">-(0", paste("+",first.sel[grep(SEL,second.sel)],
                                      sep='',collapse=''),")",sep=''),Constrain.List1)
      }
    }
    Constrain.List1<-as.character(Constrain.List1)
  }else{
    Constrain.List1<-NULL
  }
  
  itemParmName<-c(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==2],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==3],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==4],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==5],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==6],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==7],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==8],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==9],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==10],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==11],OUTPUT[[5]])
  numMainEffect<-length(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1])
  Constrain.List<-paste('  real <lower=0>',itemParmName[1:numMainEffect],';\n ')
  Unconstrain.List<-paste('  real',itemParmName[-(1:numMainEffect)],';\n ')
  Reparm<-as.data.frame(matrix(0,nr,nclass))
  
  
  
  
  #############################################################
  #######061719update: Attribute Match Table###################
  #############################################################
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  #for main effect
  Attr4Main<-cbind(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],
                   paste("A",substrRight(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],1),sep=''))
  Attr4Main<-cbind(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],
                   apply(Attr4Main,1,function(x){paste(x,collapse = "*")}))
  #for interaction effect
  Attr4Interaction<-NULL
  if(length(name.inter)!=0){
    for (inter in 1: length(name.inter)){
      temp.nw<-numway.inter[inter]
      temp.nm<-name.inter[inter]
      temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
      Attr4Interaction<-rbind(Attr4Interaction,c(temp.nm,
                                                 paste(temp.nm,paste(paste("A",temp.subnm,sep=""),collapse = '*'),sep= '*')))
    }}else{
      Attr4Interaction<-NULL
    }
  
  replaceWithAttribute<-rbind(Attr4Main,Attr4Interaction)
  replaceWithAttributeLoop<-replaceWithAttribute
  replaceWithAttributeLoop[,1]<-paste(replaceWithAttribute[,1],"[itero]",sep='')
  replaceAttributeLoop<-cbind(paste("A",1:nc,sep=''),
                              paste(paste("A",1:nc,sep=''),"[iterp,itero]",sep=''))
  
  for(j in 1:nrow(replaceWithAttributeLoop)){
    replaceWithAttributeLoop[j,2]<-str_replace_all(replaceWithAttributeLoop[j,2],replaceWithAttribute[j,1],replaceWithAttributeLoop[j,1])
  }
  for(j in 1:nrow(replaceWithAttributeLoop)){
    for(z in 1:nc){
      replaceWithAttributeLoop[j,2]<-str_replace_all(replaceWithAttributeLoop[j,2],replaceAttributeLoop[z,1],replaceAttributeLoop[z,2])
      
    }
  }
  
  befAttributeKernel<-Kernel.exp.LCDM[,Nc]
  befAttributeKernel<-paste(befAttributeKernel,"[itero]",sep='')
  befAttributeKernel<-str_replace_all(befAttributeKernel,"\\+","[itero]+")
  ############################################################
  ############################################################
  ####0617update: move after Nt####### 
  
  
  #############################################################
  #######060719update: Multiple Group##########################
  #############################################################
  aftAttributeKernel<-befAttributeKernel
  for(j in 1:length(befAttributeKernel)){
    for(z in 1:nrow(replaceWithAttribute)){
      aftAttributeKernel[j]<-str_replace_all(aftAttributeKernel[j],fixed(replaceWithAttributeLoop[z,1]),replaceWithAttributeLoop[z,2])
    }
  }
  
  for(j in 1:nr){
    aftAttributeKernel[j]<-paste('Eta[',j,",itero][iterp]=",aftAttributeKernel[j],";",sep='')
  }
  
  aftAttributeKernel<-str_replace_all(aftAttributeKernel,"A","Alpha")
  ############################################################
  #061819update:Parameter definition
  ############################################################
  itemParmStan<-c(paste('real <lower=0> ', Attr4Main[,1],'[No];',sep=''),
                  paste('real ', intercept,'[No];',sep=''),
                  paste('real ', Attr4Interaction[,1],'[No];',sep=''),
                  'real Beta0[Na];',
                  'real Beta1[Na];')
  #cat(paste(itemParmStan,"\n"))
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',c(itemParmName,"Beta0" ,"Beta1" ),'~normal(0,5)',';\n',sep='')))
  
  
  
  
  data.block<-
    'data{
  
  int Np;
  int Nr;
  int Nz;
  int No;
  int Na;
  int Ni;
  int TimeVec[No];
  int OccasionID[Nr];
  int SubjectID[Nr];
  int Y[Nr, Ni];
}'

  parameter.block<-paste("\nparameters{\n
                         matrix[Np,Nz] Zeta;
                         matrix[Np,No] Theta;",
                         '\n',
                         paste(defineAThreshold,collapse='\n'),
                         '\n',
                         paste(itemParmStan,collapse='\n'),
                         '\n',
                         '\ncov_matrix[Nz] Sigma;',
                         '\nvector[Nz] Mu;',
                         #R commend
                         if(variance.equal){
                           '\nreal <lower=0,upper=100> ResVec;'
                         }else{
                           '\nreal <lower=0,upper=100> ResMat[No];'
                         },"}")
  
  transformedParm.block<-'
  transformed parameters{
  vector[Nz] zeros;
  matrix[Nz,Nz] identity_Sigma;
  matrix[Nz,Nz] identity_Mu;
  
  zeros = rep_vector(0, Nz);
  identity_Sigma =diag_matrix(rep_vector(1.0,Nz)); 
  identity_Mu =diag_matrix(rep_vector(10,Nz)); 
  }
  '  
  ######################################################################################
  
  ######################################################################################
  if(variance.equal){
    ResPrior<-'\nResVec~exponential(.1);'
  }else{
    ResPrior<-'\nResMat~exponential(.1);'
  }
  
  AnAlpha_model<-c(paste(paste("A",1:Na,sep=""),  "[iterp,itero]= inv_logit(Beta0[1]+Beta1[1]*Theta[iterp,itero]) ;\n",sep='')
                   ,paste(paste("Alpha",1:Na,sep=""),  "[iterp,itero]=ceil(fdim(A1[iterp,itero],AThres1[itero]));\n",sep=''))
  
  
  
  
  model.block<-paste(
    "model{\n",
    paste(defineAlpha,collapse='\n'),
    "\n",
    paste(defineA,collapse='\n'),
    "\n",
    "vector[Np] Eta[Ni,No];",
    "\n",
    paste(Parmprior,collapse =''),
    "Sigma~inv_wishart(Nz,identity_Sigma);
    Mu~multi_normal(zeros,identity_Mu);",
    '\n',
    paste(AThresholdPrior,collapse ='\n'),
    ResPrior,
    '\n',
    "for (iterp in 1:Np){
    Zeta[iterp,] ~ multi_normal(Mu, Sigma);
    }",
"\n",


"for (iterp in 1:Np){
for(itero in 1:No){\n
",
paste("  ",AnAlpha_model,collapse ='  '),
"  }\n
}\n
",
"for (iterp in 1:Np){\n
for( itero in 1:No){\n ",
paste(aftAttributeKernel,collapse ='\n'),
" \n  }\n
}\n ",

if(!quad.structure){
  if(variance.equal){
    'for (iterp in 1:Np){
    for(itero in 1:No){  
    Theta[iterp,itero]~normal(Zeta[iterp,1]+Zeta[iterp,2]*TimeVec[itero],ResVec);
    }
  }'
  }else{
    'for (iterp in 1:Np){
    for(itero in 1:No){
    Theta[iterp,itero]~normal(Zeta[iterp,1]+Zeta[iterp,2]*TimeVec[itero],ResMat[itero]);
    }
}'
    
    }
}else{
  if(variance.equal){
    'for (iterp in 1:Np){
    for(itero in 1:No){  
    Theta[iterp,itero]~normal(Zeta[iterp,1]+Zeta[iterp,2]*TimeVec[itero]+Zeta[iterp,3]*TimeQuadVec[itero],ResVec);
    }
  }'
  }else{
    'for (iterp in 1:Np){
    for(itero in 1:No){  
    Theta[iterp,itero]~normal(Zeta[iterp,1]+Zeta[iterp,2]*TimeVec[itero]+Zeta[iterp,3]*TimeQuadVec[itero],ResMat[itero]);
    }
  }'
  }
  
},
"\n",
"
for (iterr in 1:Nr){
for (iteri in 1:Ni){
Y[iterr,iteri] ~ bernoulli_logit(Eta[iteri,OccasionID[iterr]][SubjectID[iterr]]);
}
}
}
",sep='')
  ###########################################################################
  ####################For posterior#############################################
  if(variance.equal){
    ResPosterior<-'\nreal <lower=0,upper=100> ResVec;'
    ResPosteriorSampling<-'\nResVec<-exponential_rng(.1);'
  }else{
    ResPosterior<-'\nreal <lower=0,upper=100> ResMat[No];'
    ResPosteriorSampling<-'\nResMat<-exponential_rng(.1);'
  }
  
  if(!quad.structure){
    if(variance.equal){
      ThetaPosteriorSampling<-'for (iterp in 1:Np){
      for(itero in 1:No){  
      Theta[iterp,itero]<-normal_rng(Zeta[iterp,1]+Zeta[iterp,2]*TimeVec[itero],ResVec);
      }
    }'
  }else{
    ThetaPosteriorSampling<-'for (iterp in 1:Np){
    for(itero in 1:No){
    Theta[iterp,itero]<-normal_rng(Zeta[iterp,1]+Zeta[iterp,2]*TimeVec[itero],ResMat[itero]);
    }
    }'
    
    }
  }else{
    if(variance.equal){
      ThetaPosteriorSampling<-'for (iterp in 1:Np){
      for(itero in 1:No){  
      Theta[iterp,itero]<-normal_rng(Zeta[iterp,1]+Zeta[iterp,2]*TimeVec[itero]+Zeta[iterp,3]*TimeQuadVec[itero],ResVec);
      }
    }'
  }else{
    ThetaPosteriorSampling<-'for (iterp in 1:Np){
    for(itero in 1:No){  
    Theta[iterp,itero]<-normal_rng(Zeta[iterp,1]+Zeta[iterp,2]*TimeVec[itero]+Zeta[iterp,3]*TimeQuadVec[itero],ResMat[itero]);
    }
    }'
  }}
  
  #####################End: For posterior###################################################
  ##########################################################################################
  
  generatedQuant.block <- paste("\n
                                generated quantities {",
                                "\n",
                                paste(defineAlpha,collapse='\n'),
                                '\n',
                                paste(defineA,collapse='\n'),
                                "\n",
                                "vector[Np] Eta[Ni,No];
                                real log_lik[Nr,Ni];
                                ",
                                "\n",
                                
                                "for (iterp in 1:Np){
                                for(itero in 1:No){\n
                                ",
                                paste("  ",AnAlpha_model,collapse ='  '),
                                "  }\n
                                }\n
                                ",                               
                                '\n',
                                "for (iterp in 1:Np){\n
                                for( itero in 1:No){\n ",
                                paste(aftAttributeKernel,collapse ='\n'),
                                " \n  }\n
                                }\n ",

                                "for (iterr in 1:Nr){
                                for (iteri in 1:Ni){
                                log_lik[iterr,iteri] = bernoulli_log(Y[iterr,iteri],inv_logit(Eta[iteri,OccasionID[iterr]][SubjectID[iterr]]));
                                }
                                }
                                }",sep="")
  
  if (.Platform$OS.type == "unix") {
    filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
  }else{
    filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
  }
  
  sink(file= filename,append=FALSE)
  cat(
    paste(c('   ',
            data.block,parameter.block,transformedParm.block,model.block,generatedQuant.block)
    ))
  sink(NULL)
  
  }

#SEPERATION#
#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanLCDM_mG.run<-function(Qmatrix,
                          response.matrix,
                          GroupID,
                          fixeditem.vector=NA,
                          class.equal=T,
                          script.path=NA,save.path=getwd(),save.name="LCDM_uninf_multiG",
                          iter=1000,warmup = 0,
                          chain.num=3,init.list='random',control.list=NA){
  group.num<-length(unique(GroupID))
  rstan.detect<-tryCatch(library("rstan"),error=function(e){"rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details."})
  if(length(rstan.detect)==1){
    stop()
  }
  Cdm.init<-F
  if(init.list=='cdm'){
    Cdm.init<-T
    Install.package(c("CDM","stringr"))
    trueParmName<-Parm.name(Qmatrix=Qmatrix)$parm.name
    Classp.exp1<-Parm.name(Qmatrix=Qmatrix)$class.expression
    mod1<-gLCDM( data =response.matrix, q.matrix = Qmatrix , maxit=700,link = "logit",progress=F)
    CDMresult<-as.data.frame(coef(mod1))
    library(stringr)
    CDM.parm.name<-paste(paste(paste('l',CDMresult[,3],sep=''),'_',sep=''),str_count(CDMresult$partype.attr,"Attr"),sep='')
    CDM.parm.name<-paste(CDM.parm.name,
                         unlist(lapply(strsplit(unlist(lapply(strsplit(CDMresult$partype.attr, 'Attr', fixed=FALSE),function(x){paste(x,collapse="")})),'-'),function(x){paste(x,collapse="")})),
                         sep='')
    CDM.parm.est<-CDMresult$est
    parm.ini<-round(CDM.parm.est[match(trueParmName,CDM.parm.name)],4)
    CDM.prop.est<-mod1$attribute.patt
    prop.ini<-CDM.prop.est[match(Classp.exp1,rownames(CDM.prop.est)),1]
    inilist1<-paste('list(',paste(noquote(paste(noquote(unlist(list(
      paste(paste('Vc=c(',paste((prop.ini),collapse=','),')',collapse=','))))
    ))),collapse=',')  ,')',collapse='')

    IniList1<-NULL
    if(!class.equal){
      temp.inilist1<-eval(parse(n =2000000 ,text=inilist1))
      eval(parse(n =2000000 ,text=paste(paste('IniList1$Vc_g',1:group.num,sep=''),"<-temp.inilist1$Vc",sep='') ))
      inilist1<-IniList1
    }else{
      inilist1<-eval(parse(n =2000000 ,text=inilist1))}


    for( i in 2:chain.num){
      temp.text<-paste('inilist',i,"<-inilist1",sep='')
      eval(parse(text=(temp.text)))
    }
    temp.text<-paste('init.list<-list(',paste(paste('inilist',1:chain.num,sep=''),collapse = ","),')',sep='')
    eval(parse(text=(temp.text)))
  }
  data.list<-Generate.datalist(Qmatrix,response.matrix,GroupID)

  if(is.na(control.list)){control.list<-list(adapt_delta=0.82)}
  if(is.na(script.path)==T){
    options(warn=-1)
    #Need to update script
    StanLCDM_mG.script(Qmatrix=Qmatrix,
                       group.num=group.num,
                       fixeditem.vector=fixeditem.vector,
                       class.equal=class.equal,
                       save.path=save.path,save.name=save.name)
    script.path<-paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
    options(warn=0)
    compiled_model<-stan_model(script.path)}
  else{
    compiled_model<-stan_model(script.path)
  }
  if(Cdm.init==T){
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})
  }else{
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})

  }

  estimated_model
}

#SEPERATION#
#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanLCDM_mG.script<-function(Qmatrix,
                             group.num,
                             fixeditem.vector=NA,
                             class.equal=T,
                             save.path=getwd(),save.name="LCDM_uninf_multiG"){

  #Load packages
  Install.package("plyr")
  Install.package('stringr')

  nc<-ncol(Qmatrix)
  nr<-nrow(Qmatrix)
  temp.table.col<-unique(apply(combn(rep(c(0,1),nc),nc),2,function(x){paste(x,collapse = "")}))
  temp.table.col<-temp.table.col[order(temp.table.col)]
  temp.table<-matrix(0,nr,length(temp.table.col))
  colnames(temp.table)<-temp.table.col
  rownames(temp.table)<-paste('item',c(1:nr),sep='')
  temp.table<-as.data.frame(temp.table)
  for (i in 1:nr){
    temp.table[i,]<-paste('l',i,'_0',sep='')
  }
  intercept<-temp.table[,1]

  #Generate attribute combinations
  comb.generator<-function(x.vector){
    if(length(x.vector)>1){
      temp.attr<-x.vector
      temp.attr.sav<-NULL
      for(i in 1:length(temp.attr)){
        temp.1<-combn(temp.attr,i)
        temp.2<-apply(temp.1,2,function(x){paste(x,collapse = "")})
        temp.attr.sav<-c(temp.attr.sav,temp.2)
      }
    }
    if(length(x.vector)==1){temp.attr.sav<-x.vector}
    temp.attr.sav
  }
  #vectors needed for combination.generator
  Item.load.id<-list()
  for ( i in 1:nr){
    Item.load.id[[i]]<-grep('1',Qmatrix[i,])}

  Attr.load.id<-list()
  attr.load.id<-matrix(0,length(temp.table.col),nc)
  for ( i in 1:length(temp.table.col)){
    attr.load.id[i,]<-unlist(strsplit(temp.table.col[i],split=''))
    Attr.load.id[[i]]<-grep('1',attr.load.id[i,])
  }

  #Generate Combination for both Item.load and Attr.load
  Item.Comb<-list()
  for ( i in 1:nr){
    Item.Comb[[i]]<-comb.generator(Item.load.id[[i]])
  }
  Attr.Comb<-list()
  for ( i in 2:length(temp.table.col)){
    Attr.Comb[[1]]<-0
    Attr.Comb[[i]]<-comb.generator(Attr.load.id[[i]])
  }
  constraints.list<-list()
  nway.inter.list<-list()
  for(i in 1:nr){
    for(a in 2:length(temp.table.col)){
      ifzero<-as.numeric(paste(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],collapse=''))
      if((!is.na(ifzero))){
        temp.table[i,a]<-paste(c(temp.table[i,a],
                                 paste("S","l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='',collapse='')
        ),collapse='')
        if(a==length(temp.table.col)){
          nway.inter.list[[i]]<-nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])])
          constraints.list[[i]]<-paste("l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='')
        }
      }
    }
  }

  #Create Lambda Table
  Lamda.Table<-temp.table
  for(i in 1:nr){
    for(a in 1:length(Lamda.Table)){
      t.ref<-unique(as.character(Lamda.Table[i,]))
      pos<-c(1:length(t.ref))[Lamda.Table[i,a]==t.ref]
      temp.table[i,a]<-paste("t",i,"_",pos,sep='')}}

  #Generate LCDM specification
  out<-list()
  out[[1]]<-Lamda.Table
  out[[2]]<-temp.table
  out[[3]]<-constraints.list
  out[[4]]<-nway.inter.list
  out[[5]]<-intercept
  OUTPUT<-out
  nclass<-ncol(OUTPUT[[1]]);Nc<-nclass

  #Produce kernel expressions across items and attributes
  Kernel.exp<-OUTPUT[[1]]
  Kernel.exp.detect<-OUTPUT[[1]] #052719updates
  Kernel.exp.LCDM<-OUTPUT[[1]] #052719updates
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])
      Kernel.exp.detect[i,j]<-NA} #052719updates
    }
  }
  for (i in 1:nrow(OUTPUT[[1]])){ #052719updates
    theClosestEffect<-which(is.na(Kernel.exp.detect[i,]))[1] #052719updates
    useToReplaceLonger<-Kernel.exp[i,theClosestEffect] #052719updates
    Kernel.exp.LCDM[i,is.na(Kernel.exp.detect[i,])]<-useToReplaceLonger #052719updates
  } #052719updates

  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                        nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))


  if(length(name.inter)!=0){
    for (inter in 1: length(name.inter)){
      temp.nw<-numway.inter[inter]
      temp.nm<-name.inter[inter]
      temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
      temp.sel<-paste(unlist(strsplit(temp.nm,split = '_'))[1],"_",(1:(temp.nw-1)),sep='')
      first.sel<-unlist(OUTPUT[[3]])[grep(paste((temp.sel),collapse="|"),unlist(OUTPUT[[3]]))]
      second.sel<-sub(".*_.", "", first.sel)
      for (sel in 1:length(temp.subnm)){
        SEL<-second.sel[sel]
        Constrain.List1<-rbind(
          paste(temp.nm,">-(0", paste("+",first.sel[grep(SEL,second.sel)],
                                      sep='',collapse=''),")",sep=''),Constrain.List1)
      }
    }
    Constrain.List1<-as.character(Constrain.List1)
  }else{
    Constrain.List1<-NULL
  }

  itemParmName<-c(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==2],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==3],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==4],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==5],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==6],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==7],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==8],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==9],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==10],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==11],OUTPUT[[5]])
  numMainEffect<-length(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1])
  Constrain.List<-paste('  real<lower=0>',itemParmName[1:numMainEffect],';\n ')
  Unconstrain.List<-paste('  real',itemParmName[-(1:numMainEffect)],';\n ')
  Reparm<-as.data.frame(matrix(0,nr,nclass))

  #############################################################
  #######060319update: Multiple Group##########################
  #############################################################
  if(!is.na(fixeditem.vector)[1]){
    fixedItem.vector<-c(1:nr)[-fixeditem.vector]
  }else{fixedItem.vector<-c(1:nr)}

  #############################################################
  #######060719update: Multiple Group##########################
  #############################################################
  Kernel.exp.LCDM<-Kernel.exp
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      Reparm[loopi,loopc]<-paste('  PImat[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp[loopi,loopc]),');\n',sep='')
    }
  }

  #############################################################
  #######060719update: Multiple Group End######################
  #############################################################

  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))  #############################################################


  Kernel.exp.LCDM.groupName<-paste("Kernel.exp.LCDM_g",c(1:group.num),sep='')
  for(i in 1:group.num){
    tempfill.Kernel.exp.LCDM<-Kernel.exp.LCDM
    temp.Kernel.exp.LCDM<-Kernel.exp.LCDM[fixedItem.vector,]
    for(j in 1:nrow(temp.Kernel.exp.LCDM)){
      for(z in 1:ncol(temp.Kernel.exp.LCDM)){
        temp.Kernel.exp.LCDM[j,z]<-paste(temp.Kernel.exp.LCDM[j,z],'_g',i,sep='')
      }
    }
    for(j in 1:nrow(temp.Kernel.exp.LCDM)){
      for(z in 1:ncol(temp.Kernel.exp.LCDM)){
        temp.Kernel.exp.LCDM[j,z]<-str_replace_all(temp.Kernel.exp.LCDM[j,z],"\\+",paste("_g",i,"+",sep=''))
      }
    }
    tempfill.Kernel.exp.LCDM[fixedItem.vector,]<-temp.Kernel.exp.LCDM
    assign(Kernel.exp.LCDM.groupName[i],tempfill.Kernel.exp.LCDM)
  }
  Kernel.exp.LCDM.list<-list()
  for(i in 1:group.num){Kernel.exp.LCDM.list[[i]]<-eval(parse(text=paste("Kernel.exp.LCDM_g",i,sep='')))}
  #############################################################
  ##########060719update: Multiple Group End###################
  #############################################################

  #############################################################
  ##########060719update: Multiple Group#######################
  #############################################################

  PImat.groupName<-paste("PImat_g",c(1:group.num),sep='')
  Reparm.multigroup<-array(0,dim = c(nr,nclass,group.num))
  #Produce Stan code for PImat parameter
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      for (loopg in 1:group.num){
        Reparm.multigroup[loopi,loopc,loopg]<-paste('  PImat_g',loopg,'[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp.LCDM.list[[loopg]][loopi,loopc]),');\n',sep='')

      }
    }
  }

  #############################################################
  ##########060719update: Multiple Group ######################
  #############################################################

  intercept.multigroup<-array(0,dim = c(nr,1,group.num))
  mainEff.multigroup<-array(0,dim = c(numMainEffect,1,group.num))
  interaction.multigroup<-array(0,dim = c(length(name.inter),1,group.num))
  #Group Invariant Parameter Name
  fixedParmName<-NULL
  if(!is.na(fixeditem.vector)[1]){
    for(i in fixeditem.vector){
      fixedParmName<-c(fixedParmName,out[[3]][[i]])
    }
    fixedParmName<-c(fixedParmName,out[[5]][fixeditem.vector])
  }
  #Group Variant Parameter Name
  freeParmName<-c(out[[5]],itemParmName)[!c(out[[5]],itemParmName)%in%fixedParmName]

  for(i in 1:group.num){
    tempfill.intercept<-out[[5]]
    tempfill.mainEff<-itemParmName[1:numMainEffect]
    tempfill.interaction<-name.inter


    temp.intercept<-tempfill.intercept[fixedItem.vector]
    temp.mainEff<-tempfill.mainEff[!tempfill.mainEff%in%fixedParmName]
    temp.interaction<-tempfill.interaction[!tempfill.interaction%in%fixedParmName]

    for(j in 1:length(temp.intercept)){
      temp.intercept[j]<-paste(temp.intercept[j],"_g",i,sep='')
    }

    for(j in 1:length(temp.mainEff)){
      temp.mainEff[j]<-paste(temp.mainEff[j],"_g",i,sep='')
    }

    for(j in 1:length(temp.interaction)){
      temp.interaction[j]<-paste(temp.interaction[j],"_g",i,sep='')
    }


    tempfill.intercept[fixedItem.vector]<-temp.intercept
    tempfill.mainEff[!tempfill.mainEff%in%fixedParmName]<-temp.mainEff
    tempfill.interaction[!tempfill.interaction%in%fixedParmName]<-temp.interaction

    intercept.multigroup[,,i]<-tempfill.intercept
    mainEff.multigroup[,,i]<-tempfill.mainEff
    interaction.multigroup[,,i]<-tempfill.interaction
  }

  Constrain.List<-paste('  real<lower=0>',unique(mainEff.multigroup),';\n ')
  Unconstrain.List<-paste('  real',unique(c(intercept.multigroup,interaction.multigroup)),';\n ')
  #############################################################
  ##########060319update: Multiple Group End###################
  #############################################################
  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))
  update.Parmprior<-Parmprior
  fix.Parmprior<-NULL

  update.Parmprior<-update.Parmprior[update.Parmprior!='']

  #############################################################
  #####060719update:Change from update.Parmprior&fix     ######
  #############################################################
  update.Parmprior.multiGroup<-NULL
  for(i in 1:length(update.Parmprior)){
    for (j in 1:length(freeParmName)){
      for (z in 1:group.num){
        if(sum(grepl(freeParmName[j],update.Parmprior[i]))>=1){
          temp.update.Parmprior<-str_replace_all(update.Parmprior[i],freeParmName[j],paste(freeParmName[j],"_g",z,sep=''))
          update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,
                                         temp.update.Parmprior)
        }
      }
    }
  }
  update.Parmprior.multiGroup<-unique(update.Parmprior.multiGroup)
  update.Parmprior.multiGroup<-c("   //Prior\n",paste('   ',fixedParmName,'~normal(0,5)',';\n',sep=''),update.Parmprior.multiGroup )
  if(class.equal){
    update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='') )
  }else{
    for(i in 1:group.num){
      update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,
                                     paste('   Vc_g',i,'~dirichlet(rep_vector(2.0, Nc));\n',sep='') )
    }
  }

  ##therefore we can use: fix.Parmprior,update.Parmprior
  #############################################################
  #############################################################

  #############################################################
  #####060419update:Likelihood Add PImat_g#####################
  #############################################################
  PImat.likelihood1<-NULL
  PImat.likelihood0<-NULL
  Vc.likelihood<-NULL
  for(loopg in 1:group.num){
    temp.PImat.likelihood1<-paste(paste('          if (GroupID[iterp]==',loopg,')',sep=''),
                                  paste('            contributionsI[iteri]=bernoulli_lpmf(1|PImat_g',loopg,'[iteri,iterc]);',sep=''),
                                  sep='\n')
    PImat.likelihood1<-paste(PImat.likelihood1,temp.PImat.likelihood1,sep='\n')
    temp.PImat.likelihood0<-paste(paste('          if (GroupID[iterp]==',loopg,')',sep=''),
                                  paste('            contributionsI[iteri]=bernoulli_lpmf(0|PImat_g',loopg,'[iteri,iterc]);',sep=''),
                                  sep='\n')
    PImat.likelihood0<-paste(PImat.likelihood0,temp.PImat.likelihood0,sep='\n')
  }
  if(!class.equal){
    for(loopg in 1:group.num){
      temp.Vc.likelihood<-paste(paste('       if (GroupID[iterp]==',loopg,')',sep=''),
                                paste('         contributionsC[iterc]=log(Vc_g',loopg,'[iterc])+sum(contributionsI);',sep=''),
                                sep='\n')
      Vc.likelihood<-paste(Vc.likelihood,temp.Vc.likelihood,sep='\n')
    }
  }

  #Likelihood Stan code
  if(class.equal){
    Likelihood<-paste('
                      \n
                      //Likelihood
                      for (iterp in 1:Np){
                      for (iterc in 1:Nc){
                      for (iteri in 1:Ni){
                      if (Y[iterp,iteri] == 1)'
                      ,PImat.likelihood1,'\n',
                      '        else'
                      ,PImat.likelihood0,
                      '}
                      contributionsC[iterc]=log(Vc[iterc])+sum(contributionsI);
                      }
                      target+=log_sum_exp(contributionsC);
                      }
                      ',sep='')}else{
                        Likelihood<-paste('
                                          \n
                                          //Likelihood
                                          for (iterp in 1:Np){
                                          for (iterc in 1:Nc){
                                          for (iteri in 1:Ni){
                                          if (Y[iterp,iteri] == 1)'
                                          ,PImat.likelihood1,'\n',
                                          '        else'
                                          ,PImat.likelihood0,
                                          '}\n',
                                          Vc.likelihood
                                          ,'


                                          }
                                          target+=log_sum_exp(contributionsC);
                                          }
                                          ',sep='')

                      }
  #############################################################
  #####060419update:Likelihood Add PImat_g  end################
  #############################################################

  #Data Specification
  data.spec<-'
  data{
  int Np;
  int Ni;
  int Nc;
  matrix[Np, Ni] Y;
  vector[Np] GroupID;
  }
  '

  Constrain.List<-unique(Constrain.List);Unconstrain.List<-unique(Unconstrain.List)
  #############################################################
  #####060419update: Stan script##############################
  #############################################################

  #Parameter Specification
  if(class.equal){parm.spec<-paste(c('
                                     parameters{
                                     simplex[Nc] Vc;\n ',paste0(Constrain.List),paste0(Unconstrain.List),
                                     '}\n'),collapse='')}else{
                                       parm.spec<-paste(c('
                                                          parameters{\n ',
                                                          paste(paste('   simplex[Nc] Vc_g',1:group.num, ";",sep=''),"\n"),
                                                          paste0(Constrain.List),paste0(Unconstrain.List),
                                                          '}\n'),collapse='')


                                     }

  #Reparameter Specification
  transparm.spec<-paste(c('
                          transformed parameters{

                          ',
                          paste('  matrix[Ni, Nc] PImat_g',1:group.num,';\n',sep=''),

                          paste0(c(unique(Reparm.multigroup) )),'}\n'),collapse='')


  update.Parmprior.multiGroup<-update.Parmprior.multiGroup[!startsWith(str_remove_all(update.Parmprior.multiGroup," "),"~")]
  #Model Specification update052619
  model.spec<-paste(c('\nmodel {\n',paste(c(Modelcontainer,update.Parmprior.multiGroup,Likelihood),sep=''),'\n}',sep=''))
  model.spec<-model.spec[!startsWith(str_remove_all(model.spec," "),"~")]



  #Generated Quantities Specification
  IC.generatedquantities<-NULL
  if(!class.equal){
    for(loopg in 1:group.num){
      temp.IC.generatedquantities<-paste(paste('       if (GroupID[iterp]==',loopg,')',sep=''),
                                         paste('         contributionsIC[iteri,iterc]=log(Vc_g',loopg,'[iterc])+contributionsI[iteri];',sep=''),
                                         sep='\n')
      IC.generatedquantities<-paste(IC.generatedquantities,temp.IC.generatedquantities,sep='\n')
    }
  }


  if(class.equal){
    generatedQuantities.spec<-paste('
  \n
  generated quantities {
  vector[Ni] log_lik[Np];
  vector[Ni] contributionsI;
  matrix[Ni,Nc] contributionsIC;

  matrix[Ni,Nc] posteriorIC;
  matrix[Np,Nc] posteriorPC;


  //Posterior
  for (iterp in 1:Np){
    for (iteri in 1:Ni){
      for (iterc in 1:Nc){
        if (Y[iterp,iteri] == 1)'
                                    ,PImat.likelihood1,'\n',
                                    '        else'
                                    ,PImat.likelihood0,
                                    '
          contributionsIC[iteri,iterc]=log(Vc[iterc])+contributionsI[iteri];
          posteriorIC[iteri,iterc]=contributionsI[iteri];
        }
      log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
    }
    for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
  }
  }
  ',sep='')}else{
    generatedQuantities.spec <- paste( '
    \n
    generated quantities {
    vector[Ni] log_lik[Np];
    vector[Ni] contributionsI;
    matrix[Ni,Nc] contributionsIC;

    matrix[Ni,Nc] posteriorIC;
    matrix[Np,Nc] posteriorPC;


    //Posterior
    for (iterp in 1:Np){
      for (iteri in 1:Ni){
        for (iterc in 1:Nc){
          if (Y[iterp,iteri] == 1)',
                                       PImat.likelihood1,'\n',
                                       '        else',
                                       PImat.likelihood0,'\n',
                                       "   ",IC.generatedquantities,'\n',
                             '          posteriorIC[iteri,iterc]=contributionsI[iteri];','\n
      }\n',
                                       '     log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
     }
     for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
   }
   }
      ',
                                       sep = ''
    )
  }





  if (.Platform$OS.type == "unix") {
    filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
  }else{
    filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
  }

  sink(file=filename, append=FALSE)
  cat(
    paste(c('   ',
            data.spec,parm.spec,transparm.spec,model.spec,generatedQuantities.spec)
    ))
  sink(NULL)

}

#SEPERATION#
#' @title Generate Stan code and Run the estimation for NCRUM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanNCRUM.run<-function(Qmatrix,response.matrix,script.path=NA,save.path=getwd(),save.name="NCRUM_uninf",iter=1000,warmup = 0,
                        chain.num=3, init.list='random',control.list=NA){
  rstan.detect<-tryCatch(library("rstan"),error=function(e){"rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details."})
  if(length(rstan.detect)==1){
    stop()
  }
  Cdm.init = F
  if(init.list=='cdm'){
    Cdm.init = T
    Install.package(c("CDM","stringr"))
    trueParmName<-Parm.name(Qmatrix=Qmatrix)$parm.name
    Classp.exp1<-Parm.name(Qmatrix=Qmatrix)$class.expression
    mod1<-gdina( data =response.matrix, q.matrix = Qmatrix , maxit=700,link = "logit",progress=F)
    CDMresult<-as.data.frame(coef(mod1))
    library(stringr)
    CDM.parm.name<-paste(paste(paste('l',CDMresult[,3],sep=''),'_',sep=''),str_count(CDMresult$partype.attr,"Attr"),sep='')
    CDM.parm.name<-paste(CDM.parm.name,
                         unlist(lapply(strsplit(unlist(lapply(strsplit(CDMresult$partype.attr, 'Attr', fixed=FALSE),function(x){paste(x,collapse="")})),'-'),function(x){paste(x,collapse="")})),
                         sep='')
    CDM.parm.est<-CDMresult$est
    # remove nagetive values in initial values
    CDM.parm.est[CDM.parm.est < 0]  <-  0.01
    parm.ini<-round(CDM.parm.est[match(trueParmName,CDM.parm.name)],4)
    parm.ini[abs(parm.ini)>10] <- 10*sign(parm.ini[abs(parm.ini)>10])
    parm.ini[parm.ini==0] <- 0.05

    CDM.prop.est<-mod1$attribute.patt
    prop.ini<-CDM.prop.est[match(Classp.exp1,rownames(CDM.prop.est)),1]
    inilist1<-paste('list(',paste(noquote(paste(noquote(unlist(list(paste(trueParmName,'=',round(parm.ini,2)),
                                                                    paste(paste('Vc=c(',paste((prop.ini),collapse=','),')',collapse=','))))
    ))),collapse=',')  ,')',collapse='')

    inilist1<-eval(parse(n =2000000 ,text=inilist1))
    for( i in 2:chain.num){
      temp.text<-paste('inilist',i,"<-inilist1",sep='')
      eval(parse(text=(temp.text)))
    }
    temp.text<-paste('init.list<-list(',paste(paste('inilist',1:chain.num,sep=''),collapse = ","),')',sep='')
    eval(parse(text=(temp.text)))
  }
  data.list<-Generate.datalist(Qmatrix,response.matrix)

  if(is.na(control.list)){control.list<-list(adapt_delta=0.82)}
  if(is.na(script.path)==T){
    options(warn=-1)
    StanNCRUM.script(Qmatrix,save.path=save.path,save.name=save.name)
    if (.Platform$OS.type == "unix") {
      filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
    }else{
      filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
    }
    script.path<-filename
    options(warn=0)
    compiled_model<-stan_model(script.path)
  }else{
    compiled_model<-stan_model(script.path)
  }
  if(Cdm.init == T){
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains= chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})
  }else{
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains= chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})

  }

  estimated_model
}


#SEPERATION#
#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanNCRUM.script<-function(Qmatrix,save.path=getwd(),save.name="NCRUM_uninf"){

  #Load packages
  Install.package("plyr")
  Install.package('stringr')

  nc<-ncol(Qmatrix)
  nr<-nrow(Qmatrix)
  temp.table.col<-unique(apply(combn(rep(c(0,1),nc),nc),2,function(x){paste(x,collapse = "")}))
  temp.table.col<-temp.table.col[order(temp.table.col)]
  temp.table<-matrix(0,nr,length(temp.table.col))
  colnames(temp.table)<-temp.table.col
  rownames(temp.table)<-paste('item',c(1:nr),sep='')
  temp.table<-as.data.frame(temp.table)
  for (i in 1:nr){
    temp.table[i,]<-paste('l',i,'_0',sep='')
  }
  intercept<-temp.table[,1]

  #Generate attribute combinations
  comb.generator<-function(x.vector){
    if(length(x.vector)>1){
      temp.attr<-x.vector
      temp.attr.sav<-NULL
      for(i in 1:length(temp.attr)){
        temp.1<-combn(temp.attr,i)
        temp.2<-apply(temp.1,2,function(x){paste(x,collapse = "")})
        temp.attr.sav<-c(temp.attr.sav,temp.2)
      }
    }
    if(length(x.vector)==1){temp.attr.sav<-x.vector}
    temp.attr.sav
  }
  #vectors needed for combination.generator
  Item.load.id<-list()
  for ( i in 1:nr){
    Item.load.id[[i]]<-grep('1',Qmatrix[i,])}

  Attr.load.id<-list()
  attr.load.id<-matrix(0,length(temp.table.col),nc)
  for ( i in 1:length(temp.table.col)){
    attr.load.id[i,]<-unlist(strsplit(temp.table.col[i],split=''))
    Attr.load.id[[i]]<-grep('1',attr.load.id[i,])
  }

  #Generate Combination for both Item.load and Attr.load
  Item.Comb<-list()
  for ( i in 1:nr){
    Item.Comb[[i]]<-comb.generator(Item.load.id[[i]])
  }
  Attr.Comb<-list()
  for ( i in 2:length(temp.table.col)){
    Attr.Comb[[1]]<-0
    Attr.Comb[[i]]<-comb.generator(Attr.load.id[[i]])
  }
  constraints.list<-list()
  nway.inter.list<-list()
  for(i in 1:nr){
    for(a in 2:length(temp.table.col)){
      ifzero<-as.numeric(paste(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],collapse=''))
      if((!is.na(ifzero))){
        temp.table[i,a]<-paste(c(temp.table[i,a],
                                 paste("S","l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='',collapse='')
        ),collapse='')
        if(a==length(temp.table.col)){
          nway.inter.list[[i]]<-nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])])
          constraints.list[[i]]<-paste("l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='')
        }
      }
    }
  }

  #Create Lambda Table
  Lamda.Table<-temp.table
  for(i in 1:nr){
    for(a in 1:length(Lamda.Table)){
      t.ref<-unique(as.character(Lamda.Table[i,]))
      pos<-c(1:length(t.ref))[Lamda.Table[i,a]==t.ref]
      temp.table[i,a]<-paste("t",i,"_",pos,sep='')}}

  #Generate LCDM specification
  out<-list()
  out[[1]]<-Lamda.Table
  out[[2]]<-temp.table
  out[[3]]<-constraints.list
  out[[4]]<-nway.inter.list
  out[[5]]<-intercept
  OUTPUT<-out
  nclass<-ncol(OUTPUT[[1]]);Nc<-nclass

  #Produce kernel expressions across items and attributes
  #Produce kernel expressions across items and attributes
  Kernel.exp<-OUTPUT[[1]]
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])}
    }
  }


  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                        nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))

  if(length(name.inter)!=0){
    for (inter in 1: length(name.inter)){
      temp.nw<-numway.inter[inter]
      temp.nm<-name.inter[inter]
      temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
      temp.sel<-paste(unlist(strsplit(temp.nm,split = '_'))[1],"_",(1:(temp.nw-1)),sep='')
      first.sel<-unlist(OUTPUT[[3]])[grep(paste((temp.sel),collapse="|"),unlist(OUTPUT[[3]]))]
      second.sel<-sub(".*_.", "", first.sel)
      for (sel in 1:length(temp.subnm)){
        SEL<-second.sel[sel]
        Constrain.List1<-rbind(
          paste(temp.nm,">-(0", paste("+",first.sel[grep(SEL,second.sel)],
                                      sep='',collapse=''),")",sep=''),Constrain.List1)
      }
    }
    Constrain.List1<-as.character(Constrain.List1)
  }else{
    Constrain.List1<-NULL
  }

  itemParmName<-c(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==2],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==3],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==4],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==5],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==6],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==7],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==8],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==9],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==10],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==11],OUTPUT[[5]])
  numMainEffect<-length(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1])
  Constrain.List<-paste('  real<lower=0>',c(itemParmName[1:numMainEffect],name.inter),';\n ')
  Unconstrain.List<-paste('  real',intercept,';\n ')
  Reparm<-as.data.frame(matrix(0,nr,nclass))


  #Produce Stan code for PImat parameter
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      Reparm[loopi,loopc]<-paste('  PImat[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp[loopi,loopc]),');\n',sep='')
    }
  }

  #########053019update:create pistar and rstar parameters###########################
  pistarParm<-rep(0,nr)
  for(loopi in 1:nr){
    pistarParm[loopi]<-paste('  piParm[',loopi,']=inv_logit(',paste(Kernel.exp[loopi,nclass]),');\n',sep='')
  }

  piParm<-str_replace_all(intercept,"l","pi");piParm<-str_replace_all(piParm,"_0","")
  rParm<-c(itemParmName[1:numMainEffect],name.inter)
  rParm<-str_replace_all(rParm,"l","r")
  Kernel.exp.NCRUM<-Kernel.exp
  for(i in 1:Nc){Kernel.exp.NCRUM[,i]<-piParm}
  #for(i in 1:nr){
  #  for (j in 1:Nc){
  #   Kernel.exp.NCRUM.inv[i,j]<-str_replace_all(Kernel.exp[i,j],
  #                                             intercept[i],
  #                                            piParm[i])
  #}
  #}

  allrParm<-NULL
  for(i in 1:nc){
    allrParm<-cbind(allrParm,paste(paste("r",1:nr,"_1",sep=''),i,sep='')
    )
  }
  allrParm[!t(apply(allrParm ,1,function(x){x%in%rParm}))]<-1
  noneOneallrParm<-allrParm[allrParm!='1']
  noneOneAllrParm<-noneOneallrParm
  noneOneAllrParm<-paste(' real',noneOneAllrParm, ';\n',sep=' ')

  for(i in 1:Nc){
    for(j in 1:nr){
      Kernel.exp.NCRUM[j,i]<-paste(paste(allrParm[j,which(strsplit(colnames(Kernel.exp)[i],'')[[1]]=='0')],"*",collapse ='',sep=''),
                                   Kernel.exp.NCRUM[j,i],sep='')
    }
  }
  Kernel.exp.NCRUM[,Nc]<-str_replace_all(c(Kernel.exp.NCRUM[,Nc]), "[*]", '')

  piTransf<-Kernel.exp.NCRUM[,Nc]
  piTransfwithoutName<-piTransf
  for(loopi in 1:nr){
    piTransf[loopi]<-paste('  piParm[',loopi,']=inv_logit(',paste(Kernel.exp[loopi,Nc]),');\n',sep='')
    piTransfwithoutName[loopi]<-paste('inv_logit(',paste(Kernel.exp[loopi,Nc]),')',sep='')
  }

  rTransf<-NULL
  for(loopi in 1:nr){
    for(loopc in 1:Nc){
      if(str_count(Kernel.exp.NCRUM[loopi,loopc],"\\*")==1){
        if(strsplit(Kernel.exp.NCRUM[loopi,loopc],"\\*")[[1]][1]=='1'){next}
        temp_r<-noneOneallrParm[(grep(strsplit(Kernel.exp.NCRUM[loopi,loopc],"\\*")[[1]][1],noneOneallrParm))]
        if(temp_r%in%noneOneallrParm){
          rTransf<-c(rTransf,paste(" ",temp_r,"=",paste( paste('inv_logit(',Kernel.exp[loopi,loopc],")",sep=''),
                                                         "/",
                                                         piTransfwithoutName[loopi])
          ))
        }
      }
    }
  }
  rTransf<-paste(rTransf,');\n',sep='')


  #########053019update:create pistar and rstar parameters END#########


  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))

  #Likelihood Stan code
  Likelihood<-'
  \n
  //Likelihood
  for (iterp in 1:Np){
    for (iterc in 1:Nc){
      for (iteri in 1:Ni){
        if (Y[iterp,iteri] == 1)
          contributionsI[iteri]=bernoulli_lpmf(1|PImat[iteri,iterc]);
        else
          contributionsI[iteri]=bernoulli_lpmf(0|PImat[iteri,iterc]);
      }
      contributionsC[iterc]=log(Vc[iterc])+sum(contributionsI);
    }
    target+=log_sum_exp(contributionsC);
  }
  '


  #Data Specification
  data.spec<-'
data{
  int Np;
  int Ni;
  int Nc;
  matrix[Np, Ni] Y;
  }
  '
  #Parameter Specification
  parm.spec<-paste(c('parameters{
  simplex[Nc] Vc;\n ',paste0(Constrain.List),paste0(Unconstrain.List),
                     '}\n'),collapse='')

  #Reparameter Specification
  transparm.spec<-paste(c('
  transformed parameters{
  matrix[Ni, Nc] PImat;
  vector[Ni] piParm;\n',
                          noneOneAllrParm,
                          piTransf,#053019update
                          rTransf,#053019update
                          paste0(unlist(Reparm)),'}\n'),collapse='')

  #Model Specification update052619
  model.spec<-paste(c('\nmodel {\n',paste(c(Modelcontainer,Parmprior,Likelihood),sep=''),'\n}',sep=''))
  model.spec<-model.spec[!startsWith(str_remove_all(model.spec," "),"~")]
  #Generated Quantities Specification
  generatedQuantities.spec<-'
  \n
generated quantities {

 vector[Ni] log_lik[Np];
 vector[Ni] contributionsI;
 matrix[Ni,Nc] contributionsIC;
 
 matrix[Ni,Nc] posteriorIC;
 matrix[Np,Nc] posteriorPC;



 //Posterior
 for (iterp in 1:Np){
   for (iteri in 1:Ni){
     for (iterc in 1:Nc){
       if (Y[iterp,iteri] == 1)
          contributionsI[iteri]=bernoulli_lpmf(1|PImat[iteri,iterc]);
       else
           contributionsI[iteri]=bernoulli_lpmf(0|PImat[iteri,iterc]);
       contributionsIC[iteri,iterc]=log(Vc[iterc])+contributionsI[iteri];
       posteriorIC[iteri,iterc]=contributionsI[iteri];
      }
      log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
    }
   for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
  }
}
'

  if (.Platform$OS.type == "unix") {
    filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
  }else{
    filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
  }

  sink(file=filename, append=FALSE)
  cat(
    paste(c('   ',
            data.spec,parm.spec,transparm.spec,model.spec,generatedQuantities.spec)
    ))
  sink(NULL)

}

#SEPERATION#
#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanNCRUM_mG.run<-function(Qmatrix,
                           response.matrix,
                           GroupID,
                           fixeditem.vector=NA,
                           class.equal=T,
                           script.path=NA,save.path=getwd(),save.name="NCRUM_uninf_multiG",
                           iter=1000,warmup = 0,
                           chain.num=3,init.list='random',control.list=NA){
  group.num<-length(unique(GroupID))
  rstan.detect<-tryCatch(library("rstan"),error=function(e){"rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details."})
  if(length(rstan.detect)==1){
    stop()
  }
  Cdm.init<-F
  if(init.list=='cdm'){
    Cdm.init<-T
    Install.package(c("CDM","stringr"))
    trueParmName<-Parm.name(Qmatrix=Qmatrix)$parm.name
    Classp.exp1<-Parm.name(Qmatrix=Qmatrix)$class.expression
    mod1<-gdina( data =response.matrix, q.matrix = Qmatrix , maxit=700,link = "logit",progress=F)
    CDMresult<-as.data.frame(coef(mod1))
    library(stringr)
    CDM.parm.name<-paste(paste(paste('l',CDMresult[,3],sep=''),'_',sep=''),str_count(CDMresult$partype.attr,"Attr"),sep='')
    CDM.parm.name<-paste(CDM.parm.name,
                         unlist(lapply(strsplit(unlist(lapply(strsplit(CDMresult$partype.attr, 'Attr', fixed=FALSE),function(x){paste(x,collapse="")})),'-'),function(x){paste(x,collapse="")})),
                         sep='')
    CDM.parm.est<-CDMresult$est
    parm.ini<-round(CDM.parm.est[match(trueParmName,CDM.parm.name)],4)
    CDM.prop.est<-mod1$attribute.patt
    prop.ini<-CDM.prop.est[match(Classp.exp1,rownames(CDM.prop.est)),1]
    inilist1<-paste('list(',paste(noquote(paste(noquote(unlist(list(
      paste(paste('Vc=c(',paste((prop.ini),collapse=','),')',collapse=','))))
    ))),collapse=',')  ,')',collapse='')

    IniList1<-NULL
    if(!class.equal){
      temp.inilist1<-eval(parse(n =2000000 ,text=inilist1))
      eval(parse(n =2000000 ,text=paste(paste('IniList1$Vc_g',1:group.num,sep=''),"<-temp.inilist1$Vc",sep='') ))
      inilist1<-IniList1
    }else{
      inilist1<-eval(parse(n =2000000 ,text=inilist1))}


    for( i in 2:chain.num){
      temp.text<-paste('inilist',i,"<-inilist1",sep='')
      eval(parse(text=(temp.text)))
    }
    temp.text<-paste('init.list<-list(',paste(paste('inilist',1:chain.num,sep=''),collapse = ","),')',sep='')
    eval(parse(text=(temp.text)))
  }
  data.list<-Generate.datalist(Qmatrix,response.matrix,GroupID)

  if(is.na(control.list)){control.list<-list(adapt_delta=0.82)}
  if(is.na(script.path)==T){
    options(warn=-1)
    #Need to update script
    StanNCRUM_mG.script(Qmatrix=Qmatrix,
                        group.num=group.num,
                        fixeditem.vector=fixeditem.vector,
                        class.equal=class.equal,
                        save.path=save.path,save.name=save.name)
    script.path<-paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
    options(warn=0)
    compiled_model<-stan_model(script.path)}
  else{
    compiled_model<-stan_model(script.path)
  }
  if(Cdm.init==T){
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})
  }else{
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})

  }

  estimated_model
}



#SEPERATION#
#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanNCRUM_mG.script<-function(Qmatrix,
                              group.num,
                              fixeditem.vector=NA,
                              class.equal=T,
                              save.path=getwd(),save.name="NCRUM_uninf_multiG"){

  #Load packages
  Install.package("plyr")
  Install.package('stringr')

  nc<-ncol(Qmatrix)
  nr<-nrow(Qmatrix)
  temp.table.col<-unique(apply(combn(rep(c(0,1),nc),nc),2,function(x){paste(x,collapse = "")}))
  temp.table.col<-temp.table.col[order(temp.table.col)]
  temp.table<-matrix(0,nr,length(temp.table.col))
  colnames(temp.table)<-temp.table.col
  rownames(temp.table)<-paste('item',c(1:nr),sep='')
  temp.table<-as.data.frame(temp.table)
  for (i in 1:nr){
    temp.table[i,]<-paste('l',i,'_0',sep='')
  }
  intercept<-temp.table[,1]

  #Generate attribute combinations
  comb.generator<-function(x.vector){
    if(length(x.vector)>1){
      temp.attr<-x.vector
      temp.attr.sav<-NULL
      for(i in 1:length(temp.attr)){
        temp.1<-combn(temp.attr,i)
        temp.2<-apply(temp.1,2,function(x){paste(x,collapse = "")})
        temp.attr.sav<-c(temp.attr.sav,temp.2)
      }
    }
    if(length(x.vector)==1){temp.attr.sav<-x.vector}
    temp.attr.sav
  }
  #vectors needed for combination.generator
  Item.load.id<-list()
  for ( i in 1:nr){
    Item.load.id[[i]]<-grep('1',Qmatrix[i,])}

  Attr.load.id<-list()
  attr.load.id<-matrix(0,length(temp.table.col),nc)
  for ( i in 1:length(temp.table.col)){
    attr.load.id[i,]<-unlist(strsplit(temp.table.col[i],split=''))
    Attr.load.id[[i]]<-grep('1',attr.load.id[i,])
  }

  #Generate Combination for both Item.load and Attr.load
  Item.Comb<-list()
  for ( i in 1:nr){
    Item.Comb[[i]]<-comb.generator(Item.load.id[[i]])
  }
  Attr.Comb<-list()
  for ( i in 2:length(temp.table.col)){
    Attr.Comb[[1]]<-0
    Attr.Comb[[i]]<-comb.generator(Attr.load.id[[i]])
  }
  constraints.list<-list()
  nway.inter.list<-list()
  for(i in 1:nr){
    for(a in 2:length(temp.table.col)){
      ifzero<-as.numeric(paste(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],collapse=''))
      if((!is.na(ifzero))){
        temp.table[i,a]<-paste(c(temp.table[i,a],
                                 paste("S","l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='',collapse='')
        ),collapse='')
        if(a==length(temp.table.col)){
          nway.inter.list[[i]]<-nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])])
          constraints.list[[i]]<-paste("l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='')
        }
      }
    }
  }

  #Create Lambda Table
  Lamda.Table<-temp.table
  for(i in 1:nr){
    for(a in 1:length(Lamda.Table)){
      t.ref<-unique(as.character(Lamda.Table[i,]))
      pos<-c(1:length(t.ref))[Lamda.Table[i,a]==t.ref]
      temp.table[i,a]<-paste("t",i,"_",pos,sep='')}}

  #Generate LCDM specification
  out<-list()
  out[[1]]<-Lamda.Table
  out[[2]]<-temp.table
  out[[3]]<-constraints.list
  out[[4]]<-nway.inter.list
  out[[5]]<-intercept
  OUTPUT<-out
  nclass<-ncol(OUTPUT[[1]]);Nc<-nclass

  #Produce kernel expressions across items and attributes
  Kernel.exp<-OUTPUT[[1]]
  Kernel.exp.detect<-OUTPUT[[1]] #052719updates
  Kernel.exp.NCRUM<-OUTPUT[[1]] #052719updates
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])
      Kernel.exp.detect[i,j]<-NA} #052719updates
    }
  }
  for (i in 1:nrow(OUTPUT[[1]])){ #052719updates
    theClosestEffect<-which(is.na(Kernel.exp.detect[i,]))[1] #052719updates
    useToReplaceLonger<-Kernel.exp[i,theClosestEffect] #052719updates
    Kernel.exp.NCRUM[i,is.na(Kernel.exp.detect[i,])]<-useToReplaceLonger #052719updates
  } #052719updates

  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                        nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))

  if(length(name.inter)!=0){
    for (inter in 1: length(name.inter)){
      temp.nw<-numway.inter[inter]
      temp.nm<-name.inter[inter]
      temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
      temp.sel<-paste(unlist(strsplit(temp.nm,split = '_'))[1],"_",(1:(temp.nw-1)),sep='')
      first.sel<-unlist(OUTPUT[[3]])[grep(paste((temp.sel),collapse="|"),unlist(OUTPUT[[3]]))]
      second.sel<-sub(".*_.", "", first.sel)
      for (sel in 1:length(temp.subnm)){
        SEL<-second.sel[sel]
        Constrain.List1<-rbind(
          paste(temp.nm,">-(0", paste("+",first.sel[grep(SEL,second.sel)],
                                      sep='',collapse=''),")",sep=''),Constrain.List1)
      }
    }
    Constrain.List1<-as.character(Constrain.List1)
  }else{
    Constrain.List1<-NULL
  }

  itemParmName<-c(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==2],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==3],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==4],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==5],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==6],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==7],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==8],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==9],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==10],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==11],OUTPUT[[5]])
  numMainEffect<-length(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1])
  Constrain.List<-paste('  real<lower=0>',itemParmName[1:numMainEffect],';\n ')
  Unconstrain.List<-paste('  real',itemParmName[-(1:numMainEffect)],';\n ')
  Reparm<-as.data.frame(matrix(0,nr,nclass))


  #########053019update:create pistar and rstar parameters###########################
  pistarParm<-rep(0,nr)
  for(loopi in 1:nr){
    pistarParm[loopi]<-paste('  piParm[',loopi,']=inv_logit(',paste(Kernel.exp[loopi,nclass]),');\n',sep='')
  }

  piParm<-str_replace_all(intercept,"l","pi");piParm<-str_replace_all(piParm,"_0","")
  rParm<-c(itemParmName[1:numMainEffect],name.inter)
  rParm<-str_replace_all(rParm,"l","r")
  Kernel.exp.NCRUM<-Kernel.exp
  for(i in 1:Nc){Kernel.exp.NCRUM[,i]<-piParm}
  #for(i in 1:nr){
  #  for (j in 1:Nc){
  #   Kernel.exp.NCRUM.inv[i,j]<-str_replace_all(Kernel.exp[i,j],
  #                                             intercept[i],
  #                                            piParm[i])
  #}
  #}
  #############################################################
  ###########060719update:The rParm update#####################

  allrParm<-NULL
  for(i in 1:nc){
    allrParm<-cbind(allrParm,paste(paste("r",1:nr,"_1",sep=''),i,sep='')
    )
  }
  allrParm[!t(apply(allrParm ,1,function(x){x%in%rParm}))]<-1
  noneOneallrParm<-allrParm[allrParm!='1']
  noneOneAllrParm<-noneOneallrParm
  noneOneAllrParm<-paste(' real',noneOneAllrParm, ';\n',sep=' ')

  for(i in 1:Nc){
    for(j in 1:nr){
      Kernel.exp.NCRUM[j,i]<-paste(paste(allrParm[j,which(strsplit(colnames(Kernel.exp)[i],'')[[1]]=='0')],"*",collapse ='',sep=''),
                                   Kernel.exp.NCRUM[j,i],sep='')
    }
  }
  Kernel.exp.NCRUM[,Nc]<-str_replace_all(c(Kernel.exp.NCRUM[,Nc]), "[*]", '')

  piTransf<-Kernel.exp.NCRUM[,Nc]
  piTransfwithoutName<-piTransf
  for(loopi in 1:nr){
    piTransf[loopi]<-paste('  piParm[',loopi,']=inv_logit(',paste(Kernel.exp[loopi,Nc]),');\n',sep='')
    piTransfwithoutName[loopi]<-paste('inv_logit(',paste(Kernel.exp[loopi,Nc]),')',sep='')
  }

  rTransf<-NULL
  for(loopi in 1:nr){
    for(loopc in 1:Nc){
      if(str_count(Kernel.exp.NCRUM[loopi,loopc],"\\*")==1){
        if(strsplit(Kernel.exp.NCRUM[loopi,loopc],"\\*")[[1]][1]=='1'){next}
        temp_r<-noneOneallrParm[(grep(strsplit(Kernel.exp.NCRUM[loopi,loopc],"\\*")[[1]][1],noneOneallrParm))]
        if(temp_r%in%noneOneallrParm){
          rTransf<-c(rTransf,paste(" ",temp_r,"=",paste( paste('inv_logit(',Kernel.exp[loopi,loopc],")",sep=''),
                                                         "/",
                                                         piTransfwithoutName[loopi])
          ))
        }
      }
    }
  }
  rTransf<-paste(rTransf,');\n',sep='')
  rTransf<-str_replace_all(rTransf,"\\)\\);","\\);")

  #############################################################
  #######060319update: Multiple Group##########################
  #############################################################
  if(!is.na(fixeditem.vector)[1]){
    fixedItem.vector<-c(1:nr)[-fixeditem.vector]
  }else{fixedItem.vector<-c(1:nr)}

  #############################################################
  #######060719update: Multiple Group##########################
  #############################################################
  Kernel.exp.NCRUM<-Kernel.exp
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      Reparm[loopi,loopc]<-paste('  PImat[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp[loopi,loopc]),');\n',sep='')
    }
  }

  #############################################################
  #######060719update: Multiple Group End######################
  #############################################################

  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))  #############################################################


  Kernel.exp.NCRUM.groupName<-paste("Kernel.exp.NCRUM_g",c(1:group.num),sep='')
  for(i in 1:group.num){
    tempfill.Kernel.exp.NCRUM<-Kernel.exp.NCRUM
    temp.Kernel.exp.NCRUM<-Kernel.exp.NCRUM[fixedItem.vector,]
    for(j in 1:nrow(temp.Kernel.exp.NCRUM)){
      for(z in 1:ncol(temp.Kernel.exp.NCRUM)){
        temp.Kernel.exp.NCRUM[j,z]<-paste(temp.Kernel.exp.NCRUM[j,z],'_g',i,sep='')
      }
    }
    for(j in 1:nrow(temp.Kernel.exp.NCRUM)){
      for(z in 1:ncol(temp.Kernel.exp.NCRUM)){
        temp.Kernel.exp.NCRUM[j,z]<-str_replace_all(temp.Kernel.exp.NCRUM[j,z],"\\+",paste("_g",i,"+",sep=''))
      }
    }
    tempfill.Kernel.exp.NCRUM[fixedItem.vector,]<-temp.Kernel.exp.NCRUM
    assign(Kernel.exp.NCRUM.groupName[i],tempfill.Kernel.exp.NCRUM)
  }
  Kernel.exp.NCRUM.list<-list()
  for(i in 1:group.num){Kernel.exp.NCRUM.list[[i]]<-eval(parse(text=paste("Kernel.exp.NCRUM_g",i,sep='')))}
  #############################################################
  ##########060719update: Multiple Group End###################
  #############################################################

  #############################################################
  ##########060719update: Multiple Group#######################
  #############################################################

  PImat.groupName<-paste("PImat_g",c(1:group.num),sep='')
  Reparm.multigroup<-array(0,dim = c(nr,nclass,group.num))
  #Produce Stan code for PImat parameter
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      for (loopg in 1:group.num){
        Reparm.multigroup[loopi,loopc,loopg]<-paste('  PImat_g',loopg,'[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp.NCRUM.list[[loopg]][loopi,loopc]),');\n',sep='')

      }
    }
  }

  #############################################################
  ##########060319update: Multiple Group ######################
  #############################################################
  piTransf.multigroup<-array(0,dim = c(nr,1,group.num))
  rTransf.multigroup<-array(0,dim = c(length(rTransf),1,group.num))
  for(i in 1:group.num){
    tempfill.piTransf<-piTransf
    tempfill.rTransf<-rTransf

    tempfill.piTransf<-str_replace_all(tempfill.piTransf,"piTransf",paste("piParm_g",i,sep=''))
    #tempfill.rTransf<-str_replace_all(tempfill.rTransf,"rTransf",paste("rTransf_g",i,sep=''))
    for (j in 1:length(noneOneallrParm)){
      for( z in 1:length(tempfill.rTransf)){
        tempfill.rTransf[z]<-str_replace_all(tempfill.rTransf[z],noneOneallrParm[j],paste(noneOneallrParm[j],"_g",i,sep=''))

      }
    }

    temp.piTransf<-piTransf
    fixedItem.vector_rTransf<-NULL;
    for(j in fixedItem.vector){
      fixedItem.vector_rTransf<-c(fixedItem.vector_rTransf,
                                  grep(j,unlist(lapply(strsplit(unlist(lapply(strsplit(rTransf,"="),function(x){x[[1]]})),"_"),function(x){x[[1]]}))))
    }

    temp.rTransf<-rTransf[fixedItem.vector_rTransf]
    #060819update:a crtical step to handle piTransf lx_0 doesnt have group but some main effects have grouped difference
    Kernel.exp.fill<-Kernel.exp
    Kernel.exp.fill<-unique(Kernel.exp.fill)
    for (w in 1:nr){
      for( j in 1:nclass){
        for (z in 1:nr){
          #Kernel.exp.fill[w,j]<-str_replace_all(Kernel.exp.fill[w,j],paste(intercept[z],"\\+",sep=''),"")
          Kernel.exp.fill[w,j]<-str_replace_all(Kernel.exp.fill[w,j],intercept[z],"")

        }
      }
    }
    Kernel.exp.fill<-strsplit(paste(unlist(c(Kernel.exp.fill[fixedItem.vector,])),collapse="+"),"\\+")[[1]]
    Kernel.exp.fill<-Kernel.exp.fill[Kernel.exp.fill!='']
    Kernel.exp.fill<-unique(Kernel.exp.fill)
    Kernel.exp.fill<-c(Kernel.exp.fill,intercept[fixedItem.vector])
    for(j in 1:length(temp.piTransf)){
      temp.piTransf[j]<-str_replace_all(temp.piTransf[j],"piParm",paste("piParm_g",i,"",sep=''))
    }
    for(z in 1:length(Kernel.exp.fill)){
      for(j in 1:length(temp.piTransf)){
        temp.piTransf[j]<-str_replace_all(temp.piTransf[j],Kernel.exp.fill[z],paste(Kernel.exp.fill[z]
                                                                                    ,"_g",i,sep=''))
      }
    }



    for(j in 1:length(temp.rTransf)){
      #    temp.piTransf<-piTransf
      temp.rTransf[j]<-str_replace_all(temp.rTransf[j],"\\)",paste("_g",i,")",sep=''))
      temp.rTransf[j]<-str_replace_all(temp.rTransf[j],"\\+",paste("_g",i,"+",sep=''))
    }


    for (j in 1:length(noneOneallrParm)){
      for( z in 1:length(temp.rTransf)){
        temp.rTransf[z]<-str_replace_all(temp.rTransf[z],noneOneallrParm[j],paste(noneOneallrParm[j],"_g",i,sep=''))

      }
    }


    tempfill.piTransf[]<-temp.piTransf
    tempfill.rTransf[fixedItem.vector_rTransf]<-temp.rTransf

    piTransf.multigroup[,,i]<-tempfill.piTransf
    rTransf.multigroup[,,i]<-tempfill.rTransf

  }

  #############################################################
  ##########060319update: Multiple Group ######################
  #############################################################

  intercept.multigroup<-array(0,dim = c(nr,1,group.num))
  mainEff.multigroup<-array(0,dim = c(numMainEffect,1,group.num))
  interaction.multigroup<-array(0,dim = c(length(name.inter),1,group.num))
  #Group Invariant Parameter Name
  fixedParmName<-NULL
  if(!is.na(fixeditem.vector)[1]){
    for(i in fixeditem.vector){
      fixedParmName<-c(fixedParmName,out[[3]][[i]])
    }
    fixedParmName<-c(fixedParmName,out[[5]][fixeditem.vector])
  }
  #Group Variant Parameter Name
  freeParmName<-c(out[[5]],itemParmName)[!c(out[[5]],itemParmName)%in%fixedParmName]

  for(i in 1:group.num){
    tempfill.intercept<-out[[5]]
    tempfill.mainEff<-itemParmName[1:numMainEffect]
    tempfill.interaction<-name.inter

    temp.intercept<-tempfill.intercept[fixedItem.vector]
    temp.mainEff<-tempfill.mainEff[!tempfill.mainEff%in%fixedParmName]
    temp.interaction<-tempfill.interaction[!tempfill.interaction%in%fixedParmName]

    for(j in 1:length(temp.intercept)){
      temp.intercept[j]<-paste(temp.intercept[j],"_g",i,sep='')
    }

    for(j in 1:length(temp.mainEff)){
      temp.mainEff[j]<-paste(temp.mainEff[j],"_g",i,sep='')
    }

    for(j in 1:length(temp.interaction)){
      temp.interaction[j]<-paste(temp.interaction[j],"_g",i,sep='')
    }

    tempfill.intercept[fixedItem.vector]<-temp.intercept
    tempfill.mainEff[!tempfill.mainEff%in%fixedParmName]<-temp.mainEff
    tempfill.interaction[!tempfill.interaction%in%fixedParmName]<-temp.interaction

    intercept.multigroup[,,i]<-tempfill.intercept
    mainEff.multigroup[,,i]<-tempfill.mainEff
    interaction.multigroup[,,i]<-tempfill.interaction
  }

  Constrain.List<-paste('  real<lower=0>',unique(mainEff.multigroup),';\n ')
  Constrain.List<-c(Constrain.List,paste('  real<lower=0>',unique(interaction.multigroup),';\n '))
  Unconstrain.List<-paste('  real',unique(c(intercept.multigroup)),';\n ')
  #############################################################
  ##########060319update: Multiple Group End###################
  #############################################################


  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))
  update.Parmprior<-Parmprior
  fix.Parmprior<-NULL
  update.Parmprior<-update.Parmprior[update.Parmprior!='']

  #############################################################
  #####060319update:Change from update.Parmprior&fix     ######
  #############################################################
  update.Parmprior.multiGroup<-NULL
  for(i in 1:length(update.Parmprior)){
    for (j in 1:length(freeParmName)){
      for (z in 1:group.num){
        if(sum(grepl(freeParmName[j],update.Parmprior[i]))>=1){
          temp.update.Parmprior<-str_replace_all(update.Parmprior[i],freeParmName[j],paste(freeParmName[j],"_g",z,sep=''))
          update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,
                                         temp.update.Parmprior)
        }
      }
    }
  }
  update.Parmprior.multiGroup<-unique(update.Parmprior.multiGroup)
  update.Parmprior.multiGroup<-c("   //Prior\n",paste('   ',fixedParmName,'~normal(0,5)',';\n',sep=''),update.Parmprior.multiGroup )
  if(class.equal){
    update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='') )
  }else{
    for(i in 1:group.num){
      update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,
                                     paste('   Vc_g',i,'~dirichlet(rep_vector(2.0, Nc));\n',sep='') )
    }
  }

  ##therefore we can use: fix.Parmprior,update.Parmprior
  #############################################################
  #############################################################

  #############################################################
  #####060419update:Likelihood Add PImat_g#####################
  #############################################################
  PImat.likelihood1<-NULL
  PImat.likelihood0<-NULL
  Vc.likelihood<-NULL
  for(loopg in 1:group.num){
    temp.PImat.likelihood1<-paste(paste('          if (GroupID[iterp]==',loopg,')',sep=''),
                                  paste('            contributionsI[iteri]=bernoulli_lpmf(1|PImat_g',loopg,'[iteri,iterc]);',sep=''),
                                  sep='\n')
    PImat.likelihood1<-paste(PImat.likelihood1,temp.PImat.likelihood1,sep='\n')
    temp.PImat.likelihood0<-paste(paste('          if (GroupID[iterp]==',loopg,')',sep=''),
                                  paste('            contributionsI[iteri]=bernoulli_lpmf(0|PImat_g',loopg,'[iteri,iterc]);',sep=''),
                                  sep='\n')
    PImat.likelihood0<-paste(PImat.likelihood0,temp.PImat.likelihood0,sep='\n')
  }
  if(!class.equal){
    for(loopg in 1:group.num){
      temp.Vc.likelihood<-paste(paste('       if (GroupID[iterp]==',loopg,')',sep=''),
                                paste('         contributionsC[iterc]=log(Vc_g',loopg,'[iterc])+sum(contributionsI);',sep=''),
                                sep='\n')
      Vc.likelihood<-paste(Vc.likelihood,temp.Vc.likelihood,sep='\n')
    }
  }

  #Likelihood Stan code
  if(class.equal){
    Likelihood<-paste('
                      \n
                      //Likelihood
                      for (iterp in 1:Np){
                      for (iterc in 1:Nc){
                      for (iteri in 1:Ni){
                      if (Y[iterp,iteri] == 1)'
                      ,PImat.likelihood1,'\n',
                      '        else'
                      ,PImat.likelihood0,
                      '}
                      contributionsC[iterc]=log(Vc[iterc])+sum(contributionsI);
                      }
                      target+=log_sum_exp(contributionsC);
                      }
                      ',sep='')}else{
                        Likelihood<-paste('
                                          \n
                                          //Likelihood
                                          for (iterp in 1:Np){
                                          for (iterc in 1:Nc){
                                          for (iteri in 1:Ni){
                                          if (Y[iterp,iteri] == 1)'
                                          ,PImat.likelihood1,'\n',
                                          '        else'
                                          ,PImat.likelihood0,
                                          '}\n',
                                          Vc.likelihood
                                          ,'


                                          }
                                          target+=log_sum_exp(contributionsC);
                                          }
                                          ',sep='')

                      }
  #############################################################
  #####060419update:Likelihood Add PImat_g  end################
  #############################################################

  #Data Specification
  data.spec<-'
  data{
  int Np;
  int Ni;
  int Nc;
  matrix[Np, Ni] Y;
  vector[Np] GroupID;
  }
  '


  #############################################################
  #####060419update: Stan script##############################
  #############################################################
  Constrain.List<-unique(Constrain.List);Unconstrain.List<-unique(Unconstrain.List)
  #Parameter Specification
  if(class.equal){parm.spec<-paste(c('
                                     parameters{
                                     simplex[Nc] Vc;\n ',paste0(Constrain.List),paste0(Unconstrain.List),
                                     '}\n'),collapse='')}else{
                                       parm.spec<-paste(c('
                                                          parameters{\n ',
                                                          paste(paste('   simplex[Nc] Vc_g',1:group.num, ";",sep=''),"\n"),
                                                          paste0(Constrain.List),paste0(Unconstrain.List),
                                                          '}\n'),collapse='')


                                     }
  #To set rParm to real value
  rParm.real<-
    paste('real', str_replace_all(unique(unlist(lapply(strsplit(rTransf.multigroup,"="),function(x){x[[1]]})))," ",""),';\n',sep=' ')

  #Reparameter Specification
  tranrTransf.spec<-paste(c('
                            transformed parameters{

                            ',
                            paste('  matrix[Ni, Nc] PImat_g',1:group.num,';\n',sep=''),
                            paste('  vector[Ni] piParm_g',1:group.num,';\n',sep=''),
                            rParm.real,

                            unique(c(piTransf.multigroup)), #060419update
                            unique(c(rTransf.multigroup)), #060419update
                            paste0(c(Reparm.multigroup)),'}\n'),collapse='')

  #Model Specification update052619
  model.spec<-paste(c('\nmodel {\n',paste(c(Modelcontainer,update.Parmprior.multiGroup,Likelihood),sep=''),'\n}',sep=''))
  model.spec<-model.spec[!startsWith(str_remove_all(model.spec," "),"~")]

  #Generated Quantities Specification
  IC.generatedquantities<-NULL
  if(!class.equal){
    for(loopg in 1:group.num){
      temp.IC.generatedquantities<-paste(paste('       if (GroupID[iterp]==',loopg,')',sep=''),
                                         paste('         contributionsIC[iteri,iterc]=log(Vc_g',loopg,'[iterc])+contributionsI[iteri];',sep=''),
                                         sep='\n')
      IC.generatedquantities<-paste(IC.generatedquantities,temp.IC.generatedquantities,sep='\n')
    }
  }


 if(class.equal){
    generatedQuantities.spec<-paste('
  \n
  generated quantities {
  vector[Ni] log_lik[Np];
  vector[Ni] contributionsI;
  matrix[Ni,Nc] contributionsIC;

  matrix[Ni,Nc] posteriorIC;
  matrix[Np,Nc] posteriorPC;


  //Posterior
  for (iterp in 1:Np){
    for (iteri in 1:Ni){
      for (iterc in 1:Nc){
        if (Y[iterp,iteri] == 1)'
                                    ,PImat.likelihood1,'\n',
                                    '        else'
                                    ,PImat.likelihood0,
                                    '
          contributionsIC[iteri,iterc]=log(Vc[iterc])+contributionsI[iteri];
          posteriorIC[iteri,iterc]=contributionsI[iteri];
        }
      log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
    }
    for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
  }
  }
  ',sep='')}else{
    generatedQuantities.spec <- paste( '
    \n
    generated quantities {
    vector[Ni] log_lik[Np];
    vector[Ni] contributionsI;
    matrix[Ni,Nc] contributionsIC;

    matrix[Ni,Nc] posteriorIC;
    matrix[Np,Nc] posteriorPC;


    //Posterior
    for (iterp in 1:Np){
      for (iteri in 1:Ni){
        for (iterc in 1:Nc){
          if (Y[iterp,iteri] == 1)',
                                       PImat.likelihood1,'\n',
                                       '        else',
                                       PImat.likelihood0,'\n',
                                       "   ",IC.generatedquantities,'\n',
                             '          posteriorIC[iteri,iterc]=contributionsI[iteri];','\n
      }\n',
                                       '     log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
     }
     for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
   }
   }
      ',
                                       sep = ''
    )
  }




  if (.Platform$OS.type == "unix") {
    filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
  }else{
    filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
  }

  sink(file=filename, append=FALSE)
  cat(
    paste(c('   ',
            data.spec,parm.spec,tranrTransf.spec,model.spec,generatedQuantities.spec)
    ))
  sink(NULL)

}


#SEPERATION#
#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\data.RData")


StanORDM.run<-function(Qmatrix,response.matrix,
                       script.path=NA,save.path=getwd(),save.name="ORDM_uninf",
                       iter=1000,warmup = 0,
                       chain.num=3,init.list='random',control.list=NA){
  rstan.detect<-tryCatch(library("rstan"),error=function(e){"rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details."})
  if(length(rstan.detect)==1){
    stop()
  }
  Cdm.init<-F
  if(init.list=='cdm'){
    Cdm.init<-T
    Install.package(c("CDM","stringr"))
    trueParmName<-Parm.name(Qmatrix=Qmatrix)$parm.name
    Classp.exp1<-Parm.name(Qmatrix=Qmatrix)$class.expression
    mod1<-gdina( data =respMatrix, q.matrix = Qmatrix , maxit=700,link = "logit",progress=F)
    CDMresult<-as.data.frame(coef(mod1))
    library(stringr)
    CDM.parm.name<-paste(paste(paste('l',CDMresult[,3],sep=''),'_',sep=''),str_count(CDMresult$partype.attr,"Attr"),sep='')
    CDM.parm.name<-paste(CDM.parm.name,
                         unlist(lapply(strsplit(unlist(lapply(strsplit(CDMresult$partype.attr, 'Attr', fixed=FALSE),function(x){paste(x,collapse="")})),'-'),function(x){paste(x,collapse="")})),
                         sep='')
    CDM.parm.est<-CDMresult$est
    parm.ini<-round(CDM.parm.est[match(trueParmName,CDM.parm.name)],4)
    CDM.prop.est<-mod1$attribute.patt
    prop.ini<-CDM.prop.est[match(Classp.exp1,rownames(CDM.prop.est)),1]
    inilist1<-paste('list(',paste(noquote(paste(noquote(unlist(list(
      paste(paste('Vc=c(',paste((prop.ini),collapse=','),')',collapse=','))))
    ))),collapse=',')  ,')',collapse='')

    inilist1<-eval(parse(n =2000000 ,text=inilist1))
    for( i in 2:chain.num){
      temp.text<-paste('inilist',i,"<-inilist1",sep='')
      eval(parse(text=(temp.text)))
    }
    temp.text<-paste('init.list<-list(',paste(paste('inilist',1:chain.num,sep=''),collapse = ","),')',sep='')
    eval(parse(text=(temp.text)))
  }
  data.list<-Generate.datalist(Qmatrix,response.matrix)

  if(is.na(control.list)){control.list<-list(adapt_delta=0.82)}
  if(is.na(script.path)==T){
    options(warn=-1)
    StanORDM.script(Qmatrix,save.path=save.path,save.name=save.name)
    script.path<-paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
    options(warn=0)
    compiled_model<-stan_model(script.path)
  }else{
    compiled_model<-stan_model(script.path)
  }
  if(Cdm.init==T){
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})
  }else{
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})

  }

  estimated_model
}



#SEPERATION#
#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanORDM.script<-function(Qmatrix,scale.num,save.path=getwd(),save.name="ORDM_uninf"){
  nc<-ncol(Qmatrix)
  nr<-nrow(Qmatrix)
  temp.table.col<-unique(apply(combn(rep(c(0,1),nc),nc),2,function(x){paste(x,collapse = "")}))
  temp.table.col<-temp.table.col[order(temp.table.col)]
  temp.table<-matrix(0,nr,length(temp.table.col))
  colnames(temp.table)<-temp.table.col
  rownames(temp.table)<-paste('item',c(1:nr),sep='')
  temp.table<-as.data.frame(temp.table)
  for (i in 1:nr){
    temp.table[i,]<-paste('l',i,'_0',sep='')
  }
  intercept<-temp.table[,1]

  #Generate attribute combinations
  comb.generator<-function(x.vector){
    if(length(x.vector)>1){
      temp.attr<-x.vector
      temp.attr.sav<-NULL
      for(i in 1:length(temp.attr)){
        temp.1<-combn(temp.attr,i)
        temp.2<-apply(temp.1,2,function(x){paste(x,collapse = "")})
        temp.attr.sav<-c(temp.attr.sav,temp.2)
      }
    }
    if(length(x.vector)==1){temp.attr.sav<-x.vector}
    temp.attr.sav
  }
  #vectors needed for combination.generator
  Item.load.id<-list()
  for ( i in 1:nr){
    Item.load.id[[i]]<-grep('1',Qmatrix[i,])}

  Attr.load.id<-list()
  attr.load.id<-matrix(0,length(temp.table.col),nc)
  for ( i in 1:length(temp.table.col)){
    attr.load.id[i,]<-unlist(strsplit(temp.table.col[i],split=''))
    Attr.load.id[[i]]<-grep('1',attr.load.id[i,])
  }

  #Generate Combination for both Item.load and Attr.load
  Item.Comb<-list()
  for ( i in 1:nr){
    Item.Comb[[i]]<-comb.generator(Item.load.id[[i]])
  }
  Attr.Comb<-list()
  for ( i in 2:length(temp.table.col)){
    Attr.Comb[[1]]<-0
    Attr.Comb[[i]]<-comb.generator(Attr.load.id[[i]])
  }
  constraints.list<-list()
  nway.inter.list<-list()
  for(i in 1:nr){
    for(a in 2:length(temp.table.col)){
      ifzero<-as.numeric(paste(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],collapse=''))
      if((!is.na(ifzero))){
        temp.table[i,a]<-paste(c(temp.table[i,a],
                                 paste("S","l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='',collapse='')
        ),collapse='')
        if(a==length(temp.table.col)){
          nway.inter.list[[i]]<-nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])])
          constraints.list[[i]]<-paste("l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='')
        }
      }
    }
  }

  #Create Lambda Table
  Lamda.Table<-temp.table
  for(i in 1:nr){
    for(a in 1:length(Lamda.Table)){
      t.ref<-unique(as.character(Lamda.Table[i,]))
      pos<-c(1:length(t.ref))[Lamda.Table[i,a]==t.ref]
      temp.table[i,a]<-paste("t",i,"_",pos,sep='')}}

  #Generate LCDM specification
  out<-list()
  out[[1]]<-Lamda.Table
  out[[2]]<-temp.table
  out[[3]]<-constraints.list
  out[[4]]<-nway.inter.list
  out[[5]]<-intercept
  OUTPUT<-out
  nclass<-ncol(OUTPUT[[1]]);Nc<-nclass

  #Produce kernel expressions across items and attributes
  Kernel.exp<-OUTPUT[[1]]
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])}
    }
  }


  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                        nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))

  if(length(name.inter)!=0){
    for (inter in 1: length(name.inter)){
      temp.nw<-numway.inter[inter]
      temp.nm<-name.inter[inter]
      temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
      temp.sel<-paste(unlist(strsplit(temp.nm,split = '_'))[1],"_",(1:(temp.nw-1)),sep='')
      first.sel<-unlist(OUTPUT[[3]])[grep(paste((temp.sel),collapse="|"),unlist(OUTPUT[[3]]))]
      second.sel<-sub(".*_.", "", first.sel)
      for (sel in 1:length(temp.subnm)){
        SEL<-second.sel[sel]
        Constrain.List1<-rbind(
          paste(temp.nm,">-(0", paste("+",first.sel[grep(SEL,second.sel)],
                                      sep='',collapse=''),")",sep=''),Constrain.List1)
      }
    }
    Constrain.List1<-as.character(Constrain.List1)
  }else{
    Constrain.List1<-NULL
  }

  itemParmName<-c(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==2],OUTPUT[[5]])
  numMainEffect<-length(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1])
  Constrain.List<-paste('  real<lower=0>',itemParmName[1:numMainEffect],';\n ')
  Unconstrain.List<-paste('  real',itemParmName[-(1:numMainEffect)],';\n ')

  Reparm<-array(rep(0,nr*nclass*(scale.num)),dim = c(nr,nclass,(scale.num)))

  #step parameters
  StepParmName<-paste('step',1:scale.num-1,sep='')
  StepParmName[1]<-0
  CulStepParmName<-StepParmName

  ##########0529update Kernal.exp is updated to a list so that list[[?]] is a PImat for a cate
  Kernel.exp.list<-list(scale.num)
  for(i in 1:scale.num){
    Kernel.exp.list[[i]]<-Kernel.exp
  }
  Kernel.exp.list.copy<-Kernel.exp.list
  for(i in 1:scale.num){
    for (j in 1:nr){
      for( z in 1:nclass){
        Kernel.exp.list[[i]][j,z]<-str_replace_all(Kernel.exp.list[[i]][j,z],Kernel.exp.list.copy[[i]][j,1],
                                                   paste(Kernel.exp.list.copy[[i]][j,1],"r",i-1,sep=''))

      }
    }
  }
  for(i in 1:scale.num){
    for (j in 1:nr){
      for( z in 2:nclass){
        Kernel.exp.list[[i]][j,z]<-str_replace_all(Kernel.exp.list[[i]][j,z],"\\+",
                                                   paste("+",i-1,"*",sep=''))

      }
    }
  }
  Kernel.exp.list[[1]]<-0
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      for( loops in 1:1){Reparm[loopi,loopc,loops]<-paste('  PImat[',loopi,',',loopc,'][',loops,']=0;\n',sep='')}
      for( loops in 2:scale.num){
        Reparm[loopi,loopc,loops]<-paste('  PImat[',loopi,',',loopc,'][',loops,']=',
                                         paste(Kernel.exp.list[[loops]][loopi,loopc]),';\n',sep='')
      }
    }
  }
  itemParmName<-c(unlist(OUTPUT[[3]]))
  for( i in 2:scale.num){
    itemParmName<-c(itemParmName,c(Kernel.exp.list[[i]][,1]))
  }
  Unconstrain.List<-paste('  real',itemParmName[-(1:numMainEffect)],';\n ')

  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,5)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))
  #Likelihood Stan code
  Likelihood<-'
  \n
  //Likelihood
  for (iterp in 1:Np){
    for (iterc in 1:Nc){
      for (iteri in 1:Ni){
        contributionsI[iteri]= categorical_lpmf(Y[iterp,iteri]| softmax(((PImat[iteri,iterc]))));
      }
      contributionsC[iterc]=log(Vc[iterc])+sum(contributionsI);
    }
    target+=log_sum_exp(contributionsC);
  }
  '
  data.spec<-'
  data{
  int Np;
  int Ni;
  int Nc;
  int Ns;
  matrix[Np, Ni] Y;
  }
  '
  parm.spec<-paste(c('
                     parameters{
                     simplex[Nc] Vc;\n ',paste0(Constrain.List),paste0(Unconstrain.List),
                     '}\n'),collapse='')

  #Parameter Specification
  parm.spec<-paste(c('parameters{
                     simplex[Nc] Vc;\n ',paste0(Constrain.List),paste0(Unconstrain.List),'}\n'),collapse='')

  #Reparameter Specification

  transparm.spec<-paste(c('transformed parameters{
                          vector[Ns] PImat[Ni, Nc];\n',
                          paste0(unlist(Reparm)),'}\n'),collapse='')

  #Model Specification
  model.spec<-paste(c('\nmodel {\n',paste(c(Modelcontainer,Parmprior,Likelihood),sep=''),'\n}',sep=''))
  model.spec<-model.spec[!startsWith(str_remove_all(model.spec," "),"~")]

  generatedQuantities.spec<-'
  \n
  generated quantities {
  vector[Ni] log_lik[Np];
  vector[Ni] contributionsI;
  matrix[Ni,Nc] contributionsIC;
 matrix[Ni,Nc] posteriorIC;
 matrix[Np,Nc] posteriorPC;



  //Posterior
  for (iterp in 1:Np){
    for (iteri in 1:Ni){
      for (iterc in 1:Nc){
        contributionsI[iteri]= categorical_lpmf(Y[iterp,iteri]| softmax(((PImat[iteri,iterc]))));
        contributionsIC[iteri,iterc]=log(Vc[iterc])+contributionsI[iteri];
        posteriorIC[iteri,iterc]=contributionsI[iteri];
      }
      log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
    }
    for (iterc in 1:Nc){posteriorPC[iterp,iterc]=prod(exp(posteriorIC[,iterc]));}
  }
  }
  '
  if (.Platform$OS.type == "unix") {
    filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
  }else{
    filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
  }

  sink(file= filename,append=FALSE)
  cat(
    paste(c('   ',
            data.spec,parm.spec,transparm.spec,model.spec,generatedQuantities.spec)
    ))
  sink(NULL)

}

#SEPERATION#
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

#SEPERATION#
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

#SEPERATION#
#' @title Generate Stan code and Run the estimation for LCDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param savepath save the .stan file to somewhere; the default path is getwd()
#' @param savename name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\data.RData") ;Qmatrix<-cbind(Qmatrix,rep(1,9));Qmatrix[1,1]<-0

StanPrior.show<-function(script.path){
  Install.package("plyr")
  Install.package('stringr')
  readStan.script<- readLines(script.path)
  readStan.script<-readStan.script[grep("~",readStan.script)]
  readStan.script<-str_remove_all(readStan.script," ")
  readStan.script<-str_remove_all(readStan.script,";")
  readStan.script
}



#SEPERATION#
#' @title Generate Stan code and Run the estimation for LCDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param savepath save the .stan file to somewhere; the default path is getwd()
#' @param savename name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\data.RData") ;Qmatrix<-cbind(Qmatrix,rep(1,9));Qmatrix[1,1]<-0
StanPrior.update<-function(priorUpdate.matrix,script.path,save.path=getwd(),save.name=NA){
  Install.package("plyr")
  Install.package('stringr')
  options(warn=-1)
  if(is.na(save.name)){
    saveUpdate.name=paste(format(Sys.time(), "%a%b%Y"),'_priorUpdate.stan',sep='')}else{
      saveUpdate.name=paste(paste(save.name,sep='/'),'.stan',sep='')
    }
  readStan.script<- readLines(script.path)
for(i in 1:nrow(priorUpdate.matrix)){
  beReplaced<-readStan.script[str_detect(str_remove_all(str_remove_all(readStan.script," "),";"),priorUpdate.matrix[i,1])][grep("~",readStan.script[str_detect(str_remove_all(str_remove_all(readStan.script," "),";"),priorUpdate.matrix[i,1])])
  ]
  beReplaced<-str_remove_all(str_remove_all(beReplaced," "),";")
  
  readStan.script[which(str_remove_all(str_remove_all(readStan.script," "),";")==beReplaced)]<-paste("   ",
                                                                                                    paste(priorUpdate.matrix[i,1],priorUpdate.matrix[i,2],sep="~")  ,
                                                                                                                  ';//UpdatedPrior ')
}
  options(warn=0)
  filename = paste(paste(save.path,saveUpdate.name,sep='/'),sep='')
  writeLines(readStan.script,filename)
}