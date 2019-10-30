
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
