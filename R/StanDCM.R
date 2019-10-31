#targetPath<-"\\\\libfs1\\users$\\home\\zjiang17\\My Documents\\GitHub\\StanDCM\\R\\";setwd(targetPath)
#AllFiles<-list.files(targetPath)
#AllFiles<-AllFiles[AllFiles!="StanDCM.R"]

#for(i in 1:length(AllFiles)){
  
#  line <- readLines(AllFiles[i])
#  write(line, "StanDCM.txt",append = T)
#  write("\n#SEPERATION#", "StanDCM.txt",append = T)
#}



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


