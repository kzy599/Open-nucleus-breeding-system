if(g==0){
#collect data of the founder parents

#recored id of them
parents <- c(unique(pop@mother),unique(pop@father))

#abstract the data form the founder population
pop_parents <- pop_founder[parents]

foundergeno = pullSnpGeno(pop_parents)
    
rownames(foundergeno)<- paste("Animal_",sprintf("%0*d",9,as.numeric(rownames(foundergeno))),sep = "")


#collect data
collcetdata_founder <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop_parents@id)),sep = ""),
                                           rep("0",pop_parents@nInd),
                                           rep("0",pop_parents@nInd),
                                           pheno(pop_parents),
                                           pop_parents@sex,
                                           pop_parents@gv))
#collect current generation                                           
colnames(collcetdata_founder) <- c("AnimalID","SireID","DamID","M2BW","SexID","GV")

#record generation
collcetdata_founder[,c("Generation"):=rep(paste("G","founder",sep = ""),pop_parents@nInd)]

#make the type of the pheno be numeric
collcetdata_founder$M2BW <- as.numeric(collcetdata_founder$M2BW)

#add col for family
collcetdata_founder[,c("FamilyID"):=NA]

#add col for EBV
collcetdata_founder[,c("ebv"):=NA]
collcetdata_founder[,c("ebvib"):=NA]
if(gem=="ssGBLUP"){

  if(!exists("snpfre_v")){
  snp_012_dt = pullSnpGeno(pop_founder)
  snpfre_v <- apply(snp_012_dt, 2, function(x){
    single_snpfre_s <- sum(x, na.rm = TRUE)/(2*length(x))
    return(single_snpfre_s)
  })
}

  collcetdata_founder[,c("ib"):=calinbgeno(poplist = pop_parents,snpfre = snpfre_v)]
}else{
  collcetdata_founder[,c("ib"):= 0]
}

collcetdata_founder[,c("poptype"):= "founder"]

}

if(g>0 & !exists("mix_FemalePart") ){
    parents <- c(unique(pop@mother),unique(pop@father))
    
    pop_parents <- pop_nucleus_candidate[parents]

}else if(exists("mix_FemalePart")){
    parents_MP <- c(unique(mix_FemalePart@mother),unique(mix_MalePart@father))
    
    parents_nucleus <- setdiff(c(unique(pop@mother),unique(pop@father)),parents_MP)
    
    pop_parents_nucleus <- pop_nucleus_candidate[parents_nucleus]
    
    pop_parents_MP <- pop_MP_candidate[parents_MP]
}

if(exists("pop_parents_MP")){
   collcetdata_MP_parents <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop_parents_MP@id)),sep = ""),
                                                      rep("0",pop_parents_MP@nInd),
                                                      rep("0",pop_parents_MP@nInd),
                                                      pheno(pop_parents_MP),
                                                      pop_parents_MP@sex,
                                                      pop_parents_MP@gv))
      colnames(collcetdata_MP_parents) <- c("AnimalID","SireID","DamID","M2BW","SexID","GV")
      
      collcetdata_MP_parents[,c("Generation"):=rep(paste("G",(g-1),sep = ""),pop_parents_MP@nInd)]
      
      collcetdata_MP_parents$M2BW <- as.numeric(collcetdata_MP_parents$M2BW)
      
      collcetdata_MP_parents[,c("FamilyID"):=NA]
      
      collcetdata_MP_parents[,c("ebv"):=NA]
      collcetdata_MP_parents[,c("ebvib"):=NA]

      if(gem == "ssGBLUP"){
        collcetdata_MP_parents[,c("ib"):=data_allMP$ib[match(collcetdata_MP_parents$AnimalID,data_allMP$AnimalID)]]
      }else{
        collcetdata_MP_parents[,c("ib"):=data_allMP$ib[match(collcetdata_MP_parents$AnimalID,data_allMP$AnimalID)]]
      }
      collcetdata_MP_parents[,c("poptype"):="MP"]
      ped_MP_parents_full<- collcetdata_MP_parents[,1:3]
      
      ped_MP_parents_full[,SireID:=paste("Animal_",sprintf("%0*d",9,as.numeric(pop_parents_MP@father)),sep = "")]
      
      ped_MP_parents_full[,DamID:=paste("Animal_",sprintf("%0*d",9,as.numeric(pop_parents_MP@mother)),sep = "")]

      ped_MP_parents <- collcetdata_MP_parents[,1:3]
      
      if(any(c(ped_MP_parents$DamID,ped_MP_parents$SireID)!="0")){
       stop("data.table is altered")
      }

}

#collect current generation
CollectData <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop@id)),sep = ""),
                                   paste("Animal_",sprintf("%0*d",9,as.numeric(pop@father)),sep = ""),
                                   paste("Animal_",sprintf("%0*d",9,as.numeric(pop@mother)),sep = ""),
                                   pheno(pop),
                                   pop@sex,
                                   pop@gv))

#format colnames
colnames(CollectData) <- c("AnimalID","SireID","DamID","M2BW","SexID","GV")

#record generation
CollectData[,c("Generation"):=rep(paste("G",g,sep = ""),nIndPerGeneration)]

#record family ID
CollectData[,c("FamilyID"):=paste(CollectData$SireID,"_",CollectData$DamID,sep="")]

#add col for ebv
CollectData[,c("ebv"):=NA]
CollectData[,c("ebvib"):=NA]

if(gem == "ssGBLUP"){
  CollectData[,c("ib"):=calinbgeno(poplist = pop,snpfre = snpfre_v)]
}else{
  
  CollectData[,c("ib"):= NA ]
}
CollectData[,c("poptype"):= "NP" ]
#make the type of M2BW be numeric
CollectData$M2BW <- as.numeric(CollectData$M2BW)

#pedigree of current generation
ped_thisGeneration <- CollectData[,1:3]

#record animalID
Aid [(g+1)]<- list(CollectData$AnimalID)

#record parents
Pid[(g+1)] <- list(c(unique(CollectData$SireID),unique(CollectData$DamID)))

#数据整合
if(g==0){
  alldata <- rbind(collcetdata_founder,CollectData)
}else if(g>0&!exists("collcetdata_MP_parents")){
  alldata <- rbind(alldata,CollectData)
}else if(exists("collcetdata_MP_parents")){
  alldata <- rbind(alldata,collcetdata_MP_parents,CollectData)
}

if(gem != "ssGBLUP") {alldata[AnimalID%in%ped_thisGeneration$AnimalID,ib:=calinbped(alldata = alldata,ped_thisGeneration = ped_thisGeneration)]}

if(g==0){
  ped <- ped_thisGeneration
  ped_full<- ped_thisGeneration
}else if(g>0&!exists("ped_MP_parents")){
  ped <- rbind(ped,ped_thisGeneration)
  ped_full <- rbind(ped_full,ped_thisGeneration)
}else if(exists("ped_MP_parents")){
  ped <- rbind(ped,ped_MP_parents,ped_thisGeneration)
  ped_full <- rbind(ped_full,ped_MP_parents_full,ped_thisGeneration)
}

#output pedigree that multipliers no parents' information. 
ped_tidy <- visPedigree::tidyped(ped,cand = alldata$AnimalID)

#output pedigree with no error
ped_full_tidy <- visPedigree::tidyped(ped_full,cand = alldata$AnimalID)


#the order process is not necessary as the blupf90 software would do it automatic
setorder(alldata,AnimalID)
setorder(ped_tidy,Ind)
setorder(ped_full_tidy,Ind)

fwrite(alldata,file = "pheno.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "0")
fwrite(ped_tidy[,1:3],file = "ped.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "0")
fwrite(ped_full_tidy[,1:3],file = "ped_full.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "0")

gc()

#if(g == 0){
 # tempid= paste("Animal_",sprintf("%0*d",9,as.numeric(pop@id)),sep = "")

#popdt= data.table(gen = g, gvalue = pop@gv,bvalue = bv(pop),phe = pop@pheno,ib = alldata$ibmatch(tempid,alldata$AnimalID))

  #colnames(popdt) = c("gen","gvalue","bvalue","phe","ib")
#}