require(data.table)
program = "ONBS_CNBS"
require(AlphaSimR)
require(nadiv)
require(sampling)
require(optiSel)
require(visPedigree)
require(AlphaMME)
require(AlphaLearn)
source("ocs_SC.R")
source("ocsBLUP_SC.R")

abstractebv_blupf90 = function(pop,dt){
  reqid = paste("Animal_",sprintf("%0*d",9,as.numeric(pop@id)),sep = "")
  popebv = dt$ebv[match(reqid,dt$AnimalID)]
  popebv = matrix(popebv,ncol = 1)
  return(popebv)
}

abstractebv = function(pop,dt){
  reqid = paste("Animal_",sprintf("%0*d",9,as.numeric(pop@id)),sep = "")
  popebv = dt$ebv[match(reqid,dt$AnimalID)]
  return(popebv)
}

calinbya2 <- function(snp_012_dt,snpfre){
  calformula <- function(x, snpfre){
    y <- mean((x^2 - ((1 + 2 * snpfre) * x) + 2 * (snpfre^2)) / (2 * snpfre * (1 - snpfre)), na.rm = TRUE)
    return(y)
  }
  inb <- apply(snp_012_dt, 1, calformula,snpfre = snpfre)
  return(inb)
}

makeped = function(pop){
    
    z=pullSnpGeno(pop)
    
    z[z==2]=22
    z[z==1]=12
    z[z==0]=11

    #z = as.data.table(cbind(pop@id,z))
    z = as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop@id)),sep = ""),z))
    return(z)
  }

makemap = function(...){
  mapck = getSnpMap()
  mapck = mapck[colnames(peddt)[-1],]
  mapck$id = mapck$chr
  mapck$chr = rownames(mapck)
  mapck$site = mapck$pos
  mapck$pos = rep(0,nrow(mapck))
  return(mapck)
  }

#ooutout ped and pheno file with MP candidates to work director
#output ped_tidy with MP candidates in environment as varible
outputdataforMPcand = function(pop_MP_candidate,g,alldata,ped,wkdir){

collcetdata_MP_candidate <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop_MP_candidate@id)),sep = ""),
                                                      rep("0",pop_MP_candidate@nInd),
                                                      rep("0",pop_MP_candidate@nInd),
                                                      pheno(pop_MP_candidate),
                                                      pop_MP_candidate@sex,
                                                      pop_MP_candidate@gv))
      colnames(collcetdata_MP_candidate) <- c("AnimalID","SireID","DamID","M2BW","SexID","GV")
      collcetdata_MP_candidate[,c("Generation"):=rep(paste("G",g,sep = ""),pop_MP_candidate@nInd)]
      collcetdata_MP_candidate$M2BW <- as.numeric(collcetdata_MP_candidate$M2BW)
      collcetdata_MP_candidate[,c("FamilyID"):=NA]
      collcetdata_MP_candidate[,c("ebv"):=NA]
      data_allMP<- rbind(alldata,collcetdata_MP_candidate)

ped_MP_current = collcetdata_MP_candidate[,1:3]

ped_allMP = rbind(ped,ped_MP_current)

ped_tidy <- visPedigree::tidyped(ped_allMP,cand = data_allMP$AnimalID)
setorder(data_allMP,AnimalID)
setorder(ped_tidy,Ind)

fwrite(data_allMP,file = paste(wkdir,"/","pheno.csv",sep=""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "0")

fwrite(ped_tidy[,1:3],file = paste(wkdir,"/","ped.csv",sep=""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "0")

#return ped_tidy for ocsBLUP
return(ped_tidy)
}

creatmpdata = function(pop_MP_candidate,g,alldata){

collcetdata_MP_candidate <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop_MP_candidate@id)),sep = ""),
                                                      rep("0",pop_MP_candidate@nInd),
                                                      rep("0",pop_MP_candidate@nInd),
                                                      pheno(pop_MP_candidate),
                                                      pop_MP_candidate@sex,
                                                      pop_MP_candidate@gv))
      colnames(collcetdata_MP_candidate) <- c("AnimalID","SireID","DamID","M2BW","SexID","GV")
      collcetdata_MP_candidate[,c("Generation"):=rep(paste("G",g,sep = ""),pop_MP_candidate@nInd)]
      collcetdata_MP_candidate$M2BW <- as.numeric(collcetdata_MP_candidate$M2BW)
      collcetdata_MP_candidate[,c("FamilyID"):=NA]
      collcetdata_MP_candidate[,c("ebv"):=NA]
      data_allMP<- rbind(alldata,collcetdata_MP_candidate)
      setorder(data_allMP,AnimalID)
#return ped_tidy for ocsBLUP
return(data_allMP)
}

calibd = function(pop_founder){
 pop_inb_founder= makeDH(pop_founder)
 x = round(((meanG(pop_founder) - meanG(pop_inb_founder))/meanG(pop_founder))*100,digits =3)
 return(x)
}

calmeanD = function(pop){
  mean_D = meanG(pop) - meanG(makeDH(pop))
  return(mean_D)
}