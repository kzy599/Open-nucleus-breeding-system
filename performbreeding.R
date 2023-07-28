for(g in 0:nGeneration){

    source("tidydtformati.R")

    source("blupf90format.R")

   if(gem=="ssGBLUP"&!exists("nucleusgeno")){
    nucleusgeno = foundergeno
   }
   
   if(gem=="ssGBLUP"&exists("parents_MP")){
     MPpgeno = MPgeno[rownames(MPgeno)%in%paste("Animal_",sprintf("%0*d",9,as.numeric(parents_MP)),sep = ""),]
  
     nucleusgeno = rbind(nucleusgeno,MPpgeno)
   }
   
   if(gem=="ssGBLUP"){
     genodt<- as.data.table(nucleusgeno,keep.rownames = TRUE)

     rn<- genodt$rn

     genodt<- genodt[,-1]

     genodt <- apply(genodt,1,function(x){
     return(paste(x,sep="",collapse = ""))
     })

     genodt <-data.table(rn,genodt)

     fwrite(genodt,file = "geno_selectedparents.txt",sep = " ",col.names = FALSE, row.names = FALSE,quote = FALSE)
     source("ssGBLUP.R")
   }else {
    source("pblup.R")
   }


    pop@ebv = abstractebv_blupf90(pop = pop, dt = alldata )

    initialFam<- selectFam(pop,nFam = ncf,trait = 1,use = "ebv",simParam = SP)

    initialInd_F <- selectWithinFam(initialFam,nInd = ngi*(2/3),trait = 1,use = "ebv",simParam = SP,sex  ="F")

    initialInd_M <- selectWithinFam(initialFam,nInd = ngi*(1/3),trait = 1,use = "ebv",simParam = SP,sex  ="M")

    initialInd <- c(initialInd_F,initialInd_M)

    if(exists("pop_MP")){

      pop_MP_candidate_F <- selectWithinFam(pop_MP,nInd = 4,trait = 1,use = "pheno",simParam = SP,sex = "F")

      pop_MP_candidate_M <- selectWithinFam(pop_MP,nInd = 2,trait = 1,use = "pheno",simParam = SP,sex = "M")
      
      pop_MP_candidate <- c(pop_MP_candidate_M,pop_MP_candidate_F)
    }


#collect genomic information

if(exists("nucleusgeno")){
    nucleusgeno1 = pullSnpGeno(initialInd)
    
    rownames(nucleusgeno1)<- paste("Animal_",sprintf("%0*d",9,as.numeric(rownames(nucleusgeno1))),sep = "")
    
    nucleusgeno = rbind(nucleusgeno,nucleusgeno1) 
}

if(gem=="ssGBLUP"&exists("pop_MP_candidate")){
    MPgeno= pullSnpGeno(pop_MP_candidate)
    rownames(MPgeno)<- paste("Animal_",sprintf("%0*d",9,as.numeric(rownames(MPgeno))),sep = "") 
    
    genodt= rbind(nucleusgeno,MPgeno)

    genodt<- as.data.table(genodt,keep.rownames = TRUE)
}else if(gem=="ssGBLUP"){
    genodt<- as.data.table(nucleusgeno,keep.rownames = TRUE)
}




if(exists("pop_MP_candidate")){
pedmore = outputdataforMPcand(pop_MP_candidate = pop_MP_candidate, ped = ped,g = g, alldata = alldata, wkdir = getwd())
data_allMP = creatmpdata(pop_MP_candidate = pop_MP_candidate,g = g,alldata = alldata)
source("blupf90format.R")
}

if(gem=="ssGBLUP"){
  rn<- genodt$rn

  genodt<- genodt[,-1]

  genodt <- apply(genodt,1,function(x){
  return(paste(x,sep="",collapse = ""))
  })

  genodt <-data.table(rn,genodt)

  fwrite(genodt,file = "geno_selectedparents.txt",sep = " ",col.names = FALSE, row.names = FALSE,quote = FALSE)
  source("ssGBLUP.R")


  if(exists("pop_MP_candidate")){
  initialInd@ebv = abstractebv_blupf90(pop = initialInd, dt = data_allMP )
  
  pop_MP_candidate@ebv = abstractebv_blupf90(pop = pop_MP_candidate, dt = data_allMP)

  pop_MP_candidate_M@ebv<- abstractebv_blupf90(pop = pop_MP_candidate_M, dt = data_allMP)

  pop_MP_candidate_F@ebv<- abstractebv_blupf90(pop = pop_MP_candidate_F, dt = data_allMP)
}else{
  
  initialInd@ebv = abstractebv_blupf90(pop = initialInd, dt = alldata )

}

}else{
  if(exists("pop_MP_candidate")){
  
  pop_MP_candidate@ebv = pop_MP_candidate@pheno

  pop_MP_candidate_M@ebv<- pop_MP_candidate_M@pheno

  pop_MP_candidate_F@ebv<- pop_MP_candidate_F@pheno
}
}


#select families and individuals according to EBV or GEBV
#select the top 8 families as the candidate families to consturct MP in next generation
#select the top 9 individuals (6 females and 3 males) of each those families 
if(pname%in%c("op","os")){
MP_fam <- selectFam(initialInd,nFam = ncf_MP,trait = 1,use = "ebv",simParam = SP)
MP_fam_ind_F <- selectWithinFam(MP_fam,nInd = nci_MP*(2/3),trait = 1,use = "ebv",simParam = SP,sex = "F")
MP_fam_ind_M <- selectWithinFam(MP_fam,nInd = nci_MP*(1/3),trait = 1,use = "ebv",simParam = SP,sex = "M")
MP_fam_ind <- c(MP_fam_ind_M,MP_fam_ind_F)
}

#seleced top 6 individuals (4 females and 2 males) of genotyped individuals in each top 50 families
nucleus_fam_ind_F <- selectWithinFam(initialInd,nInd = nci*(2/3),trait = 1,use = "ebv",simParam = SP,sex = "F")
nucleus_fam_ind_M <- selectWithinFam(initialInd,nInd = nci*(1/3),trait = 1,use = "ebv",simParam = SP,sex = "M")
pop_nucleus_candidate <- c(nucleus_fam_ind_M,nucleus_fam_ind_F)



source("calculation_population_parameters.R")

#
source("creatNPandMP.R")

}
if(pname == "os"){
  if(nooverlap==TRUE){
     fwrite(output,file = paste("ONBS","realg",0,r,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)
  }else{
    fwrite(output,file = paste("ONBS","realg",(nMP/150)*10,r,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)
  }
}else if(pname == "op"){
  if(nooverlap==TRUE){
    fwrite(output,file = paste("ONBS","realgp",0,r,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)
  }else{
    fwrite(output,file = paste("ONBS","realgp",(nMP/150)*10,r,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)
  }
}else if(pname == "cs"){
  fwrite(output,file = paste("CNBS","real1",r,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)
}else if(pname == "cp"){
  fwrite(output,file = paste("CNBS","real2",r,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)
}

