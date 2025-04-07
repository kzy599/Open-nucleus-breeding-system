rm(list=ls())
gc()
library(MASS)
setwd("/home/GTDisk1/kangziyi/onbsresultac/")
getgg<- function(W,idc,Finitial=0,Va,Finb){
  getpervalue = function(x){
    pervalue = c(numeric(length = length(x)))
    for(i in 1:length(x)){
      if(i == 1){
        pervalue[i] = 0  
      }else{
        pervalue[i] = x[i] - x[i-1]
      }
    }
    return(pervalue)
  }
  ggper = getpervalue(W)
  inbper = getpervalue(Finb)
  idc = idc*(Va/5.39)
  idc[1] = 0.325
  idc = idc[1:20]
  idc = c(0.325,idc)
  ggperinb = ggper*(1 - (idc*inbper))
  gg = c(numeric(length = length(W)))
  for(i in 1:length(ggperinb)){
    gg[i] = sum(ggperinb[1:i])
  }
  return(gg)
}

getprate = function(x){
  prate = numeric(length = (length(x)-1))
  
  for(i in 2:length(x)){
    prate[i-1] = (x[i]-x[i-1])/(1-x[i-1])
  }
  prate = rep(mean(prate),(length(x)))
  return(prate)
}

abstractdata= function(ONBSfinalg11,pname,generation){
  gv<- ONBSfinalg11$mean_gv-ONBSfinalg11$mean_gv[1]
  Avalue<- ONBSfinalg11$mean_A-ONBSfinalg11$mean_A[1]
  Dvalue<- ONBSfinalg11$mean_D-ONBSfinalg11$mean_D[1]
  gg = c(0:generation)
  ld <- ONBSfinalg11$LD
  inb_plink <- ONBSfinalg11$inbreeding_plink
  inb <- ONBSfinalg11$inbreeding
  inbprate = getprate(x = inb_plink) 
  inbprates = rep(unique(getprate(x = inb_plink[1:11])),length(gg))
  inbrate = getprate(x = inb)
  Ne <- ONBSfinalg11$Ne
  Va <- ONBSfinalg11$Va
  Vg <- ONBSfinalg11$Vg
  Vd <- ONBSfinalg11$Vd
  genicVa <- ONBSfinalg11$genicVa
  #gv = gv/sqrt(Vg[1])
  #losvg= (sqrt(Vg[1]) - sqrt(Vg))/sqrt(Vg[1])
  calconrate = function(gva,agv,ngen){
    losvg= sqrt((genicVa[1] - genicVa[1:(ngen+1)])/genicVa[1])
    
    conrate= data.table(Avalue[1:(ngen+1)]/sqrt(genicVa[1]),losvg)
    
    colnames(conrate) = c("savalue","losvg")
    
    conrate= rlm(formula =savalue~losvg ,data = conrate)
    
    conrate= as.numeric(conrate$coefficients[2])
    
    conrate= rep(conrate,generation+1)
    return(conrate)
  }
  calrgv = function(gv,ngen){
    if(ngen ==10){
     rag = (gv[(ngen+1)] - gv[1])/10
    }else{
      rag = (gv[(ngen+1)] - gv[11])/10
    }
    return(rep(rag,generation+1))
  }
  
  conrate = calconrate(gva = genicVa,agv = Avalue,ngen = 20)
  
  conrate10 = calconrate(gva = genicVa,agv = Avalue,ngen = 10)
  #水平值（绝对值）的变化可以这样来分段，
  #如果看不同世代区间的增长率（相对值），需要确保基础一致，
  #要用与第0世代的差值作为进展，再用方程同时拟合所有数据再去分段
  #线性模型的系数需要相对一个基准来解释，所以分段比较时需要在一个统一基准下
  rgv10 = calrgv(gv = gv, ngen = 10)
  rgv20 = calrgv(gv = gv, ngen = 20)
  rav10 = calrgv(gv = Avalue, ngen = 10)
  rav20 = calrgv(gv = Avalue, ngen = 20)
  rdv10 = calrgv(gv = Dvalue, ngen = 10)
  rdv20 = calrgv(gv = Dvalue, ngen = 20)
  pgv_NP<- ONBSfinalg11$pgv_NP
  pgv_MP<- ONBSfinalg11$pgv_MP
  pgvrat= ((pgv_MP - pgv_NP)/pgv_NP)*100
  prel<- ONBSfinalg11$prel_mpvsnp
  nCoancestor<- ONBSfinalg11$nCoancestor
  identicalp<- ONBSfinalg11$identicalp
  
  sameparents = ONBSfinalg11$sameparents/150
  samecandidates = ONBSfinalg11$samecandidates
  if(!pname%in%c("CS","CP")){
    program<- paste(pname,g,sep = "") 
  }else{
    program<- pname
  }
  dt<- as.data.table(cbind(program,gg,gv,inb,inb_plink,inbprate,inbprates,inbrate,ld,Ne,Va,Vg,Vd,pgv_MP,pgv_NP,prel,nCoancestor,identicalp,pgvrat,conrate,conrate10,Avalue,Dvalue,sameparents,samecandidates,
                           rgv10,rgv20,rdv10,rdv20,rav10,rav20))
  return(dt)
}
nCycle=10
nGeneration=20
#nGeneration=10
#idc <- 0.325
idc = 0
for(g in c(1,3,5)){
  for (r in 1:nCycle) {
    ONBSfinalg11<- fread(paste("ONBSrealg",g,r,".csv",sep = ""),sep = ",")
    if(r==1){
      dt1<- abstractdata(ONBSfinalg11 = ONBSfinalg11,generation = 20,pname = "OS")
    }else{
      dt1<- rbind(dt1,abstractdata(ONBSfinalg11 = ONBSfinalg11,generation = 20,pname = "OS"))  
    }
    
  }
  if(g==1){
    dtfinalos<- dt1
  }else{
    dtfinalos<- rbind(dtfinalos,dt1)
  }
}

mymodel<- aov(gv~program,data=dtfinalos[program%in%c("OS1","OS3","OS5")&gg==20,])
summary(mymodel)
###########op

for(g in c(1,3,5)){
  for (r in 1:nCycle) {
    ONBSfinalg11<- fread(paste("ONBSrealgp",g,r,".csv",sep = ""),sep = ",")
    if(r==1){
      dt1<- abstractdata(ONBSfinalg11 = ONBSfinalg11,generation = 20,pname = "OP")
    }else{
      dt1<- rbind(dt1,abstractdata(ONBSfinalg11 = ONBSfinalg11,generation = 20,pname = "OP"))
    }
  }
  if(g==1){
    dtfinalop<- dt1
  }else{
    dtfinalop<- rbind(dtfinalop,dt1)
  }
}
mymodel<- aov(gv~program,data=dtfinalop[program%in%c("OP1","OP3","OP5")&gg==20,])
summary(mymodel)


#cs
for (r in 1:nCycle) {
  ONBSfinalg11<- fread(paste("CNBSreal",1,r,".csv",sep = ""),sep = ",")
  if(r==1){
    dtfinalcs = abstractdata(ONBSfinalg11 = ONBSfinalg11,generation = 20,pname = "CS")
  }else{
    dtfinalcs<- rbind(dtfinalcs,abstractdata(ONBSfinalg11 = ONBSfinalg11,generation = 20,pname = "CS")) 
  }
  
}


#cp

for (r in 1:nCycle) {
  
  ONBSfinalg11<- fread(paste("CNBSreal",2,r,".csv",sep = ""),sep = ",")
  if(r==1){
    dtfinalcp<- abstractdata(ONBSfinalg11 = ONBSfinalg11,generation = 20,pname = "CP")
  } else{
    dtfinalcp <- rbind(dtfinalcp,abstractdata(ONBSfinalg11 = ONBSfinalg11,generation = 20,pname = "CP"))
  }   
  
  
}

mymodel<- aov(gv~program,data=dtfinalos)
summary(mymodel)
dtfinal<- rbind(dtfinalos,dtfinalop,dtfinalcs,dtfinalcp)
fwrite(dtfinal,"dtallfinalac.csv",sep = ",")

