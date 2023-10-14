#!~/bin/Rscript
for(r in 1:10){

#clear environment
#consever r which account repeats
rm(list=setdiff(ls(),"r"))
gc()
#source parameters
source("parameter.R")

#source package
source("package.R")

nchar("Open_nucleus")
if(subwdyw)

setwd(substr(wdyw,1,nchar(wdyw)-1))

scrdir= getwd()

setwd(paste(wdyw,"FG0",sep=""))

#load repeat G0 population
load(paste("G0",r,".rda",sep=""))


if(r==1){
dir.creat(paste(wdyw,pname,nMP/150*100,sep=""))

scr= list.files(scrdir)
scr = scr[scr%flike%".R"]

file.copy(from = paste(scrdir,scr,sep = "/"),to = paste(wdyw,pname,nMP/150*100,sep=""))

}

setwd(paste(wdyw,pname,nMP/150*100,sep=""))

output <- data.frame(generation=0:(nGeneration),
                     mean_gv=numeric(nGeneration+1),
                     mean_gv_MP = numeric(nGeneration+1),
                     mean_pgv = numeric(nGeneration+1),
                     inbreeding = numeric(nGeneration+1),
                     inbreeding_plink = numeric(nGeneration+1),
                     LD=numeric(nGeneration+1),
                     LDscore=numeric(nGeneration+1),
                     Ne = numeric(nGeneration+1),
                     Va = numeric(nGeneration+1),
                     Vd = numeric(nGeneration+1),
                     Vp = numeric(nGeneration+1),
                     Vg = numeric(nGeneration+1),
                     Vg_MP = numeric(nGeneration+1),
                     genicVg = numeric(nGeneration+1),
                     genicVg_MP = numeric(nGeneration+1),
                     pheno=numeric(nGeneration+1),
                     ibdle = rep(calibd(pop_founder),(nGeneration+1)),
                     mean_gv_founder = rep(meanG(pop_founder),(nGeneration+1)),
                     mean_D_founder = rep(mean(needpara_founder$gv_d),(nGeneration+1)),
                     mean_A_founder = rep(mean(needpara_founder$gv_a),(nGeneration+1)),
                     mean_gvu_founder = rep(mean(needpara_founder$gv_mu),(nGeneration+1)),
                     d2infounder = rep(varD(pop_founder)/varP(pop_founder),(nGeneration+1)),
                     d2 = numeric(nGeneration+1),
                     d2_MP = numeric(nGeneration+1),
                     Va_MP = numeric(nGeneration+1),
                     genicVa_MP = numeric(nGeneration+1),
                     Vd_MP = numeric(nGeneration+1),
                     genicVd_MP = numeric(nGeneration+1),
                     Vp_MP = numeric(nGeneration+1),
                     h2_MP = numeric(nGeneration+1),
                     pheno_MP = numeric(nGeneration+1),
                     pgv = numeric(nGeneration+1),
                     pgv_NP = numeric(nGeneration+1),
                     pgv_MP = numeric(nGeneration+1),
                     accuracy_mp = numeric(nGeneration+1),
                     prel_mpvsnp = numeric(nGeneration+1),
                     genicVa = numeric(nGeneration+1),
                     genicVd = numeric(nGeneration+1),
                     h2 = numeric(nGeneration+1),
                     nCoancestor = numeric(nGeneration+1),
                     identicalp = numeric(nGeneration+1),
                     accuracy = numeric(nGeneration+1),
                     accuracy_con = numeric(nGeneration+1),
                     accuracy_mp_con = numeric(nGeneration+1),
                     sameparents = numeric(nGeneration+1),
                     samecandidates = numeric(nGeneration+1),
                     mean_D = numeric(nGeneration+1),
                     mean_A = numeric(nGeneration+1),
                     mean_gvu = numeric(nGeneration+1))
Geno <-vector(mode="list", length=21)
Aid <-vector(mode="list", length=21)
Pid <- vector(mode="list", length=21)
rt_ebv <- vector(mode="list",length = 21)
source("performbreeding.R") 
}
