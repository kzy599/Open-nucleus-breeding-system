#!~/bin/Rscript
#globe script
for(r in c(1:10)){
rm(list=setdiff(ls(),"r"))
gc()
#prepare for simulation
#load R package
source("package.R")

#globe parameters
source("parameter.R")

#creat dir for founder and G0 population
if(r==1){
dir.creat(paste(wdyw,"FG0",sep=""))
}

founderdir= paste(wdyw,"FG0",sep="")

#creat base and G0 population
source("population.R")
}