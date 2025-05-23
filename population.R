#creadt base poopulation and G0 population
#Historical effective population size 
#histNe=c(305,447,509,588,676,737,842,965,1115,1317,1554,1887,2195,2651,3266,3837,4495,4978,5619,6427,7129,7551,7568,8517,8253)
histNe=c(1115,1317,1554,1887,2195,2651,3266,3837,4495,4978,5619,6427,7129,7551,7568,8517,8253)
#histGen=c(13,15,17,20,23,27,32,38,45,54,65,80,98,122,151,188,238,294,357,454,555,708,784,877,952)
histGen=c(45,54,65,80,98,122,151,188,238,294,357,454,555,708,784,877,952)
BaseNe=1000
MaCSeNFlags = ""
if(length(histNe)>0){
  histNe = histNe/BaseNe
  histGen = histGen/(4*BaseNe)
  for(i in 1:length(histNe)){
    MaCSeNFlags = paste(MaCSeNFlags,"-eN",
                        histGen[i],histNe[i])
  }
}
ChrSize = (2.6 * 10^9) / 44
MutRate = 2.5E-7#突变率
RecRate = 1.67E-8#重组率，RecRate=GenLen(遗传距离)/ChrSize
founderPop = runMacs(nInd = 1000,
                     nChr = 44,
                     segSites = 2700,
                     manualCommand = paste(as.integer(ChrSize),
                                           "-t", MutRate * 4 * BaseNe,
                                           "-r", RecRate * 4 * BaseNe,
                                           MaCSeNFlags),
                     manualGenLen = RecRate * ChrSize)

#easy population
SP = SimParam$new(founderPop)#
SP$restrSegSites(minQtlPerChr = 100,minSnpPerChr = 1250,overlap = FALSE)
SP$addTraitAD(nQtlPerChr = 100,mean = 20,var = 5.39,meanDD = 0.12,varDD = 0.2)
SP$setVarE(h2=0.41)#遗传力和偏差
SP$addSnpChip(nSnpPerChr = 1250,minSnpFreq = 0.05)#55kSNP芯片
SP$setSexes("yes_sys")#按照一个雌一个雄来分配个体的性别
pop_founder = newPop(founderPop, simParam=SP)
pop <- selectCross(pop_founder,
                   nFemale = nFemale,nMale = nMale,
                   nCrosses = nCrosses,nProgeny = nProgeny,
                   use = "rand",
                   simParam = SP)

rm(.Random.seed)

save.image(file = paste(founderdir,"/","G0",r,".rda",sep=""))
