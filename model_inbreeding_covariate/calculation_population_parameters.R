
pop_select <- selectWithinFam(pop = pop,nInd = 10,use = "rand",simParam = SP)

keep<- paste("Animal_",sprintf("%0*d",9,as.numeric(pop_select@id)),sep = "")

#pop_LD<- pop_select
pop_Inb <-mergePops(list(pop_founder,pop_select)) 

if(!exists("snpfre_v")){
  snp_012_dt = pullSnpGeno(pop_founder)
  snpfre_v <- apply(snp_012_dt, 2, function(x){
  single_snpfre_s <- sum(x, na.rm = TRUE)/(2*length(x))
  return(single_snpfre_s)
})
  #locincon= names(snpfre_v[snpfre_v > 0.5])
  #snpfre_v[snpfre_v > 0.5] = 1 - snpfre_v[snpfre_v > 0.5] 
}

snp_012_dt = pullSnpGeno(pop_select)
#snp_012_dt[,colnames(snp_012_dt)%in%locincon] = 2 - snp_012_dt[,colnames(snp_012_dt)%in%locincon]

inbreeding_plink = mean(calinbya2(snp_012_dt = snp_012_dt, snpfre = snpfre_v))

#writePlink(pop_Inb,baseName = "inb",chromLength = ChrSize,simParam = SP)

#system("plink --file inb  --ibc --chr-set 44 --out inb")

#plinkibc_dt <- fread("inb.ibc",header = TRUE,sep = "\t")

#inb= plinkibc_dt$Fhat3[match(pop_select@id,plinkibc_dt$IID)]

#inbreeding_plink <- mean(inb)

ped_calparameters <-  fread(
  input = "ped_full.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)
Pedig <- prePed(ped_calparameters, keep=keep)
Res   <- pedInbreeding(Pedig)
inbreeding <- mean(Res$Inbr[Res$Indiv %in% keep])
#

peddt = makeped(pop = pop_select)

mapdt = makemap()

fwrite(peddt,file="popLD.ped",col.names = FALSE,sep = "\t")

fwrite(mapdt, file = "popLD.map",col.names = FALSE,sep = " ")

system("plink --file popLD  --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out popLD")

system("gcta64 --bfile popLD --autosome-num 44 --ld-score --ld-wind 10000 --ld-rsq-cutoff 0.01 --out popLD")

LDfile<- fread(input="popLD.score.ld",sep = " ")
LD <- mean(LDfile$mean_rsq)
LDscore <- mean(LDfile$ldscore)

if (g == 0) {
  Ne <- 150
} else
  if (g > 0) {
    ped_calparameters  <-  fread(
      input = "ped_full.csv",
      sep = ",",
      header = TRUE,
      stringsAsFactors = FALSE
    )
    
    pedig <- prePed(ped_calparameters)
    pKin   <- pedIBD(pedig, keep.only = keep)
    Summary <- summary(pedig)
    id     <- keep
    x      <- Summary[Indiv %in% id]$equiGen
    N      <- length(x)
    n      <- (matrix(x, N, N, byrow = TRUE) + matrix(x, N, N, byrow = FALSE)) / 2
    deltaC <- 1 - (1 - pKin[id, id]) ^ (1 / n)
    Ne   <- 1 / (2 * mean(deltaC))
  }
   ped_calparameters<- visPedigree::tidyped(ped_calparameters,cand = keep)
   output$nCoancestor[g+1]<- length(ped_calparameters[Gen==1,Ind])
  if(exists("pop_parents_MP")){
      pop_Inb <-c(pop_founder,pop_parents_nucleus,pop_parents_MP)
      writePlink(pop_Inb,baseName = "inb",chromLength = ChrSize,simParam = SP)
      system("plink --file inb --geno 0.05 --mind 0.1 --maf 0.05  --make-bed --chr-set 44 --out inb")
      system("plink --bfile inb --chr-set 44 --make-rel square")
      Gmat<- fread(input = "plink.rel")
      Gmatname<- fread(input = "plink.rel.id")
      Gmat<- as.matrix(Gmat)
      rownames(Gmat) <- Gmatname$V2
      colnames(Gmat) <- Gmatname$V2
      output$prel_mpvsnp[g+1] <- mean(Gmat[rownames(Gmat)%in%pop_parents_nucleus@id,colnames(Gmat)%in%pop_parents_MP@id])
  }

    output$accuracy[g+1] <- cor(bv(pop_nucleus_candidate),pop_nucleus_candidate@ebv)
    if(exists("pop_MP_candidate")){
      output$accuracy_mp[g+1] <- cor(bv(pop_MP_candidate),pop_MP_candidate@ebv)
    }
    output$inbreeding[g+1] <- inbreeding
    output$inbreeding_plink[g+1] <- inbreeding_plink
    output$LD[g+1] <- LD
    output$LDscore[g+1] <- LDscore
    output$Ne[g+1] <- Ne
    output$d2[g+1] <- varD(pop,simParam = SP)/varP(pop)
    output$Va[g+1] <- varA(pop,simParam = SP)
    output$Vg[g+1] <- varG(pop)
    output$Vd[g+1] <- varD(pop,simParam = SP)
    output$Vp[g+1] <- varP(pop)
    output$pheno[g+1] <- mean(pop@pheno)
    output$mean_gv[g+1] <- mean(pop@gv)
    output$genicVa[g+1] <- genicVarA(pop)
    output$genicVd[g+1] <- genicVarD(pop)
    output$genicVg[g+1] <- genicVarG(pop)
    if(exists("pop_MP")){
    output$d2_MP[g+1] <- varD(pop_MP,simParam = SP)/varP(pop_MP)
    output$Va_MP[g+1] <- varA(pop_MP,simParam = SP)
    output$genicVa_MP[g+1] <- genicVarA(pop_MP)
    output$Vd_MP[g+1] <- varD(pop_MP,simParam = SP)
    output$Vp_MP[g+1] <- varP(pop_MP)
    output$Vg_MP[g+1] <- varG(pop_MP)
    output$genicVd_MP[g+1] <- genicVarD(pop_MP)
    output$genicVg_MP[g+1] <- genicVarG(pop_MP)
    output$h2_MP[g+1] <- varA(pop_MP)/varP(pop_MP)
    output$pheno_MP[g+1] <- mean(pop_MP@pheno)
    output$mean_gv_MP[g+1] <- mean(pop_MP@gv)
    output$identicalp[g+1] <- sum(unique(c(pop@mother,pop@father))%in%unique(c(pop_MP@mother,pop_MP@father)))
    }
    output$h2[g+1] = varA(pop)/varP(pop)
    if(!exists("mix_FemalePart")){
      output$pgv[g+1] = output$pgv_NP[g+1] = mean(pop_parents@gv)
    }else{
     output$pgv_NP[g+1] = mean(pop_parents_nucleus@gv)
     output$pgv_MP[g+1] = mean(pop_parents_MP@gv)
     output$pgv[g+1] = mean(c(pop_parents_MP@gv,pop_parents_nucleus@gv))
    }

    if(contrast){
      if(nrow(trueorder)==nrow(alldata)){
     
       output$corac[g+1] = cor(abstractebv(pop = pop, dt = alldata),abstractebvib(pop = pop, dt = alldata))

       output$coracand[g+1] = cor(abstractebv(pop = pop_nucleus_candidate, dt = alldata),abstractebvib(pop = pop_nucleus_candidate, dt = alldata))

       if(exists("pop_MP_candidate")){
        output$coracmp[g+1] = cor(abstractebv(pop = pop_MP_candidate, dt = alldata),abstractebvib(pop = pop_MP_candidate, dt = alldata))
        output$accuracyib_mp[g+1] = cor(abstractebvib(pop = pop_MP_candidate, dt = alldata),bv(pop_MP_candidate))
        }
       output$accuracyib[g+1] = cor(abstractebvib(pop = pop_nucleus_candidate, dt = alldata),bv(pop_nucleus_candidate))
     
      }else if(nrow(trueorder)==nrow(data_allMP)){

        output$corac[g+1] = cor(abstractebv(pop = pop, dt = data_allMP),abstractebvib(pop = pop, dt = data_allMP))

        output$accuracyib[g+1] = cor(abstractebvib(pop = pop_nucleus_candidate, dt = data_allMP),bv(pop_nucleus_candidate))

        output$coracand[g+1] = cor(abstractebv(pop = pop_nucleus_candidate, dt = data_allMP),abstractebvib(pop = pop_nucleus_candidate, dt = data_allMP))

       if(exists("pop_MP_candidate")){
        output$coracmp[g+1] = cor(abstractebv(pop = pop_MP_candidate, dt = data_allMP),abstractebvib(pop = pop_MP_candidate, dt = data_allMP))
        output$accuracyib_mp[g+1] = cor(abstractebvib(pop = pop_MP_candidate, dt = data_allMP),bv(pop_MP_candidate))
        }

        }
    }
