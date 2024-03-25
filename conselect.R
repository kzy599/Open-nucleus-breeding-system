    pop@ebv = abstractebv_blupf90(pop = pop, dt = alldata )

    initialFam<- selectFam(pop,nFam = ncf,trait = 1,use = "ebv",simParam = SP)

    initialInd_F <- selectWithinFam(initialFam,nInd = ngi*(2/3),trait = 1,use = "ebv",simParam = SP,sex  ="F")

    initialInd_M <- selectWithinFam(initialFam,nInd = ngi*(1/3),trait = 1,use = "ebv",simParam = SP,sex  ="M")

    initialInd <- c(initialInd_F,initialInd_M)

    if(exists("pop_MP")){

      pop_MP_candidate_F_con <- selectWithinFam(pop_MP,nInd = 4,trait = 1,use = "pheno",simParam = SP,sex = "F")

      pop_MP_candidate_M_con <- selectWithinFam(pop_MP,nInd = 2,trait = 1,use = "pheno",simParam = SP,sex = "M")
      
      pop_MP_candidate_con <- c(pop_MP_candidate_M_con,pop_MP_candidate_F_con)
    }

  if(exists("pop_MP_candidate_con")){
  
  pop_MP_candidate_con@ebv = pop_MP_candidate_con@pheno

  pop_MP_candidate_M_con@ebv<- pop_MP_candidate_M_con@pheno

  pop_MP_candidate_F_con@ebv<- pop_MP_candidate_F_con@pheno
}


#seleced top 6 individuals (4 females and 2 males) of genotyped individuals in each top 50 families
nucleus_fam_ind_F_con <- selectWithinFam(initialInd,nInd = nci*(2/3),trait = 1,use = "ebv",simParam = SP,sex = "F")
nucleus_fam_ind_M_con <- selectWithinFam(initialInd,nInd = nci*(1/3),trait = 1,use = "ebv",simParam = SP,sex = "M")
pop_nucleus_candidate_con <- c(nucleus_fam_ind_M_con,nucleus_fam_ind_F_con)

    
if(exists("pop_MP_candidate_con")){
      output$accuracy_mp_con[g+1] <- cor(pop_MP_candidate_con@gv,pop_MP_candidate_con@ebv)
    }

