if(pname%in%c("os","op")&!exists("pop_MP")){
  if(gem=="ssGBLUP"){
  pop <- ocs(
  pop = pop_nucleus_candidate,
  nCrosses = nCrosses,
  nProgenyPerCross = nProgeny,
  nFemalesMax = nFemale,
  nMalesMax = nMale,
  equalizeFemaleContributions =TRUE,
  equalizeMaleContributions = TRUE,
  targetDegree = 45,
  use = "ebv"
)
  if(nooverlap==TRUE){
  MP_fam_ind = MP_fam_ind[setdiff(MP_fam_ind@id,unique(c(pop@mother,pop@father)))]
  }
  pop_MP <- ocs(
  pop=MP_fam_ind,
  nCrosses = nCrosses_MP,
  nProgenyPerCross = nProgeny_MP,
  nFemalesMax = nFemale_MP,
  nMalesMax = nMale_MP,
  equalizeFemaleContributions = TRUE,
  equalizeMaleContributions = TRUE,
  targetDegree = 45,
  use = "ebv"
)
  }else{
    pop <- ocsBLUP(
  pop = pop_nucleus_candidate,
  pedmore = ped,
  number = 9,
  pedkf = FALSE,
  nCrosses = nCrosses,
  nProgenyPerCross = nProgeny,
  nFemalesMax = nFemale,
  nMalesMax = nMale,
  equalizeFemaleContributions =TRUE,
  equalizeMaleContributions = TRUE,
  targetDegree = 45,
  use = "ebv"
)
  if(nooverlap==TRUE){
  MP_fam_ind = MP_fam_ind[setdiff(MP_fam_ind@id,unique(c(pop@mother,pop@father)))]
  }
  pop_MP <- ocsBLUP(
  pop = MP_fam_ind,
  pedmore = ped,
  number = 9,
  pedkf = FALSE,
  nCrosses = nCrosses_MP,
  nProgenyPerCross = nProgeny_MP,
  nFemalesMax = nFemale_MP,
  nMalesMax = nMale_MP,
  equalizeFemaleContributions = TRUE,
  equalizeMaleContributions = TRUE,
  targetDegree = 45,
  use = "ebv"
)

  }
 
}else if(pname%in%c("os","op") & g<nGeneration){
    pop_nucleus_candidate_mixF <- c(pop_MP_candidate_F,nucleus_fam_ind_M)

    pop_nucleus_candidate_mixM <-c(pop_MP_candidate_M,nucleus_fam_ind_F)

    if(gem=="ssGBLUP"){
     mix_FemalePart <- ocs(
        pop = pop_nucleus_candidate_mixF,
        nCrosses = nMP*(2/3),
        nProgenyPerCross = 200,
        nFemalesMax = nMP*(2/3),
        nMalesMax = nMP*(1/3),
        equalizeFemaleContributions = TRUE,
        equalizeMaleContributions = TRUE,
        targetDegree = 45,
        use = "ebv"
      )
      #建立扩繁父本x核心母本的混合家系
      mix_MalePart <- ocs(
        pop = pop_nucleus_candidate_mixM,
        nCrosses = nMP*(2/3),
        nProgenyPerCross = 200,
        nFemalesMax = nMP*(2/3),
        nMalesMax = nMP*(1/3),
        equalizeFemaleContributions = TRUE,
        equalizeMaleContributions = TRUE,
        targetDegree = 45,
        use = "ebv"
      )

      if(nMP!=75){
        mix_parents_nucleus <- c(unique(mix_FemalePart@father),unique(mix_MalePart@mother))

        pop_nucleus_candidate_removeMix <- pop_nucleus_candidate[setdiff(pop_nucleus_candidate@id,mix_parents_nucleus)]

        pop_nucleus <- ocs(
        pop = pop_nucleus_candidate_removeMix,
        nCrosses = (100-(nMP*(4/3))),
        nProgenyPerCross = 200,
        nFemalesMax = (100-(nMP*(4/3))),
        nMalesMax = ((100-(nMP*(4/3)))/2),
        equalizeFemaleContributions = TRUE,
        equalizeMaleContributions = TRUE,
        targetDegree = 45,
        use = "ebv"
      )
      pop <- c(mix_FemalePart,mix_MalePart,pop_nucleus)
      }else{
      pop <- c(mix_FemalePart,mix_MalePart)
      }

      if(nooverlap==TRUE){
  MP_fam_ind = MP_fam_ind[setdiff(MP_fam_ind@id,unique(c(pop@mother,pop@father)))]
}
pop_MP <- ocs(
  pop=MP_fam_ind,
  nCrosses = nCrosses_MP,
  nProgenyPerCross = nProgeny_MP,
  nFemalesMax = nFemale_MP,
  nMalesMax = nMale_MP,
  equalizeFemaleContributions = TRUE,
  equalizeMaleContributions = TRUE,
  targetDegree = 45,
  use = "ebv"
)
  }else{
     mix_FemalePart <- ocsBLUP(
        pop = pop_nucleus_candidate_mixF,
        pedmore = pedmore,
        number = 9,
        pedkf = TRUE,
        nCrosses = nMP*(2/3),
        nProgenyPerCross = 200,
        nFemalesMax = nMP*(2/3),
        nMalesMax = nMP*(1/3),
        equalizeFemaleContributions = TRUE,
        equalizeMaleContributions = TRUE,
        targetDegree = 45,
        use = "ebv"
      )
      #建立扩繁父本x核心母本的混合家系
      mix_MalePart <- ocsBLUP(
        pop = pop_nucleus_candidate_mixM,
        pedmore = pedmore,
        number = 9,
        pedkf = TRUE,
        nCrosses = nMP*(2/3),
        nProgenyPerCross = 200,
        nFemalesMax = nMP*(2/3),
        nMalesMax = nMP*(1/3),
        equalizeFemaleContributions = TRUE,
        equalizeMaleContributions = TRUE,
        targetDegree = 45,
        use = "ebv"
      )

      if(nMP!=75){
        mix_parents_nucleus <- c(unique(mix_FemalePart@father),unique(mix_MalePart@mother))

        pop_nucleus_candidate_removeMix <- pop_nucleus_candidate[setdiff(pop_nucleus_candidate@id,mix_parents_nucleus)]

        pop_nucleus <- ocsBLUP(
        pop = pop_nucleus_candidate_removeMix,
        pedmore = ped,
        number = 9,
        pedkf = FALSE,
        nCrosses = (100-(nMP*(4/3))),
        nProgenyPerCross = 200,
        nFemalesMax = (100-(nMP*(4/3))),
        nMalesMax = ((100-(nMP*(4/3)))/2),
        equalizeFemaleContributions = TRUE,
        equalizeMaleContributions = TRUE,
        targetDegree = 45,
        use = "ebv"
      )
      pop <- c(mix_FemalePart,mix_MalePart,pop_nucleus)
      }else{
      pop <- c(mix_FemalePart,mix_MalePart)
      }

      if(nooverlap==TRUE){
  MP_fam_ind = MP_fam_ind[setdiff(MP_fam_ind@id,unique(c(pop@mother,pop@father)))]
}
pop_MP <- ocsBLUP(
  pop=MP_fam_ind,
  pedmore = ped,
  pedkf = FALSE,
  number = 9,
  nCrosses = nCrosses_MP,
  nProgenyPerCross = nProgeny_MP,
  nFemalesMax = nFemale_MP,
  nMalesMax = nMale_MP,
  equalizeFemaleContributions = TRUE,
  equalizeMaleContributions = TRUE,
  targetDegree = 45,
  use = "ebv"
)
  }

}else if(pname%in%c("cs","cp") & g<nGeneration){
  if(gem=="ssGBLUP"){
        pop <- ocs(
  pop = pop_nucleus_candidate,
  nCrosses = nCrosses,
  nProgenyPerCross = nProgeny,
  nFemalesMax = nFemale,
  nMalesMax = nMale,
  equalizeFemaleContributions =TRUE,
  equalizeMaleContributions = TRUE,
  targetDegree = 45,
  use = "ebv"
)
  }else{
pop <- ocsBLUP(
  pop = pop_nucleus_candidate,
  pedmore = ped,
  number = 9,
  pedkf = FALSE,
  nCrosses = nCrosses,
  nProgenyPerCross = nProgeny,
  nFemalesMax = nFemale,
  nMalesMax = nMale,
  equalizeFemaleContributions =TRUE,
  equalizeMaleContributions = TRUE,
  targetDegree = 45,
  use = "ebv"
)
  }
}
