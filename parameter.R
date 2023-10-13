#program name
pname = "os" #os cp op 

#genetic evolution method
gem = "ssGBLUP"#PBLUP

#gradients 10% 30% 50%
#the number of parents from MP account for all parents
#15/150 = 10% 45/150 = 30%  75/150 = 50% 
nMP = 75 # 45 75 

#numeber of female
nFemale = 100
nFemale_MP = 20

#number of male
nMale = 50 
nMale_MP = 10

#the number of family
nCrosses = 100 
nCrosses_MP = 20

#the number of perogeny of per family
nProgeny = 200 
nProgeny_MP = 10000

#breeding length
nGeneration = 20

#the number of individuals in each generation
nIndPerGeneration = nProgeny*nCrosses 

#the number of candidate family
ncf = 50

#the number of genotyped individuals of each candidate family per generation
ngi = 21

#the number of candidate individuals of each candidate family per generation
nci = 6

#the number of candidate family for constructing MP
ncf_MP = 8

#the number of candidate individuals in each candidate family for construting MP
nci_MP = 12

#the overlap between NP parents for construcing the next generation MP and NP
nooverlap = FALSE


ChrSize = (2.6 * 10^9) / 44

MutRate = 2.5E-7

RecRate = 1.67E-8

#You can give any value for "contrast",
# because the program would check whether the "contrast" exists in the globe environment
# and to take the corresponding precedure
#if you don't want to contrast the accuracy between genomic and pedigree-based selection
#you should remove the parameters
contrast = TRUE