# Open-nucleus-breeding-system

Author:
Ziyi Kang kangziyi1998@163.com;
Sheng Luan luansheng@ysfri.ac.cn

Description

This project presents a simulation study published in 2023. The R script provided here facilitates the simulation of breeding schemes for ONBS and CNBS, while also calculating essential population parameters.

Requirements

Before running the simulation program, ensure that the following software and R packages are installed:

AlphaSimR (Version 1.4.2)
BLUPF90
Plink 1.9
AlphaMate (v0.2.0)
GCTA64
Additional R packages listed in package.R
Please make sure you have the "devtools" package installed to fetch any required packages from GitHub or other code repositories.

Getting Started

Place the "Open-nucleus-breeding-system-main" in the work directory which you want ("wdyw"). Please do not change the name of "Open-nucleus-breeding-system-main"
Set the desired parameters in the parameter.R file.
Execute the loadfile.R script to generate founder populations. The number of founder populations is determined by the repeat times for each breeding scheme. Adjust the loops accordingly in runme.R and loadfile.R.
Finally, load runme.R to produce the required breeding schemes. Be prepared for longer execution times, and consider distributing the load by running the scripts on multiple files simultaneously.
If you want to contrat the accuracy of EBV for model include inbreeding covariate versus not include, please "cd "wdyw"/Open-nucleus-breeding-system-main/model_inbreeding_covariate" and then load the runme.R.

Note

Two files, ocs_SC.R and ocsBLUP_sc.R, have been modified to change "gender" to "sex". Additionally, ocsBLUP_sc.R includes a function to generate the A matrix, replacing the G matrix. ocs_SC.R and ocsBLUP_sc.R are sourced from the code of AlphaLearn. These scripts have been included in this project to aid in specific functionalities. Please refer to the original AlphaLearn documentation or source code for additional information regarding these scripts.
