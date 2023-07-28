parameter_txt_v <- c(
  # 指定数据文件名称，默认名称，不能修改
  "DATAFILE",
  "pheno.txt",
  # *待分析性状所在列位置
  "TRAITS",
  4,
  # 追加到输出文件中的列位置，可以为空
  "FIELDS_PASSED TO OUTPUT",
  "",
  # 性状加权，可以为空  
  "WEIGHT(S)",
  "",
  # 残差方差
  "RESIDUAL_VARIANCE",
  (varP(pop_founder)-varG(pop_founder)),
  "EFFECT",
  # *个体编号列位置，crossclassified, 字符类型
  "1 cross alpha",
  # 随机效应
  "RANDOM",
  # 个体动物模型
  "animal",
  # 定义系谱文件名称，默认名称，不能修改
  "FILE",
  "ped.txt",
  # *定义个体、父本和母本在系谱文件中的列位置
  "FILE_POS",
  "1 2 3 0 0",
  # 系谱深度，0表示追溯到奠基者世代
  "PED_DEPTH",
  0,
  # 加性遗传方差
  "(CO)VARIANCES",
  varG(pop_founder)
)
# 生成参数文件renumf90.par，utf-8格式
parameter_con <- file("renumf90.par")
writeLines(text = paste(parameter_txt_v, collapse = "\n"),
           con = parameter_con)
close(parameter_con)

# 调用renumf90程序
system2("renumf90", args = "renumf90.par", stdout="renumf90.log", wait = TRUE)

system2("blupf90", args = "renf90.par", stdout="blupf90.log",  wait = TRUE)
solution<- fread(input="solutions",sep = ' ')
setorder(solution,level)
rename<- fread(input = "renadd01.ped",sep = " ")
trueoder<- rename[,.(V1,V10)]
setorder(trueoder,V1)
ebv<- solution[`trait/effect`==1,solution]
trueorder<- cbind(trueoder,ebv)
setorder(trueorder,V10)
trueorder<- as.data.table(trueorder)
if(nrow(trueorder)==nrow(alldata)){
 alldata[,ebv:=trueorder[,ebv]]
}else if(nrow(trueorder)==nrow(data_allMP)){
 data_allMP[,ebv:=trueorder[,ebv]]
}else{
  stop("no data.table to assign ebv in equal length")
}
