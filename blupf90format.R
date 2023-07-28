pheno_file_s <- "pheno.csv"

ped_file_s <- "ped.csv"
snp_file_s <- "geno_selectedparents.txt"

pheno_dt <-fread(
  pheno_file_s,
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE,
  na.strings = "NA"
)
ped_dt <-
  fread(
    ped_file_s,
    sep = ",",
    header = TRUE,
    stringsAsFactors = FALSE,
    na.strings = "NA"
  )

fwrite(
  pheno_dt,
  file = "pheno.txt",
  na = "0",
  append = FALSE,
  sep = " ",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
fwrite(
  ped_dt[,1:3],
  file = "ped.txt",
  append = FALSE,
  sep = " ",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  na="0"
)
