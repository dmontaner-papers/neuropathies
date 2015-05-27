##s020_read_cel_files.r
##2011-03-11 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Leemos los CEL files para hacer un control de calidad:
##Vemos como son las distribuciones de los ficheros CEL para revisar por que la calidad de los datos es tan mala

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 2.12.1 (2010-12-16)"
library (affy); packageDescription ("affy", fields = "Version") # "1.28.0"
#library (kinship); packageDescription ("kinship", fields = "Version") #"1.1.0-23"

try (source (".job.r")); try (.job)


###DATOS
ficheros.carmen <- dir (file.path (.job$datadir, "data_raw", "snp"), full.names = TRUE)
ficheros.hapmap <- dir (file.path (.job$datadir, "data_raw", "hapmap"), full.names = TRUE)
ficheros.hapmap <- sample (ficheros.hapmap, size = 15)

ficheros.carmen 
ficheros.hapmap


system.time (datos <- ReadAffy (filenames = c (ficheros.carmen, ficheros.hapmap)))


datos <- exprs (datos)
dim (datos)
colnames (datos)


###SALVAMOS
save (list = "datos", file = file.path (.job$datadir, "data_processed", "raw_data.RData"))


###EXIT
warnings ()
sessionInfo ()
q ("no")
