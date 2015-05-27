##r010_data_checking.r
##2011-02-03 dmontaner@cipf.es
##2011-03-11 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Comprobamos los datos del pedigree

### 14 individuos:
##                   11 antiguos
##                    3 nuevos - Con uno de los arrays por DUPLICADO: el 907
## 15 arrays

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 2.12.1 (2010-12-16)"
library (affy); packageDescription ("affy", fields = "Version") #"1.28.0"
library (affyio); packageDescription ("affyio", fields = "Version") #"1.18.0"
try (source (".job.r")); try (.job)


###DATOS
setwd (file.path (.job$datadir, "data_raw", "snp"))
dir ()

ficheros <- dir ()
ficheros ## ok son unicos


array.type <- unique (sapply (ficheros, function (x) read.celfile.header (x)$cdfName))
array.type

################################################################################

patient.id <- sub (pattern = ".CEL", replacement = "", ficheros)
patient.id <- sub (pattern = "PP ", replacement = "", patient.id)
patient.id <- sub (pattern = "PP", replacement = "", patient.id)
patient.id <- sub (pattern = "DIC", replacement = "", patient.id)

sinfo <- as.data.frame (list (filename = ficheros, samplename = patient.id), stringsAsFactors = FALSE)
sinfo

sapply (sinfo, class)

###EXIT
warnings ()
sessionInfo ()
q ("no")
