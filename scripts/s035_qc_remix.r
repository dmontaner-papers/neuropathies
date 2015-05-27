##s035_qc_remix.r
##2011-03-11 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Chequeamos la calidad de los resultados del APT

##HAY MUCHISIMOS MISSING CALLS

##Como era de esperar cuanto mas baja es la señal mas missing calls hay

### 14 individuos:
##                   11 antiguos
##                    3 nuevos - Con uno de los arrays por DUPLICADO
## 15 arrays

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 2.12.1 (2010-12-16)"
#library (kinship); packageDescription ("kinship", fields = "Version") #"1.1.0-23"

try (source (".job.r")); try (.job)


###PED
load (file = file.path (.job$datadir, "data_processed", "ped.RData"))
ped


###DATOS
load (file = file.path (.job$datadir, "data_processed", "quality_control.RData"))
load (file = file.path (.job$datadir, "data_processed", "raw_data_medians.RData"))
ls ()

medianas <- medianas[1:15]

################################################################################

a.id <- sub (pattern = ".CEL", replacement = "", rownames (array.qc))
a.id <- sub (pattern = "\\.", replacement = "", a.id)

b.id <- sub (pattern = ".CEL", replacement = "", names (medianas))
b.id <- sub (pattern = " ", replacement = "", b.id)

table (a.id == b.id) ##OK COINCIDEN

datos.qc <- cbind (array.qc, medianas)
rownames (datos.qc) <- b.id

datos.qc

plot (datos.qc[,"medianas"], datos.qc[,"na.rate"])
points (datos.qc[,"medianas"], datos.qc[,"na.rate.hapmap"], col = "red")

###EXIT
dev.off ()
warnings ()
sessionInfo ()
q ("no")
