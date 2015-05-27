##s020_estudio_cel_files.r
##2011-03-11 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
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
system.time (load (file = file.path (.job$datadir, "data_processed", "raw_data.RData")))

dim (datos)

summary (datos)

medianas <- apply (datos, 2, median)
medianas

###SALVAMOS
save (list = "medianas", file = file.path (.job$datadir, "data_processed", "raw_data_medians.RData"))

################################################################################

misfilas <- sample (1:nrow (datos), size = 10000)

datos <- datos[misfilas,]

gc ()
gc ()

colores <- rep (c("red", "blue"), each = 15)
colores

nombres <- colnames (datos)
nombres <- sub (".CEL", "", nombres)
nombres

pdf (file = file.path (.job$plotsdir, "boxplots.pdf"))
boxplot (datos, log = "y", las = 3, col = colores, names = nombres)
dev.off ()
     

###EXIT
warnings ()
sessionInfo ()
q ("no")
