##u010_estudio_eroresm_merlin.r
##2011-03-12 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Vemos los ERRORES Mendelianos detectados por el pedstats del Merlin

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 2.12.1 (2010-12-16)"
#library (kinship); packageDescription ("kinship", fields = "Version") #"1.1.0-23"

try (source (".job.r")); try (.job)

###DATOS
load (file = file.path (.job$datadir, "data_processed", "ped.RData"))
load (file = file.path (.job$datadir, "data_processed", "datos_filtrados.RData"))
load (file = file.path (.job$datadir, "data_processed", "errores_mendel_plink.RData"))

dim (ped)
dim (datos)

table (colnames (datos) %in% rownames (ped)) ##OK

################################################################################


##Datos Errores Mendelianos - pedstats Merlin

lineas <- readLines (con = file.path (.job$datadir, "data_processed", "merlin", "salida_pedstats_0.txt"))
length (lineas)
lineas[100:110]

lineas  <- grep  (pattern = "- Fam", lineas, ignore.case = FALSE, value = TRUE)
length (lineas)
lineas[1:3]

men <- matrix (unlist (strsplit (lineas, split = " - ")), ncol = 2, byrow = TRUE)
men[1:10,]
dim (men)

men.merlin <- unique (men[,1])
length (men.merlin) #44377


men[1,]
datos["rs1705415",]

men[100,]
datos["rs11153413",]

################################################################################

##comparamos con el Plink

length (men.plink)
length (men.merlin)

length (intersect (men.plink, men.merlin))

setdiff (men.plink, men.merlin) #OK todos los detectados por PLINK estan tambien en MERLIN

solo.merlin <- setdiff (men.merlin, men.plink) #OK todos los detectados por PLINK estan tambien en MERLIN
length (solo.merlin)

table (annot[solo.merlin,"Chromosome"]) ##OBS: todos los no detectados estan en el X
table (annot[men.plink,"Chromosome"])   ##OBS: pero el plink SI detecta algunos en el X


touse <- men[,1] %in% solo.merlin

men[touse,][1:10,]
datos["rs5985612",]
datos["rs845445",]
datos["rs5974328",]
annot["rs845445",]


###SAL
save (list = "men.merlin", file = file.path (.job$datadir, "data_processed", "errores_mendel_merlin.RData"))


###EXIT
warnings ()
sessionInfo ()
q ("no")
