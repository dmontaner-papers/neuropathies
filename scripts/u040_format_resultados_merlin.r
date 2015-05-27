##u040_format_resultados_merlin.r
##2011-03-14 dmontaner@cipf.es
##2011-05-03 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Formateamos los resultados de los analisis del Merlin.

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 2.12.1 (2010-12-16)"
#library (kinship); packageDescription ("kinship", fields = "Version") #"1.1.0-23"

try (source (".job.r")); try (.job)

###DATOS
load (file = file.path (.job$datadir, "data_processed", "datos_por_cromosoma.RData"))
ls ()

################################################################################


##LEEMOS: afectados MISSING COMO MISSING
resm0 <- list()
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  
  res <- read.table (file = file.path (.job$datadir, "data_processed", "merlin_0", cro, "merlin-parametric.tbl"),
                     header = TRUE, sep = "\t", quote = "", as.is = TRUE)
  
  map <- read.table (file = file.path (.job$datadir, "data_processed", "merlin_0", cro, "datos.map"),
                     header = FALSE, sep = "\t", quote = "", as.is = TRUE)
  orden <- order (map[,3])
  map <- map[orden,]
  
  if (any (res[,"LABEL"] != round (map[,3], 3))) stop ("PROBLEMA")
  
  res[,c("snp", "GeneticMap")] <- map[,2:3]
  res[,c("deCODE.genetic.map", "PhysicalPosition")] <- lanot[[cro]][res[,"snp"], c ("deCODE.genetic.map", "Physical.Position")]

  resm0[[cro]] <- res

  print (resm0[[cro]][1:3,])
}

t (sapply (resm0, dim))
t (sapply (resm0, function (x) summary (x[,"LOD"])))
t (sapply (resm0, function (x) sapply (x, class)))

save (list = "cromosomas", "resm0", file = file.path (.job$datadir, "data_processed", "res_merlin0.RData"))
################################################################################


##LEEMOS: afectados MISSING COMO SANOS
resm1 <- list()
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  
  res <- read.table (file = file.path (.job$datadir, "data_processed", "merlin_1", cro, "merlin-parametric.tbl"),
                     header = TRUE, sep = "\t", quote = "", as.is = TRUE)
  
  map <- read.table (file = file.path (.job$datadir, "data_processed", "merlin_1", cro, "datos.map"),
                     header = FALSE, sep = "\t", quote = "", as.is = TRUE)
  orden <- order (map[,3])
  map <- map[orden,]
  
  if (any (res[,"LABEL"] != round (map[,3], 3))) stop ("PROBLEMA")
  
  res[,c("snp", "GeneticMap")] <- map[,2:3]
  res[,c("deCODE.genetic.map", "PhysicalPosition")] <- lanot[[cro]][res[,"snp"], c ("deCODE.genetic.map", "Physical.Position")]

  resm1[[cro]] <- res

  print (resm1[[cro]][1:3,])
}

t (sapply (resm1, dim))
t (sapply (resm1, function (x) summary (x[,"LOD"])))
t (sapply (resm1, function (x) sapply (x, class)))

save (list = "cromosomas", "resm1", file = file.path (.job$datadir, "data_processed", "res_merlin1.RData"))
################################################################################


##LEEMOS: afectados MISSING COMO ENFERMOS
resm2 <- list()
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  
  res <- read.table (file = file.path (.job$datadir, "data_processed", "merlin_2", cro, "merlin-parametric.tbl"),
                     header = TRUE, sep = "\t", quote = "", as.is = TRUE)
  
  map <- read.table (file = file.path (.job$datadir, "data_processed", "merlin_2", cro, "datos.map"),
                     header = FALSE, sep = "\t", quote = "", as.is = TRUE)
  orden <- order (map[,3])
  map <- map[orden,]
  
  if (any (res[,"LABEL"] != round (map[,3], 3))) stop ("PROBLEMA")
  
  res[,c("snp", "GeneticMap")] <- map[,2:3]
  res[,c("deCODE.genetic.map", "PhysicalPosition")] <- lanot[[cro]][res[,"snp"], c ("deCODE.genetic.map", "Physical.Position")]

  resm2[[cro]] <- res

  print (resm2[[cro]][1:3,])
}

t (sapply (resm2, dim))
t (sapply (resm2, function (x) summary (x[,"LOD"])))
t (sapply (resm2, function (x) sapply (x, class)))

save (list = "cromosomas", "resm2", file = file.path (.job$datadir, "data_processed", "res_merlin2.RData"))
################################################################################


##LEEMOS: afectados MISSING corregidos segun las indicaciones de Eduardo
resm3 <- list()
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  
  res <- read.table (file = file.path (.job$datadir, "data_processed", "merlin_3", cro, "merlin-parametric.tbl"),
                     header = TRUE, sep = "\t", quote = "", as.is = TRUE)
  
  map <- read.table (file = file.path (.job$datadir, "data_processed", "merlin_3", cro, "datos.map"),
                     header = FALSE, sep = "\t", quote = "", as.is = TRUE)
  orden <- order (map[,3])
  map <- map[orden,]
  
  if (any (res[,"LABEL"] != round (map[,3], 3))) stop ("PROBLEMA")
  
  res[,c("snp", "GeneticMap")] <- map[,2:3]
  res[,c("deCODE.genetic.map", "PhysicalPosition")] <- lanot[[cro]][res[,"snp"], c ("deCODE.genetic.map", "Physical.Position")]

  resm3[[cro]] <- res

  print (resm3[[cro]][1:3,])
}

t (sapply (resm3, dim))
t (sapply (resm3, function (x) summary (x[,"LOD"])))
t (sapply (resm3, function (x) sapply (x, class)))

save (list = "cromosomas", "resm3", file = file.path (.job$datadir, "data_processed", "res_merlin3.RData"))
################################################################################


###EXIT
warnings ()
sessionInfo ()
q ("no")
