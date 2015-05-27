##s050_clean_data.r
##2011-03-11 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Filtramos los SNPs que vamos a usar en el Merlin

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 2.12.1 (2010-12-16)"
#library (kinship); packageDescription ("kinship", fields = "Version") #"1.1.0-23"

try (source (".job.r")); try (.job)

###DATOS
load (file = file.path (.job$datadir, "data_processed", "ped.RData"))
load (file = file.path (.job$datadir, "data_processed", "all_calls_limpias.RData"))
load (file = file.path (.job$datadir, "data_processed", "anotacion.RData"))
load (file = file.path (.job$datadir, "data_processed", "quality_control.RData"))

dim (ped)
dim (annot)
dim (datos)
dim (snp.qc)

annot[1:3,]
snp.qc[1:10,]

table (rownames (annot) == rownames (datos)) ##OK
table (rownames (datos) %in% rownames (snp.qc)) ##OK

snp.qc <- snp.qc[rownames (datos),]

################################################################################

###Filtramos por calidades

table (snp.qc[,"con.call"] > 10, snp.qc[,"como.hap.map"] > 10)
100 * table (snp.qc[,"con.call"] > 10, snp.qc[,"como.hap.map"] > 10) / nrow (snp.qc)

touse.qc <- snp.qc[,"con.call"] > 10 & snp.qc[,"como.hap.map"] > 10
table (touse.qc)

table (annot[,"Chromosome"] %in% c("---", "Y"))
touse.chr <- !annot[,"Chromosome"] %in% c("---", "Y")
table (touse.chr)

table (touse.qc, touse.chr)
100 * table (touse.qc, touse.chr) / length (touse.qc)


touse <- touse.qc & touse.chr
table (touse)  

annot <- annot[touse,]
datos <- datos[touse,]
snp.qc <- snp.qc[touse,]


###quitamos el snp repetido

dup <- duplicated (annot[,"dbSNP.RS.ID"])
table (dup)

annot <- annot[!dup,]
datos <- datos[!dup,]
snp.qc <- snp.qc[!dup,]

dim (annot)
dim (datos)
dim (snp.qc)



###derivamos distancias genomicas y eliminamos los SNPs que no las tienen
fun <- function (x) {
  y <- strsplit (x, split = " ")[[1]][1]
  return (y)
}
system.time (deCODE.genetic.map <- sapply (annot[,"Genetic.Map"], fun, USE.NAMES = FALSE))

deCODE.genetic.map[1:10]

annot[,"deCODE.genetic.map"] <- as.numeric (deCODE.genetic.map)

esna <- is.na (annot[,"deCODE.genetic.map"])
table (esna) #63

deCODE.genetic.map[esna]

table (annot[,"Chromosome"], esna) #todos en el cromosoma X

annot <- annot[!esna,]
datos <- datos[!esna,]
snp.qc <- snp.qc[!esna,]

dim (annot)
dim (datos)
dim (snp.qc)


##formateamos a numerica la posicion cromosomica
annot[,"posicion"] <- as.numeric (annot[,"Physical.Position"])
table (is.na (annot[,"posicion"]))

################################################################################

table (rownames (annot) == rownames (datos)) ##OK
table (rownames (datos) == rownames (snp.qc)) ##OK


###nombres de las filas
table (duplicated (annot[,"dbSNP.RS.ID"]))

rownames (annot) <- annot[,"dbSNP.RS.ID"]
rownames (datos) <- annot[,"dbSNP.RS.ID"]

################################################################################

datos <- datos[,-15] ##el ultimo array esta repetido

###nombres de las columnas
rownames (ped)
colnames (datos) ##no coinciden

array.id <- sub (pattern = ".CEL", replacement = "", colnames (datos))
array.id <- sub (pattern = "PP",   replacement = "id", array.id)
array.id <- sub (pattern = "\\.",  replacement = "", array.id)
array.id

colnames (datos) <- array.id

table (colnames (datos) %in% rownames (ped))

###SALVAMOS
save (list = c("datos", "annot"), file = file.path (.job$datadir, "data_processed", "datos_filtrados.RData"))


###EXIT
warnings ()
sessionInfo ()
q ("no")

