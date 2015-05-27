##u020_data_to_list.r
##2011-03-13 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Eliminamos errores mendelianos.
##Formateamos los datos como una lista segun los cromosomas.

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
load (file = file.path (.job$datadir, "data_processed", "errores_mendel_merlin.RData"))

ls ()

dim (ped)
dim (datos)
dim (annot)
length (unique (men.merlin))
length (men.merlin)

table (colnames (datos) %in% rownames (ped)) ##OK
table (rownames (datos)  ==  rownames (annot)) ##OK

table (annot[,"Chromosome"])

################################################################################


###Eliminamos los Errores Mendelianos

error <- rownames (datos) %in% men.merlin

datos <- datos[!error,]
annot <- annot[!error,]

table (colnames (datos) %in% rownames (ped)) ##OK
table (rownames (datos)  ==  rownames (annot)) ##OK

################################################################################


###Incluimos los missing y ordenamos segun el pedigree.

faltan <- setdiff (rownames (ped), colnames (datos))
faltan

zeros <- matrix ("0\t0", nrow = nrow (datos), ncol = length (faltan))
colnames (zeros) <- faltan
dim (zeros)
zeros[1:3,]

datos <- cbind (datos, zeros)
datos <- datos[,rownames (ped)]

dim (datos)
dim (ped)

table (colnames (datos) == rownames (ped)) ##OK

################################################################################


###Pasamos a formato de lista

cromosomas <- sort (unique (annot[,"Chromosome"]))
cromosomas

latos <- list ()
lanot <- list ()
##
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  ##
  touse <- annot[,"Chromosome"] == cro
  print (sum (touse))
  ##
  latos[[cro]] <- datos[touse,]
  lanot[[cro]] <- annot[touse,]
}

names (latos)
names (lanot)

latos.ids <- lapply (latos, rownames)
lanot.ids <- lapply (lanot, rownames)

identical (latos.ids, lanot.ids) ##OK

latos.col <- t (sapply (latos, colnames))
latos.col

table (latos.col[1,] == rownames (ped))

latos[["2"]][1:3,]
lanot[["2"]][1:3,]

################################################################################

###SALVAMOS
save (list = c("latos", "lanot", "ped", "cromosomas"), file = file.path (.job$datadir, "data_processed", "datos_por_cromosoma.RData"))


###EXIT
warnings ()
sessionInfo ()
q ("no")
