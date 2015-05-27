##u030_analisis_merlin_0.r
##2011-03-13 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Formateamos los datos y hacemos el analisis con Merlin cromosoma a coromosoma

##TRABAJAMOS CON LOS DOS afectados MISSING COMO MISSING
##TRABAJAMOS CON LOS DOS afectados MISSING COMO SANOS
##TRABAJAMOS CON LOS DOS afectados MISSING COMO ENFERMOS

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 2.12.1 (2010-12-16)"
#library (kinship); packageDescription ("kinship", fields = "Version") #"1.1.0-23"

try (source (".job.r")); try (.job)

###PARAMETROS
NMM <- 2500  #numero MAXIMO de marcadores para incluir en el Merlin

###DATOS
load (file = file.path (.job$datadir, "data_processed", "datos_por_cromosoma.RData"))
ls ()

cromosomas <- cromosomas[12:13]
cromosomas

##cambio missing por nada
#setwd (file.path (.job$datadir, "data_processed", "merlin_0"))
##cambio missing por sanos
#setwd (file.path (.job$datadir, "data_processed", "merlin_1")); ped[ped$afec == 0, "afec"] <- 1
##cambio missing por enfermos
setwd (file.path (.job$datadir, "data_processed", "merlin_2")); ped[ped$afec == 0, "afec"] <- 2 

ped

ped.t <- t (ped[,1:6])
ped.t
################################################################################


### preparamos los datos para Merlin

dir ()

for (cro in cromosomas) {
  cat ("\n")
  print (cro)

  ## DIRECTORIO
  dir.create (cro)
  file.copy ("modelo.txt", to = file.path (cro, "modelo.txt"))
  
  ## RECORTE
  N <- nrow (latos[[cro]])
  print (N)
  if (N < NMM){
    indices <- 1:N
  } else {
    indices <- round (seq (from = 1, to = N, length.out = NMM))
  }
  
  ## ### PED ###
  chr.ped <- t(rbind (ped.t, latos[[cro]][indices,]))
  write.table (chr.ped, file = file.path (cro, "datos.ped"),
               append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  ## ### MAP ###
  write.table (lanot[[cro]][indices, c("Chromosome", "dbSNP.RS.ID", "deCODE.genetic.map")],
               file = file.path (cro, "datos.map"),
               append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  ## ### DAT ###
  dat <- rbind (c("A", "disease"),                                  ##   A - affection status
                cbind ("M", lanot[[cro]][indices, "dbSNP.RS.ID"]))  ##   M - marker
  write.table (dat, file.path (cro, "datos.dat"),
               append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
}
################################################################################


#rm (list = setdiff (ls (), "cromosomas"))
ls ()
gc ()
gc ()

### Merlin
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  
  setwd (cro)

  print  ("pedstats")
  system ("pedstats -d datos.dat -p datos.ped                                            > salida_pedstats.txt")
  print  ("merlin --error")
  system ("merlin   -d datos.dat -p datos.ped -m datos.map --error                       > salida_merlin_error.txt")
  print  ("pedwipe")
  system ("pedwipe  -d datos.dat -p datos.ped                                            > salida_pedwipe.txt")
  print  ("merlin")
  system ("merlin   -d wiped.dat -p wiped.ped -m datos.map --model modelo.txt --tabulate > salida_merlin.txt")

  setwd ("..")
}  

###EXIT
warnings ()
sessionInfo ()
q ("no")
