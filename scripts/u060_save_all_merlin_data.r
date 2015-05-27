##u060_save_all_merlin_data.r
##2011-05-03 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Recopilamos todos los datos en formato de Merlin para envairselos a Eduardo

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

##cambio missing segun las indicaciones de Eduardo
ped["id724", "afec"] <- 2 ##El 724 es ENFERMO
ped["id726", "afec"] <- 1 ##el 726 es SANO

ped

################################################################################

setwd (file.path (.job$datadir, "data_results", "from_script_results", "all_merlin_data"))

### preparamos los datos para Merlin

all.ped <- ped[,1:6]
all.map <- NULL
all.dat <- c("A", "disease")

all.ped

for (cro in cromosomas) {
  cat ("\n")
  print (cro)

  ## ### PED ###
  all.ped <- cbind (all.ped, t (latos[[cro]]))
  
  ## ### MAP ###
  all.map <- rbind (all.map, lanot[[cro]][,c("Chromosome", "dbSNP.RS.ID", "deCODE.genetic.map")])
  
  ## ### DAT ###
  all.dat <- rbind (all.dat, cbind ("M", lanot[[cro]][,"dbSNP.RS.ID"]))
}

###SALVAMOS

write.table (all.ped, file = "all_data.ped",
             append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table (all.map, file = "all_data.map",
             append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table (all.dat, file = "all_data.dat",
             append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

###EXIT
warnings ()
sessionInfo ()
q ("no")
