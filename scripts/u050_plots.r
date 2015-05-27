##u050_plots.r
##2011-03-14 dmontaner@cipf.es
##2011-05-03 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Hacemos plots de los resultados del Merlin

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 2.12.1 (2010-12-16)"
#library (kinship); packageDescription ("kinship", fields = "Version") #"1.1.0-23"

try (source (".job.r")); try (.job)

###DATOS
load (file = file.path (.job$datadir, "data_processed", "res_merlin0.RData"))
load (file = file.path (.job$datadir, "data_processed", "res_merlin1.RData"))
load (file = file.path (.job$datadir, "data_processed", "res_merlin2.RData"))
load (file = file.path (.job$datadir, "data_processed", "res_merlin3.RData"))
ls ()

graphics.off ()

################################################################################

###merlin0
pdf (file = file.path (.job$plotsdir, "res_merlin0_GeneticMap.pdf"), onefile = TRUE)
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  ##
  plot (resm0[[cro]][,c("GeneticMap", "LOD")], type = "l", main = cro)
  abline (h = 0, col = "blue")
  abline (h = 2, col = "red")
  axis (side = 2, at = 2, labels = 2, col = "red")
}
dev.off ()
###
pdf (file = file.path (.job$plotsdir, "res_merlin0_PhysicalPosition.pdf"), onefile = TRUE)
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  ##
  plot (resm0[[cro]][,c("PhysicalPosition", "LOD")], type = "l", main = cro)
  abline (h = 0, col = "blue")
  abline (h = 2, col = "red")
  axis (side = 2, at = 2, labels = 2, col = "red")
}
dev.off ()
################################################################################

###merlin1
pdf (file = file.path (.job$plotsdir, "res_merlin1_GeneticMap.pdf"), onefile = TRUE)
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  ##
  plot (resm1[[cro]][,c("GeneticMap", "LOD")], type = "l", main = cro)
  abline (h = 0, col = "blue")
  abline (h = 2, col = "red")
  axis (side = 2, at = 2, labels = 2, col = "red")
}
dev.off ()
###
pdf (file = file.path (.job$plotsdir, "res_merlin1_PhysicalPosition.pdf"), onefile = TRUE)
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  ##
  plot (resm1[[cro]][,c("PhysicalPosition", "LOD")], type = "l", main = cro)
  abline (h = 0, col = "blue")
  abline (h = 2, col = "red")
  axis (side = 2, at = 2, labels = 2, col = "red")
}
dev.off ()
################################################################################

###merlin2
pdf (file = file.path (.job$plotsdir, "res_merlin2_GeneticMap.pdf"), onefile = TRUE)
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  ##
  plot (resm2[[cro]][,c("GeneticMap", "LOD")], type = "l", main = cro)
  abline (h = 0, col = "blue")
  abline (h = 2, col = "red")
  axis (side = 2, at = 2, labels = 2, col = "red")
}
dev.off ()
###
pdf (file = file.path (.job$plotsdir, "res_merlin2_PhysicalPosition.pdf"), onefile = TRUE)
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  ##
  plot (resm2[[cro]][,c("PhysicalPosition", "LOD")], type = "l", main = cro)
  abline (h = 0, col = "blue")
  abline (h = 2, col = "red")
  axis (side = 2, at = 2, labels = 2, col = "red")
}
dev.off ()
################################################################################
################################################################################


###merlin012 incluimos los tres analisis (0-1-2) en la misma
pdf (file = file.path (.job$plotsdir, "res_merlin012_GeneticMap.pdf"), onefile = TRUE)
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  ##
  plot (resm0[[cro]][,c("GeneticMap", "LOD")], type = "l", main = cro)
  abline (h = 0, col = "blue")
  abline (h = 2, col = "red")
  axis (side = 2, at = 2, labels = 2, col = "red")
  ##
  lines (resm1[[cro]][,c("GeneticMap", "LOD")], type = "l", col = "blue", lty = 2)
  lines (resm2[[cro]][,c("GeneticMap", "LOD")], type = "l", col = "red",  lty = 3)
}
dev.off ()
###
pdf (file = file.path (.job$plotsdir, "res_merlin012_PhysicalPosition.pdf"), onefile = TRUE)
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  ##
  plot (resm0[[cro]][,c("PhysicalPosition", "LOD")], type = "l", main = cro)
  abline (h = 0, col = "blue")
  abline (h = 2, col = "red")
  axis (side = 2, at = 2, labels = 2, col = "red")
  ##
  lines (resm1[[cro]][,c("PhysicalPosition", "LOD")], type = "l", col = "blue", lty = 2)
  lines (resm2[[cro]][,c("PhysicalPosition", "LOD")], type = "l", col = "red",  lty = 3)
}
dev.off ()

################################################################################
################################################################################


###merlin3 (analisisi con el pedigree bien definido segun Eduardo)
pdf (file = file.path (.job$plotsdir, "res_merlin3_GeneticMap.pdf"), onefile = TRUE)
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  ##
  plot (resm3[[cro]][,c("GeneticMap", "LOD")], type = "l", main = cro)
  abline (h = 0, col = "blue")
  abline (h = 2, col = "red")
  axis (side = 2, at = 2, labels = 2, col = "red")
}
dev.off ()
###
pdf (file = file.path (.job$plotsdir, "res_merlin3_PhysicalPosition.pdf"), onefile = TRUE)
for (cro in cromosomas) {
  cat ("\n")
  print (cro)
  ##
  plot (resm3[[cro]][,c("PhysicalPosition", "LOD")], type = "l", main = cro)
  abline (h = 0, col = "blue")
  abline (h = 2, col = "red")
  axis (side = 2, at = 2, labels = 2, col = "red")
}
dev.off ()


###EXIT
warnings ()
sessionInfo ()
q ("no")
