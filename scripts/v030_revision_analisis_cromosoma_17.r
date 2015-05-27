##v030_revision_analisis_cromosoma_17.r
##2011-03-14 dmontaner@cipf.es
##2011-05-03 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Hacemos reports de las regiones significativas del Merlin

###SNPs con posiciones duplicadas en la nueva version de dbSNP "rs2610383"  "rs4055285"  "rs41348644"

###SNP que parece qeu se ha movido al cromosoma 10 "rs41480146"
## HuRef	36.3	10	83363119	NW_001838005.2	10756703	+	A	-	view	blast

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 2.12.1 (2010-12-16)"
#library (kinship); packageDescription ("kinship", fields = "Version") #"1.1.0-23"

try (source (".job.r")); try (.job)

###CALLS
load (file = file.path (.job$datadir, "data_processed", "datos_filtrados.RData"))
calls <- datos
rm (datos)
calls <- sub ("\t", "/", calls)
calls <- as.data.frame (calls, stringsAsFactors = FALSE)
dim (calls)
calls[1:3,]



###ANOTACION
load (file = file.path (.job$datadir, "data_processed", "anotacion.RData"))
touse <- annot$Chromosome == "17"
table (touse)
annot <- annot[touse,]

table (duplicated (annot[,"Probe.Set.ID"]))
table (duplicated (annot[,"dbSNP.RS.ID"]))

annot[,"Physical.Position"] <- as.numeric (annot[,"Physical.Position"])
orden <- order (annot[,"Physical.Position"])
annot <- annot[orden,]

annot[1:3,]



###DATOS
load (file = file.path (.job$datadir, "data_processed", "res_merlin1.RData"))
datos <- resm1[["17"]]
dim (datos)

datos[1:3,]

datos <- datos[,c ("snp", "CHR", "PhysicalPosition", "GeneticMap", "LOD")]
sapply (datos, class)

datos[,"PhysicalPosition"] <- as.numeric (datos[,"PhysicalPosition"])

table (datos[,"PhysicalPosition"] == sort (datos[,"PhysicalPosition"]))

################################################################################


###BUSCAMOS LA REGION
## plot (datos[,"LOD"], type = "l")
## plot (datos[,"PhysicalPosition"], datos[,"LOD"], type = "l")
## ##
## plot  (datos[,"GeneticMap"], datos[, "LOD"], type = "l")
## #lines (datos[,"GeneticMap"], datos[,"HLOD"], col = "red")
## abline (h = 0, col = "blue")
## abline (h = 2, col = "red")


which (datos[,"LOD"] > 0)
datos[1:100,]

corte <- max (datos[1:100, "PhysicalPosition"])
corte

################################################################################
################################################################################

datosc <- datos[1:100,]
rownames (datosc) <- datosc$snp

c("rs2610383", "rs4055285", "rs41348644") %in% datosc$snp  ##OK los SNPs con varias posiciones no estan aqui

"rs41480146" %in% datosc$snp  ##OK el SNP desaparedico no esta aqui

################################################################################

## plot  (datosc[, "GeneticMap"], datosc[, "LOD"], type = "l")
## abline (h = 0, col = "blue")
## abline (h = 2, col = "red")
## abline (h = 3, col = "grey")

## plot  (datosc[, "PhysicalPosition"], datosc[, "LOD"], type = "l")
## abline (h = 0, col = "blue")
## abline (h = 2, col = "red")
## abline (h = 3, col = "grey")
## abline (h = 2.5, col = "grey")

################################################################################

touse <- annot[,"Physical.Position"] <= corte
table (touse)

annot <- annot[touse,]
rownames (annot) <- annot[,"dbSNP.RS.ID"]

###Plot
posiciones <- annot[,"Physical.Position"]
names (posiciones) <- annot[,"dbSNP.RS.ID"]

posin <- names (posiciones) %in% datos[,"snp"]
table (posin)

altos <- datosc[, "LOD"] > 2.5
datosc[altos,]
##
datosc["rs7223663", "PhysicalPosition"]
datosc["rs8067685", "PhysicalPosition"]


pdf (file = file.path (.job$plotsdir, "cromosoma17_primera_parte.pdf"), width = 7*2, height = 7, onefile = TRUE)
##
plot  (datosc[, "PhysicalPosition"], datosc[, "LOD"], type = "l", main = "chr 17; primera parte...")
abline (h = 0, col = "blue")
abline (h = 2, col = "red")
abline (h = 3, col = "grey")
abline (h = 2.5, col = "grey")
##
points (cbind (posiciones[!posin],   0), col = "red",   pch = "|")
points (cbind (posiciones[posin], -0.5), col = "green", pch = "|")
#abline (v = datosc["rs902966", "PhysicalPosition"], col = "grey", lty = 2)
abline (v = datosc["rs7223663", "PhysicalPosition"], col = "grey", lty = 2)
abline (v = datosc["rs8067685", "PhysicalPosition"], col = "grey", lty = 2)
##
dev.off ()

################################################################################

###FORMATEAMOS
add <- matrix (NA, nrow = sum (!posin), ncol = ncol (datosc))
colnames (add) <- colnames (datosc)
add <- as.data.frame (add)
dim (add)

tocomplete <- setdiff (annot$dbSNP.RS.ID, datosc$snp)
length (tocomplete)

add[,"snp"] <- tocomplete 
add[,"CHR"] <- 17
add[,"PhysicalPosition"] <- annot[tocomplete, "Physical.Position"]

datos2 <- rbind (datosc, add)
sapply (datos2, class)

orden <- order (datos2$PhysicalPosition)
datos2 <- datos2[orden,]

table (sort (datosc$PhysicalPosition) ==  (datosc$PhysicalPosition))

################################################################################

###incluimos las CALLS
table (datos2$snp %in% rownames (calls))
table (datosc$snp %in% rownames (calls))

datosc <- cbind (datosc, calls[datosc$snp,])
sapply (datosc, class)

datos2 <- cbind (datos2, calls[datos2$snp,])
sapply (datos2, class)

save (list = c("datosc", "datos2"), file = file.path (.job$datadir, "data_processed", "datos_principio_chr17.RData"))


###SALVAMOS
setwd (file.path (.job$datadir, "data_results", "from_script_results", "estudio_cromosoma_17"))

write.table (datosc, file = "principio_chr17_conLOD.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table (datos2, file = "principio_chr17_array.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



###EXIT
warnings ()
sessionInfo ()
q ("no")
