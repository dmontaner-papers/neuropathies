##r020_pedigree_checking.r
##2011-02-03 dmontaner@cipf.es
##2011-03-11 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Comprobamos los datos del pedigree

### 14 individuos:
##                   11 antiguos
##                    3 nuevos - Con uno de los arrays por DUPLICADO: el 907
## 15 arrays

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 2.12.1 (2010-12-16)"
library (kinship); packageDescription ("kinship", fields = "Version") #"1.1.0-23"
## library (LDheatmap)
## library (mapLD)
## library (LDtests)

try (source (".job.r")); try (.job)


###DATOS: NOMBRES DE LOS ARRAYS
ficheros <- dir (file.path (.job$datadir, "data_raw", "snp"))
ficheros ## ok son unicos

patient.id <- sub (pattern = ".CEL", replacement = "", ficheros)
patient.id <- sub (pattern = "PP ", replacement = "", patient.id)
patient.id <- sub (pattern = "PP", replacement = "", patient.id)
patient.id <- sub (pattern = "DIC", replacement = "", patient.id)

sinfo <- as.data.frame (list (filename = ficheros, samplename = patient.id), stringsAsFactors = FALSE)
sinfo

sinfo <- sinfo[-15,] ##eliminamos el array PP907DIC.CEL #que parece que tiene peor calidad
sinfo

rownames (sinfo) <- paste ("id", sinfo$samplename, sep = "")
sinfo

################################################################################


###PEDIGREE
## sexo 1:hombre 2:m
## status 0:desconocido; 1:sano; 2:enfermo

ped <- read.table (file = file.path (.job$sinfodir, "ped.txt"), header = FALSE, sep = "\t", quote = "", as.is = TRUE)
colnames (ped) <- c("fam", "id", "padre", "madre", "sexo", "afec")
ped[,"sinarray"] <- 1 * ped[,"id"] %in% c(300, 400, 921, 923)
sapply (ped, class)
ped

a <- rep (NA, times = nrow (ped))
a <- rep (0, times = nrow (ped)) #MISSING

a <- rep (1, times = nrow (ped)) #UNAFECTED - VACIO
a <- rep (2, times = nrow (ped)) #AFECTED   - LLENO

a <- ped[,"afec"]
a[a==0] <- 1
a

color <- rep ("black", times = nrow (ped)) #UNAFECTED - VACIO
color[ped[,"afec"] == 0] <- "red"
color
ped

myped <- pedigree (id = ped[,"id"], dadid = ped[,"padre"], momid = ped[,"madre"], sex = ped[,"sexo"], affected = a, status = ped[,"sinarray"])
myped


###PLOT
graphics.off ()
pdf(file = file.path (.job$plotsdir, "pedigree.pdf"))
##
par (xpd = TRUE)
plot (myped, main = "familia 266", col = color, angle = 0)
##
dev.off ()
    
##plot.pedigree (myped, main = "familia 266", col = color, angle = 90, density=50)

################################################################################


##ADD cel names to ped
rownames (ped) <- paste ("id", ped$id, sep = "")
ped

setdiff (rownames (sinfo), rownames (ped))
setdiff (rownames (ped), rownames (sinfo)) ##OK

ped[,"cel"] <- NA
ped[rownames (sinfo),"cel"] <- sinfo$filename
ped

###SALVAMOS EL PEDIGREE
save (list = "ped", file = file.path (.job$datadir, "data_processed", "ped.RData"))


###EXIT
warnings ()
sessionInfo ()
q ("no")
