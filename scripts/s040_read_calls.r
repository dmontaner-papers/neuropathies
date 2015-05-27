##s040_read_calls.r
##2011-03-11 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Leemos los calls para convertirlos a sus alelos

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 2.12.1 (2010-12-16)"
#library (kinship); packageDescription ("kinship", fields = "Version") #"1.1.0-23"

try (source (".job.r")); try (.job)


###ANOTACION
system.time (annot <- read.csv (file = file.path (.job$datadir, "data_raw", "annotation", "Mapping250K_Nsp.na31.annot.csv",
                                  "Mapping250K_Nsp.na31.annot.csv"), header = TRUE, sep = ",", quote="\"", dec=".",
                                fill = TRUE, comment.char="#", as.is = TRUE))
dim (annot) #262264     27
clases <- sapply (annot, class)
clases

colnames (annot)

table (annot[,"Allele.A"], exclude = NULL)
table (annot[,"Allele.B"], exclude = NULL)

table (annot[,"Chromosome"], exclude = NULL)

table (duplicated (annot[,"Probe.Set.ID"])) #OK no hay duplicados
rownames (annot) <- annot[,"Probe.Set.ID"]

table (duplicated (annot[,"dbSNP.RS.ID"]))  #PROBLEMA HAY UN DUPLICADO

annot[duplicated (annot[,"dbSNP.RS.ID"]),"dbSNP.RS.ID"]

table (annot[,"Strand"], exclude = NULL)

annot[1000:1010, "OMIM"]
annot[1000:1010, "Copy.Number.Variation"]

################################################################################


###DATOS sin HapMap
system.time (datos.sin <- read.table (file = file.path (.job$datadir, "data_processed", "apt", "probeset-genotype-chp",
                                        "brlmm.calls.txt"),
                                      header = TRUE, sep = "\t", quote = "", as.is = TRUE, comment.char = "#", nrows = 100))
mis.clases <- sapply (datos.sin, class)
mis.clases

system.time (datos.sin <- read.table (file = file.path (.job$datadir, "data_processed", "apt", "probeset-genotype-chp",
                                        "brlmm.calls.txt"),
                                      header = TRUE, sep = "\t", quote = "", row.names = 1, colClasses = mis.clases,
                                      comment.char = "#"))

datos.sin[1:10,]
dim (datos.sin) # 262314     15

################################################################################
################################################################################

table (rownames (datos.sin) %in% annot[,"Probe.Set.ID"]) #50 CONTROLES

table (annot[,"Probe.Set.ID"] %in% rownames (datos.sin)) #OK estan todos

datos.sin <- datos.sin[annot[,"Probe.Set.ID"],]
datos.sin[1:3,]

##datos <- matrix (NA, nrow = nrow (datos.sin), ncol = ncol (datos.sin))
datos <- matrix ("0\t0", nrow = nrow (datos.sin), ncol = ncol (datos.sin))
dimnames (datos) <- dimnames (datos.sin)
datos[1:3,]

table (rownames (datos) == annot[,"Probe.Set.ID"])

aleloA <- matrix (annot[,"Allele.A"], nrow = nrow (datos.sin), ncol = ncol (datos.sin))
aleloB <- matrix (annot[,"Allele.B"], nrow = nrow (datos.sin), ncol = ncol (datos.sin))

alelo0 <- matrix (paste (aleloA, aleloA, sep = "\t"), nrow = nrow (datos.sin), ncol = ncol (datos.sin))
alelo1 <- matrix (paste (aleloA, aleloB, sep = "\t"), nrow = nrow (datos.sin), ncol = ncol (datos.sin))
alelo2 <- matrix (paste (aleloB, aleloB, sep = "\t"), nrow = nrow (datos.sin), ncol = ncol (datos.sin))

annot[1:5, "Allele.A"]
annot[1:5, "Allele.B"]

aleloA[1:5,]
aleloB[1:5,]

alelo0[1:5,]
alelo1[1:5,]
alelo2[1:5,]

datos[datos.sin == 0] <- alelo0[datos.sin == 0]
datos[datos.sin == 1] <- alelo1[datos.sin == 1]
datos[datos.sin == 2] <- alelo2[datos.sin == 2]


datos.sin[1000:1003, 1:5]
annot[1000:1003, c("Probe.Set.ID", "Allele.A", "Allele.B")]
datos[1000:1003, 1:5]

################################################################################

###SALVAMOS
save (list = "datos", file = file.path (.job$datadir, "data_processed", "all_calls_limpias.RData"))
save (list = "annot", file = file.path (.job$datadir, "data_processed", "anotacion.RData"))


###EXIT
warnings ()
sessionInfo ()
q ("no")
