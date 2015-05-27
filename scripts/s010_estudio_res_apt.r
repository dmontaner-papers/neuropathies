##s010_estudio_res_apt.r
##2011-03-11 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Chequeamos la calidad de los resultados del APT

##HAY MUCHISIMOS MISSING CALLS

### 14 individuos:
##                   11 antiguos
##                    3 nuevos - Con uno de los arrays por DUPLICADO
## 15 arrays

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 2.12.1 (2010-12-16)"
#library (kinship); packageDescription ("kinship", fields = "Version") #"1.1.0-23"

try (source (".job.r")); try (.job)

###PED
load (file = file.path (.job$datadir, "data_processed", "ped.RData"))
ped

################################################################################

pdf (file = file.path (.job$plotsdir, "plots_calidad.pdf"),  onefile = TRUE)

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
dim (datos.sin)
summary (datos.sin)

conteos.sin <- t (apply (datos.sin, 2, table))
conteos.sin

malos.rate.sin <- 100 * conteos.sin[,1] / rowSums (conteos.sin[,2:4])
malos.rate.sin

cbind (conteos.sin, malos.rate.sin)

###

misarrays <- colnames (datos.sin)
misarrays

################################################################################


###DATOS con HapMap

system.time (datos.con <- read.table (file = file.path (.job$datadir, "data_processed", "apt_hapmap",
                                        "probeset-genotype-chp", "brlmm.calls.txt"),
                                      header = TRUE, sep = "\t", quote = "", as.is = TRUE, comment.char = "#", nrows = 100))
mis.clases <- sapply (datos.con, class)
mis.clases

##15 segundos
system.time (datos.con <- read.table (file = file.path (.job$datadir, "data_processed", "apt_hapmap",
                                        "probeset-genotype-chp", "brlmm.calls.txt"),
                                      header = TRUE, sep = "\t", quote = "", row.names = 1, colClasses = mis.clases,
                                      comment.char = "#"))

datos.con[1:10,]
dim (datos.con)
summary (datos.con[,1:10])

system.time (conteos.con <- t (apply (datos.con, 2, table))) #25 segundos
conteos.con[misarrays,]

malos.rate.con <- 100 * conteos.con[,1] / rowSums (conteos.con[,2:4])

cbind (conteos.con, malos.rate.con)[misarrays,]

################################################################################

cbind (malos.rate.sin, malos.rate.con[misarrays])

mio <- colnames (datos.con) %in% misarrays
table (mio)

summary (malos.rate.sin)
summary (malos.rate.con[mio])
summary (malos.rate.con[!mio])


plot (cbind (malos.rate.sin, malos.rate.con[misarrays]))
abline (0,1, col = "red")

boxplot (malos.rate.con ~ mio)
boxplot (malos.rate.con[!mio])
boxplot (malos.rate.con[mio])

################################################################################

table (rownames (datos.sin) == rownames (datos.con[,misarrays]))
table (colnames (datos.sin) == colnames (datos.con[,misarrays]))


iguales <- datos.sin == datos.con[,misarrays]
dim (iguales)

100 * sum (iguales) / length (iguales) #80 % coincidencias

conteos.iguales <- apply (iguales, 2, sum)
conteos.iguales.pct <- 100 * conteos.iguales / nrow (datos.sin)

cbind (conteos.iguales, conteos.iguales.pct)

cbind (malos.rate.sin, malos.rate.con[misarrays], conteos.iguales.pct)

plot (malos.rate.sin/malos.rate.con[misarrays], conteos.iguales.pct) #parece aleatorio

plot (malos.rate.sin, malos.rate.sin/malos.rate.con[misarrays])

###

dim (iguales)
dim (datos.sin)

##OBS: de los calls que coinciden entre las dos normalizaciones casi todos son no missing calls
100 * table (datos.sin[iguales])  / sum (iguales)
100 * table (datos.sin[!iguales]) / sum (!iguales)


###Analizamos calidad de los SNP en general: esto es por filas
dim (iguales)
iguales[1:3,]

ig.r <- apply (iguales, 1, all)
table (ig.r, exclude = NULL)

100 * sum (ig.r) / length (ig.r) ##solo un 12% coinciden exactamente

ig.s <- apply (iguales, 1, sum)
table (ig.s, exclude = NULL)

100 * table (ig.s > 5) / length (ig.s)
100 * table (ig.s > 10) / length (ig.s) #79% no esta del todo mal


###filas con missing CALL
con.call <- datos.sin != -1 
dim (con.call)

100 * table (con.call) / length (con.call)

100 * table (apply (con.call, 1, all), exclude = NULL) / nrow (con.call) ##21% de snps con call para todos los arrays

con.call[1:3,]

cc.s <- apply (con.call, 1, sum)
table (cc.s, exclude = NULL)

100 * table (cc.s > 5) / length (cc.s)
100 * table (cc.s > 10) / length (cc.s)


table (con.call=cc.s > 10, iguales=ig.s > 10)
100 * table (con.call=cc.s > 10, iguales=ig.s > 10) / length (cc.s) ##DEBERIA DE USAR ESTE

100 * table (con.call=cc.s > 11, iguales=ig.s > 10) / length (cc.s)

################################################################################


###REPORT
snp.qc <- cbind (con.call = cc.s, como.hap.map = ig.s)
snp.qc[1:100,]

array.qc <- cbind (na.rate = malos.rate.sin,
                   na.rate.hapmap = malos.rate.con[misarrays],
                   same.when.hapmap = conteos.iguales.pct)
array.qc


###SALVAMOS
save (list = c("snp.qc", "array.qc"), file = file.path (.job$datadir, "data_processed", "quality_control.RData"))

array.qc <- cbind (array = rownames (array.qc), array.qc)
write.table (array.qc, file = file.path (.job$datadir, "data_results", "from_script_results", "array_qc.xls"), 
             append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


###EXIT
dev.off ()
warnings ()
sessionInfo ()
q ("no")
