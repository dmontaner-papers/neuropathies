##t010_estudio_eroresm_plink.r
##2011-03-12 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Vemos los ERRORES Mendelianos detectados por el plink


## To generate a list of Mendel errors for SNPs and families, use the option:
  
##   plink --file plinkfile --mendel

## which will create files:
##      plink.mendel
##      plink.imendel
##      plink.fmendel
##      plink.lmendel

## The *.mendel file contains all Mendel errors (i.e. one line per error);
## the *.imendel file contains a summary of per-individual error rates;
## the *.fmendel file contains a summary of per-family error rates;
## the *.lmendel file contains a summary of per-SNP error rates.


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

dim (ped)
dim (datos)

table (colnames (datos) %in% rownames (ped)) ##OK

################################################################################


##Datos Errores Mendelianos - Plink

imen <- read.table (file = file.path (.job$datadir, "data_processed", "plink", "plink.imendel"),
                    header = TRUE, as.is = TRUE)

sapply (imen, class)

imen[,'pct'] <- 100 * imen$N / nrow (datos)

imen

################################################################################


hea <- read.table (file = file.path (.job$datadir, "data_processed", "plink", "plink.mendel"),
                   header = FALSE, as.is = TRUE, nrow = 1)
hea <- unlist (hea)
hea

men <- read.table (file = file.path (.job$datadir, "data_processed", "plink", "plink.mendel"),
                   header = FALSE, as.is = TRUE, skip = 1)

dim (men) #50233    10

sapply (men, class)
colnames (men)[1:5] <- hea[1:5]

length (unique (men$SNP)) #42334
100 *length (unique (men$SNP)) / nrow (datos)

men[1:10,]
datos['rs2272908',]
datos['rs12073797',]


###SALVAMOS errores mendelianos segun plink
men.plink <- unique (men$SNP)
length (men.plink)
men.plink[1:10]

save (list = "men.plink", file = file.path (.job$datadir, "data_processed", "errores_mendel_plink.RData"))



###EXIT
warnings ()
sessionInfo ()
q ("no")
