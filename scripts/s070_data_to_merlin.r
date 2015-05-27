##s070_data_to_merlin.r
##2011-03-13 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Filtramos los SNPs que vamos a usar en el Merlin


################################################################################
###                                Merlin PED                                ###
################################################################################

## all the information that is necessary to reconstruct individual relationships in a 
## pedigree file can be summarized in five items:

##   family identifier
##   individual identifier
##   link to each parent (if available) 
##   indicator of each individual's sex: recoding sexes as 2 (female) and 1 (male)

## FAMILY     PERSON   FATHER   MOTHER   SEX

## text identifiers are usually replaced by unique numeric values.

## Usually the five standard columns are followed by various types of genetic data, 
## including phenotypes for discrete and quantitative traits and marker genotypes.

## Disease status is usually encoded in a single column as

##    U or 1 for unaffecteds, 
##    A or 2 for affecteds, and 
##    X or 0 for missing phenotypes.

## 0 como missing vale para el estatus y para el genotipo


################################################################################
###                                Merlin DAT                                ###
################################################################################

## Since each pedigree file has a unique structure (apart from the first five columns),
## its contents must be described in a companion data file.

## The data file includes one row per data item in the pedigree file,
## indicating the data type encoded as:

##   M - marker
##   A - affection status
##   T - Quantitative Trait
##   C - Covariate

## and providing a one-word label for each item.


################################################################################
###                                Merlin MAP                                ###
################################################################################

## To analyse genetic markers, MERLIN requires information on their chromosomal location.
## This is usually provided in a map file.

## this file has one line per marker with three columns, indicating:
##   chromosome  (parece que vale la notacion X para el cromosoma X)
##   marker name
##   position (in centiMorgans)

## CHROMOSOME     MARKER      POSITION

## If you are using sex-specific maps...

## The data file and map file can include different sets of markers,
## but markers that are absent from the map file will be ignored by MERLIN.

## ... allows MERLIN to analyse multiple chromosomes in a single run.


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

setwd (file.path (.job$datadir, "data_processed", "merlin"))
dir ()


### PED - como el de plink

zeros <- rep (0, times = 2 * nrow (datos))
length (zeros)

system.time ({
  d <- NULL
  for (pt in rownames (ped)) {
    if (pt %in% colnames (datos)) {
      print (pt)
      d[pt] <- paste (c(ped[pt, 1:6], datos[,pt]), collapse = "\t")
    } else {
      d[pt] <- paste (c(ped[pt, 1:6], zeros), collapse = "\t")
    }
  }
})
length (d)

writeLines (d, con = "merlinfile.ped", sep = "\n")

################################################################################


### MAP
annot[1:10, c("Physical.Position", "posicion")]

annot[10, c("Genetic.Map", "deCODE.genetic.map")]
annot[10000, c("Genetic.Map", "deCODE.genetic.map")]

##write.table (annot[,c("Chromosome", "dbSNP.RS.ID", "deCODE.genetic.map", "posicion")], ##plink 
write.table (annot[,c("Chromosome", "dbSNP.RS.ID", "deCODE.genetic.map")],               ##merlin
             file = "merlinfile.map", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



### DAT

dat <- rbind (c("A", "diseasecolumn"),             ##   A - affection status
             cbind ("M", annot[,"dbSNP.RS.ID"]))  ##   M - marker
dat[1:10,]

write.table (dat, file = "merlinfile.dat", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



###EXIT
warnings ()
sessionInfo ()
q ("no")
