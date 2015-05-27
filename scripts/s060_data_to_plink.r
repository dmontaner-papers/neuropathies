##s060_data_to_plink.e
##2011-03-12 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Formateamos los datos para enviar al Plink

################################################################################
###                         Fichero ped de Plink                             ###
################################################################################

## The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:
##      Family ID
##      Individual ID
##      Paternal ID
##      Maternal ID
##      Sex (1=male; 2=female; other=unknown)
##      Phenotype

## The phenotype can be either a quantitative trait or an affection status column:
## PLINK will automatically detect which type (i.e. based on whether a value other than 0, 1, 2
## or the missing genotype code is observed).

## Affection status, by default, should be coded:
##     -9 missing 
##      0 missing
##      1 unaffected
##      2 affected

## If your file is coded 0/1 to represent unaffected/affected, then use the --1 flag:

## plink --file mydata --1

## which will specify a disease phenotype coded:
##      -9 missing
##       0 unaffected
##       1 affected

## Genotypes (column 7 onwards) should also be white-space delimited;
## they can be any character (e.g. 1,2,3,4 or A,C,G,T or anything else)
## except 0 which is, by default, the missing genotype character.

## All markers should be biallelic.
## All SNPs (whether haploid or not) must have two alleles specified.
## Either Both alleles should be missing (i.e. 0) or neither.
## No header row should be given.


################################################################################
###                         Fichero map de Plink                             ###
################################################################################

## By default, each line of the MAP file describes a single marker and must contain exactly 4 columns:
##      chromosome (1-22, X, Y or 0 if unplaced)
##      rs# or snp identifier
##      Genetic distance (morgans)
##      Base-pair position (bp units)

## Genetic distance can be specified in centimorgans with the --cm flag.
## Alternatively, you can use a MAP file with the genetic distance excluded by adding the flag --map3, i.e.
## plink --file mydata --map3

## In this case, the three columns are expected to be
##      chromosome (1-22, X, Y or 0 if unplaced)
##      rs# or snp identifier
##      Base-pair position (bp units)

## Base-pair positions are expected to correspond to positive integers
## within the range of typical human chromosome sizes.

## The autosomes should be coded 1 through 22.
## The following other codes can be used to specify other chromosome types:
##      X    X chromosome                    -> 23
##      Y    Y chromosome                    -> 24
##      XY   Pseudo-autosomal region of X    -> 25
##      MT   Mitochondrial                   -> 26

################################################################################
###                          Distancia Genetica                              ###
################################################################################

##DE la anotacion de affymetrix

##    3. Genetic maps.

##    This annotation is to provide a rough estimate on SNP genetic 
## distances to the p-telomere. Those estimates may be used as seed input 
## for linkage analysis programs like MERLIN 
## (http://www.sph.umich.edu/csg/abecasis/Merlin/). As a requirement of 
## MERLIN, every SNP has to have a unique genetic distance. SNP genetic 
## distances were extrapolated from three experimentally obtained genetic 
## maps: deCODE map, Marshfield map, and SLM1 map.

## 	1). deCODE genetic map was built by genotyping 5,136 
## microsatellite markers for 146 families (1), and is available through 
## Nature Genetics (see reference).

##   	2). Marshfield genetic map (2) was based on CEPH family genotypes 
## for 7,740 microsatellite markers, and is available from 
## http://research.marshfieldclinic.org/genetics/Map_Markers/maps/IndexMap
## Frames.html.

##   	3). SLM1 (SNP Linkage Map) map was generated from unpublished 
## data from Affymetrix and Dr. Aravinda Chakravarti group at Johns 
## Hopkins University. It was based on genotypes for 2,022 microsatellite 
## markers and 6,205 SNPs.

##    Physical locations of markers used in each genetic map were obtained 
## from the UCSC database (ftp://genome.ucsc.edu). Markers are removed 
## when their genetic order is opposite of their physical order. When 
## physically neighboring markers share the same genetic distance, only 
## the one with the largest physical position was kept in order that no 
## two SNPs have exactly the same genetic distance. Physical positions of 
## SNPs and physical locations of the markers in the cleaned genetic maps 
## were compared to infer genetic distances for SNPs. We assume that 
## genetic distance changes linearly with physical distance between any 
## two neighboring markers in each genetic map. 


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

setwd (file.path (.job$datadir, "data_processed", "plink"))
dir ()


### PED

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

writeLines (d, con = "plinkfile.ped", sep = "\n")

################################################################################


### MAP
annot[1:10, c("Physical.Position", "posicion")]

annot[10, c("Genetic.Map", "deCODE.genetic.map")]
annot[10000, c("Genetic.Map", "deCODE.genetic.map")]


write.table (annot[,c("Chromosome", "dbSNP.RS.ID", "deCODE.genetic.map", "posicion")],
             file = "plinkfile.map", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


###EXIT
warnings ()
sessionInfo ()
q ("no")
