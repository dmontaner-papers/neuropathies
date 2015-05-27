##v010_datos_dbSNP_cromosoma_17.r
##2012-02-16 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Revisamos los resultados del cromosoma 17
##preparo la anotacion descargada de dbSNP

## FORMAT DESCRIPTIONS CONTINUED
## --------------------
## CHROMOSOME REPORTS
## --------------------

## The Chromosome Reports format provides an ordered list of RefSNPs in approximate
## chromosome coordinates (the same coordinate system used for the
## NCBI genome MapViewer); it is a small file to download, and contains a great
## deal of information about each SNP.

## The lines of the Chromosome Report format give the following information
## for a single RefSNP in tab-delimited columns:

## Column   Data
##   1      RefSNP id (rs#)
##   2      mapweight where
##             1 = Unmapped
##             2 = Mapped to single position in genome
##             3 = Mapped to 2 positions on a single chromosome
##             4 = Mapped to 3-10 positions in genome (possible paralog hits)
##             5 = Mapped to >10 positions in genome.  
## 		Please note that the number used to code mapweight for the Chromosome Report Format
## 		is different from the database tables.  Please see the online FAQ at:
## 		http://www.ncbi.nlm.nih.gov/books/bv.fcgi?mapweight&rid=helpsnpfaq.
## 		section.Build.The_dbSNP_Mapping_Pr#Build.Mapping_Weight
##   3      snp_type where
##             0 = Not withdrawn.
##             1 = Withdrawn. There are several reasons for withdrawn, the
##                 withdrawn status is fully defined in the asn1, flatfile,
##                 and XML descriptions of the RefSNP. See /specs/docsum_3.0.asn
##                 for a full definition of snp-type values.
##   4      Total number of chromosomes hit by this RefSNP during mapping
##   5      Total number of contigs hit by this RefSNP during mapping
##   6      Total number of hits to genome by this RefSNP during mapping
##   7      Chromosome for this hit to genome
##   8      Contig accession for this hit to genome
##   9      Version number of contig accession for this hit to genome
##  10      Contig ID for this hit to genome
##  11      Position of RefSNP in contig coordinates
##  12      Position of RefSNP in chromosome coordinates (used to order report)
##             Locations are specified in NCBI sequence location convention where:
##                    x, a single number, indicates a feature at base position x
##                    x..y, denotes a feature that spans from x to y inclusive.
##                    x^y, denotes a feature that is inserted between bases x and y
##  13      Genes at this same position on the chromosome
##  14      Average heterozygosity of this RefSNP
##  15      Standard error of average heterozygosity
##  16      Maximum reported probability that RefSNP is real. (For computationally-
##              predicted submissions)
##  17      Validated status
##              0 = No validation information
##              1 = Cluster has 2+ submissions, with 1+ submission assayed 
##                  with a non-computational method
##              2 = At least one subsnp in cluster has frequency data submitted
##              3 = Non-computational method in cluster and frequency data present
##              4 = At lease one subsnp in cluster has been experimentally 
##                  validated by submitter
##                  for other validation status values, please see:
##                  <a href="ftp://ftp.ncbi.nih.gov/snp/database/organism_shared_data/SnpValidationCode.bcp.gz">ftp://ftp.ncbi.nih.gov/snp/database/organism_shared_data
##                  /SnpValidationCode.bcp.gz</a>
##  18      Genotypes available in dbSNP for this RefSNP
##              1 = yes
##              0 = no
##  19      Linkout available to submitter website for further data on the RefSNP
##              1 = yes
##              0 = no
##  20      dbSNP build ID when the refSNP was first created (i.e. the create date)
##  21      dbSNP build ID of the most recent change to the refSNP cluster. The
##          date of the change is represented by the build ID which has an
##          approximate date/time associated with it. (see:
##          http://www.ncbi.nlm.nih.gov/projects/SNP/buildhistory.cgi)
## 22       Mapped to a reference or alternate (e.g. Celera) assembly 

## Also included within the chr_rpt file are two additional files:
##               multi/ contains a list (in chromosome report format) of SNPs that
##               hit multiple chromosomes in the genome

##               noton/ contains a list (in chromosome report format) of SNPs that
##               didn't hit any chromosomes in the genome

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 2.12.1 (2010-12-16)"
#library (kinship); packageDescription ("kinship", fields = "Version") #"1.1.0-23"

try (source (".job.r")); try (.job)


###DATOS
setwd (file.path (.job$datadir, "data_raw", "dbSNP"))
system.time (dbsnp <- read.table (file = "chr_17.txt", header = FALSE, sep = "\t", quote = "", as.is = TRUE, skip = 7))
dim (dbsnp) #2143929      26

dbsnp <- dbsnp[,c(1, 2, 4, 7, 12, 13, 16, 22)]
colnames (dbsnp) <- c("id", "mapweigh", "NtotCHR", "chr", "pos", "genes", "prob", "ref")
dbsnp[,"id"] <- paste ("rs", dbsnp$id, sep = "")
sapply (dbsnp, class)

table (duplicated (dbsnp))

dbsnp <- unique (dbsnp)
dim (dbsnp)

summary (dbsnp)

table (duplicated (dbsnp[,c("id", "chr", "pos", "ref")]))
table (duplicated (dbsnp[,c("id")]))


table (dbsnp[,"chr"])
table (dbsnp[,"NtotCHR"]) ##???

table (dbsnp[,"ref"])


###eliminamos los que tienen la pos NA
dbsnp <- dbsnp[!is.na (dbsnp$pos),]


###SALVAMOS
save (list = "dbsnp", file = file.path (.job$datadir, "data_processed", "anotacion_dbSNP_cromosoma17.RData"))

###EXIT
warnings ()
sessionInfo ()
q ("no")
