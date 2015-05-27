##v020_revision_anotacion_cromosoma17.r
##2011-03-14 dmontaner@cipf.es
##2011-05-03 dmontaner@cipf.es
##analisis de datos para Carmen Espinos y Eduardo Calpena.
##Hacemos reports de las regiones significativas del Merlin

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 2.12.1 (2010-12-16)"
#library (kinship); packageDescription ("kinship", fields = "Version") #"1.1.0-23"

try (source (".job.r")); try (.job)


###ANOTACION dbSNP actualizada
load (file = file.path (.job$datadir, "data_processed", "anotacion_dbSNP_cromosoma17.RData"))
dim (dbsnp)


###ANOTACION antigua
load (file = file.path (.job$datadir, "data_processed", "anotacion.RData"))
annot[,"Physical.Position"] <- as.numeric (annot[,"Physical.Position"])


###un snp cualquiera
dbsnp[dbsnp$id == "rs11713",]
annot[annot$dbSNP.RS.ID == "rs11713", "Physical.Position"]  ##OK MIS posiciones coinciden con las "GRCh37.p5"


################################################################################

length (setdiff (dbsnp[,"id"], annot[,"dbSNP.RS.ID"]))
length (setdiff (annot[,"dbSNP.RS.ID"], dbsnp[,"id"]))

touse17 <- annot$Chromosome == "17"
table (touse17)
setdiff (annot[touse17, "dbSNP.RS.ID"], dbsnp[,"id"])  ##ok SOLO UNO ha salido del 17

comunes <- intersect (dbsnp[,"id"], annot[,"dbSNP.RS.ID"])
length (comunes)

table (comunes %in% annot[touse17, "dbSNP.RS.ID"])

annot <- annot[touse17,]

dbsnp <- dbsnp[dbsnp$id %in% comunes,]

################################################################################

dim (annot)
dim (dbsnp)

table (dbsnp[,"ref"], exclude = NULL)

dbsnp.gr <- dbsnp[dbsnp$ref == "GRCh37.p5",]
dbsnp.hr <- dbsnp[dbsnp$ref == "HuRef",]
dbsnp.ta <- dbsnp[dbsnp$ref == "CRA_TCAGchr7v2",]

table (duplicated (dbsnp.gr[,"id"]))
table (duplicated (dbsnp.hr[,"id"]))
table (duplicated (dbsnp.ta[,"id"]))

###

length (intersect (dbsnp.gr[,"id"], dbsnp.hr[,"id"]))
setdiff (dbsnp.gr[,"id"], dbsnp.hr[,"id"])
setdiff (dbsnp.hr[,"id"], dbsnp.gr[,"id"])   ### ok: HR < GR

faltan <- setdiff (annot[,"dbSNP.RS.ID"], dbsnp.gr[,"id"])
faltan ##"rs41480146"

###

dup <- duplicated (dbsnp.gr[,"id"])
dups <- unique (dbsnp.gr[dup, "id"])
dups %in% annot[,"dbSNP.RS.ID"]

dbs.sin.dup <- dbsnp.gr[!dbsnp.gr$id %in% dups,]
dbs.con.dup <- dbsnp.gr[dbsnp.gr$id %in% dups,]
dbs.con.dup

table (annot[,"dbSNP.RS.ID"] %in% dbs.sin.dup[!dup, "id"])

rownames (dbs.sin.dup) <- dbs.sin.dup[,"id"]
annot[,"pos"] <- dbs.sin.dup[annot$dbSNP.RS.ID, "pos"]

table  (annot[,"pos"] == annot[,"Physical.Position"])


#### los MALOS tambien coinciden
dbs.con.dup
annot[annot$dbSNP.RS.ID %in% dups, c("dbSNP.RS.ID", "Physical.Position")]

dups


###EXIT
warnings ()
sessionInfo ()
q ("no")
