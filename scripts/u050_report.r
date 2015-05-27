##u050_report.r
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

###DATOS
load (file = file.path (.job$datadir, "data_processed", "res_merlin0.RData"))
load (file = file.path (.job$datadir, "data_processed", "res_merlin1.RData"))
load (file = file.path (.job$datadir, "data_processed", "res_merlin2.RData"))
load (file = file.path (.job$datadir, "data_processed", "res_merlin3.RData"))
ls ()

## ################################################################################
## for (corte in c(0, 2)) {
##   for (cro in cromosomas) {
##     cat ("\n")
##     print (cro)
    
##     ##detectamos runs
##     signf <- resm0[[cro]][,"LOD"] > corte       #significative
##     prece <- c (FALSE, signf[-length (signf)])  #preceding
##     start <- signf & !prece                     #starting poin
##     cumul <- cumsum (start)                     #cumulative distribution
##     myrun <- cumul * signf                      #run of significative snps 
    
##     runs <- setdiff (unique (myrun), 0)
    
##     ##report
##     repfile <- file.path (.job$datadir, "data_results", "from_script_results", paste ("merlin0_lor", corte, "_cromosoma_", cro, ".txt", sep = ""))
    
##     cat (paste ("########## Regiones seleccionadas cromosoma", cro, "##########"),
##          file = repfile, append = FALSE, fill = TRUE)
    
##     for (run in runs) {
##       print (run)
##       ##
##       cat (paste ("##### Region:", run, "#####"),
##            file = repfile,
##            append = TRUE, fill = TRUE)
##       ##
##       write.table (resm0[[cro]][myrun == run, c("snp", "deCODE.genetic.map", "PhysicalPosition", "LOD", "ALPHA", "HLOD")],
##                    file = repfile,
##                    append = TRUE, quote = FALSE, sep = "\t",
##                    row.names = FALSE, col.names = TRUE)
##     }
##   }
## }
## ################################################################################

###merlin0
for (corte in c(0, 2)) {
  
  repfile <- file.path (.job$datadir, "data_results", "from_script_results", paste ("merlin0_lor", corte, ".txt", sep = ""))
  
  cat (paste ("########## Regiones seleccionadas con LOR:", corte, "##########"),
       file = repfile, append = TRUE, fill = TRUE)
  
  for (cro in cromosomas) {
    cat ("\n")
    print (cro)
    
    ##detectamos runs
    signf <- resm0[[cro]][,"LOD"] > corte       #significative
    prece <- c (FALSE, signf[-length (signf)])  #preceding
    start <- signf & !prece                     #starting poin
    cumul <- cumsum (start)                     #cumulative distribution
    myrun <- cumul * signf                      #run of significative snps 
    
    runs <- setdiff (unique (myrun), 0)
    
    ##report
    cat (paste ("########## Regiones seleccionadas cromosoma", cro, "##########"),
         file = repfile, append = TRUE, fill = TRUE)
    
    for (run in runs) {
      print (run)
      ##
      cat (paste ("##### Region:", run, "#####"),
           file = repfile,
           append = TRUE, fill = TRUE)
      ##
      write.table (resm0[[cro]][myrun == run, c("snp", "deCODE.genetic.map", "PhysicalPosition", "LOD", "ALPHA", "HLOD")],
                   file = repfile,
                   append = TRUE, quote = FALSE, sep = "\t",
                   row.names = FALSE, col.names = TRUE)
    }
  }
}
################################################################################

###merlin1
for (corte in c(0, 2)) {
  
  repfile <- file.path (.job$datadir, "data_results", "from_script_results", paste ("merlin1_lor", corte, ".txt", sep = ""))
  
  cat (paste ("########## Regiones seleccionadas con LOR:", corte, "##########"),
       file = repfile, append = TRUE, fill = TRUE)
  
  for (cro in cromosomas) {
    cat ("\n")
    print (cro)
    
    ##detectamos runs
    signf <- resm1[[cro]][,"LOD"] > corte       #significative
    prece <- c (FALSE, signf[-length (signf)])  #preceding
    start <- signf & !prece                     #starting poin
    cumul <- cumsum (start)                     #cumulative distribution
    myrun <- cumul * signf                      #run of significative snps 
    
    runs <- setdiff (unique (myrun), 0)
    
    ##report
    cat (paste ("########## Regiones seleccionadas cromosoma", cro, "##########"),
         file = repfile, append = TRUE, fill = TRUE)
    
    for (run in runs) {
      print (run)
      ##
      cat (paste ("##### Region:", run, "#####"),
           file = repfile,
           append = TRUE, fill = TRUE)
      ##
      write.table (resm1[[cro]][myrun == run, c("snp", "deCODE.genetic.map", "PhysicalPosition", "LOD", "ALPHA", "HLOD")],
                   file = repfile,
                   append = TRUE, quote = FALSE, sep = "\t",
                   row.names = FALSE, col.names = TRUE)
    }
  }
}
################################################################################

###merlin2
for (corte in c(0, 2)) {
  
  repfile <- file.path (.job$datadir, "data_results", "from_script_results", paste ("merlin2_lor", corte, ".txt", sep = ""))
  
  cat (paste ("########## Regiones seleccionadas con LOR:", corte, "##########"),
       file = repfile, append = TRUE, fill = TRUE)
  
  for (cro in cromosomas) {
    cat ("\n")
    print (cro)
    
    ##detectamos runs
    signf <- resm2[[cro]][,"LOD"] > corte       #significative
    prece <- c (FALSE, signf[-length (signf)])  #preceding
    start <- signf & !prece                     #starting poin
    cumul <- cumsum (start)                     #cumulative distribution
    myrun <- cumul * signf                      #run of significative snps 
    
    runs <- setdiff (unique (myrun), 0)
    
    ##report
    cat (paste ("########## Regiones seleccionadas cromosoma", cro, "##########"),
         file = repfile, append = TRUE, fill = TRUE)
    
    for (run in runs) {
      print (run)
      ##
      cat (paste ("##### Region:", run, "#####"),
           file = repfile,
           append = TRUE, fill = TRUE)
      ##
      write.table (resm2[[cro]][myrun == run, c("snp", "deCODE.genetic.map", "PhysicalPosition", "LOD", "ALPHA", "HLOD")],
                   file = repfile,
                   append = TRUE, quote = FALSE, sep = "\t",
                   row.names = FALSE, col.names = TRUE)
    }
  }
}
################################################################################

###merlin3
for (corte in c(0, 2)) {
  
  repfile <- file.path (.job$datadir, "data_results", "from_script_results", paste ("merlin3_lor", corte, ".txt", sep = ""))
  
  cat (paste ("########## Regiones seleccionadas con LOR:", corte, "##########"),
       file = repfile, append = TRUE, fill = TRUE)
  
  for (cro in cromosomas) {
    cat ("\n")
    print (cro)
    
    ##detectamos runs
    signf <- resm3[[cro]][,"LOD"] > corte       #significative
    prece <- c (FALSE, signf[-length (signf)])  #preceding
    start <- signf & !prece                     #starting poin
    cumul <- cumsum (start)                     #cumulative distribution
    myrun <- cumul * signf                      #run of significative snps 
    
    runs <- setdiff (unique (myrun), 0)
    
    ##report
    cat (paste ("########## Regiones seleccionadas cromosoma", cro, "##########"),
         file = repfile, append = TRUE, fill = TRUE)
    
    for (run in runs) {
      print (run)
      ##
      cat (paste ("##### Region:", run, "#####"),
           file = repfile,
           append = TRUE, fill = TRUE)
      ##
      write.table (resm3[[cro]][myrun == run, c("snp", "deCODE.genetic.map", "PhysicalPosition", "LOD", "ALPHA", "HLOD")],
                   file = repfile,
                   append = TRUE, quote = FALSE, sep = "\t",
                   row.names = FALSE, col.names = TRUE)
    }
  }
}
################################################################################

###EXIT
warnings ()
sessionInfo ()
q ("no")
