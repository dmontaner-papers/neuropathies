## .job.r
## 2010-03-09 dmontaner@cipf.es
## script that keeps the settings for each job
## bear in mind that it is dependent on the data folder structure

##  .OWNER: name of the people for whom the analysis is done.
##          It will be used as the JOB_NAME (no spaces, lower case, etc...)

## datadir: path to the directory where all data (raw, generated and results) are stored.
##          The last directory in the path has to be JOB_NAME.

## docsdir: path to the directory where scripts, sample information and documentation are stored.
##          The last directory in the path has to be JOB_NAME.
##          It can be the same or different than the datadir.

################################################################################

.OWNER = "eduardo_calpena" #nombre de quien encarga el analisis (sin espacios)

################################################################################

.job <- list ()
.job$owner <- .OWNER

## ### rootDir Location in SHARED DISK analisis1
## .job$datadir <- file.path ("", "mnt", "analisis", "analisis1", "trabajos-2010", .job$owner) #starting with ("") to get an absolute path
## .job$docsdir <- file.path ("", "mnt", "analisis", "analisis1", "trabajos-2010", .job$owner) #starting by /mnt

### rootDir Location in MY COMPUTER
.job$datadir <- file.path ("~", "datos",    "2011", .job$owner) #starting with ("~") if working in your home directory
.job$docsdir <- file.path ("~", "trabajos", "2011", .job$owner) #or ("") if working in the root directory

###MORE
.job$scriptdir <- file.path (.job$docsdir, "scripts") 
.job$sinfodir  <- file.path (.job$docsdir, "sampleinfo")

.job$plotsdir  <- file.path (.job$datadir, "data_results", "from_script_plots")

.job$testmode <- FALSE ##para pruebas. No siempre se usa (descomenta)
.job$dec <- "."

################################################################################

rm (list = ".OWNER")

##MESSAGE
cat ("\n")
message (".job.r has been sourced")
cat ("\n")
