options(repos=c(CRAN="https://cloud.r-project.org/"))
install.packages(c("pals","fields"))
tryCatch(install.packages("locfdr"), error=function(e) cat("locfdr CRAN fail
"))
if (!require("locfdr", quietly=TRUE)) {
  tryCatch(install.packages("https://cran.r-project.org/src/contrib/Archive/locfdr/locfdr_1.1-16.tar.gz", repos=NULL, type="source"), error=function(e) cat("archive fail
"))
}
install.packages("C:/Users/HONOR/AppData/Local/Temp/spacexr-master", repos=NULL, type="source")
cat(ifelse(suppressMessages(require("spacexr", quietly=TRUE)),"[OK] spacexr","[FAIL] spacexr"),"
")
