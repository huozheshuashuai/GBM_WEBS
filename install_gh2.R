options(timeout=1800, repos=c(CRAN="https://cloud.r-project.org/"))
tryCatch(install.packages("C:/Users/HONOR/AppData/Local/Temp/spacexr.tar.gz", repos=NULL, type="source"), error=function(e) cat("spacexr ERR:", e$message, "
"))
tryCatch(install.packages("C:/Users/HONOR/AppData/Local/Temp/CellChat.tar.gz", repos=NULL, type="source"), error=function(e) cat("CellChat ERR:", e$message, "
"))
for (p in c("spacexr","CellChat")) cat(sprintf("[%s] %s
", ifelse(suppressMessages(require(p, character.only=TRUE, quietly=TRUE)),"OK","FAIL"), p))
