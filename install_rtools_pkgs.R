Sys.setenv(PATH=paste("E:/Rtools45/x86_64-w64-mingw32.static.posix/bin","E:/Rtools45/usr/bin",Sys.getenv("PATH"),sep=";"))
options(repos=c(CRAN="https://cloud.r-project.org/"))
install.packages("C:/Users/HONOR/AppData/Local/Temp/CellChat-main", repos=NULL, type="source")
install.packages("C:/Users/HONOR/AppData/Local/Temp/spacexr-master", repos=NULL, type="source")
for (p in c("CellChat","spacexr")) cat(sprintf("[%s] %s
", ifelse(suppressMessages(require(p, character.only=TRUE, quietly=TRUE)),"OK","FAIL"), p))
