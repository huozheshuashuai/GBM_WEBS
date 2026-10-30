options(timeout=1800, repos=c(CRAN="https://cloud.r-project.org/"))
deps <- c("parallelDist","dlm","MCMCpack","transport","NMF","ggalluvial","sna","ggnetwork","collapse")
install.packages(deps)
for (f in c("C:/Users/HONOR/AppData/Local/Temp/copykat.tar.gz", "C:/Users/HONOR/AppData/Local/Temp/spacexr.tar.gz", "C:/Users/HONOR/AppData/Local/Temp/CellChat.tar.gz")) {
  cat(">>> installing", f, "
")
  tryCatch(install.packages(f, repos=NULL, type="source"), error=function(e) cat("ERR:", e$message, "
"))
}
for (p in c("copykat","spacexr","CellChat","fastshap")) cat(sprintf("[%s] %s
", ifelse(suppressMessages(require(p, character.only=TRUE, quietly=TRUE)),"OK","FAIL"), p))
