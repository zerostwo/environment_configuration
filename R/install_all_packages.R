r_packages <- read.table("./R/R_packages_list.txt")[, 1]
root <- list.files(.libPaths()[1])
print(root)
need2install <- r_packages[!(r_packages %in% root)]
print(need2install)
for (i in need2install) {
  print(i)
  BiocManager::install(i)
}
