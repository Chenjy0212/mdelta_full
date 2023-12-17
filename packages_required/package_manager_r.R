# 定义要使用的包列表
required_packages <- c("rjson", "BiocManager")

# 检查每个包是否已经安装
missing_packages <- setdiff(required_packages, installed.packages())

# 如果有缺失的包，则进行自动安装
if (length(missing_packages) > 0) {
  install.packages(missing_packages, , repos = "http://cran.us.r-project.org")
}

# 定义要使用的包列表
required_packages_bioc <- c("ggtree", "ggtreeExtra")

# 检查每个包是否已经安装
missing_packages_bioc <- setdiff(required_packages_bioc, installed.packages())

# 如果有缺失的包，则进行自动安装
if (length(missing_packages_bioc) > 0) {
  BiocManager::install(missing_packages_bioc)
}
# 加载所需的包
# lapply(required_packages, library, character.only = TRUE)

# if(!require("rjson", quietly=TRUE)){
# 	install.packages("rjson", repos = "http://cran.us.r-project.org")
# }
# if (!require("BiocManager", quietly = TRUE)){
#     install.packages("BiocManager", repos = "http://cran.us.r-project.org")
# }
# if(!require("ggtree", quietly=TRUE)){
# 	BiocManager::install("ggtree")
# }
# if(!require("ggtreeExtra", quietly=TRUE)){
# 	BiocManager::install("ggtreeExtra")
# }

