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
