# 用于安装R包们

install.packages("httpgd")
install.packages("harmony")
install.packages("pals")
install.packages("argparse")
# install.packages("freetype-devel")
# install.packages("libpng-devel")
# install.packages("libtiff-devel")
# install.packages("libjpeg-turbo-devel")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Seurat")
# install.packages("Seurat")

install.packages("glue")
install.packages("shiny")

install.packages("sitmo")

install.packages("spam")
install.packages("RcppEigen")
install.packages("curl")
install.packages("RcppTOML")
install.packages("reticulate")

install.packages("igraph")
install.packages("httr")
install.packages("RcppArmadillo")
install.packages("RcppAnnoy")
install.packages("RSpectra")
install.packages("dqrng")

install.packages("uwot")

BiocManager::install("Seurat")




# 
# library(shiny)
