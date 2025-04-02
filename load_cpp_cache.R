# load cpp code to cache
library(Rcpp)
library(RcppArmadillo)
sourceCpp("cpp_func_cls.cpp", cacheDir = "/scratch/ylong5/rcpp-cache/robust-fda")
quit()