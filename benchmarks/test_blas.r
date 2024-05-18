#library(RhpcBLASctl)
#threads <- as.integer(Sys.getenv("NOMPT"))
#omp_set_num_threads(threads)
#blas_set_num_threads(threads)
n=3000
p=1000
q=5000
#cat(sprintf("Using %d threads:\n",threads))
for (iter in 1:3) {
  cat(sprintf("iteration %d\n",iter))
  A = matrix(runif(n*p),nrow=n, ncol=p)
  B = matrix(runif(p*q),nrow=p, ncol=q)
  C = A %*% B # multi-threaded matrix product
  D=svd(C)
}
quit(save='no')
