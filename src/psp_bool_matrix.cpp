#include "PspMatrix.h"

#include <Rcpp.h>
using namespace Rcpp;


/*** R
ii <- c(3,5,3,4,5,2,3,6,3,4,5)-1L
jj <- c(1,1,2,3,3,4,4,5,6,6,6)-1L
x <- 1:11

# m <- new(phtools:::PspIntMatrix, 6, 6)
# m$construct(ii,jj,x)
# phtools:::test_spmatrix(ii, jj, x, m = 6, n = 6)
# Matrix::sparseMatrix(i=ii,j=jj,x=x,dims=c(6,6),index1=FALSE)
m <- phtools::psp_matrix(c(6,6), ii, jj, x)
m$add_columns(0,1)

p <- c(2,1,3,5,4,6)-1L
m$permute_rows(p)
*/
