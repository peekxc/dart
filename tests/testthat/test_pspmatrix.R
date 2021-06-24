testthat::test_that("can construct a matrix", {
	testthat::expect_silent(dart::psp_matrix(dims = c(10,10)))
	A <- dart::psp_matrix(dims = c(10,10))
	testthat::expect_true("PspMatrix" %in% class(A))
})

testthat::test_that("Matrix is valid", {
	A <- dart::psp_matrix(dims = c(10,10))
	all(A$otc[A$cto+1] == seq(A$n_rows)-1L)
})

testthat::test_that("Can remove zeroed entries", {
	A <- dart::psp_matrix(dims = c(3,3))
	A$matrix$insert(0,0, TRUE)
	A$matrix$insert(1,0, TRUE)
	A$matrix$add_columns(1,0)
	A$matrix$add_columns(1,0)
	A$matrix$remove(1,1)
})

testthat::test_that("Can clean zeroed entries", {
	A <- dart::psp_matrix(dims = c(5,5))
	A$matrix$insert(3,3,TRUE)
	A$matrix$insert(2,3,TRUE)
	A$matrix$add_columns(4,3)
	A$matrix$add_columns(4,3)
	testthat::expect_equal(A$matrix$nnz, 4L)
	testthat::expect_silent(A$matrix$clean(0))
	testthat::expect_equal(A$matrix$nnz, 2L)
})

A <- as(Matrix::rsparsematrix(10,15,density = 0.35), "lgCMatrix")
P <- psp_matrix(x = A) 
for (i in seq(10)){
	pr <- sample(seq(nrow(A)))
	P$matrix$permute_rows(pr-1L)
	print(all(A[pr,] == P$as.Matrix()))
	A <- A[pr,]
}

low_A <- apply(A, 2, dart:::low_entry)
low_P <- sapply(seq(ncol(A))-1L, function(j){ P$matrix$low_entry(j) }) + 1L
low_A == low_P


testthat::test_that("Can do random operations", {
	M <- Matrix::rsparsematrix(10,10,density = 0.35)
	ij <- which(M != 0, arr.ind = TRUE)
	A <- phtools::psp_matrix(dims = dim(M), i = ij[,1]-1L, j = ij[,2]-1L, x = rep(1, nrow(ij)))
	M <- M != 0

	testthat::expect_true(all(M == A$as.Matrix()))
	
	
	swaps <- combn(nrow(M), 2)[, sample(x = seq(choose(nrow(M), 2)), size = 15, replace = TRUE)]
	for ( i in seq(ncol(swaps)) ){
		if (i %% 2 == 0){
			s <- swaps[,i]
			A$swap_rows(s[1]-1L, s[2]-1L)
			A$swap_cols(s[1]-1L, s[2]-1L)
			M[s,] <- M[rev(s),]
			M[,s] <- M[,rev(s)]
			testthat::expect_true(all((M != 0) == A$as.Matrix()))
			
			A_prev <- A$as.Matrix()
			M_prev <- M
			A$add_columns(s[2]-1L, s[1]-1L)
			M[,s[2]] <- xor(M[,s[2]], M[,s[1]])
			testthat::expect_true(all((M != 0) == A$as.Matrix()))
		}
	}
	
	
	A$insert(3,3,TRUE)
	A$insert(2,3,TRUE)
	A$add_columns(4,3)
	A$add_columns(4,3)
	A$as.Matrix()
	
	testthat::expect_equal(A$nnz, 2L)
	A$print()
	A$clean(FALSE)
	
})

