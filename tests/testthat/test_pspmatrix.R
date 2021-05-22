A <- phtools::psp_matrix(dims = c(10,10))


A$insert(5,5, 1)

testthat::test_that("Matrix is valid", {
	all(A$otc[A$cto+1] == seq(A$n_rows)-1L)
})

testthat::test_that("Can remove zeroed entries", {
	A <- phtools::psp_matrix(dims = c(3,3))
	A$insert(0,0, TRUE)
	A$insert(1,0, TRUE)
	A$add_columns(1,0)
	A$add_columns(1,0)
	A$remove(1,1)
})

testthat::test_that("Can clean zeroed entries", {
	A <- phtools::psp_matrix(dims = c(5,5))
	A$insert(3,3,TRUE)
	A$insert(2,3,TRUE)
	A$add_columns(4,3)
	A$add_columns(4,3)
	A$as.Matrix()
	
	testthat::expect_equal(A$nnz, 2L)
	A$print()
	A$clean(FALSE)
	
})

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

