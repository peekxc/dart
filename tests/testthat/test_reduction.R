library(testthat)

library(simplextree)
library(dart)

## Test bare bones reduction 
set.seed(1234)
st <- simplextree::simplex_tree(combn(3,2)) %>% simplextree::expand(k=2)
D <- boundary_matrix(st)
D_psp <- psp_matrix(x=D$matrix)
V_psp <- psp_matrix(x=Matrix::Diagonal(ncol(D$matrix)))
dart:::reduce_pspbool(D_psp$matrix$as_XPtr(), V_psp$matrix$as_XPtr(), show_progress = FALSE)
D_psp
V_psp

## cohomology 
# db <- boundary_matrix(st, labeled=TRUE)$matrix
# dc <- t(db[rev(seq(nrow(db))), rev(seq(ncol(db)))])

## Test dimension specific reduction
set.seed(1234)
st <- simplextree::simplex_tree(combn(3,2)) %>% simplextree::expand(k=2)
D <- boundary_matrix(st, dims = c(1,2))
D1_psp <- psp_matrix(x=D$matrix[[1]])
D2_psp <- psp_matrix(x=D$matrix[[2]])
V1_psp <- psp_matrix(x=Matrix::Diagonal(ncol(D$matrix[[1]])))
V2_psp <- psp_matrix(x=Matrix::Diagonal(ncol(D$matrix[[2]])))
dart:::reduce_local_pspbool(
	D1_psp$matrix$as_XPtr(), V1_psp$matrix$as_XPtr(), 
	D2_psp$matrix$as_XPtr(), V2_psp$matrix$as_XPtr(),
	clearing = TRUE, show_progress = FALSE
)
D1_psp
D2_psp
V1_psp
V2_psp


set.seed(1234)
st <- simplextree::simplex_tree(combn(3,2)) %>% simplextree::expand(k=2)
D <- boundary_matrix(st, dims = c(1,2))
D1_psp <- psp_matrix(x=D$matrix[[1]])
D2_psp <- psp_matrix(x=D$matrix[[2]])
V1_psp <- psp_matrix(x=Matrix::Diagonal(ncol(D$matrix[[1]])))
V2_psp <- psp_matrix(x=Matrix::Diagonal(ncol(D$matrix[[2]])))
dart:::reduce_local_pspbool(
	D1_psp$matrix$as_XPtr(), V1_psp$matrix$as_XPtr(), 
	D2_psp$matrix$as_XPtr(), V2_psp$matrix$as_XPtr(),
	clearing = TRUE, show_progress = FALSE
)
D1_psp
D2_psp
V1_psp
V2_psp

set.seed(1234)
st <- simplextree::simplex_tree(combn(3,2)) %>% simplextree::expand(k=2)
D <- boundary_matrix(st)
rv <- reduce(D, validate = FALSE)
R <- rv$R$as.Matrix()
dv <- ((D$matrix %% 2) %*% rv$V$as.Matrix()) %% 2
testthat::expect_true(all(R == dv))

rv <- reduce(D, options = c(clearing = FALSE), validate = FALSE)
R <- rv$R$as.Matrix()
dv <- as(((D$matrix %% 2) %*% rv$V$as.Matrix()) %% 2, "lgCMatrix")
testthat::expect_true(all(R == dv))

# all(rank_combn(rf$simplices, n) == rf$ranks)

set.seed(1234)
x <- replicate(2, runif(8))
rf <- rips_filtration(x, diameter = 0.30, dim = 2)
D <- boundary_matrix(rf$complex, dim = c(1,2))
rv <- reduce(D, options = c(clearing=FALSE), validate = FALSE)

r1 <- as(rv$R[[1]], "sparseMatrix")
r2 <- as(rv$R[[2]], "sparseMatrix")
d1 <- as(D$matrix[[1]], "lgCMatrix")
d2 <- as(D$matrix[[2]], "lgCMatrix")
v1 <- as(rv$V[[1]], "sparseMatrix")
v2 <- as(rv$V[[2]], "sparseMatrix")
all(r1 == as((d1 %*% v1) %% 2, "lgCMatrix"))

D <- boundary_matrix(rf, dim = c(1,2))
rv <- reduce(D, options = c(clearing=FALSE), validate = FALSE)
r1 <- as(rv$R[[1]], "sparseMatrix")
r2 <- as(rv$R[[2]], "sparseMatrix")
d1 <- as(D$matrix[[1]], "lgCMatrix")
d2 <- as(D$matrix[[2]], "lgCMatrix")
v1 <- as(rv$V[[1]], "sparseMatrix")
v2 <- as(rv$V[[2]], "sparseMatrix")
all(r1 == as((d1 %*% v1) %% 2, "lgCMatrix"))


D <- boundary_matrix(rf$complex, dim = c(1,2))
rv <- reduce(D, options = c(clearing=FALSE), validate = FALSE)
dart::validate_decomp(rv, D = D)

D <- boundary_matrix(rf$complex, dim = c(1,2))
rv <- reduce(D, options = c(clearing=TRUE), validate = FALSE)
dart::validate_decomp(rv, D = D)

D <- boundary_matrix(rf, dim = c(1,2))
rv <- reduce(D, options = c(clearing=FALSE), validate = FALSE)
dart::validate_decomp(rv, D = D)

D <- boundary_matrix(rf)
rv <- reduce(D, options = c(clearing=TRUE), validate = FALSE)
dart::validate_decomp(rv, D = D)

## Test larger example 
set.seed(1234)
x <- replicate(2, runif(16))
rf <- rips_filtration(x, diameter = 0.30, dim = 2)
D <- boundary_matrix(rf, c(1,2))
rv <- reduce(D, options = c(clearing=FALSE), validate = TRUE)
all(validate_decomp(rv, D))

## Test larger example still
set.seed(1234)
test_12 <- sapply(seq(100), function(i){
	x <- replicate(2, runif(32))
	rf <- rips_filtration(x, diameter = 0.30, dim = 2)
	D <- boundary_matrix(rf, c(1,2))
	rv <- reduce(D, options = c(clearing=TRUE), validate = TRUE)
	all(validate_decomp(rv, D))
})
all(test_12)


x <- replicate(2, runif(100))
rf <- rips_filtration(x, diameter = 0.30, dim = 2)
D <- boundary_matrix(rf, c(1,2))
rv <- reduce(D, options = c(clearing=TRUE), validate = TRUE, show_progress = TRUE)
all(validate_decomp(rv, D))
	
# dart::pivots(rv)


## Test moves 
library(simplextree)
st <- simplex_tree(combn(3,2)) %>%  expand(k = 2)
D <- boundary_matrix(st, dims = c(1, 2), labeled = "both")
rv <- reduce(D, validate = TRUE)
as(rv$R[[1]], "sparseMatrix")
as(rv$R[[2]], "sparseMatrix")
as(rv$V[[1]], "sparseMatrix")
as(rv$V[[2]], "sparseMatrix")
dart:::move_schedule_local(r1 = rv$R[[1]]$matrix$as_XPtr(), v1 = rv$V[[1]]$matrix$as_XPtr(), 
													 r2 = rv$R[[2]]$matrix$as_XPtr(), v2 = rv$V[[2]]$matrix$as_XPtr(), 
													 schedule = c(1,2), f = NULL)
as(rv$R[[1]], "sparseMatrix")
as(rv$R[[2]], "sparseMatrix")
as(rv$V[[1]], "sparseMatrix")
as(rv$V[[2]], "sparseMatrix")

## Test moves 
library(simplextree)
st <- simplex_tree(combn(3,2)) %>%  expand(k = 2)
D <- boundary_matrix(st, dims = c(1, 2), labeled = "both")
rv <- reduce(D, validate = TRUE)
as(rv$R[[1]], "sparseMatrix")
as(rv$R[[2]], "sparseMatrix")
as(rv$V[[1]], "sparseMatrix")
as(rv$V[[2]], "sparseMatrix")
dart:::move_schedule_local(r1 = rv$R[[1]]$matrix$as_XPtr(), v1 = rv$V[[1]]$matrix$as_XPtr(), 
													 r2 = rv$R[[2]]$matrix$as_XPtr(), v2 = rv$V[[2]]$matrix$as_XPtr(), 
													 schedule = c(0,1), f = NULL)
as(rv$R[[1]], "sparseMatrix")
as(rv$R[[2]], "sparseMatrix")
as(rv$V[[1]], "sparseMatrix")
as(rv$V[[2]], "sparseMatrix")


## Test all possible move rights within a simplex
library(simplextree)
st <- simplex_tree(combn(5,2)) %>% simplextree::expand(k = 2)
moves <- combn(st$n_simplices[2],2)
for (i in seq(ncol(moves))){
	D <- boundary_matrix(st, dims = c(1, 2), labeled = "both")
	rv <- reduce(D, validate = FALSE)
	original_valid <- all(validate_decomp(rv, D))
	
	R1 <- as(rv$R[[1]], "sparseMatrix")
	R2 <- as(rv$R[[2]], "sparseMatrix")
	V1 <- as(rv$V[[1]], "sparseMatrix")
	V2 <- as(rv$V[[2]], "sparseMatrix")
	
	dart:::move_schedule_local(
		r1 = rv$R[[1]]$matrix$as_XPtr(), v1 = rv$V[[1]]$matrix$as_XPtr(),
		r2 = rv$R[[2]]$matrix$as_XPtr(), v2 = rv$V[[2]]$matrix$as_XPtr(), 
		schedule = moves[,i]-1L, f = NULL
	)
	
	rv_correct <- move_right(list(R1, R2), list(V1, V2), i = moves[1,i], j = moves[2,i])
	# all(as(rv$V[[1]], "sparseMatrix") == rv_correct$V[[1]])
	
	r1_new <- as(rv$R[[1]], "sparseMatrix")
	r2_new <- as(rv$R[[2]], "sparseMatrix")
	v1_new <- as(rv$V[[1]], "sparseMatrix")
	v2_new <- as(rv$V[[2]], "sparseMatrix")
	
	D$matrix[[1]] <- permute_move(D$matrix[[1]], i = moves[1,i], j = moves[2,i], dims = "cols")
	D$matrix[[2]] <- permute_move(D$matrix[[2]], i = moves[1,i], j = moves[2,i], dims = "rows")
	new_valid <- all(validate_decomp(rv, D))
	print(sprintf("move right (%d, %d), original: %s, new: %s", moves[1,i], moves[2,i], as.character(original_valid), as.character(new_valid)))
	if (!new_valid){
		break
		# (D$matrix[[1]] %*% as(rv$V[[1]], "sparseMatrix")) %% 2
		# (D$matrix[[1]][1:2,] %% 2) %*% as(rv$V[[1]], "sparseMatrix")[,6,drop=FALSE]
		# as(rv$R[[1]], "sparseMatrix")
	}
	R1 <- as(rv$R[[1]], "sparseMatrix")
	D1 <- as(D$matrix[[1]], "lgCMatrix")
	V1 <- as(rv$V[[1]], "sparseMatrix")
	to_check <- as((D1 %*% V1) %% 2, "lgCMatrix")
}



as(rv$R[[1]], "sparseMatrix")
as(rv$R[[2]], "sparseMatrix")
as(rv$V[[1]], "sparseMatrix")
as(rv$V[[2]], "sparseMatrix")

as(rv$R[[1]], "sparseMatrix")
as(rv$R[[2]], "sparseMatrix")
as(rv$V[[1]], "sparseMatrix")
as(rv$V[[2]], "sparseMatrix")


set.seed(1234)
R <- dart::r_rips_complex(75, radius = 0.35, dim = 2, coords = TRUE)
rf <- filtration(R)
D <- boundary_matrix(rf, labeled = FALSE)
rv <- reduce(D, options = c(clearing=TRUE), validate = TRUE)

set.seed(1234)
R <- dart::r_rips_complex(16, radius = 0.25, dim = 2, coords = TRUE)

## Option 1 
{
	D <- dart::boundary_matrix(as.list(level_order(R)))
}

## Option 2
{
	dx <- parallelDist::parallelDist(attr(R, "coords"))
	W <- dx[rankr::rank_comb(t(R$edges))]
	FI <- new(dart:::ImplicitFiltration, R$as_XPtr(), W)
	D <- dart:::boundary_matrix_fi_full(FI$as_XPtr())
}

## Test move operations
st <- simplextree::simplex_tree(combn(3,2)) %>% expand(k=2)
RV <- dart::reduce(D)

dart:::simulate_vineyard_cpp(RV$R$matrix$as_XPtr(), V_ptr = RV$V$matrix$as_XPtr(), schedule = matrix(c(16,17,18,19), nrow = 2))


D_psp <- dart::psp_matrix(x = D$matrix)
V_psp <- dart::psp_matrix(x = Matrix::Diagonal(ncol(D$matrix)))
dart:::reduce_pspbool(D_psp$matrix$as_XPtr(), V_psp$matrix$as_XPtr())
dart::validate_decomp(RV)




set.seed(1234)
R <- pbgrad::r_geometric_complex(n = 30, radius = 0.25, dim = 2)

dart:::boundary_matrix_fi()

D <- dart::boundary_matrix(as.list(level_order(R)))
RV <- dart::reduce(D)

D_psp <- dart::psp_matrix(x = D$matrix)
V_psp <- dart::psp_matrix(x = Matrix::Diagonal(ncol(D$matrix)))
dart:::reduce_pspbool(D_psp$matrix$as_XPtr(), V_psp$matrix$as_XPtr())