library(testthat)

library(simplextree)
library(dart)

set.seed(1234)
R <- pbgrad::r_geometric_complex(16, radius = 0.25, dim = 2, coords = TRUE)

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