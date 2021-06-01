
## Wide matrix column major 
M <- matrix(0, nrow = 5, ncol = 8)
indices <- t(phtools:::unrank_gridR(0:((5*8)-1), nr = 5, nc = 8, column_major = TRUE))+1L
M[indices] <- seq(5*8)

subscripts <- sapply(seq(5*8), function(i){ which(M == i, arr.ind = TRUE) })
phtools:::rank_gridR(subscripts, nc = 5, column_major = TRUE) # truth

## Wide matrix row major 
M <- matrix(0, nrow = 5, ncol = 8)
indices <- t(phtools:::unrank_gridR(0:((5*8)-1), nr = 5, nc = 8, column_major = FALSE))+1L
M[indices] <- seq(5*8)

subscripts <- sapply(seq(5*8), function(i){ which(M == i, arr.ind = TRUE) })
phtools:::rank_gridR(subscripts, nc = 8, column_major = FALSE) # truth

## Tall matrix column major 
M <- matrix(0, nrow = 8, ncol = 5)
indices <- t(phtools:::unrank_gridR(0:((5*8)-1), nr = 8, nc = 5, column_major = TRUE))+1L
M[indices] <- seq(5*8)

subscripts <- sapply(seq(5*8), function(i){ which(M == i, arr.ind = TRUE) })
phtools:::rank_gridR(subscripts, nc = 8, column_major = TRUE) # truth

## Tall matrix row major 
M <- matrix(0, nrow = 8, ncol = 5)
indices <- t(phtools:::unrank_gridR(0:((5*8)-1), nr = 8, nc = 5, column_major = FALSE))+1L
M[indices] <- seq(5*8)

subscripts <- sapply(seq(5*8), function(i){ which(M == i, arr.ind = TRUE) })
phtools:::rank_gridR(subscripts, nr = 8, nc = 5, column_major = FALSE) # truth


all(dart::unrank_combn(seq(choose(30,3)), n = 30, k = 3) == combn(30,3))


all(unrank_combn(rank_combn(combn(10,2), n = 10), n = 10, k = 2) == combn(10,2))
rank_combn(unrank_combn(seq(choose(10,2)), n = 10, k = 2), n = 10)