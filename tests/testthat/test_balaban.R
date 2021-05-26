

yb <- sort(runif(250))
ye <- runif(250)


dart::kendall_dist(seq_along(yb), order(ye))


plot(rbind(cbind(0, yb), cbind(1,ye)), pch = 20)
for (i in seq_along(yb)){
	segments(x0 = 0, y0 = yb[i], x1 = 1, y1 = ye[i])
}
points(t(int_points), col = "red", pch = 20, cex = 0.2)



microbenchmark::microbenchmark({
	int_points <- sf::st_intersection(L, L)
}, setup = {
	L <- sf::st_sfc(lapply(seq_along(yb), function(i){
		sf::st_linestring(matrix(c(0,yb[i],1,ye[i]), nrow = 2, byrow = TRUE))
	}))
}, times = 5L)
# int_points <- sf::st_collection_extract(intersections, "POINT")
# int_points <- unique(do.call(rbind, lapply(int_points, function(pt){ as.numeric(pt) })))

microbenchmark::microbenchmark({
  base_intersections <- combn(length(yb), 2, function(x){ xor(yb[x[1]] < yb[x[2]], ye[x[1]] < ye[x[2]]) })
	int_pairs <- simplextree::nat_to_sub(which(base_intersections), n = length(yb), k = 2)
	int_points <- dart:::span_intersections(int_pairs, yb, ye, 0, 1)
})
microbenchmark::microbenchmark({
	dart:::SearchInStrip(yb, ye, 0, 1)
})

microbenchmark::microbenchmark({
	indices <- dart:::SearchInStrip(yb, ye, 0, 1)
	int_points <- dart:::span_intersections(indices, yb, ye, 0, 1)
})
	
indices <- dart:::SearchInStrip(yb, ye, 0, 1)
int_points <- dart:::span_intersections(indices, yb, ye, 0, 1)

