#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#define INF 100000000
#define MAXN 1024

typedef vector<int> VI;
typedef pair<int,int> PII;
typedef vector<PII> VPII;

char used[MAXN+1];
int dp[MAXN+1][MAXN+1]; 
int previous[MAXN+1][MAXN+1];

VPII simulate(char *used, VI &p) {
	int n = p.size();
	int cost = 0;
	VPII moves;

	for (int i=n-1; i>=0; i--) {
		if (used[i])
			continue;

		// find element i
		VI::iterator it = find(p.begin(),p.end(),i);
		assert(it != p.end());

		int from_pos = distance(p.begin(),it) + 1;
		p.erase(it);

		// find element i+1
		it = find(p.begin(),p.end(),i+1);

		// reinsert i before element i+1
		p.insert(it, i);

		int to_pos = distance(p.begin(),it) + 1;
		cost += from_pos + to_pos;
		moves.push_back(PII(from_pos,to_pos));
	}

	for (int i=0; i<n; i++)
		assert(p[i] == i);

	// print the result
	Rprintf("%d\n",(int)moves.size());
	for (VPII::iterator it=moves.begin(); it!=moves.end(); ++it)
		Rprintf("%d %d\n",it->first,it->second);
	// return cost;
	return(moves);
}

VPII solve(VI &p) {
	int n = p.size();
	VI pos(n+1);

	p.push_back(n);
	for (int i=0; i<=n; i++)
		pos[p[i]] = i;
	for (int i=0; i<n; i++)
		dp[n][i] = INF;
	dp[n][n] = 0;

	for (int i=n-1; i>=0; --i) {
		int v = 1;

		// determine the current position v, assuming all bigger elements already moved
		// we can only make this assumption because we adjust the costs of not moving an
		// element accordingly (see below)
		for (int j=0; j<pos[i]; ++j)
			if (p[j]<i)
				++v;

		int v2 = 1;
		// move i before element i+1
		// h is the old index position where it moves to
		// v2 is the "real" position where it moves to
		for (int h=0; h<=n; h++) {
			if (p[h] < i)
				++v2;
			dp[i][h] = dp[i+1][h]+v+v2;
			previous[i][h] = h;
		}

		// try to let i at its current position
		// add specifies the additional costs for elements which are skipped
		int add = 0;
		for (int k=pos[i]+1; k<=n; k++) {
			int tcost = dp[i+1][k]+add;
			if (dp[i][pos[i]] > tcost) {
				dp[i][pos[i]] = tcost;
				previous[i][pos[i]] = k;
			}
			if (p[k]<i) {
				// there are i - p[k] elements which will be moved before it
				// since they are all bigger than p[k] they wouldn't be counted
				// when determining the current position
				add += (i-p[k]);
			}
		}
	}

	int mincost = INF,ind = -1;
	for (int i=0; i<=n; i++)
		if (dp[0][i] < mincost) {
			mincost = dp[0][i];
			ind = i;
		}
	assert(ind >= 0);

	// determine which numbers kept their position
	for (int i=0; i<n; ++i) {
		used[i] = (ind == pos[i])?1:0;
		ind = previous[i][ind];
	}
	assert(ind == n);

	used[n] = 1;
	// int tcost = simulate(used,p);
	//assert(tcost == mincost);
	// Rcout << "cost: " << tcost << std::endl;
	return simulate(used,p);
}

// [[Rcpp::export]]
IntegerMatrix test_moves(IntegerVector scores){
	const int n = scores.size(); 
	map<int,int> byscore;
	for (int i=0; i<n; ++i) {
		byscore[scores[i]] = i;
	}

	int pos = 0; // position counter
	VI p(n); // desired final positions
	for (map<int,int>::reverse_iterator it=byscore.rbegin(); it!=byscore.rend(); ++it,++pos)
		p[it->second] = pos;
	VPII moves = solve(p);
	
	// for (auto pp: p){
	// 	Rcout << pp << ","; 
	// 	// Rcout << p.first << ", " << p.second << std::endl;
	// }
	
	IntegerMatrix result(2, moves.size());
	for (size_t i = 0; i < moves.size(); ++i){
		result(0,i) = moves[i].first; 
		result(1,i) = moves[i].second; 
	}
	
	return result; 
}

/*** R

x <- sample(1:15)
M <- phtools:::test_moves(x)

nm <- length(x) - length(phtools:::longest_subseq(x))
num_crossing(M)
y <- x
for (i in seq(ncol(M))){
	 y <- permute_move(y, M[1,i], M[2,i]) 
}
src <- x
tgt <- rev(seq(length(x)))
lcs <- tgt[phtools:::longest_subseq(match(tgt, src))]
M_strats <- move_sequence(symbols = setdiff(x, lcs), s = src, t = tgt, lis = lcs, ordered = TRUE)

sapply(M_strats, function(x){ num_crossing(x$moves) })



x_dist <- combn(120,2, function(x){ kernrank::kendall_weight(P[x[1],], P[x[2],], u =c(1,2,3,4,5), method = "mult") })
D[t(combn(120,2))] <- x_dist 
D <- as.dist(kernrank::AllKendall(P, P))


xy <- cmdscale(D)
text(xy, labels = apply(P, 1, function(s){ paste0(s, collapse = ",")}), pos = 3, cex = 0.25)


bubble_sort <- function(A){
	P <- matrix(0, ncol = 0, nrow = 2)
	n <- length(A)
	while(!all(order(A) == seq(n))){
		for (i in seq(n-1)){
			if (A[i] > A[i+1]){
				A[c(i, i+1)] <- A[c(i+1, i)]
				P <- cbind(P, c(i, i+1))
			}
		}
	}
	return(P)
}

## Kendall-weighted distance determining transposition order
P <- rankr::permn(4)
apply(P, 2, permutations::as.word)

D <- as.dist(kernrank::AllKendall(P, P))
D <- as.matrix(D)

# D[D > 1] <- 0
G <- igraph::graph_from_adjacency_matrix(as.matrix(D), weighted = TRUE, mode = "undirected")
S <- P[sample(1:factorial(4), size = 2),]
i <- rankr::rank_permutation(S[1,])
j <- rankr::rank_permutation(S[2,])
sp <- igraph::shortest_paths(G, from = i, to = j, weights = NULL)




ij <- t(combn(24,2)[,combn(24,2, function(x){ (x[1] < x[2]) & abs(diff(x)) == 1 })])

spearman_dist <- function(p, q=seq_along(p)){
	stopifnot(is.vector(p), all(seq_along(p) %in% p))	
	I <- seq_along(p)
	return(sum(abs(match(I, p) - match(I, q))))
}

move_permutation <- function(i,j,n){ return(permute_move(seq(n), i, j)) }

x <- seq(8)
y <- sample(x)

lis <- y[phtools:::longest_subseq(match(y, x))]
MS <- move_sequence(symbols = setdiff(x,lis), s = x, t = y, lis = lis, ordered = FALSE, rule = "all")

## Applies a sequence of permutations given by the columns of 'P' to 'p' successively
compose <- function(p, M){
	for (i in seq(ncol(M))){ p <- permute_move(p, i = M[1,i], j = M[2,i]) }
	return(p)
}
move_p <- apply(MS[[1]]$moves, 2, function(p){
	move_permutation(i = p[1], j = p[2], n = length(x))
})

schedule_cost <- function(p, M, cost=c("spearman", "kendall")){
	if (ncol(M) == 0){ return(0) }
	cost_f <- ifelse(missing(cost) || cost == "spearman", spearman_dist, kendall_dist)
	crossing_cost <- vector("numeric", length= ncol(M))
	# s1 <- p
	# s2 <- permute_move(p, i = M[1,1], j = M[2,1])
	# crossing_cost[1] <- cost_f(s1, s2)
	# for (i in 1:(ncol(M)-1)){
	# 	s1 <- permute_move(s1, M[1,i], M[2,i])
	# 	s2 <- permute_move(s2, M[1,i+1], M[2,i+1])
	# 	crossing_cost[i+1] <- cost_f(s1, s2)
	# }
	crossing_cost <- sapply(1:(ncol(M) - 1), function(i){
		s1 <- compose(p, M[,1:i])
		s2 <- compose(p, M[,1:(i+1)])
		cost_f(s1, s2)
	})
	return(sum(crossing_cost))
}

## Footrule distance *not* invariant to the order of symbols to move
MS <- move_sequence(symbols = setdiff(x,lis), s = x, t = y, lis = lis, ordered = FALSE, rule = "all")
sapply(MS, function(m){ schedule_cost(x, m$moves) })

## Footrule distance + Kendall distance w/ fixed heuristic *not* invariant to ordering
MS <- move_sequence(symbols = setdiff(x,lis), s = x, t = y, lis = lis, ordered = FALSE, rule = "closest")
sapply(MS, function(m){ schedule_cost(x, m$moves) })

## Footrule distance + Kendall distance w/ fixed order *not* invariant to heuristic 
to_move <- setdiff(x,lis)
TM <- rankr::permn(to_move)
apply(TM, 2, function(tm){
	MS <- move_sequence(symbols = tm, s = x, t = y, lis = lis, ordered = TRUE, rule = "all")
	length(unique(sapply(MS, function(m){ schedule_cost(x, m$moves) }))) == 1
})

x <- seq(54)
y <- sample(x)
lis <- y[phtools:::longest_subseq(match(y, x))]

## Footrule distance w/ fixed order --- 'latest' heuristic is always better
to_move <- setdiff(x,lis)
# TM <- rankr::permn(to_move)
TM <- apply(t(as.matrix(permutations::rperm(150, r = length(to_move)))), 2, function(p){ to_move[p] })
heuristic_costs <- apply(TM, 2, function(tm){
	rules <- c("earliest", "latest", "closest", "furthest")
	costs <- sapply(rules, function(rr){
		MS <- move_sequence(symbols = tm, s = x, t = y, lis = lis, ordered = TRUE, rule = rr)
		schedule_cost(x, MS[[1]]$moves)
	})
	costs
})

sapply(seq(4), function(i){
	other <- setdiff(seq(4), i)
	all(sapply(other, function(o){
		all(heuristic_costs[i,] <= heuristic_costs[o,])
	}))
})

all(heuristic_costs[2,] <= heuristic_costs[4,])
all(heuristic_costs[2,] <= heuristic_costs[3,])


# tm <- TM[,sample(seq(ncol(TM)), size = 1)]
tm <- to_move[sample(seq(length(to_move)))]

## Plot heuristic strategies 
plot_moves <- function(p, M, ...){
	plot.default(NULL, NULL, xlim = c(0, length(p)+1), ylim = c(0, ncol(M)+1), ...)
	for (i in seq(ncol(M))){
		segments(x0 = M[1,i], y0 = i, x1 = M[2,i], y1 = i)
	}
	for (i in seq(ncol(M))){
		if (M[1,i] < M[2,i]){
			points(cbind(M[1,i], i), pch = 4)
			points(cbind(M[2,i], i), pch = 20)
		} else {
			points(cbind(M[1,i], i), pch = 4)
			points(cbind(M[2,i], i), pch = 20)
		}
	}
}
closest_M <- move_sequence(symbols = tm, s = x, t = y, lis = lis, ordered = TRUE, rule = "closest")[[1]]$moves
earliest_M <- move_sequence(symbols = tm, s = x, t = y, lis = lis, ordered = TRUE, rule = "earliest")[[1]]$moves
latest_M <- move_sequence(symbols = tm, s = x, t = y, lis = lis, ordered = TRUE, rule = "latest")[[1]]$moves
furthest_M <- move_sequence(symbols = tm, s = x, t = y, lis = lis, ordered = TRUE, rule = "furthest")[[1]]$moves

layout(matrix(1:4, ncol = 2))
plot_moves(x, closest_M, main = sprintf("closest (%d)", schedule_cost(x, closest_M)))
plot_moves(x, earliest_M, main = sprintf("earliest (%d)", schedule_cost(x, earliest_M)))
plot_moves(x, latest_M, main = sprintf("latest (%d)", schedule_cost(x, latest_M)))
plot_moves(x, furthest_M, main = sprintf("furthest (%d)", schedule_cost(x, furthest_M)))

## Minimium weight Perfect Bipartite matching 
x <- seq(12)
y <- sample(x)
lis <- y[phtools:::longest_subseq(match(y, x))]

to_move <- setdiff(x, lis)
mp <- move_pair(symbol = to_move[1], source = x, target = y, lis = lis, rule = "all")

match(permute_move(x, 1, 3), x)

C <- matrix(Inf, nrow = length(x), ncol = length(x))
# diag(C) <- 0
for (sym in to_move){
	mp <- move_pair(symbol = sym, source = x, target = y, lis = lis, rule = "all")
	for (i in seq(nrow(mp))){
		C[mp[i,1],mp[i,2]] <- abs(diff(mp[i,]))
	}
}
RcppHungarian::HungarianSolver(C)



# apply(TM, 2, function(tm){
# 	MS <- move_sequence(symbols = tm, s = x, t = y, lis = lis, ordered = TRUE, rule = "all")
# 	length(unique(sapply(MS, function(m){ schedule_cost(x, m$moves) }))) == 1
# })


a <- 1:54
b <- sample(a)
lis <- b[phtools:::longest_subseq(match(b, a))]
P <- t(as.matrix(permutations::rperm(n = 1500, r = length(a) - length(lis))))

strat_costs <- lapply(c("earliest", "latest", "closest", "furthest"), function(strat){
	MS <- move_sequence(symbols = setdiff(a, lis), s = a, t = b, lis = lis, ordered = FALSE, rule = strat)
	sapply(MS, function(x){ schedule_cost(a, M = x$moves, cost = "kendall") })
})

strat_costs <- do.call(cbind, strat_costs)
all(strat_costs[,1] <= strat_costs[,2]) ## False
all(strat_costs[,2] <= strat_costs[,1]) ## True
all(strat_costs[,2] <= strat_costs[,3]) ## True
all(strat_costs[,2] <= strat_costs[,4]) ## True 

random_order_costs <- apply(P, 2, function(p){
	MS <- move_sequence(symbols = setdiff(a, lis)[p], s = a, t = b, lis = lis, ordered = TRUE, rule = "latest")
	schedule_cost(a, MS[[1]]$moves)
})


plot(table(random_order_costs))
schedule_cost(a, greedy_min_cross(a = a, b = b, lcs = lis, use_lcs = TRUE))
schedule_cost(a, greedy_min_cross(a = a, b = b, lcs = lis, use_lcs = TRUE, opt = "maximize"))

## Motivation: Order of transpositions can dramatically affect the # of column operations
R <- pbgrad::r_geometric_complex(n = 15, radius = 0.20, dim = 2, filtered = TRUE)
ds <- sapply(R$simplices, length)
D <- pbgrad::boundary_matrix(R)


sample_filtration <- function(S){
	stopifnot(is.character(S))
	ds <- sapply(str_to_simplex(S), length)
	P <- seq(length(ds))
	for (u in unique(ds)){
		P[ds == u] <- sample(which(ds == u))
	}
	S <- S[P]
	for (i in seq(length(S))){
		sigma <- as.vector(str_to_simplex(S[i]))
		d <- length(sigma)
		face_idx <- match(simplex_to_str(combn(sigma, d-1)), S)
		if (!is.na(face_idx) && any(face_idx > i)){
			j <- max(face_idx)
			tmp <- S[j]
			S[j] <- S[i]
			S[i] <- tmp
		}
	}
	return(S)
}

P <- colnames(D)
Q <- sample_filtration(P)
kendall_dist(seq_along(P), match(Q, P))
spearman_dist(seq_along(P), match(Q, P))

RV <- pbgrad::reduce(D)
f0 <- structure(seq_along(colnames(D)), names = colnames(D))
f1 <- structure(match(Q, P), names = Q)
res <- pbgrad::update_RV(RV = RV, f0 = f0, f1 = f1)

# simplices <- str_to_simplex(S)


n  <- 10
f0 <- seq(n)/n
f1 <- sample(f0)
xy0 <- unname(cbind(0, f0))
xy1 <- unname(cbind(1, f1))


plot.default(NULL, xlim = c(-0.25, 1.25), ylim = range(f0), asp = 1)
points(xy0)
points(xy1)
segments(y0 = xy0[1,2], y1 = xy1[1,2], x0 = 0, x1 = 1)
v <- (xy1[1,] - xy0[1,])
u <- v / norm(matrix(v, ncol = 2), "F")
up <- -c(u[2], -u[1])

interp_pt <- xy0[1,] + 0.5*v
points(matrix(interp_pt, ncol = 2), col = "red")
points(matrix(interp_pt + 0.15*up, ncol = 2), col = "purple")

## generates n points randomly about a segment (x0,y0) - (x1, y1) by sampling 
## points some percentage of the lines distance about its normal
points_around_line <- function(pt1, pt2, n, p){
	v <- (pt2 - pt1)
	u <- v / norm(matrix(v, ncol = 2), "F")
	up <- -c(u[2], -u[1])
	alpha <- runif(n)
	pts_along_line <- matrix(pt1 + rep(alpha, each = 2)*v, ncol = 2, byrow = TRUE)
	signs <- sample(x = c(TRUE, FALSE), size = n, replace = TRUE)
	line_dist <- norm(matrix(rbind(pt1, pt2), ncol = 2), "F")
	r <- line_dist*p
	out <- as.vector(t(pts_along_line) + rep(ifelse(signs, r*up, -r*up), each = 2))
	return(matrix(out, ncol = 2, byrow = TRUE))
	# out <- pts_along_line + ifelse(signs, r*up, -r*up)
	# return(matrix(out, ncol = 2, byrow = TRUE))
}

pt1 <- xy0[1,]
pt2 <- xy1[1,]
points(points_around_line(pt1, pt2, n = 10, p = 0.005), col = "blue")

s1 <- points_around_line(pt1, pt2, n = 10, p = 0.005)
s1 <- s1[order(s1[,1]),]

line1 <- sf::st_linestring(s1, dim = 2)
segments_lst <- lapply(seq(nrow(xy1)), function(i){
	pt1 <- xy0[i,]
	pt2 <- xy1[i,]
	s1 <- points_around_line(pt1, pt2, n = 10, p = 0.05)
	to_keep <- (s1[,1] > pt1[1]) & (s1[,1] < pt2[1])
	s1 <- s1[to_keep,]
	s1 <- s1[order(s1[,1]),]
	s1 <- do.call(rbind, list(pt1, s1, pt2))
	return(s1)
})

lines <- lapply(segments_lst, function(s){ sf::st_linestring(s, dim = 2) })
smoothed_lines <- lapply(lines, smoothr::smooth)

plot(sf::st_multilinestring(lines))
plot(sf::st_multilinestring(smoothed_lines))


sf::st_intersection(lines[[1]], lines[[5]])

segments2 <- lapply(seq(nrow(xy1)), function(i){
	pt1 <- xy0[i,]
	pt2 <- xy1[i,]
	s1 <- points_around_line(pt1, pt2, n = 10, p = 0.05)
	to_keep <- (s1[,1] > pt1[1]) & (s1[,1] < pt2[1])
	s1 <- s1[to_keep,]
	s1 <- s1[order(s1[,1]),]
	s1 <- do.call(rbind, list(pt1, s1, pt2))
	return(s1)
})

lines2 <- lapply(segments2, function(s){ sf::st_linestring(s, dim = 2) })

sl1 <- lapply(lines, smoothr::smooth)
sl2 <- lapply(lines2, smoothr::smooth)

l1 <- data.frame(geo = sf::st_sfc(sl1))
l2 <- data.frame(geo = sf::st_sfc(sl2))
interp_lines <- transformr::tween_sf(l1, l2, 'linear', 50)

animation::saveGIF({
	frames <- unique(interp_lines$.frame)
	for (frame in frames){
		plot(interp_lines$geometry[interp_lines$.frame == frame])
	}
}, movie.name = "homotopy1.gif", interval = 0.1)




*/




