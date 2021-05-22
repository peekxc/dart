// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// reduce_pspbool
void reduce_pspbool(SEXP Ds, SEXP Vs);
RcppExport SEXP _dart_reduce_pspbool(SEXP DsSEXP, SEXP VsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Ds(DsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Vs(VsSEXP);
    reduce_pspbool(Ds, Vs);
    return R_NilValue;
END_RCPP
}
// reduce_local_pspbool
void reduce_local_pspbool(SEXP D1s, SEXP V1s, SEXP D2s, SEXP V2s, bool clearing);
RcppExport SEXP _dart_reduce_local_pspbool(SEXP D1sSEXP, SEXP V1sSEXP, SEXP D2sSEXP, SEXP V2sSEXP, SEXP clearingSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type D1s(D1sSEXP);
    Rcpp::traits::input_parameter< SEXP >::type V1s(V1sSEXP);
    Rcpp::traits::input_parameter< SEXP >::type D2s(D2sSEXP);
    Rcpp::traits::input_parameter< SEXP >::type V2s(V2sSEXP);
    Rcpp::traits::input_parameter< bool >::type clearing(clearingSEXP);
    reduce_local_pspbool(D1s, V1s, D2s, V2s, clearing);
    return R_NilValue;
END_RCPP
}
// simulate_vineyard_cpp
void simulate_vineyard_cpp(SEXP Rs, SEXP Vs, IntegerVector schedule, Nullable< Function > f);
RcppExport SEXP _dart_simulate_vineyard_cpp(SEXP RsSEXP, SEXP VsSEXP, SEXP scheduleSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Rs(RsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Vs(VsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type schedule(scheduleSEXP);
    Rcpp::traits::input_parameter< Nullable< Function > >::type f(fSEXP);
    simulate_vineyard_cpp(Rs, Vs, schedule, f);
    return R_NilValue;
END_RCPP
}
// reduce_arma
Rcpp::List reduce_arma(arma::sp_mat& D, arma::sp_mat& v);
RcppExport SEXP _dart_reduce_arma(SEXP DSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(reduce_arma(D, v));
    return rcpp_result_gen;
END_RCPP
}
// reduce_local_arma
Rcpp::List reduce_local_arma(arma::sp_mat& D1, arma::sp_mat& v1, arma::sp_mat& D2, arma::sp_mat& v2, bool clearing);
RcppExport SEXP _dart_reduce_local_arma(SEXP D1SEXP, SEXP v1SEXP, SEXP D2SEXP, SEXP v2SEXP, SEXP clearingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type D1(D1SEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type D2(D2SEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type v2(v2SEXP);
    Rcpp::traits::input_parameter< bool >::type clearing(clearingSEXP);
    rcpp_result_gen = Rcpp::wrap(reduce_local_arma(D1, v1, D2, v2, clearing));
    return rcpp_result_gen;
END_RCPP
}
// boundary_matrix_st
S4 boundary_matrix_st(SEXP stree, const size_t k);
RcppExport SEXP _dart_boundary_matrix_st(SEXP streeSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type stree(streeSEXP);
    Rcpp::traits::input_parameter< const size_t >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(boundary_matrix_st(stree, k));
    return rcpp_result_gen;
END_RCPP
}
// all_lcs
ListOf< IntegerVector > all_lcs(IntegerVector a, IntegerVector b);
RcppExport SEXP _dart_all_lcs(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(all_lcs(a, b));
    return rcpp_result_gen;
END_RCPP
}
// LIS
List LIS(const IntegerVector& a);
RcppExport SEXP _dart_LIS(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(LIS(a));
    return rcpp_result_gen;
END_RCPP
}
// countNumLIS
int countNumLIS(std::vector<int> x);
RcppExport SEXP _dart_countNumLIS(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(countNumLIS(x));
    return rcpp_result_gen;
END_RCPP
}
// inversion_count
int inversion_count(IntegerVector iv);
RcppExport SEXP _dart_inversion_count(SEXP ivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type iv(ivSEXP);
    rcpp_result_gen = Rcpp::wrap(inversion_count(iv));
    return rcpp_result_gen;
END_RCPP
}
// perm_dist_mat
NumericVector perm_dist_mat(IntegerMatrix P, const bool kendall, const bool normalize);
RcppExport SEXP _dart_perm_dist_mat(SEXP PSEXP, SEXP kendallSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< const bool >::type kendall(kendallSEXP);
    Rcpp::traits::input_parameter< const bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(perm_dist_mat(P, kendall, normalize));
    return rcpp_result_gen;
END_RCPP
}
// fast_choose
size_t fast_choose(const size_t n, const size_t k);
RcppExport SEXP _dart_fast_choose(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const size_t >::type n(nSEXP);
    Rcpp::traits::input_parameter< const size_t >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_choose(n, k));
    return rcpp_result_gen;
END_RCPP
}
// spearman_perm
size_t spearman_perm(const IntegerVector& x, const IntegerVector& y);
RcppExport SEXP _dart_spearman_perm(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(spearman_perm(x, y));
    return rcpp_result_gen;
END_RCPP
}
// reduce_local_dense
List reduce_local_dense(const NumericMatrix& D1, const NumericMatrix& v1, const NumericMatrix& D2, const NumericMatrix& v2);
RcppExport SEXP _dart_reduce_local_dense(SEXP D1SEXP, SEXP v1SEXP, SEXP D2SEXP, SEXP v2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type D1(D1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type D2(D2SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type v2(v2SEXP);
    rcpp_result_gen = Rcpp::wrap(reduce_local_dense(D1, v1, D2, v2));
    return rcpp_result_gen;
END_RCPP
}
// reduce_dense
List reduce_dense(const NumericMatrix& D, const NumericMatrix& v);
RcppExport SEXP _dart_reduce_dense(SEXP DSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(reduce_dense(D, v));
    return rcpp_result_gen;
END_RCPP
}
// interval_cost_rcpp
int interval_cost_rcpp(IntegerVector s, IntegerVector O);
RcppExport SEXP _dart_interval_cost_rcpp(SEXP sSEXP, SEXP OSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type O(OSEXP);
    rcpp_result_gen = Rcpp::wrap(interval_cost_rcpp(s, O));
    return rcpp_result_gen;
END_RCPP
}
// inverse_permutation
IntegerVector inverse_permutation(IntegerVector p);
RcppExport SEXP _dart_inverse_permutation(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(inverse_permutation(p));
    return rcpp_result_gen;
END_RCPP
}
// inverse_permutation2
IntegerVector inverse_permutation2(IntegerVector p);
RcppExport SEXP _dart_inverse_permutation2(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(inverse_permutation2(p));
    return rcpp_result_gen;
END_RCPP
}
// test_unrank
IntegerVector test_unrank(const size_t r, const size_t n, const size_t k);
RcppExport SEXP _dart_test_unrank(SEXP rSEXP, SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const size_t >::type r(rSEXP);
    Rcpp::traits::input_parameter< const size_t >::type n(nSEXP);
    Rcpp::traits::input_parameter< const size_t >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(test_unrank(r, n, k));
    return rcpp_result_gen;
END_RCPP
}
// order_simplices
IntegerVector order_simplices(List simplices, NumericVector weights);
RcppExport SEXP _dart_order_simplices(SEXP simplicesSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type simplices(simplicesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(order_simplices(simplices, weights));
    return rcpp_result_gen;
END_RCPP
}
// unique_numeric
NumericVector unique_numeric(std::vector< double > x, const double eps);
RcppExport SEXP _dart_unique_numeric(SEXP xSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector< double > >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(unique_numeric(x, eps));
    return rcpp_result_gen;
END_RCPP
}
// longest_inc_subseq
Rcpp::IntegerVector longest_inc_subseq(SEXP X);
RcppExport SEXP _dart_longest_inc_subseq(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(longest_inc_subseq(X));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_PspBoolMatrix();
RcppExport SEXP _rcpp_module_boot_implicit_filtration_module();

static const R_CallMethodDef CallEntries[] = {
    {"_dart_reduce_pspbool", (DL_FUNC) &_dart_reduce_pspbool, 2},
    {"_dart_reduce_local_pspbool", (DL_FUNC) &_dart_reduce_local_pspbool, 5},
    {"_dart_simulate_vineyard_cpp", (DL_FUNC) &_dart_simulate_vineyard_cpp, 4},
    {"_dart_reduce_arma", (DL_FUNC) &_dart_reduce_arma, 2},
    {"_dart_reduce_local_arma", (DL_FUNC) &_dart_reduce_local_arma, 5},
    {"_dart_boundary_matrix_st", (DL_FUNC) &_dart_boundary_matrix_st, 2},
    {"_dart_all_lcs", (DL_FUNC) &_dart_all_lcs, 2},
    {"_dart_LIS", (DL_FUNC) &_dart_LIS, 1},
    {"_dart_countNumLIS", (DL_FUNC) &_dart_countNumLIS, 1},
    {"_dart_inversion_count", (DL_FUNC) &_dart_inversion_count, 1},
    {"_dart_perm_dist_mat", (DL_FUNC) &_dart_perm_dist_mat, 3},
    {"_dart_fast_choose", (DL_FUNC) &_dart_fast_choose, 2},
    {"_dart_spearman_perm", (DL_FUNC) &_dart_spearman_perm, 2},
    {"_dart_reduce_local_dense", (DL_FUNC) &_dart_reduce_local_dense, 4},
    {"_dart_reduce_dense", (DL_FUNC) &_dart_reduce_dense, 2},
    {"_dart_interval_cost_rcpp", (DL_FUNC) &_dart_interval_cost_rcpp, 2},
    {"_dart_inverse_permutation", (DL_FUNC) &_dart_inverse_permutation, 1},
    {"_dart_inverse_permutation2", (DL_FUNC) &_dart_inverse_permutation2, 1},
    {"_dart_test_unrank", (DL_FUNC) &_dart_test_unrank, 3},
    {"_dart_order_simplices", (DL_FUNC) &_dart_order_simplices, 2},
    {"_dart_unique_numeric", (DL_FUNC) &_dart_unique_numeric, 2},
    {"_dart_longest_inc_subseq", (DL_FUNC) &_dart_longest_inc_subseq, 1},
    {"_rcpp_module_boot_PspBoolMatrix", (DL_FUNC) &_rcpp_module_boot_PspBoolMatrix, 0},
    {"_rcpp_module_boot_implicit_filtration_module", (DL_FUNC) &_rcpp_module_boot_implicit_filtration_module, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_dart(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}