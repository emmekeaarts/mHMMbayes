#include <Rcpp.h>
using namespace Rcpp;

// This is a c++ version of the the forward part of the forward backward algorithm
// for a categorical emission distribution including timevarying covariate(s) to predict gamma
//' @keywords internal
// [[Rcpp::export(rng = false)]]
List cat_tv_mult_fw_cpp(NumericMatrix allprobs, NumericMatrix tgamma, int m, int n, NumericVector delta) {

  int i, t, k;
  NumericVector foo(m);
  NumericVector foo1(m);
  double sumfoo;
  NumericMatrix alpha_prob(m,n);
  NumericMatrix gamma(m,m);
  double lscale;
  NumericMatrix lalpha(m, n);

  foo = delta * allprobs(0,_);
  sumfoo = 0;
  for (i = 0; i < m; i++){
    sumfoo += foo[i];
  }
  for (i = 0; i < m; i++){
    alpha_prob(i,0) = foo[i]/sumfoo;
  }
  lscale = log(sumfoo);
  for(i = 0; i < m; i++){
    lalpha(i, 0) = log(alpha_prob(i, 0)) + lscale;
  }

  for (t = 1; t < n; t++){
    NumericVector gamma_vec = tgamma(t,_);
    gamma_vec.attr("dim") = Dimension(m, m);
    gamma = as<NumericMatrix>(gamma_vec);
    gamma = transpose(gamma);

    for (i = 0; i < m; i++){
      foo1[i] = 0;
      for (k = 0; k < m; k++){
        foo1[i] += alpha_prob(k, (t - 1)) * gamma(k,i);
      }
    }
    foo = foo1 * allprobs(t,_);
    sumfoo = 0;
    for (i = 0; i < m; i++){
      sumfoo += foo[i];
    }
    for (i = 0; i < m; i++){
      alpha_prob(i,t) = foo[i]/sumfoo;
    }
    lscale = lscale + log(sumfoo);
    for(i = 0; i < m; i++){
      lalpha(i, t) = log(alpha_prob(i, t)) + lscale;
    }
  }

  return List::create(alpha_prob, lalpha);
}

