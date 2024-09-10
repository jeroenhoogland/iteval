#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List mbcb_cpp(NumericVector p_0,
              NumericVector p_1,
              NumericVector tau_hat,
              NumericVector sample_weight) {
  // initialize components of standard deviation
  double sumr = 0;
  double sumr2 = 0;
  double sumw = 0;
  double sumw2 = 0;
  double sumrw = 0;
  // initialize numerator and denominator of c-index
  unsigned long long int n = tau_hat.size();
  double snum = 0;
  double sden = 0;
  // define patterns for estimating P(d_j < d_t)
  NumericVector pl1(5);
  NumericVector pl2(5);
  pl1[0] = 0; pl2[0] = 0;
  pl1[1] = 1; pl2[1] = 1;
  pl1[2] = 1; pl2[2] = 0;
  pl1[3] = 1; pl2[3] = 0;
  pl1[4] = 1; pl2[4] = 0;
  // definme patterns for estimating P(d_j > d_t)
  NumericVector pg1(5);
  NumericVector pg2(5);
  pg1[0] = 0; pg2[0] = 1;
  pg1[1] = 0; pg2[1] = 1;
  pg1[2] = 0; pg2[2] = 1;
  pg1[3] = 0; pg2[3] = 0;
  pg1[4] = 1; pg2[4] = 1;
  // loop over observations
  for (unsigned long long int j = 0; j < n; j++) {
    // initialize component of sumw and sumr for observation j
    double wj = 0;
    double rj = 0;
    // compute j's component in P(d_j < d_t) for each of the five patterns
    NumericVector plj(5);
    plj[0] = (1 - p_1[j]) * p_0[j];
    plj[1] = (1 - p_1[j]) * p_0[j];
    plj[2] = (1 - p_1[j]) * p_0[j];
    plj[3] = (1 - p_1[j]) * (1 - p_0[j]);
    plj[4] = p_1[j] * p_0[j];
    // compute j's component in P(d_j > d_t) for each of the five patterns
    NumericVector pgj(5);
    pgj[0] = (1 - p_1[j]) * (1 - p_0[j]);
    pgj[1] = p_1[j] * p_0[j];
    pgj[2] = p_1[j] * (1 - p_0[j]);
    pgj[3] = p_1[j] * (1 - p_0[j]);
    pgj[4] = p_1[j] * (1 - p_0[j]);
    // inner loop over observations
    for (unsigned long long int t = 0; t < n; t++) {
      // compare only pairs of different observations
      if (t == j) continue;
      int tauless = 0;
      int taueq = 0;
      int taugrt = 0;
      // for those indices t where expected benefit of index j is lower than for index t, include the probability
      // of observing that direction of differential benefit.
      if (tau_hat[j] < tau_hat[t]) tauless++;
      if (tau_hat[j] == tau_hat[t]) taueq++;
      if (tau_hat[j] > tau_hat[t]) taugrt++;
      // compute t's component in P(d_j < d_t) for each of the five patterns. Note the contribution of the comparison
      // between observation j and observation t to the c-for-benefit is weighted according to the provided
      // sample weights.
      for (int k = 0; k < 5; k++) {
        double pltk = (pl1[k] * p_1[t] + (1 - pl1[k]) * (1 - p_1[t])) *
          (pl2[k] * p_0[t] + (1 - pl2[k]) * (1 - p_0[t]));
        double pgtk = (pg1[k] * p_1[t] + (1 - pg1[k]) * (1 - p_1[t])) *
          (pg2[k] * p_0[t] + (1 - pg2[k]) * (1 - p_0[t]));
        snum += sample_weight[j] * sample_weight[t] * (tauless * pltk * plj[k] + 0.5 * taueq * pltk * plj[k]);
        sden += sample_weight[j] * sample_weight[t] * (pltk * plj[k]);
        wj += sample_weight[j] * sample_weight[t] *
          (tauless * pltk * plj[k] + taugrt  * pgtk * pgj[k] -
          taugrt  * pltk * plj[k] - tauless * pgtk * pgj[k]);
        rj += sample_weight[j] * sample_weight[t] * (pltk * plj[k] + pgtk * pgj[k]);
      }
    }
    // add j's component to standard error sums
    sumr += rj;
    sumr2 += rj * rj;
    sumw += wj;
    sumw2 += wj * wj;
    sumrw += rj * wj;
  }
  // calculate c-for-benefit
  double cindex = snum / sden;
  // calculate standard error
  double sd = 0;
  sd = sumr2 * sumw * sumw - 2 * sumr * sumw * sumrw + sumw2 * sumr * sumr;
  sd = sqrt(sd) / sumr / sumr;
  // return list with c-for-benefit and its standard error
  return List::create(
    _["C Index"] = cindex,
    _["S.D."] = sd
  );
}
