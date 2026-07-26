#ifndef SBLR_BLR_MT_COVARIANCE_RNG_H
#define SBLR_BLR_MT_COVARIANCE_RNG_H

#include <armadillo>

#include <cmath>
#include <random>
#include <stdexcept>

namespace sblr {
namespace mt {

inline arma::mat draw_wishart(
 double df,
 const arma::mat& scale,
 std::mt19937& rng
) {
 const arma::uword p=scale.n_rows;
 if (p==0 || scale.n_cols!=p)
  throw std::invalid_argument("draw_wishart: scale must be nonempty and square.");
 if (!std::isfinite(df) || df<=static_cast<double>(p)-1.0)
  throw std::invalid_argument("draw_wishart: df must be > dimension - 1.");
 if (!scale.is_finite())
  throw std::invalid_argument("draw_wishart: scale must be finite.");

 arma::mat lower;
 if (!arma::chol(lower, 0.5*(scale+scale.t()), "lower"))
  throw std::invalid_argument("draw_wishart: scale must be symmetric positive definite.");

 arma::mat bartlett(p, p, arma::fill::zeros);
 for (arma::uword i=0; i<p; ++i) {
  std::chi_squared_distribution<double> chi_square(df-static_cast<double>(i));
  bartlett(i,i)=std::sqrt(chi_square(rng));
  std::normal_distribution<double> normal(0.0, 1.0);
  for (arma::uword j=0; j<i; ++j) bartlett(i,j)=normal(rng);
 }

 arma::mat value=lower*bartlett*bartlett.t()*lower.t();
 return 0.5*(value+value.t());
}

inline arma::mat draw_inverse_wishart(
 double df,
 const arma::mat& scale,
 std::mt19937& rng
) {
 if (scale.n_rows==0 || scale.n_cols!=scale.n_rows)
  throw std::invalid_argument(
   "draw_inverse_wishart: scale must be nonempty and square.");
 if (!scale.is_finite())
  throw std::invalid_argument("draw_inverse_wishart: scale must be finite.");
 arma::mat symmetric_scale=0.5*(scale+scale.t());
 arma::mat scale_inverse;
 if (!arma::inv_sympd(scale_inverse, symmetric_scale))
  throw std::invalid_argument(
   "draw_inverse_wishart: scale must be symmetric positive definite.");
 arma::mat wishart=draw_wishart(df, scale_inverse, rng);
 arma::mat value;
 if (!arma::inv_sympd(value, wishart))
  throw std::runtime_error("draw_inverse_wishart: Wishart draw was not SPD.");
 return 0.5*(value+value.t());
}

}  // namespace mt
}  // namespace sblr

#endif
