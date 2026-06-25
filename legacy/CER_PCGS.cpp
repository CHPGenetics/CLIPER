#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

inline bool chol_logdet(const arma::mat& A, double& logdet_out) {
  arma::mat L;
  if (!arma::chol(L, A)) return false;
  logdet_out = 2.0 * arma::sum(arma::log(L.diag()));
  return true;
}

arma::vec rdirichlet(const arma::vec& alpha) {
  int K = alpha.n_elem;
  arma::vec g(K);
  for (int k = 0; k < K; ++k) g(k) = R::rgamma(alpha(k), 1.0);
  double sumg = arma::accu(g);
  return g / sumg;
}

// [[Rcpp::export]]
Rcpp::List CER_PCGS_rcpp(const arma::mat& X, 
                         const arma::vec& y, 
                         int K, 
                         int n_iter, 
                         int burn_in, 
                         int alpha_conc,
                         double tau2, 
                         double rho0, 
                         double p1, 
                         bool add_b_mu,
                         bool predict, 
                         const arma::mat& m_init) {
  
  int n = X.n_rows;
  int p = X.n_cols;
  const double eps = 1e-6;
  
  if (n_iter <= 0) Rcpp::stop("n_iter must be positive.");
  if (burn_in < 0 || burn_in >= n_iter) Rcpp::stop("Require 0 <= burn_in < n_iter.");
  if (K < 2) Rcpp::stop("K must be at least 2: one null cluster plus at least one effect cluster.");
  if (m_init.n_rows != p || m_init.n_cols != K) Rcpp::stop("m_init must be p x K.");
  if (tau2 <= 0 || rho0 <= 0) Rcpp::stop("tau2 and rho0 must be positive.");
  if (p1 <= 0 || p1 >= 1) Rcpp::stop("p1 must be in (0,1).");
  
  arma::mat m = m_init;
  arma::vec p_cluster(K);
  p_cluster.fill((1.0 - p1) / std::max(1, K - 1));
  p_cluster(0) = p1;
  
  arma::vec b_mu0(K, arma::fill::zeros);
  
  double sigma2 = 1.0;
  
  if (K > 1) {
    arma::mat m_sub = m.cols(1, K - 1);
    arma::mat Z = X * m_sub;
    
    arma::mat At = Z.t() * Z + (1.0 / tau2) * arma::eye(K - 1, K - 1);
    arma::vec Ct = Z.t() * y;
    arma::vec beta_hat = arma::solve(At, Ct, arma::solve_opts::fast);
    
    double b_mu_init = arma::mean(beta_hat);
    if (add_b_mu) {
      for (int k = 1; k < K; ++k) b_mu0(k) = b_mu_init;
    }
    
    arma::vec yhat = Z * beta_hat;
    arma::vec resid = y - yhat;
    sigma2 = arma::as_scalar(arma::mean(arma::square(resid)));
  } else {
    arma::vec resid = y;
    sigma2 = arma::as_scalar(arma::mean(arma::square(resid)));
  }
  sigma2 = std::max(sigma2, 1e-12);
  
  arma::vec mu_k(K, arma::fill::zeros);
  arma::vec b_k(K, arma::fill::zeros);
  
  if (K > 1) {
    for (int k = 1; k < K; ++k) {
      if (arma::accu(m.col(k)) > 0) {
        arma::vec Zk(n, arma::fill::zeros);
        for (int j = 0; j < p; ++j) {
          if (m(j, k) == 1.0) Zk += X.col(j);
        }
        double num = arma::dot(Zk, y) / sigma2 + b_mu0(k) / tau2;
        double den = arma::dot(Zk, Zk) / sigma2 + (1.0 / tau2);
        mu_k(k) = num / den;
      } else {
        mu_k(k) = 0.0;
      }
    }
    b_k = mu_k;
    b_k(0) = 0.0;
  }
  
  arma::vec sigma2_samples(n_iter, arma::fill::zeros);
  arma::mat mu_samples(n_iter, K, arma::fill::zeros);
  arma::mat b_samples(n_iter, K, arma::fill::zeros);
  arma::mat y_samples;
  if (predict) { y_samples.set_size(n_iter, n); y_samples.fill(0.0); }
  
  arma::mat posterior_m_sum(p, K, arma::fill::zeros);
  arma::vec b_pred_sum(p, arma::fill::zeros);
  arma::vec b_pred_sq_sum(p, arma::fill::zeros);
  
  for (int iter = 0; iter < n_iter; ++iter) {
    
    for (int j = 0; j < p; ++j) {
      arma::vec logp(K, arma::fill::zeros);
      
      for (int k = 0; k < K; ++k) {
        arma::mat m_temp = m;
        m_temp.row(j).zeros();
        m_temp(j, k) = 1.0;
        
        if (K > 1) {
          arma::mat m_sub = m_temp.cols(1, K - 1);
          arma::mat Z = X * m_sub;
          
          arma::mat A = (1.0 / sigma2) * (Z.t() * Z)
            + (1.0 / tau2) * arma::eye(K - 1, K - 1)
            + eps * arma::eye(K - 1, K - 1);
            arma::vec C = (1.0 / sigma2) * (Z.t() * y)
              + (1.0 / tau2) * b_mu0.subvec(1, K - 1);
            
            double logdetA;
            if (!chol_logdet(A, logdetA)) {
              logp(k) = -std::numeric_limits<double>::infinity();
            } else {
              arma::vec sol = arma::solve(A, C, arma::solve_opts::fast);
              double quad = 0.5 * arma::dot(C, sol);
              logp(k) = quad - 0.5 * logdetA + std::log(std::max(p_cluster(k), 1e-300));
            }
        } else {
          logp(k) = std::log(std::max(p_cluster(k), 1e-300));
        }
      }
      
      double maxlog = logp.max();
      arma::vec probs = arma::exp(logp - maxlog);
      probs /= arma::accu(probs);
      
      double u = R::runif(0.0, 1.0);
      double acc = 0.0;
      int pick = 0;
      for (int k = 0; k < K; ++k) { acc += probs(k); if (u <= acc) { pick = k; break; } }
      
      m.row(j).zeros();
      m(j, pick) = 1.0;
    }
    
    arma::mat Zfull = X * m;
    arma::mat A = (1.0 / sigma2) * (Zfull.t() * Zfull) + (1.0 / tau2) * arma::eye(K, K);
    arma::vec C = (1.0 / sigma2) * (Zfull.t() * y)      + (1.0 / tau2) * b_mu0;
    
    if (K > 1) {
      arma::mat A_sub = A.submat(1, 1, K - 1, K - 1);
      
      arma::mat L;
      double jitter = 1e-10;
      int tries = 0;
      while (tries < 10 && !arma::chol(L, A_sub)) {
        A_sub += jitter * arma::eye(K - 1, K - 1);
        jitter *= 10.0;
        ++tries;
      }
      if (!arma::chol(L, A_sub)) {
        Rcpp::stop("A_sub not PD even after jitter.");
      }
      
      arma::mat A_inv = arma::inv_sympd(A_sub);
      arma::vec mu_full = A_inv * C.subvec(1, K - 1);
      
      mu_k.zeros();
      mu_k(0) = 0.0;
      for (int kk = 1; kk < K; ++kk) {
        mu_k(kk) = (arma::accu(m.col(kk)) > 0.0) ? mu_full(kk - 1) : 0.0;
      }
      
      arma::vec draw = arma::mvnrnd(mu_full, A_inv);
      
      b_k.zeros();
      b_k(0) = 0.0;
      for (int kk = 1; kk < K; ++kk) b_k(kk) = draw(kk - 1);
      
    } else {
      mu_k.zeros(); mu_k(0) = 0.0;
      b_k.zeros();  b_k(0)  = 0.0;
    }
    
    arma::vec b_for_pred(p);
    for (int j = 0; j < p; ++j) {
      arma::uword cid = arma::index_max(m.row(j));
      b_for_pred(j) = b_k(cid);
    }
    arma::vec resid = y - X * b_for_pred;
    double shape = rho0 + static_cast<double>(n) / 2.0;
    double rate  = rho0 + 0.5 * arma::dot(resid, resid);
    sigma2 = 1.0 / R::rgamma(shape, 1.0 / rate);
    
    arma::vec counts = arma::sum(m, 0).t();
    if (p1 <= 0.0 || p1 >= 1.0) {
      Rcpp::stop("p1 must be in (0, 1).");
    }
    if (alpha_conc <= 0.0) {
      Rcpp::stop("alpha_conc must be positive.");
    }
    
    arma::vec pi0(K);
    pi0.fill((1.0 - p1) / std::max(1, K - 1));
    pi0(0) = p1;
    
    arma::vec alpha0 = alpha_conc * pi0;
    arma::vec alpha = alpha0 + counts;
    alpha = arma::clamp(alpha, 1e-8, arma::datum::inf);
    p_cluster = rdirichlet(alpha);
    
    sigma2_samples(iter) = sigma2;
    b_samples.row(iter)  = b_k.t();
    mu_samples.row(iter) = mu_k.t();
    if (predict) y_samples.row(iter) = (X * b_for_pred).t();
    
    if (iter >= burn_in) {
      posterior_m_sum += m;
      b_pred_sum      += b_for_pred;
      b_pred_sq_sum   += arma::square(b_for_pred);
    }
  }
  
  int keep0 = std::max(burn_in, 0);
  int keepN = std::max(0, n_iter - keep0);
  
  double posterior_sigma2 = arma::mean(sigma2_samples.subvec(keep0, n_iter - 1));
  arma::rowvec posterior_b = arma::mean(b_samples.rows(keep0, n_iter - 1), 0);
  arma::rowvec posterior_mu = arma::mean(mu_samples.rows(keep0, n_iter - 1), 0);
  
  arma::mat posterior_m(p, K, arma::fill::zeros);
  arma::vec posterior_b_pred(p, arma::fill::zeros);
  arma::vec posterior_b_pred_sd(p, arma::fill::zeros);
  arma::vec posterior_b_pred_q025(p, arma::fill::zeros);
  arma::vec posterior_b_pred_q975(p, arma::fill::zeros);
  
  if (keepN > 0) {
    posterior_m      = posterior_m_sum / double(keepN);
    posterior_b_pred = b_pred_sum      / double(keepN);
    
    if (keepN > 1) {
      arma::vec var_b_pred =
        (b_pred_sq_sum - double(keepN) * arma::square(posterior_b_pred)) /
          double(keepN - 1);
      
      var_b_pred = arma::clamp(var_b_pred, 0.0, arma::datum::inf);
      posterior_b_pred_sd = arma::sqrt(var_b_pred);
      
      posterior_b_pred_q025 = posterior_b_pred - 1.96 * posterior_b_pred_sd;
      posterior_b_pred_q975 = posterior_b_pred + 1.96 * posterior_b_pred_sd;
    }
  }
  
  arma::vec posterior_y;
  if (predict) {
    posterior_y = arma::mean(y_samples.rows(keep0, n_iter - 1), 0).t();
  }
  
  return Rcpp::List::create(
    Rcpp::Named("posterior_sigma2")     = posterior_sigma2,
    Rcpp::Named("posterior_m")          = posterior_m,
    Rcpp::Named("posterior_mu")         = posterior_mu,
    Rcpp::Named("posterior_b")          = posterior_b,
    Rcpp::Named("posterior_b_pred")     = posterior_b_pred,
    Rcpp::Named("posterior_b_pred_sd")  = posterior_b_pred_sd,
    Rcpp::Named("posterior_b_pred_q025")= posterior_b_pred_q025,
    Rcpp::Named("posterior_b_pred_q975")= posterior_b_pred_q975,
    Rcpp::Named("b_samples")            = b_samples,
    Rcpp::Named("posterior_y")          = predict ? posterior_y : arma::vec()
  );
}
