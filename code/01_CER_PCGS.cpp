#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Utility: stable log-determinant via Cholesky (for PD matrices)
inline bool chol_logdet(const arma::mat& A, double& logdet_out) {
  arma::mat L;
  if (!arma::chol(L, A)) return false;
  // log|A| = 2 * sum(log(diag(L)))
  logdet_out = 2.0 * arma::sum(arma::log(L.diag()));
  return true;
}

// Dirichlet sampler
arma::vec rdirichlet(const arma::vec& alpha) {
  int K = alpha.n_elem;
  arma::vec g(K);
  for (int k = 0; k < K; ++k) g(k) = R::rgamma(alpha(k), 1.0);
  double sumg = arma::accu(g);
  return g / sumg;
}

// This function performs the PCGS for CER
// 
// X: n x p data matrix (n observations, p predictors)
// y: response vector (length n)
// K: number of clusters
// n_iter: total number of MCMC iterations
// burn_in: number of burn-in iterations
// tau2, rho0, p1: hyperparameters
// predict: if TRUE, record the predicted y in each iteration
// m_init: initial assignment matrix (p x K) with one-hot encoding
//
// [[Rcpp::export]]
Rcpp::List CER_PCGS_rcpp(const arma::mat& X, 
                         const arma::vec& y, 
                         int K, 
                         int n_iter, 
                         int burn_in, 
                         double tau2, 
                         double rho0, 
                         double p1, 
                         bool add_b_mu,
                         bool predict, 
                         const arma::mat& m_init) {
  
  int n = X.n_rows;
  int p = X.n_cols;
  const double eps = 1e-6;
  
  // Use the provided initial assignment matrix (p x K)
  arma::mat m = m_init;
  
  // Initialize cluster probabilities: first cluster gets p1; others share (1-p1)
  arma::vec p_cluster(K);
  p_cluster.fill((1.0 - p1) / std::max(1, K - 1));
  p_cluster(0) = p1;
  
  // ---- Build Z = X %*% m[ , 2:K ] and initialize b_mu0, sigma2 ----
  arma::vec b_mu0(K, arma::fill::zeros); // b_mu0[0]=0
  
  double sigma2 = 1.0; // placeholder; will be updated below
  
  if (K > 1) {
    arma::mat m_sub = m.cols(1, K - 1);      // p x (K-1)
    arma::mat Z = X * m_sub;                 // n x (K-1)
    
    // "Ridge" coefficient using Bayesian prior tau2 (no CV; consistent with prior structure)
    // beta_hat = (Z'Z + I/tau2)^(-1) Z' y
    arma::mat At = Z.t() * Z + (1.0 / tau2) * arma::eye(K - 1, K - 1);
    arma::vec Ct = Z.t() * y;
    arma::vec beta_hat = arma::solve(At, Ct, arma::solve_opts::fast);
    
    double b_mu_init = arma::mean(beta_hat);
    if (add_b_mu) {
      for (int k = 1; k < K; ++k) b_mu0(k) = b_mu_init;
    } // else remains zeros
    
    // Residuals and sigma2 init
    arma::vec yhat = Z * beta_hat;
    arma::vec resid = y - yhat;
    sigma2 = arma::as_scalar(arma::mean(arma::square(resid)));
  } else {
    // K == 1: only the null cluster b1=0; set sigma2 to y variance
    arma::vec resid = y; // y - 0
    sigma2 = arma::as_scalar(arma::mean(arma::square(resid)));
  }
  
  Rcpp::Rcout << "Init b_mu0[1]: " << b_mu0(1) << std::endl;
  Rcpp::Rcout << "Init sigma2: " << sigma2 << std::endl;
  
  // ---- Initialize mu_k, b_k (b1 = 0 fixed) ----
  arma::vec mu_k(K, arma::fill::zeros);
  arma::vec b_k(K, arma::fill::zeros);
  
  // Practical init: if cluster k non-empty, take a shrinkage estimate
  if (K > 1) {
    for (int k = 1; k < K; ++k) {
      if (arma::accu(m.col(k)) > 0) {
        // Z_k = sum of X columns assigned to cluster k
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
    b_k = mu_k; // start at mean
    b_k(0) = 0.0;
  }
  
  // ---- Storage ----
  arma::vec sigma2_samples(n_iter, arma::fill::zeros);
  arma::mat mu_samples(n_iter, K, arma::fill::zeros);
  arma::mat b_samples(n_iter, K, arma::fill::zeros);
  arma::cube m_samples;               // store only post burn-in (like your R version)
  if (n_iter > burn_in) m_samples.set_size(p, K, n_iter - burn_in);
  arma::mat y_samples;
  if (predict) { y_samples.set_size(n_iter, n); y_samples.fill(0.0); }
  
  // ===================== GIBBS =====================
  for (int iter = 0; iter < n_iter; ++iter) {
    
    // ---- (1) Update m (cluster label for each predictor) ----
    for (int j = 0; j < p; ++j) {
      arma::vec logp(K, arma::fill::zeros);
      
      for (int k = 0; k < K; ++k) {
        arma::mat m_temp = m;
        m_temp.row(j).zeros();
        m_temp(j, k) = 1.0;
        
        if (K > 1) {
          arma::mat m_sub = m_temp.cols(1, K - 1);         // p x (K-1)
          arma::mat Z = X * m_sub;                         // n x (K-1)
          
          arma::mat A = (1.0 / sigma2) * (Z.t() * Z)
            + (1.0 / tau2) * arma::eye(K - 1, K - 1)
            + eps * arma::eye(K - 1, K - 1);
            arma::vec C = (1.0 / sigma2) * (Z.t() * y)
              + (1.0 / tau2) * b_mu0.subvec(1, K - 1);
            
            double logdetA;
            if (!chol_logdet(A, logdetA)) {
              logp(k) = -std::numeric_limits<double>::infinity();
            } else {
              // 0.5 * C' A^{-1} C  -  0.5 * log|A|  + log p_cluster[k]
              arma::vec sol = arma::solve(A, C, arma::solve_opts::fast);
              double quad = 0.5 * arma::dot(C, sol);
              logp(k) = quad - 0.5 * logdetA + std::log(std::max(p_cluster(k), 1e-300));
            }
        } else {
          // K == 1
          logp(k) = std::log(std::max(p_cluster(k), 1e-300));
        }
      }
      
      // normalize (log-sum-exp)
      double maxlog = logp.max();
      arma::vec probs = arma::exp(logp - maxlog);
      probs /= arma::accu(probs);
      
      // sample new cluster
      double u = R::runif(0.0, 1.0);
      double acc = 0.0;
      int pick = 0;
      for (int k = 0; k < K; ++k) { acc += probs(k); if (u <= acc) { pick = k; break; } }
      
      m.row(j).zeros();
      m(j, pick) = 1.0;
    }
    
    // ---- (2) Jointly update mu_k and b_k (fix b1=0) ----
    arma::mat Zfull = X * m; // n x K
    arma::mat A = (1.0 / sigma2) * (Zfull.t() * Zfull)
      + (1.0 / tau2) * arma::eye(K, K);
    arma::vec C = (1.0 / sigma2) * (Zfull.t() * y)
      + (1.0 / tau2) * b_mu0;
    
    // PD safeguard
    int tries = 0;
    while (tries < 5 && arma::min(arma::eig_sym(A)) <= 0.0) {
      A += eps * arma::eye(K, K);
      ++tries;
    }
    
    if (K > 1) {
      arma::mat A_sub = A.submat(1, 1, K - 1, K - 1);
      arma::mat A_inv = arma::inv_sympd(A_sub);
      arma::vec mu_full = A_inv * C.subvec(1, K - 1);
      
      mu_k(0) = 0.0;
      for (int k = 1; k < K; ++k) {
        mu_k(k) = (arma::accu(m.col(k)) > 0.0) ? mu_full(k - 1) : 0.0;
      }
      
      // reject if all small (|b| <= 0.05); require ANY > 0.05
      arma::vec draw;
      do {
        draw = arma::mvnrnd(mu_full, A_inv);
      } while (arma::all(arma::abs(draw) <= 0.05));
      
      b_k.zeros();
      b_k(0) = 0.0;
      for (int k = 1; k < K; ++k) b_k(k) = draw(k - 1);
    } else {
      mu_k(0) = 0.0;
      b_k(0)  = 0.0;
    }
    
    // ---- (3) Update sigma2 (Inv-Gamma) ----
    arma::vec b_for_pred(p);
    for (int j = 0; j < p; ++j) {
      arma::uword cid = arma::index_max(m.row(j));
      b_for_pred(j) = b_k(cid);
    }
    arma::vec resid = y - X * b_for_pred;
    double shape = rho0 + static_cast<double>(n) / 2.0;
    double rate  = rho0 + 0.5 * arma::dot(resid, resid); // 'rate' in R::rgamma is 1/scale
    sigma2 = 1.0 / R::rgamma(shape, 1.0 / rate);
    
    // ---- (4) Update p_cluster (Dirichlet) ----
    arma::vec counts = arma::sum(m, 0).t();
    arma::vec alpha0(K);
    alpha0.fill((1.0 - p1) / std::max(1, K - 1));
    alpha0(0) = p1;
    p_cluster = rdirichlet(alpha0 + counts);
    
    // ---- Save ----
    sigma2_samples(iter) = sigma2;
    b_samples.row(iter)  = b_k.t();
    mu_samples.row(iter) = mu_k.t();
    if (predict) y_samples.row(iter) = (X * b_for_pred).t();
    if (iter >= burn_in && n_iter > burn_in) {
      m_samples.slice(iter - burn_in) = m;
    }
  }
  // =================== END GIBBS ===================
  
  // ---- Posterior summaries (burn-in removed) ----
  int keep0 = std::max(burn_in, 0);
  int keepN = std::max(0, n_iter - keep0);
  double posterior_sigma2 = arma::mean(sigma2_samples.subvec(keep0, n_iter - 1));
  arma::rowvec posterior_b = arma::mean(b_samples.rows(keep0, n_iter - 1), 0);
  arma::rowvec posterior_mu = arma::mean(mu_samples.rows(keep0, n_iter - 1), 0);
  
  arma::mat posterior_m;
  if (keepN > 0) posterior_m = arma::mean(m_samples, 2); // p x K
  
  arma::vec posterior_y;
  if (predict) {
    posterior_y = arma::mean(y_samples.rows(keep0, n_iter - 1), 0).t();
  }
  
  return Rcpp::List::create(
    Rcpp::Named("posterior_sigma2") = posterior_sigma2,
    Rcpp::Named("posterior_m")      = posterior_m,
    Rcpp::Named("posterior_mu")     = posterior_mu,
    Rcpp::Named("posterior_b")      = posterior_b,
    Rcpp::Named("b_samples")        = b_samples,
    Rcpp::Named("y_samples")        = predict ? y_samples : arma::mat(),
    Rcpp::Named("posterior_y")      = predict ? posterior_y : arma::vec(),
    Rcpp::Named("m_samples")        = m_samples
  );
}

