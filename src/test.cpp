// Penalized Reduced Rank Multi-Outcome Cox Model
// alpha: pxR matrix
// Gamma: KxR matrix
// B: pxK matrix

#include <Rcpp.h>
// #include <omp.h> // OpenMP for parallel computing
#include <vector>
#include <iostream>
#include <sys/time.h>
#include <cmath>
#include "mrcox_types.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;


// Reverse Cumulative Sum
// In cox model, we compute the risk set at each time point:
// risk_denom(i,k) = Σ_{j: Y_{jk} ≥ Y_{ik}} exp(η_{jk})
void rev_cumsum_assign(const MatrixXd &src, MatrixXd &dest)
{
  const int K = src.cols();
  const int N = src.rows();
  
  for(int k = 0; k < K; ++k){
    double current = 0;
    for (int i = 0; i < N; ++i){
      current += src(N-1-i, k);
      dest(N-1-i, k) = current;
    }
  }
}


// 1. compute cox partial likelihood for reduced rank structure
// 2. handle tied event times using Efron's method
// 3. compute gradients using chain rule
// linear predictor: η_{ik} = Σ_r Z_{ir} * γ_{kr}
// Partial log-likelihood for outcome k: ℓ_k = Σ_{i∈D_k} [η_{ik} - log(Σ_{j∈R_i(k)} exp(η_{jk}))]
class MCox_ReducedRank
{
  // constant
  const int N;  // number of observations
  const int K;  // number of outcomes
  const int p;  // number of predictors
  const int R;  // number of rank
  
  // data pointers (point to R-allocated memory, no copying)
  // use Eigen::Map to provide matrix interface to raw arrays
  MapMatd X;        // N×p
  MapMatd status;   // N×K
  MapMati rankmin;  // N×K: minimum rank for tied events
  MapMati rankmax;  // N×K: maximum rank for tied events
  
  // order[k] is a permutation matrix that sorts outcome k by event time
  std::vector<PermMat> orders;
  
  // intermediate computation storage
  MatrixXd eta;          // N×K: linear predictors
  MatrixXd exp_eta;      // N×K: exp(eta)
  MatrixXd risk_denom;   // N×K: risk set denominators Σ_{j∈R_i} exp(η_j)
  MatrixXd outer_accumu; // N×K: cumulative hazards
  MatrixXd residual;     // N×K: ∂ℓ/∂η, used for gradients
  MatrixXd Z;            // N×R: Z = X %*% αlpha
  
  // outer_accumu(i,k) = Σ_{j: Y_{jk}≥Y_{ik}} δ_{jk} / risk_denom(j,k)
  // residual(i,k) = ∂(-ℓ)/∂η_{ik} 
  //                 = Σ_{j: Y_{jk}≥Y_{ik}} [δ_{jk} * exp(η_{ik})] / [Σ_{m∈R_j(k)} exp(η_{mk})] - δ_{ik}
  // residual(i,k) = outer_accumu(i,k) * exp(η_{ik}) - δ_{ik}
  // get_val: if true, calculate and return objective value
  double get_residual_rr(const MapMatd &alpha, const MapMatd &Gamma, bool get_val=false){
    
    // (1) compute linear predictors
    Z.noalias() = X * alpha;               // Z = Xα (N×R)
    eta.noalias() = Z * Gamma.transpose(); // η = ZΓ^T (N×K)
    
    // (2) sort by event time
    //     η_sorted = P_k * η
    //     where P_k is the permutation matrix for outcome k
    // #pragma omp parallel for
    for(int k = 0; k < K; ++k){
      eta.col(k) = orders[k] * eta.col(k);
    }
    
    // (3) exponentials (for risk set denominators)
    exp_eta.noalias() = eta.array().exp().matrix();
    
    // (4) compute risk set denominators
    // risk_denom(i,k) = Σ_{j≥i} exp(η_{jk})
    // #pragma omp parallel for 
    for(int k = 0; k < K; ++k){
      double current = 0;
      for (int i = 0; i < N; ++i){
        current += exp_eta(N-1-i, k);
        risk_denom(N-1-i, k) = current;
      }
      
      // adjust for tied events
      //   If rows 5, 6, 7 have the same event time:
      //   - rankmin(5,k) = rankmin(6,k) = rankmin(7,k) = 5
      //   - All three use risk_denom(5,k) (the risk set at earliest tied time)
      for(int i = 0; i < N; ++i){
        risk_denom(i, k) = risk_denom(rankmin(i, k), k);
      }
    }
    
    // (5) compute outer_accumu
    // outer_accumu(i,k) = δ_{ik} / risk_denom(i,k)
    // δ_{ik} (1 if event, 0 if censored)
    outer_accumu.noalias() = (status.array() / risk_denom.array()).matrix();
    
    // (6) cumulative sum for outer_accumu
    // outer_accumu(i,k) = Σ_{j=i}^{N-1} [δ_{jk} / risk_denom(j,k)]
    // #pragma omp parallel for
    for(int k = 0; k < K; ++k){
      double current = 0;
      for(int i = 0; i < N; ++i){
        current += outer_accumu(i, k);
        outer_accumu(i, k) = current;
      }
      
      // adjust for tied events
      //   If rows 5, 6, 7 have the same event time:
      //   - rankmax(5,k) = rankmax(6,k) = rankmax(7,k) = 7
      //   - All three use outer_accumu(7,k) (cumsum up to last tied individual)
      for(int i = 0; i < N; ++i){
        outer_accumu(i, k) = outer_accumu(rankmax(i, k), k);
      }
    }
    
    // (7) computer residuals
    // residual(i,k) = ∂(-ℓ)/∂η_{ik} = outer_accumu(i,k) * exp(η_{ik}) - δ_{ik}
    residual.noalias() = (outer_accumu.array() * exp_eta.array() - status.array()).matrix();
    
    // (8) restore oringinal order (for gradient computation)
    // #pragma omp parallel for
    for(int k = 0; k < K; ++k){
      residual.col(k) = orders[k].transpose() * residual.col(k);
    }
    
    // (9) compute objective value
    // -ℓ = Σ_{i,k} δ_{ik} [log(risk_denom(i,k)) - η_{ik}]
    double cox_val = 0;
    if(get_val){
      cox_val = ((risk_denom.array().log() - eta.array()) * status.array()).sum();
    }
    
    return cox_val;
  }
  
public:
  // initialize the Mcox_ReducedRank object
  MCox_ReducedRank(int N,
                   int K,
                   int p,
                   int R,
                   const double *X, // Pointer to first element of X
                   const double *status,
                   const int *rankmin,
                   const int *rankmax,
                   const Rcpp::List order_list) : 
  N(N), //initialize const int N
  K(K),
  p(p),
  R(R),                    
  X(X, N, p), // Create Map: X points to external data
  status(status, N, K),
  rankmin(rankmin, N, K),
  rankmax(rankmax, N, K),
  eta(N, K), // Allocate N×K matrix
  exp_eta(N, K),
  risk_denom(N, K),
  outer_accumu(N, K),
  residual(N, K),
  Z(N, R)                  
  {
    // Rcpp::as<VectorXi>: convert R integer vector to Eigen integer vector
    // order_list[k] contains an integer vector representing the permutation
    for (int k = 0; k < K; ++k){
      orders.emplace_back(Rcpp::as<VectorXi>(order_list[k]));
    }
  }
  
  // compute gradients
  // ∇_α(-ℓ) = (X^T * residual) * Gamma  (p×K) × (K×R) = (p×R)
  // ∇_Γ(-ℓ) = (X^T * residual)^T * alpha = residual^T * X * alpha
  //         = residual^T * Z  (K×N) × (N×R) = (K×R)
  // alpha_ptr: Pointer to alpha matrix data (p×R)
  // Gamma_ptr: Pointer to Gamma matrix data (K×R)
  // grad_alpha: Output matrix for ∇_α(-ℓ) (p×R)
  // grad_Gamma: Output matrix for ∇_Γ(-ℓ) (K×R)
  // get_val: If true, also return objective function value
  double get_gradients(const double *alpha_ptr, 
                       const double *Gamma_ptr,
                       MatrixXd &grad_alpha, 
                       MatrixXd &grad_Gamma,
                       bool get_val=false){
    
    // create Map objects (view raw data as Eigen matrices)
    // no copying occurs here - just wrapping pointers
    MapMatd alpha(alpha_ptr, p, R);
    MapMatd Gamma(Gamma_ptr, K, R);
    
    // compute residuals 
    double cox_val = get_residual_rr(alpha, Gamma, get_val);
    
    grad_alpha.noalias() = X.transpose() * residual * Gamma;
    grad_Gamma.noalias() = residual.transpose() * Z;
    
    return cox_val;
  }
  
  // compute objective value
  double get_value_only(const double *alpha_ptr, const double *Gamma_ptr){
    MapMatd alpha(alpha_ptr, p, R);
    MapMatd Gamma(Gamma_ptr, K, R);
    
    // compute η
    Z.noalias() = X * alpha;
    eta.noalias() = Z * Gamma.transpose();
    
    // sort by event time
    // #pragma omp parallel for
    for(int k = 0; k < K; ++k){
      eta.col(k) = orders[k] * eta.col(k);
    }
    
    exp_eta.noalias() = eta.array().exp().matrix();
    
    // compute risk_denom
    // #pragma omp parallel for  
    for(int k = 0; k < K; ++k){
      double current = 0;
      for (int i = 0; i < N; ++i){
        current += exp_eta(N-1-i, k);
        risk_denom(N-1-i, k) = current;
      }
      
      for(int i = 0; i < N; ++i){
        risk_denom(i, k) = risk_denom(rankmin(i, k), k);
      }
    }
    
    // compute objective value
    double cox_val = ((risk_denom.array().log() - eta.array()) * status.array()).sum();
    return cox_val;
  }
  
  // get residual matrix (for R to access)
  MatrixXd get_residual_matrix(){
    return residual;
  }
};


// group lasso soft-thresholding
// Pen(M) = λ * Σ_i ||M_{i·}||_2
// prox_{ηλ}(v) = argmin_u { (1/2)||u - v||_F^2 + ηλ * Σ_i ||u_{i·}||_2 }
// u_{i·} = { (1 - ηλ/||v_{i·}||_2) * v_{i·}  if ||v_{i·}||_2 > ηλ
//          { 0   otherwise
void prox_group_lasso(MatrixXd &M_out, 
                      const MatrixXd &M_in,
                      double step_size,
                      double lambda,
                      const VectorXd &penalty_factor)
{
  int nrows = M_in.rows();
  int ncols = M_in.cols();
  
  // process each row independently
  // Could be parallelized: #pragma omp parallel for
  for(int i = 0; i < nrows; ++i) {
    // compute row norm
    double row_norm = M_in.row(i).norm();
    
    // compute threshold (threshold = η * λ * penalty_factor_i)
    double threshold = step_size * lambda * penalty_factor[i];
    
    // Group soft-thresholding
    //   If ||v|| = c*threshold (c > 1), then:
    //   ||output|| = ||v|| * (1 - 1/c) = ||v|| * (c-1)/c
    if(row_norm > threshold) {
      double shrinkage = 1.0 - threshold / row_norm;
      M_out.row(i) = shrinkage * M_in.row(i);
    } else {
      M_out.row(i).setZero();
    }
  }
}


// proximal gradient descent with Nesterov
// minimize: L(α,Γ) = f(α,Γ) + g(α) + h(Γ)
// - f(α,Γ): smooth part (negative log-likelihood)
// - g(α) = λ_α * Σ_p ||α_{p·}||_2 (non-smooth group LASSO)
// - h(Γ) = λ_Γ * Σ_k ||Γ_{k·}||_2 (non-smooth group LASSO)
// convergence: stop when |L(t) - L(t-1)| < tolerance or when maximum iterations reached
 //' @param X Covariate matrix (N x p)
 //' @param status Event indicator matrix (N x K)
 //' @param rankmin Minimum rank for ties (N x K)
 //' @param rankmax Maximum rank for ties (N x K)
 //' @param order_list List of ordering vectors
 //' @param alpha0 Initial alpha matrix (p x R)
 //' @param Gamma0 Initial Gamma matrix (K x R)
 //' @param lambda_alpha_all Penalty parameter sequence for alpha
 //' @param lambda_Gamma_all Penalty parameter sequence for Gamma
 //' @param pfac_alpha Penalty factor for alpha (length p)
 //' @param pfac_Gamma Penalty factor for Gamma (length K)
 //' @param step_size Initial step size
 //' @param niter Maximum number of iterations
 //' @param tol Convergence tolerance
 //' @param linesearch_beta Line search shrinkage factor
 //' @param verbose Print progress
 //' @export
 // [[Rcpp::export]]
 Rcpp::List fit_reduced_rank(Rcpp::NumericMatrix X,
                             Rcpp::NumericMatrix status,
                             Rcpp::IntegerMatrix rankmin,
                             Rcpp::IntegerMatrix rankmax,
                             Rcpp::List order_list,
                             Rcpp::NumericMatrix alpha0,
                             Rcpp::NumericMatrix Gamma0,
                             Rcpp::NumericVector lambda_alpha_all,
                             Rcpp::NumericVector lambda_Gamma_all,
                             Rcpp::NumericVector pfac_alpha,
                             Rcpp::NumericVector pfac_Gamma,
                             double step_size = 0.01,
                             int niter = 500,
                             double tol = 1e-4,
                             double linesearch_beta = 0.5,
                             bool verbose = true)
 {
   // extract dimensions
   int N = X.rows();
   int p = X.cols();
   int K = status.cols();
   int R = alpha0.cols();
   
   // create problem object
   MCox_ReducedRank prob(N, K, p, R,
                         &X(0,0),
                         &status(0,0),
                         &rankmin(0,0),
                         &rankmax(0,0),
                         order_list);
   
   // initialize parameters
   // copy initial values from R to C++ Eigen matrices
   MapMatd alpha0_map(&alpha0(0,0), p, R);
   MapMatd Gamma0_map(&Gamma0(0,0), K, R);
   
   MatrixXd alpha(p, R);
   MatrixXd Gamma(K, R);
   alpha = alpha0_map; // copy initial alpha
   Gamma = Gamma0_map;
   
   // Nesterov acceleration points
   // v^(t+1) = x^(t+1) + β(x^(t+1) - x^(t))
   MatrixXd v_alpha(p, R);
   MatrixXd v_Gamma(K, R);
   v_alpha = alpha;
   v_Gamma = Gamma;
   
   MatrixXd alpha_prev(p, R);
   MatrixXd Gamma_prev(K, R);
   
   // gradients computed at accelerated point v
   MatrixXd grad_alpha(p, R);
   MatrixXd grad_Gamma(K, R);
   
   MatrixXd alpha_temp(p, R);
   MatrixXd Gamma_temp(K, R);
   
   VectorXd pf_alpha = as<VectorXd>(pfac_alpha);
   VectorXd pf_Gamma = as<VectorXd>(pfac_Gamma);
   
   const int nlambda = lambda_alpha_all.size();
   Rcpp::List result(nlambda);
   Rcpp::NumericVector obj_values(nlambda);
   Rcpp::IntegerVector num_iters(nlambda);
   
   struct timeval start, end;
   
   // iterate over lambda sequence
   for(int lam_ind = 0; lam_ind < nlambda; ++lam_ind){
     
     gettimeofday(&start, NULL);
     
     double lambda_alpha = lambda_alpha_all[lam_ind];
     double lambda_Gamma = lambda_Gamma_all[lam_ind];
     
     double current_step = step_size;
     double weight_old = 1.0;
     double weight_new;
     
     double obj_prev = R_PosInf;
     
     if(verbose){
       Rcpp::Rcout << "\n=== Lambda pair " << lam_ind + 1 << "/" << nlambda << " ===" << std::endl;
       Rcpp::Rcout << "lambda_alpha = " << lambda_alpha << ", lambda_Gamma = " << lambda_Gamma << std::endl;
     }
     
     // At iteration t:
     // 1. Gradient step: x̃ = v - η∇f(v)
     // 2. Proximal step: x^(t+1) = prox_{ηg}(x̃)
     // 3. Momentum: v^(t+1) = x^(t+1) + β(x^(t+1) - x^(t))
     for(int iter = 0; iter < niter; ++iter){
       
       // check for user interrupt
       Rcpp::checkUserInterrupt();
       
       alpha_prev = alpha;
       Gamma_prev = Gamma;
       
       // (1) compute gradients at accelerated point
       double cox_val = prob.get_gradients(v_alpha.data(), v_Gamma.data(),
                                           grad_alpha, grad_Gamma, true);
       
       // (2) backtracking line search: find appropriate step size η to ensure descent
       bool ls_success = false;
       double rhs_ls;
       int ls_iter = 0;
       const int max_ls_iter = 20;
       
       while(ls_iter < max_ls_iter){
         
         alpha_temp = v_alpha - current_step * grad_alpha;
         Gamma_temp = v_Gamma - current_step * grad_Gamma;

         prox_group_lasso(alpha, alpha_temp, current_step, lambda_alpha, pf_alpha);
         prox_group_lasso(Gamma, Gamma_temp, current_step, lambda_Gamma, pf_Gamma);
         
         double cox_val_new = prob.get_value_only(alpha.data(), Gamma.data());
         
         if(!std::isfinite(cox_val_new)){
           current_step *= linesearch_beta;
           ls_iter++;
           continue;
         }
         
         double grad_inner = (grad_alpha.array() * (alpha - v_alpha).array()).sum() +
           (grad_Gamma.array() * (Gamma - v_Gamma).array()).sum();
         
         double prox_penalty = (alpha - v_alpha).squaredNorm() + (Gamma - v_Gamma).squaredNorm();
         
         rhs_ls = cox_val + grad_inner + prox_penalty / (2.0 * current_step);
         
         if(cox_val_new <= rhs_ls){
           ls_success = true;
           break;
         }

         current_step *= linesearch_beta; // multiply by 0.5 (default)
         ls_iter++;
       }
       
       // line search warning
       if(!ls_success && verbose){
         Rcpp::Rcout << "Warning: Line search failed at iteration " << iter << std::endl;
       }
       
       // compute complete objective value
       double penalty_alpha = 0.0;
       double penalty_Gamma = 0.0;
       
       for(int i = 0; i < p; ++i){
         penalty_alpha += lambda_alpha * pf_alpha[i] * alpha.row(i).norm();
       }
       
       for(int k = 0; k < K; ++k){
         penalty_Gamma += lambda_Gamma * pf_Gamma[k] * Gamma.row(k).norm();
       }
       
       double obj_current = cox_val + penalty_alpha + penalty_Gamma;
       
       // check convergence
       double obj_change = std::abs(obj_current - obj_prev);
       
       if(obj_change < tol && iter > 0){
         num_iters[lam_ind] = iter + 1;
         obj_values[lam_ind] = obj_current;
         
         if(verbose){
           Rcpp::Rcout << "Converged at iteration " << iter + 1 
                       << ", objective = " << obj_current << std::endl;
         }
         break;
       }
       
       obj_prev = obj_current;
       
       // Nesterov Momentum update
       weight_new = 0.5 * (1.0 + sqrt(1.0 + 4.0 * weight_old * weight_old));
       double momentum = (weight_old - 1.0) / weight_new;
       
       v_alpha = alpha + momentum * (alpha - alpha_prev);
       v_Gamma = Gamma + momentum * (Gamma - Gamma_prev);
       
       weight_old = weight_new;
       
       // progress output (optional)
       if(verbose && iter % 10 == 0){
         Rcpp::Rcout << "Iter " << iter << ": obj = " << obj_current 
                     << ", change = " << obj_change << std::endl;
       }

       if(iter == niter - 1){
         num_iters[lam_ind] = niter;
         obj_values[lam_ind] = obj_current;
         
         if(verbose){
           Rcpp::Rcout << "Reached max iterations (" << niter << ")" << std::endl;
         }
       }
     }
     
     // compute elapsed time
     gettimeofday(&end, NULL);
     double elapsed = ((end.tv_sec - start.tv_sec) * 1000000u + 
                       end.tv_usec - start.tv_usec) / 1.0e6;
     
     if(verbose){
       Rcpp::Rcout << "Time elapsed: " << elapsed << " seconds\n" << std::endl;
     }
     
     // store results
     result[lam_ind] = Rcpp::List::create(
       Rcpp::Named("alpha") = alpha,
       Rcpp::Named("Gamma") = Gamma,
       Rcpp::Named("B") = alpha * Gamma.transpose()
     );
   }
   
   // return results
   return Rcpp::List::create(
     Rcpp::Named("result") = result,
     Rcpp::Named("objective_values") = obj_values,
     Rcpp::Named("num_iterations") = num_iters
   );
 }


// somtimes people want to inspect residuals after fitting
//' Compute Residuals for Reduced Rank Cox Model
 //' @export
 // [[Rcpp::export]]
 MatrixXd compute_residual_rr(Rcpp::NumericMatrix X,
                              Rcpp::NumericMatrix status,
                              Rcpp::IntegerMatrix rankmin,
                              Rcpp::IntegerMatrix rankmax,
                              Rcpp::List order_list,
                              MatrixXd alpha,
                              MatrixXd Gamma)
 {
   int N = X.rows();
   int p = X.cols();
   int K = status.cols();
   int R = alpha.cols();
   
   MCox_ReducedRank prob(N, K, p, R,
                         &X(0,0),
                         &status(0,0),
                         &rankmin(0,0),
                         &rankmax(0,0),
                         order_list);
   
   MapMatd alpha_map(alpha.data(), p, R);
   MapMatd Gamma_map(Gamma.data(), K, R);
   prob.get_gradients(alpha.data(), Gamma.data(), alpha, Gamma, false);
   
   return prob.get_residual_matrix();
 }


//' Compute Dual Norm (for lambda selection)
 //' @export
 // [[Rcpp::export]]
 VectorXd compute_dual_norm(MatrixXd grad,
                            double alpha,
                            double tol)
 {
   int p = grad.rows();
   VectorXd upperbound((grad.cwiseAbs().rowwise().maxCoeff()).cwiseMin(grad.rowwise().norm()/alpha));
   VectorXd dual_norm(p);
   
   // #pragma omp parallel for schedule(dynamic, 1)  
   for (int i = 0; i < p; ++i){
     double lower = 0.0;
     double upper = upperbound[i];
     if (upper <= tol){
       dual_norm[i] = 0.0;
     } else {
       int num_iter = (int)ceil(log2(upper/tol));
       for (int j = 0; j < num_iter; ++j){
         double bound = (lower + upper)/2;
         bool less = ((grad.row(i).array().abs() - bound).max(0).matrix().norm()) <= alpha * bound;
         if (less){
           upper = bound;
         } else {
           lower = bound;
         }
       }
       dual_norm[i] = (lower + upper)/2;
     }
   }
   return dual_norm;
 }