// Penalized Reduced Rank Multi-Outcome Cox Model
// alpha: pxR matrix (lasso penalty)
// Gamma: KxR matrix (lasso penalty, NO Stiefel constraint)
// B = alpha %*% t(Gamma): pxK coefficient matrix

#include <Rcpp.h>
#include <vector>
#include <iostream>
#include <sys/time.h>
#include <cmath>
#include "lasso_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace Eigen;


// ============================================================
// Cox partial likelihood class
// ============================================================
class MCox_RR_Lasso
{
  const int N, K, p, R;
  
  MapMatd X;
  MapMatd status;
  MapMati rankmin;
  MapMati rankmax;
  std::vector<PermMat> orders;
  
  MatrixXd eta, exp_eta, risk_denom, outer_accumu, residual, Z;
  
  double get_residual_RR_lasso(const MapMatd &alpha, const MapMatd &Gamma, bool get_val = false){
    Z.noalias()   = X * alpha;
    eta.noalias() = Z * Gamma.transpose();
    
    for(int k = 0; k < K; ++k)
      eta.col(k) = orders[k] * eta.col(k);
    
    exp_eta.noalias() = eta.array().exp().matrix();
    
    for(int k = 0; k < K; ++k){
      double cur = 0;
      for(int i = 0; i < N; ++i){ cur += exp_eta(N-1-i,k); risk_denom(N-1-i,k) = cur; }
      for(int i = 0; i < N; ++i) risk_denom(i,k) = risk_denom(rankmin(i,k), k);
    }
    
    outer_accumu.noalias() = (status.array() / risk_denom.array()).matrix();
    
    for(int k = 0; k < K; ++k){
      double cur = 0;
      for(int i = 0; i < N; ++i){ cur += outer_accumu(i,k); outer_accumu(i,k) = cur; }
      for(int i = 0; i < N; ++i) outer_accumu(i,k) = outer_accumu(rankmax(i,k), k);
    }
    
    residual.noalias() = (outer_accumu.array() * exp_eta.array() - status.array()).matrix();
    
    for(int k = 0; k < K; ++k)
      residual.col(k) = orders[k].transpose() * residual.col(k);
    
    double cox_val = 0;
    if(get_val)
      cox_val = ((risk_denom.array().log() - eta.array()) * status.array()).sum();
    
    return cox_val;
  }
  
public:
  MCox_RR_Lasso(int N, int K, int p, int R,
             const double *X_, const double *status_,
             const int *rankmin_, const int *rankmax_,
             const Rcpp::List order_list)
    : N(N), K(K), p(p), R(R),
      X(X_, N, p), status(status_, N, K),
      rankmin(rankmin_, N, K), rankmax(rankmax_, N, K),
      eta(N,K), exp_eta(N,K), risk_denom(N,K),
      outer_accumu(N,K), residual(N,K), Z(N,R)
  {
    for(int k = 0; k < K; ++k)
      orders.emplace_back(Rcpp::as<VectorXi>(order_list[k]));
  }
  
  double get_gradients_RR_Lasso(const double *alpha_ptr, const double *Gamma_ptr,
                       MatrixXd &grad_alpha, MatrixXd &grad_Gamma,
                       bool get_val = false){
    MapMatd alpha(alpha_ptr, p, R);
    MapMatd Gamma(Gamma_ptr, K, R);
    double cox_val = get_residual_RR_lasso(alpha, Gamma, get_val);
    grad_alpha.noalias() = X.transpose() * residual * Gamma;
    grad_Gamma.noalias() = residual.transpose() * Z;
    return cox_val;
  }
  
  double get_value_only_RR_Lasso(const double *alpha_ptr, const double *Gamma_ptr){
    MapMatd alpha(alpha_ptr, p, R);
    MapMatd Gamma(Gamma_ptr, K, R);
    Z.noalias()   = X * alpha;
    eta.noalias() = Z * Gamma.transpose();
    for(int k = 0; k < K; ++k)
      eta.col(k) = orders[k] * eta.col(k);
    exp_eta.noalias() = eta.array().exp().matrix();
    for(int k = 0; k < K; ++k){
      double cur = 0;
      for(int i = 0; i < N; ++i){ cur += exp_eta(N-1-i,k); risk_denom(N-1-i,k) = cur; }
      for(int i = 0; i < N; ++i) risk_denom(i,k) = risk_denom(rankmin(i,k), k);
    }
    return ((risk_denom.array().log() - eta.array()) * status.array()).sum();
  }
  
  MatrixXd get_residual_matrix_RR_Lasso(){ return residual; }
};


// ============================================================
// Lasso proximal operator (element-wise soft-thresholding)
// prox_{η λ}(v)_j = sign(v_j) * max(|v_j| - ηλ, 0)
// Applied to both alpha AND Gamma (same function, no SVD needed for Gamma)
// ============================================================
void prox_RR_Lasso(MatrixXd &M_out,
                const MatrixXd &M_in,
                double step_size,
                double lambda)
{
  double threshold = step_size * lambda;
  M_out = M_in.array().sign() * (M_in.array().abs() - threshold).max(0.0);
}


// ============================================================
// Main fitting function
// ============================================================
//' @export
// [[Rcpp::export]]
 Rcpp::List fit_RR_Lasso(Rcpp::NumericMatrix X,
                                   Rcpp::NumericMatrix status,
                                   Rcpp::IntegerMatrix rankmin,
                                   Rcpp::IntegerMatrix rankmax,
                                   Rcpp::List order_list,
                                   Rcpp::NumericMatrix alpha0,
                                   Rcpp::NumericMatrix Gamma0,
                                   Rcpp::NumericVector lambda_alpha_all,
                                   Rcpp::NumericVector lambda_gamma_all,
                                   double step_size       = 0.5,
                                   int    niter           = 500,
                                   double tol             = 1e-4,
                                   double linesearch_beta = 0.5,
                                   bool   verbose         = true)
 {
   int N = X.rows(), p = X.cols(), K = status.cols(), R = alpha0.cols();
   
   if(lambda_alpha_all.size() != lambda_gamma_all.size())
     Rcpp::stop("lambda_alpha_all and lambda_gamma_all must have the same length");
   
   if(verbose)
     Rcpp::Rcout << "p=" << p << "  R=" << R << "  K=" << K
                 << "  penalty=lasso (alpha + Gamma)" << std::endl;
     
     MCox_RR_Lasso prob(N, K, p, R,
                     &X(0,0), &status(0,0),
                     &rankmin(0,0), &rankmax(0,0),
                     order_list);
     
     MapMatd alpha0_map(&alpha0(0,0), p, R);
     MapMatd Gamma0_map(&Gamma0(0,0), K, R);
     
     MatrixXd alpha(p,R), Gamma(K,R);
     alpha = alpha0_map;
     Gamma = Gamma0_map;
     
     MatrixXd v_alpha(p,R), v_Gamma(K,R);
     v_alpha = alpha;  v_Gamma = Gamma;
     
     MatrixXd alpha_prev(p,R), Gamma_prev(K,R);
     MatrixXd grad_alpha(p,R), grad_Gamma(K,R);
     MatrixXd alpha_temp(p,R), Gamma_temp(K,R);
     
     const int nlambda = lambda_alpha_all.size();
     Rcpp::List result(nlambda);
     Rcpp::NumericVector obj_values(nlambda);
     Rcpp::IntegerVector num_iters(nlambda);
     
     struct timeval t0, t1;
     
     for(int lam_ind = 0; lam_ind < nlambda; ++lam_ind){
       
       gettimeofday(&t0, NULL);
       
       double lambda_alpha  = lambda_alpha_all[lam_ind];
       double lambda_gamma  = lambda_gamma_all[lam_ind];
       double current_step  = step_size;
       double weight_old    = 1.0, weight_new;
       double obj_prev      = R_PosInf;
       double cox_val_new   = R_PosInf;
       int    ls_fail_count = 0;
       
       for(int iter = 0; iter < niter; ++iter){
         
         Rcpp::checkUserInterrupt();
         current_step = step_size;   // reset each iteration
         alpha_prev = alpha;
         Gamma_prev = Gamma;
         
         // (1) gradients at accelerated point v
         double cox_val = prob.get_gradients_RR_Lasso(v_alpha.data(), v_Gamma.data(),
                                             grad_alpha, grad_Gamma, true);
         
         // (2) backtracking line search
         // Gamma uses plain gradient step + soft-threshold (no Riemannian / SVD)
         bool ls_success = false;
         for(int ls_iter = 0; ls_iter < 20; ++ls_iter){
           
           alpha_temp = v_alpha - current_step * grad_alpha;
           Gamma_temp = v_Gamma - current_step * grad_Gamma;
           
           // alpha: lasso proximal step
           prox_RR_Lasso(alpha, alpha_temp, current_step, lambda_alpha);
           
           // Gamma: lasso proximal step (no orthogonality constraint)
           prox_RR_Lasso(Gamma, Gamma_temp, current_step, lambda_gamma);
           
           cox_val_new = prob.get_value_only_RR_Lasso(alpha.data(), Gamma.data());
           
           if(!std::isfinite(cox_val_new)){
             current_step *= linesearch_beta; continue;
           }
           
           double grad_inner = (grad_alpha.array() * (alpha - v_alpha).array()).sum()
             + (grad_Gamma.array() * (Gamma - v_Gamma).array()).sum();
           double prox_sq    = (alpha - v_alpha).squaredNorm()
             + (Gamma - v_Gamma).squaredNorm();
           double rhs_ls     = cox_val + grad_inner + prox_sq / (2.0 * current_step);
           
           if(cox_val_new <= rhs_ls){ ls_success = true; break; }
           current_step *= linesearch_beta;
         }
         
         // line search failed: roll back and skip this iteration entirely
         if(!ls_success){
           alpha = alpha_prev;
           Gamma = Gamma_prev;
           ls_fail_count++;
           weight_old = 1.0;
           v_alpha = alpha;
           v_Gamma = Gamma;
           continue;
         }
         
         // (3) complete objective: likelihood + lasso penalties on alpha and Gamma
         double penalty_alpha = lambda_alpha * alpha.cwiseAbs().sum();
         double penalty_gamma = lambda_gamma * Gamma.cwiseAbs().sum();
         double obj_current   = cox_val_new + penalty_alpha + penalty_gamma;
         
         // (4) gradient mapping norm
         double grad_map_norm = (alpha - alpha_prev).norm() / current_step
         + (Gamma - Gamma_prev).norm() / current_step;
         
         // (5) Nesterov update with restart on objective increase
         if(obj_current > obj_prev + 1e-10){
           weight_old = 1.0;
           v_alpha = alpha;
           v_Gamma = Gamma;
         } else {
           weight_new = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * weight_old * weight_old));
           double momentum = (weight_old - 1.0) / weight_new;
           v_alpha    = alpha + momentum * (alpha - alpha_prev);
           v_Gamma    = Gamma + momentum * (Gamma - Gamma_prev);
           weight_old = weight_new;
         }
         
         obj_prev = obj_current;
         
         // (6) convergence check
         if(grad_map_norm < tol && iter > 0){
           num_iters[lam_ind]  = iter + 1;
           obj_values[lam_ind] = obj_current;
           if(verbose){
             Rcpp::Rcout << "Lambda " << lam_ind+1 << "/" << nlambda
                         << "  lambda_alpha=" << lambda_alpha
                         << "  lambda_gamma=" << lambda_gamma
                         << "  converged  iter=" << iter+1
                         << "  obj=" << obj_current
                         << "  grad_map=" << grad_map_norm;
             if(ls_fail_count > 0)
               Rcpp::Rcout << "  ls_fail=" << ls_fail_count;
             Rcpp::Rcout << std::endl;
           }
           break;
         }
         
         if(iter == niter - 1){
           num_iters[lam_ind]  = niter;
           obj_values[lam_ind] = obj_current;
           if(verbose){
             Rcpp::Rcout << "Lambda " << lam_ind+1 << "/" << nlambda
                         << "  lambda_alpha=" << lambda_alpha
                         << "  lambda_gamma=" << lambda_gamma
                         << "  MAX_ITER=" << niter
                         << "  obj=" << obj_current
                         << "  grad_map=" << grad_map_norm;
             if(ls_fail_count > 0)
               Rcpp::Rcout << "  ls_fail=" << ls_fail_count;
             Rcpp::Rcout << std::endl;
           }
         }
       }  // end iter loop
       
       result[lam_ind] = Rcpp::List::create(
         Rcpp::Named("alpha")        = alpha,
         Rcpp::Named("Gamma")        = Gamma,
         Rcpp::Named("B")            = alpha * Gamma.transpose(),
         Rcpp::Named("lambda_alpha") = lambda_alpha,
         Rcpp::Named("lambda_gamma") = lambda_gamma
       );
     }
     
     return Rcpp::List::create(
       Rcpp::Named("result")            = result,
       Rcpp::Named("objective_values")  = obj_values,
       Rcpp::Named("num_iterations")    = num_iters
     );
 }


// ============================================================
// Compute residuals (diagnostic)
// ============================================================
//' @export
 // [[Rcpp::export]]
 MatrixXd compute_residual_RR_Lasso(Rcpp::NumericMatrix X,
                                 Rcpp::NumericMatrix status,
                                 Rcpp::IntegerMatrix rankmin,
                                 Rcpp::IntegerMatrix rankmax,
                                 Rcpp::List order_list,
                                 MatrixXd alpha,
                                 MatrixXd Gamma)
 {
   int N = X.rows(), p = X.cols(), K = status.cols(), R = alpha.cols();
   MCox_RR_Lasso prob(N, K, p, R,
                   &X(0,0), &status(0,0),
                   &rankmin(0,0), &rankmax(0,0),
                   order_list);
   MatrixXd grad_alpha_tmp(p, R), grad_Gamma_tmp(K, R);
   prob.get_gradients_RR_Lasso(alpha.data(), Gamma.data(),
                      grad_alpha_tmp, grad_Gamma_tmp, false);
   return prob.get_residual_matrix_RR_Lasso();
 }
