// Penalized Reduced Rank Multi-Outcome Cox Model
// Penalty: Group Lasso on alpha
// alpha: pxR matrix, Gamma: KxR matrix (Stiefel constraint)
// B = alpha %*% t(Gamma): pxK coefficient matrix

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <iostream>
#include <sys/time.h>
#include <cmath>
#include "group_lasso_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace Eigen;


// ============================================================
// Cox partial likelihood class
// ============================================================
class MCox_RR_GrpLasso
{
  const int N, K, p, R;
  
  MapMatd X;
  MapMatd status;
  MapMati rankmin;
  MapMati rankmax;
  MapMati rankmin_entry;
  
  MatrixXd entry_sorted;
  MatrixXd exit_sorted;
  
  std::vector<PermMat> orders;
  std::vector<PermMat> orders_entry;
  
  MatrixXd eta, exp_eta, risk_denom;
  MatrixXd outer_accumu, residual, Z, risk_entry;
  
  double get_residual_RR_Grplasso(const MapMatd &alpha, const MapMatd &Gamma, bool get_val = false){
    Z.noalias()   = X * alpha;
    eta.noalias() = Z * Gamma.transpose();
    
    for(int k = 0; k < K; ++k)
      eta.col(k) = orders[k] * eta.col(k);
    
    exp_eta.noalias() = eta.array().exp().matrix();
    
    for(int k = 0; k < K; ++k){
      // reverse cumsum by exit time
      double cur = 0;
      for(int i = 0; i < N; ++i){ cur += exp_eta(N-1-i,k); risk_denom(N-1-i,k) = cur; }
      for(int i = 0; i < N; ++i) risk_denom(i,k) = risk_denom(rankmin(i,k), k);
      
      // reverse cumsum by entry time
      VectorXd exp_eta_entry = orders_entry[k] * exp_eta.col(k);
      double cur_e = 0;
      for(int i = 0; i < N; ++i){ cur_e += exp_eta_entry(N-1-i); risk_entry(N-1-i,k) = cur_e; }
      for(int i = 0; i < N; ++i) risk_entry(i,k) = risk_entry(rankmin_entry(i,k), k);
      risk_entry.col(k) = orders_entry[k].transpose() * risk_entry.col(k);
      // final risk set denominator
      risk_denom.col(k) -= risk_entry.col(k);
      risk_denom.col(k) = risk_denom.col(k).cwiseMax(1e-10);
    }
    
    
    outer_accumu.noalias() = (status.array() / risk_denom.array()).matrix();
    
    for(int k = 0; k < K; ++k){
      double cur = 0;
      for(int i = 0; i < N; ++i){ cur += outer_accumu(i,k); outer_accumu(i,k) = cur; }
      for(int i = 0; i < N; ++i) outer_accumu(i,k) = outer_accumu(rankmax(i,k), k);
    }
    
    
    MatrixXd F_entry = MatrixXd::Zero(N, K);
    for (int k = 0; k < K; ++k) {
      VectorXd Fvals = outer_accumu.col(k);  
      
      const double* ex_begin = exit_sorted.col(k).data();
      const double* ex_end   = ex_begin + N;
      
      for (int i = 0; i < N; ++i) {
        double s_i = entry_sorted(i, k);
        auto it  = std::lower_bound(ex_begin, ex_end, s_i);
        int  pos = std::distance(ex_begin, it);
        F_entry(i, k) = (pos > 0) ? Fvals(pos - 1) : 0.0;
      }
      outer_accumu.col(k) -= F_entry.col(k);
    }
    
    // "residual"
    residual.noalias() = (outer_accumu.array() * exp_eta.array() - status.array()).matrix();
    
    for(int k = 0; k < K; ++k)
      residual.col(k) = orders[k].transpose() * residual.col(k);
    
    double cox_val = 0;
    if(get_val)
      cox_val = ((risk_denom.array().log() - eta.array()) * status.array()).sum();
    
    return cox_val;
  }
  
public:
  MCox_RR_GrpLasso(int N, int K, int p, int R,
                const double *X_, const double *status_,
                const int *rankmin_, const int *rankmax_,
                const int *rankmin_entry_,
                const double *entry_sorted_ptr,
                const double *exit_sorted_ptr, 
                const Rcpp::List order_list,
                const Rcpp::List order_list_entry)
    : N(N), K(K), p(p), R(R),
      X(X_, N, p), status(status_, N, K),
      rankmin(rankmin_, N, K), rankmax(rankmax_, N, K),
      rankmin_entry(rankmin_entry_, N, K),
      entry_sorted(Map<const MatrixXd>(entry_sorted_ptr, N, K)),
      exit_sorted(Map<const MatrixXd>(exit_sorted_ptr,  N, K)),
      eta(N,K), exp_eta(N,K), risk_denom(N,K),
      outer_accumu(N,K), residual(N,K), Z(N,R), risk_entry(N,K)
  {
    for(int k = 0; k < K; ++k) {
      orders.emplace_back(Rcpp::as<VectorXi>(order_list[k]));
      orders_entry.emplace_back(Rcpp::as<VectorXi>(order_list_entry[k]));
    }
  }
  
  double get_gradients_RR_GrpLasso(const double *alpha_ptr, const double *Gamma_ptr,
                                MatrixXd &grad_alpha, MatrixXd &grad_Gamma,
                                bool get_val = false){
    MapMatd alpha(alpha_ptr, p, R);
    MapMatd Gamma(Gamma_ptr, K, R);
    double cox_val = get_residual_RR_Grplasso(alpha, Gamma, get_val);
    grad_alpha.noalias() = X.transpose() * residual * Gamma;
    grad_Gamma.noalias() = residual.transpose() * Z;
    return cox_val;
  }
  
  double get_value_only_RR_GrpLasso(const double *alpha_ptr, const double *Gamma_ptr){
    MapMatd alpha(alpha_ptr, p, R);
    MapMatd Gamma(Gamma_ptr, K, R);
    Z.noalias()   = X * alpha;
    eta.noalias() = Z * Gamma.transpose();
    for(int k = 0; k < K; ++k)
      eta.col(k) = orders[k] * eta.col(k);
    exp_eta.noalias() = eta.array().exp().matrix();
    for(int k = 0; k < K; ++k){
      // reverse cumsum by exit time
      double cur = 0;
      for(int i = 0; i < N; ++i){ cur += exp_eta(N-1-i,k); risk_denom(N-1-i,k) = cur; }
      for(int i = 0; i < N; ++i) risk_denom(i,k) = risk_denom(rankmin(i,k), k);
      
      // reverse cumsum by entry time
      VectorXd exp_eta_entry = orders_entry[k] * exp_eta.col(k);
      double cur_e = 0;
      for(int i = 0; i < N; ++i){ cur_e += exp_eta_entry(N-1-i); risk_entry(N-1-i,k) = cur_e; }
      for(int i = 0; i < N; ++i) risk_entry(i,k) = risk_entry(rankmin_entry(i,k), k);
      risk_entry.col(k) = orders_entry[k].transpose() * risk_entry.col(k);
      // final risk set denominator
      risk_denom.col(k) -= risk_entry.col(k);
      risk_denom.col(k) = risk_denom.col(k).cwiseMax(1e-10);
    }
    return ((risk_denom.array().log() - eta.array()) * status.array()).sum();
  }
  
  MatrixXd get_residual_matrix_RR_GrpLasso(){ return residual; }
};

// ============================================================
// Group lasso proximal operator
// Pen(α) = λ * Σ_{g,r} ||α_{G_g, r}||_2
// prox: group soft-thresholding
// ============================================================
void prox_RR_GrpLasso(MatrixXd &M_out,
                            const MatrixXd &M_in,
                            double step_size,
                            double lambda,
                            const std::vector<std::vector<int>> &groups)
{
  int R = M_in.cols();
  int num_groups = groups.size();
  for(int g = 0; g < num_groups; ++g){
    const std::vector<int>& idx = groups[g];
    double threshold = step_size * lambda;

    double frob_sq = 0.0;
    for(int i : idx)
      for(int r = 0; r < R; ++r)
        frob_sq += M_in(i,r) * M_in(i,r);
    double frob = std::sqrt(frob_sq);
    if(frob <= threshold){
      for(int i : idx)
        for(int r = 0; r < R; ++r)
          M_out(i,r) = 0.0;
    } else {
      double scale = 1.0 - threshold / frob;
      for(int i : idx)
        for(int r = 0; r < R; ++r)
          M_out(i,r) = M_in(i,r) * scale;
    }
  }
}

std::vector<std::vector<int>> build_groups(const Rcpp::IntegerVector &group_labels,
                                           int num_groups)
{
  std::vector<std::vector<int>> groups(num_groups);
  for(int i = 0; i < group_labels.size(); ++i){
    int g = group_labels[i] - 1;
    groups[g].push_back(i);
  }
  return groups;
}

// ============================================================
// Riemannian gradient on Stiefel manifold St(K,R)
// riem_grad = grad - Gamma * (grad^T * Gamma)
// ============================================================
MatrixXd riemannian_gradient(const MatrixXd &grad, const MatrixXd &Gamma)
{
  MatrixXd A = Gamma.transpose() * grad;
  return grad - Gamma * (A + A.transpose()) * 0.5;
}


// ============================================================
// Main fitting function
// ============================================================
//' @param X Covariate matrix (N x p)
 //' @param status Event indicator matrix (N x K), normalized by number of events
 //' @param rankmin Minimum rank for ties (N x K), 0-indexed
 //' @param rankmax Maximum rank for ties (N x K), 0-indexed
 //' @param order_list List of K ordering vectors (0-indexed permutations)
 //' @param alpha0 Initial alpha matrix (p x R)
 //' @param Gamma0 Initial Gamma matrix (K x R), Gamma0^T Gamma0 = I
 //' @param lambda_alpha_all Penalty parameter sequence for alpha
 //' @param group_labels Integer vector of length p, group membership (1-indexed)
 //' @param step_size Initial step size for backtracking line search
 //' @param niter Maximum number of iterations per lambda
 //' @param tol Convergence tolerance (gradient mapping norm)
 //' @param linesearch_beta Step size shrinkage factor
 //' @param verbose Print progress
 //' @export
 // [[Rcpp::export]]
 Rcpp::List fit_RR_GrpLasso(Rcpp::NumericMatrix X,
                         Rcpp::NumericMatrix status,
                         Rcpp::IntegerMatrix rankmin,
                         Rcpp::IntegerMatrix rankmax,
                         Rcpp::IntegerMatrix rankmin_entry,
                         Rcpp::NumericMatrix entry_sorted_mat,
                         Rcpp::NumericMatrix exit_sorted_mat, 
                         Rcpp::List order_list,
                         Rcpp::List order_list_entry,
                         Rcpp::NumericMatrix alpha0,
                         Rcpp::NumericMatrix Gamma0,
                         Rcpp::NumericVector lambda_alpha_all,
                         Rcpp::IntegerVector group_labels,
                         double step_size       = 0.5,
                         int    niter           = 500,
                         double tol             = 1e-4,
                         double linesearch_beta = 0.5,
                         bool   verbose         = true)
 {
   int N = X.rows(), p = X.cols(), K = status.cols(), R = alpha0.cols();
   
   int num_groups = *std::max_element(group_labels.begin(), group_labels.end());
   std::vector<std::vector<int>> groups = build_groups(group_labels, num_groups);
   
   if(verbose)
     Rcpp::Rcout << "p=" << p << "  G=" << num_groups << "  R=" << R
                 << "  K=" << K << "  penalty=group_lasso" << std::endl;
     
     MCox_RR_GrpLasso prob(N, K, p, R,
                           &X(0,0), &status(0,0),
                           &rankmin(0,0), &rankmax(0,0), &rankmin_entry(0,0),
                           &entry_sorted_mat(0,0), &exit_sorted_mat(0,0),
                           order_list, order_list_entry);
     
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
       double current_step  = step_size;
       double weight_old    = 1.0, weight_new;
       v_alpha = alpha;
       v_Gamma = Gamma;
       double obj_prev      = R_PosInf;
       double cox_val_new   = R_PosInf;
       int    ls_fail_count = 0;
       
       for(int iter = 0; iter < niter; ++iter){
         
         Rcpp::checkUserInterrupt();
         current_step = step_size;   // reset each iteration
         alpha_prev = alpha;
         Gamma_prev = Gamma;
         
         // (1) gradients at accelerated point v
         double cox_val = prob.get_gradients_RR_GrpLasso(v_alpha.data(), v_Gamma.data(),
                                                         grad_alpha, grad_Gamma, true);
         
         // (2) Riemannian gradient for Gamma
         MatrixXd riem_grad_Gamma = riemannian_gradient(grad_Gamma, v_Gamma);
         
         // (3) backtracking line search (Armijo condition on smooth part)
         bool ls_success = false;
         for(int ls_iter = 0; ls_iter < 20; ++ls_iter){
           
           alpha_temp = v_alpha - current_step * grad_alpha;
           Gamma_temp = v_Gamma - current_step * riem_grad_Gamma;
           
           // alpha: group lasso proximal step
           prox_RR_GrpLasso(alpha, alpha_temp, current_step, lambda_alpha, groups);
           
           // Gamma: SVD retraction onto Stiefel manifold
           Eigen::JacobiSVD<MatrixXd> svd(Gamma_temp,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
           Gamma = svd.matrixU() * svd.matrixV().transpose();
           
           cox_val_new = prob.get_value_only_RR_GrpLasso(alpha.data(), Gamma.data());
           
           if(!std::isfinite(cox_val_new)){
             current_step *= linesearch_beta; continue;
           }
           
           double grad_inner = (grad_alpha.array()      * (alpha - v_alpha).array()).sum()
             + (riem_grad_Gamma.array()  * (Gamma - v_Gamma).array()).sum();
           double prox_sq    = (alpha - v_alpha).squaredNorm()
             + (Gamma - v_Gamma).squaredNorm();
           double rhs_ls     = cox_val + grad_inner + prox_sq / (2.0 * current_step);
           
           if(cox_val_new <= rhs_ls){ ls_success = true; break; }
           current_step *= linesearch_beta;
         }
         
         // line search failed: roll back and skip this iteration
         // prevents bad (alpha, Gamma) from polluting grad_map and subsequent steps
         if(!ls_success){
           alpha = alpha_prev;
           Gamma = Gamma_prev;
           ls_fail_count++;
           weight_old = 1.0;
           v_alpha = alpha;
           v_Gamma = Gamma;
           continue;
         }
         
         // (4) complete objective: likelihood + group lasso penalty
         double penalty_alpha = 0.0;
         for(int g = 0; g < num_groups; ++g){
           const std::vector<int>& idx = groups[g];
           double frob_sq = 0.0;
           for(int i : idx)
             for(int r = 0; r < R; ++r)
               frob_sq += alpha(i,r) * alpha(i,r);
           penalty_alpha += lambda_alpha * std::sqrt(frob_sq);
         }
         double obj_current = cox_val_new + penalty_alpha;
         
         // (5) gradient mapping norm — valid because line search succeeded
         double grad_map_norm = (alpha - alpha_prev).norm() / current_step
         + (Gamma - Gamma_prev).norm() / current_step;
         
         // (6) Nesterov update with restart on objective increase
         // Reference: O'Donoghue & Candes (2015)
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
         
         // (7) convergence check
         if(grad_map_norm < tol && iter > 0){
           num_iters[lam_ind]  = iter + 1;
           obj_values[lam_ind] = obj_current;
           if(verbose){
             Rcpp::Rcout << "Lambda " << lam_ind+1 << "/" << nlambda
                         << "  lambda=" << lambda_alpha
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
                         << "  lambda=" << lambda_alpha
                         << "  max_iter=" << niter
                         << "  obj=" << obj_current
                         << "  grad_map=" << grad_map_norm;
             if(ls_fail_count > 0)
               Rcpp::Rcout << "  ls_fail=" << ls_fail_count;
             Rcpp::Rcout << std::endl;
           }
         }
       }  // end iter loop
       
       result[lam_ind] = Rcpp::List::create(
         Rcpp::Named("alpha") = alpha,
         Rcpp::Named("Gamma") = Gamma,
         Rcpp::Named("B")     = alpha * Gamma.transpose()
       );
     }
     
     return Rcpp::List::create(
       Rcpp::Named("result")           = result,
       Rcpp::Named("objective_values") = obj_values,
       Rcpp::Named("num_iterations")   = num_iters
     );
 }



// ============================================================
// Compute residuals (diagnostic)
// ============================================================
//' @export
 // [[Rcpp::export]]
 MatrixXd compute_residual_RR_GrpLasso(Rcpp::NumericMatrix X,
                                    Rcpp::NumericMatrix status,
                                    Rcpp::IntegerMatrix rankmin,
                                    Rcpp::IntegerMatrix rankmax,
                                    Rcpp::IntegerMatrix rankmin_entry, 
                                    Rcpp::NumericMatrix entry_sorted_mat,
                                    Rcpp::NumericMatrix exit_sorted_mat, 
                                    Rcpp::List order_list,
                                    Rcpp::List order_list_entry,
                                    MatrixXd alpha,
                                    MatrixXd Gamma)
 {
   int N = X.rows(), p = X.cols(), K = status.cols(), R = alpha.cols();
   MCox_RR_GrpLasso prob(N, K, p, R,
                      &X(0,0), &status(0,0),
                      &rankmin(0,0), &rankmax(0,0), &rankmin_entry(0,0),
                      &entry_sorted_mat(0,0), &exit_sorted_mat(0,0),
                      order_list, order_list_entry);
   MatrixXd grad_alpha_tmp(p, R), grad_Gamma_tmp(K, R);
   prob.get_gradients_RR_GrpLasso(alpha.data(), Gamma.data(),
                               grad_alpha_tmp, grad_Gamma_tmp, false);
   return prob.get_residual_matrix_RR_GrpLasso();
 }

