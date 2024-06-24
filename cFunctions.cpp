//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;

// [[Rcpp::export]]
vec logClassificationCHelper(vec depression, double emit_nd){
	int n = depression.n_elem;
	vec log_dens_vec(n);
	log_dens_vec.zeros();

	depression.replace(datum::nan, -1);

	for (int i = 0; i < n; i++) {
		double depression_obs = depression(i);


		if (depression_obs == 0){
			log_dens_vec(i) = log( emit_nd );
		} else if (depression_obs == 1){
			log_dens_vec(i) = log( 1 - emit_nd );
		} else if (depression_obs == -1){
			log_dens_vec(i) = 0;
		}
		
	}

return log_dens_vec;
}

// [[Rcpp::export]]
vec logClassificationC(mat atbc_data, vec emit_nd){
  
  vec log_dens = logClassificationCHelper(atbc_data.col( 0 ), emit_nd[0]) + logClassificationCHelper(atbc_data.col( 1 ), emit_nd[1] ) + logClassificationCHelper(atbc_data.col( 2 ), emit_nd[2] ) + logClassificationCHelper(atbc_data.col( 3 ), emit_nd[3] )+ logClassificationCHelper(atbc_data.col( 4 ), emit_nd[4] ) + logClassificationCHelper(atbc_data.col( 5 ), emit_nd[5] )+ logClassificationCHelper(atbc_data.col( 6 ), emit_nd[6] );
  
  return log_dens;
}



/* This is from the seqHMM github*/
#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#  define EXPL expl
#else
#  define LDOUBLE double
#  define EXPL exp
#endif

// [[Rcpp::export]]
double logSumExpC(const arma::vec& x) {
  unsigned int maxi = x.index_max();
  LDOUBLE maxv = x(maxi);
  if (!(maxv > -arma::datum::inf)) {
    return -arma::datum::inf;
  }
  LDOUBLE cumsum = 0.0;
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) && (x(i) > -arma::datum::inf)) {
      cumsum += EXPL(x(i) - maxv);
    }
  }
  
  return maxv + log1p(cumsum);
}



// [[Rcpp::export]]
mat ForwardIndC(mat dep_ind, rowvec init, mat tran, vec emit_nd0_vec, vec emit_nd1_vec){

	int n = dep_ind.n_rows;
	mat alpha( n, 2 , fill::zeros);

	vec log_class_0 = logClassificationC( dep_ind, emit_nd0_vec);

	vec log_class_1 = logClassificationC( dep_ind, emit_nd1_vec);

	alpha(0,0) = log(init[0]) + log_class_0[0];
	alpha(0,1) = log(init[1]) + log_class_1[0];
	
	
	for (int i = 1; i < n; i++) {
	
	  double fp_00 = alpha(i-1,0) + log(tran(0,0)) + log_class_0[i];
	  
	  double fp_10 = alpha(i-1,1) + log(tran(1,0)) + log_class_0[i];
	  
	  double fp_01 = alpha(i-1,0) + log(tran(0,1)) + log_class_1[i];
	  
	  double fp_11 = alpha(i-1,1) + log(tran(1,1)) + log_class_1[i];
	  
	  vec fp_0 = { fp_00, fp_10 };
	  vec fp_1 = { fp_01, fp_11 };
	  
	  alpha(i,0) = logSumExpC(fp_0);
	  alpha(i,1) = logSumExpC(fp_1);
	  
	  
	}

	return(alpha);
}

// [[Rcpp::export]]
mat BackwardIndC(mat dep_ind, mat tran, vec emit_nd0_vec, vec emit_nd1_vec){
  
  int n = dep_ind.n_rows; 
  mat beta( n, 2 , fill::zeros);

  
  vec log_class_0 = logClassificationC( dep_ind, emit_nd0_vec );
  vec log_class_1 = logClassificationC( dep_ind, emit_nd1_vec );

  
  beta(n-1,0) = 0.0;
  beta(n-1,1) = 0.0;
  
  for (int i = n-2; i >= 0; i--) {
    
    double bp_00 = log(tran(0,0)) + log_class_0[i+1] + beta(i+1,0);
    
    double bp_01 = log(tran(0,1)) + log_class_1[i+1] + beta(i+1,1);
    
    double bp_10 = log(tran(1,0)) + log_class_0[i+1] + beta(i+1,0);
    
    double bp_11 = log(tran(1,1)) + log_class_1[i+1] + beta(i+1,1);
    
    vec bp_0 = { bp_00, bp_01 };
    vec bp_1 = { bp_10, bp_11 };
    
    beta(i,0) = logSumExpC(bp_0);
    beta(i,1) = logSumExpC(bp_1);

  }
  
  return beta;

}

// [[Rcpp::export]]
List ForwardC(cube depression, mat init, cube tran_array, vec emit_nd0_vec, vec emit_nd1_vec ){

	int num_people = depression.n_cols;
	int len = depression.n_rows;
	int num_re = tran_array.n_slices;
	List alpha_list(num_people);

	for (int ind = 0; ind < num_people; ind++) {
		arma::cube Cube1(len, 2, num_re);
		mat dep_ind = depression.col_as_mat(ind);

		
		for (int clust_i = 0; clust_i < num_re; clust_i++){
			mat tran = tran_array.slice(clust_i);
			rowvec init_vec = init.row(clust_i);

			Cube1.slice(clust_i) = ForwardIndC(dep_ind, init_vec, tran, emit_nd0_vec, emit_nd1_vec);

		}

		alpha_list(ind) = Cube1;
	}
	return(alpha_list);
}

// [[Rcpp::export]]
List BackwardC(cube depression, cube tran_array, vec emit_nd0_vec, vec emit_nd1_vec ){

	int num_people = depression.n_cols;
	int len = depression.n_rows;
	int num_re = tran_array.n_slices;
	List beta_list(num_people);
	
	for (int ind = 0; ind < num_people; ind++) {
		arma::cube Cube1(len, 2, num_re);
		mat dep_ind = depression.col_as_mat(ind);
		
		for (int clust_i = 0; clust_i < num_re; clust_i++){
			mat tran = tran_array.slice(clust_i);

			Cube1.slice(clust_i) = BackwardIndC(dep_ind,  tran, emit_nd0_vec, emit_nd1_vec);

		}

		beta_list(ind) = Cube1;
	}
	return(beta_list);
}  


// [[Rcpp::export]]
cube CalcTranHelperC(int init_state, int new_state, cube depression, cube tran_array, vec emit_nd_vec, vec ind_like_vec, List alpha, List beta, vec pi_l){
  int num_people = depression.n_cols;
  int len = depression.n_rows;
  int num_re = tran_array.n_slices;

  arma::cube tran_vals_re_cube( len-1, num_people, num_re );

  for (int clust_i = 0; clust_i < num_re; clust_i++){

	mat tran_vals_re_mat( len-1, num_people );
	double tran_val = tran_array(init_state,new_state,clust_i);

    for (int ind = 0; ind < num_people; ind++) {
	  
	  arma::cube alpha_ind = alpha(ind); 
	  arma::cube beta_ind = beta(ind);
	  double likelihood = ind_like_vec(ind);


	  mat depression_ind_m1 = depression.col_as_mat(ind);
	  depression_ind_m1.shed_row(0);
	
	  vec class_vec = logClassificationC( depression_ind_m1, emit_nd_vec );
	  
	  vec alpha_ind_slice = alpha_ind(span(0,len-2),span(init_state,init_state),span(clust_i,clust_i));
	  vec beta_ind_slice = beta_ind(span(1,len-1),span(new_state,new_state),span(clust_i,clust_i));

	  vec temp = alpha_ind_slice + beta_ind_slice + log(tran_val) + log(pi_l(clust_i)) + class_vec - likelihood;
	  vec tran_vals_re_ind = arma::exp(temp);

	  tran_vals_re_mat.col(ind) = tran_vals_re_ind;
	
    }

	tran_vals_re_cube.slice(clust_i) = tran_vals_re_mat;
    
  }

  return tran_vals_re_cube;
}


// // // // [[Rcpp::export]]
// // // cube CalcTranIndHelperC(int init_state, int new_state, NumericMatrix act, List tran_list_mat, NumericVector tran_ind_vec, cube emit_act, NumericVector ind_like_vec, List alpha, List beta, double lepsilon, NumericVector act_light_binom, vec pi_l){
// // //   int num_people = act.ncol();
// // //   int len = act.nrow();
// // //   int num_re = 1;
// // //   int clust_i = 0;

// // //   arma::cube tran_vals_re_cube( len-1, num_people, num_re );


// // // 	mat tran_vals_re_mat( len-1, num_people );

// // //     for (int ind = 0; ind < num_people; ind++) {
      
// // //       int tran_ind = tran_ind_vec(ind);
// // // 	  mat tran_mat = tran_list_mat(tran_ind-1);
	  

// // // 	  //0,0->0 & 1,0->1 & 0,1->2 & 1,1->3
// // // 	  int tran_vec_ind = init_state + (new_state * 2);

// // // 	  arma::cube alpha_ind = alpha(ind); 
// // // 	  arma::cube beta_ind = beta(ind);
// // // 	  double likelihood = ind_like_vec(ind);
	  



// // // 	  NumericMatrix act_ind = act( Range(1,len-1) , Range(ind,ind) );
// // // 	  NumericVector act_ind_m1 = act_ind.column(0); 
	  
	  

