//// Cpp functions 2015-3-11 Thong Pham
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export(".normalized_constant")]]
int normalized_constant(      NumericVector& norm, 
                        const NumericMatrix& degree, 
                        const NumericVector& theta,
                        const NumericVector& f, 
                        const NumericMatrix& offset_tk
                        ) {
    long T = degree.nrow();     // number of time-steps
    long N = degree.ncol();     // number of nodes
    long K = offset_tk.ncol();  // maximum degree 
    #pragma omp parallel for
    for (long i = 0; i < T; i++) {
        double total = 0;
        for (long j = 0; j < N; j++)
            if (degree(i,j) >= 0) {
                total += theta(degree(i,j))*f(j);
            }
        for (long k = 0; k < K; ++k)
            total += offset_tk(i,k)*theta(k);
        norm(i) = total;
    }
    return 0;
}

// [[Rcpp::export(".get_stats")]]
int get_stats(CharacterVector    & time_stamp,
              CharacterVector    & unique_stamp,
              const NumericVector& in_node, 
              const NumericVector& out_node,
              const NumericVector& all_node,
              const NumericVector& ok_node,
              const NumericVector& bin_vector,
              const long max_node_id,
              const int  undirected,
              const int  only_PA,
              CharacterVector& time_vector,
              NumericVector& Sum_m_k, 
              NumericMatrix& n_tk,
              NumericVector& m_tk,
              NumericVector& m_t,
              NumericMatrix& offset_tk, 
              NumericVector& z_j, 
              NumericMatrix& node_degree
              ) {
  long N     = all_node.size(); 
  long N_new = ok_node.size();
  long T     = time_vector.size();
  long K     = n_tk.ncol();
  std::vector<long> node_array(max_node_id + 1,0); //the index used when indexing all arrays whose length is N
  std::vector<int>  ok_array(max_node_id + 1,0);   //boolean check whether a node is used or not  
  std::vector<int>  is_appear(N,0);
  std::vector<int>  appear_onestep(N,0);
  std::vector<long> ok_index(max_node_id + 1,0);//the index used when indexing all arrays whose length is N_new
  std::vector<long> degree_vector(N,-1);
  std::vector<long> degree_vector_onestep(N,-1);
  std::vector<long> n_tk_vector(K,0);
  std::vector<long> m_tk_vector(K,0);
  std::vector<long> z_j_vector(N_new,0);
  std::vector<long> offset_tk_vector(K,0);
  
  for (long i = 0; i < N; ++ i) {
      node_array[all_node(i)] = i; 
  }
  for (long i = 0; i < N_new; ++ i) {
      ok_array[ok_node(i)] = 1; 
      ok_index[ok_node(i)] = i;
     
  }

  long t  = 0;
  long edge_count = 0;

  while ((t < T)) {
      checkUserInterrupt();  
      if (t > 0)  {
          if (0 == only_PA) { 
              for (long j = 0; j < ok_node.size(); ++ j) {
                  if (degree_vector.at(node_array.at(ok_node(j))) >= 0) {
                      node_degree(t - 1,ok_index.at(ok_node(j))) = bin_vector(degree_vector.at(node_array.at(ok_node(j))));
                  }else {
                      node_degree(t - 1,ok_index.at(ok_node(j))) = -1;
                  }
              }
          }
          for (long k = 0; k < K; ++k) {
              n_tk(t - 1,k)          = n_tk_vector.at(k);
              if (0 == only_PA)  
                  offset_tk(t - 1,k) = offset_tk_vector.at(k);
          }
      }    
      for (long k = 0; k < K; ++k)
          m_tk_vector.at(k) = 0;
      if (0 == only_PA)      
          for (long j = 0; j < N_new; ++j)
              z_j_vector.at(j) = 0;    
      
  while ((edge_count < in_node.size()) && (time_stamp(edge_count) == time_vector(t))) {

      long in_node_ind  = node_array.at(in_node(edge_count)); 
      long out_node_ind = node_array.at(out_node(edge_count)); 
      // consider the in-node first
      //the node has not appeared in the previous time step
      if (0 == is_appear[in_node_ind]) {
          //the node has not appeared in previous edges of the current time step
          if (0 == appear_onestep.at(in_node_ind)) { 
              appear_onestep.at(in_node_ind)  = 1;
              degree_vector.at(in_node_ind)   = 1;
              ++n_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
              if ((0 == only_PA) && (0 == ok_array.at(in_node(edge_count))))   
                  ++offset_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
              
          }
          else { //the node has already appeared in some edges of the current time step
              if ((0 == only_PA)&& (0 == ok_array.at(in_node(edge_count))))    
                  --offset_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
              --n_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
              ++degree_vector.at(in_node_ind); 
              ++n_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
              if ((0 == only_PA) && (0 == ok_array.at(in_node(edge_count))))   
                  ++offset_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
          }     
    } // the node is already appeared in the previous time-step
      else { 
           if ((0 == only_PA) && (1 == ok_array.at(in_node(edge_count))))
                   z_j(ok_index.at(in_node(edge_count)))++;
             
           ++m_tk_vector.at(bin_vector(degree_vector_onestep.at(in_node_ind))); 
           --n_tk_vector.at(bin_vector(degree_vector.at(in_node_ind)));
           if ((0 == only_PA) && (0 == ok_array.at(in_node(edge_count)))) {                  
                   -- offset_tk_vector.at(bin_vector(degree_vector.at(in_node_ind)));        
           }
           ++degree_vector.at(in_node_ind);
           ++n_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
           if ((0 == only_PA) && (0 == ok_array.at(in_node(edge_count)))) {                       
                   ++offset_tk_vector.at(bin_vector(degree_vector.at(in_node_ind)));    
          }
      }
      //consider next the Out node
      //the node has not appeared in the previous time step
      if (0 == is_appear.at(out_node_ind)) {
          //the network is undirected, so this out_node is also counted
          if (1 == undirected) {
              if (0 == appear_onestep.at(out_node_ind)) {
                  degree_vector.at(out_node_ind)   = 1;  
                   ++n_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));  
                  if ((0 == only_PA) && (0 == ok_array.at(out_node(edge_count))))
                      ++offset_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));     
                 
              }
              else {
                  if ((0 == only_PA) && (0 == ok_array.at(out_node(edge_count)))) 
                      --offset_tk_vector.at(bin_vector(degree_vector.at(out_node_ind))); 
                  --n_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));    
                  degree_vector.at(out_node_ind)   = 1;   
                  ++n_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));   
                  if ((0 == only_PA) && (0 == ok_array.at(out_node(edge_count)))) 
                      ++offset_tk_vector.at(bin_vector(degree_vector.at(out_node_ind))); 
              }
          }
          //the network is directed, so this out_node is not counted
          else {
               if (0 == appear_onestep.at(out_node_ind)) { 
                  degree_vector.at(out_node_ind)   = 0; 
                  ++n_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));
                  if ((0 == only_PA) && (0 == ok_array.at(out_node(edge_count)))) {                       
                      ++offset_tk_vector.at(bin_vector(degree_vector.at(out_node_ind))); 
                  }
              }    
          }
          appear_onestep.at(out_node_ind)       = 1;
      } 
      // the node is already appeared in the previous time-step
      else {
           //the network is undirected, so this out_node is also counted  
          if (1 == undirected) {  
              if ((0 == only_PA) && (0 == ok_array.at(out_node(edge_count)))) {                  
                  -- offset_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));        
              }  
              if ((0 == only_PA) && (1 == ok_array.at(out_node(edge_count))))
                  z_j(ok_index.at(out_node(edge_count)))++;
              ++m_tk_vector.at(bin_vector(degree_vector_onestep.at(out_node_ind))); 
              
              --n_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));
              ++degree_vector.at(out_node_ind);
              ++n_tk_vector.at(bin_vector(degree_vector.at(out_node_ind))); 
              if ((0 == only_PA) && (0 == ok_array.at(out_node(edge_count)))) {                  
                  -- offset_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));        
              }
          }
      }
      ++edge_count; 
     }

  if (t > 0)  {
      for (long k = 0; k < K; ++k) {
          m_tk(t - 1,k)     = m_tk_vector.at(k);
          Sum_m_k(k)       += m_tk_vector.at(k);
          m_t(t - 1)       += m_tk_vector.at(k);
          
      }
      if (0 == only_PA)  
          for (long i = 0; i < N_new; ++i) {
              z_j(i) += z_j_vector.at(i);
          }    
  } 
     t++;
     for (long n = 0; n < (long) degree_vector.size(); ++n)
         degree_vector_onestep.at(n) = degree_vector.at(n);
     for (long i = 0; i < (long) is_appear.size(); ++i)
         is_appear.at(i) = appear_onestep.at(i); 
  }
  return 0;
}

// [[Rcpp::export(".update_f")]]
int update_f(      NumericVector& f, 
             const NumericVector& non_zero_f,
             const NumericMatrix& degree, 
             const NumericVector& theta, 
             const NumericVector& z_j,
             const NumericVector& normalized_const, 
             const NumericVector& m_t, 
             const double         shape, 
             const double         rate
             ) {
    long T        = degree.nrow();        // number of time-steps
    long N_nozero = non_zero_f.size();   // number of nodes
    #pragma omp parallel for
    for (long j = 0; j < N_nozero; j++) {
        double total = 0;
        for (long i = 0; i < T; i++)
            if ((degree(i,non_zero_f(j) - 1) >= 0) && (normalized_const(i) != 0)) {
                total += m_t(i) / normalized_const(i) * theta(degree(i,non_zero_f(j) - 1));
            }
        f(non_zero_f(j) - 1) = (z_j(non_zero_f(j) - 1) + shape - 1)/(total + rate);
    }
    return 0;
}

// [[Rcpp::export(".coeff_theta")]]
NumericVector coeff_theta( const NumericMatrix& degree,  
                           const NumericVector& f,
                           const NumericVector& normalized_const,  
                           const NumericVector& m_t, 
                           const int            length_theta
                           ) {
    int nrow = degree.nrow();
    int ncol = degree.ncol();
    NumericVector total(length_theta);
    //#pragma omp parallel 
    //{
    NumericVector total_temp(length_theta);
    for (int j = 0; j < length_theta; ++j) {
        total_temp[j] = 0;
    } 
    //#pragma omp for
    for (int j = 0; j < ncol; j++) {
        for (int i = 0; i < nrow; i++)
           if ((degree(i,j) >= 0) && (normalized_const(i) != 0)) {
                total[degree(i, j)] += f(j)* m_t(i) / normalized_const(i) ;
            }
    }
    //#pragma omp critical
    //{
    //for (int j = 0; j < length_theta; ++j)
    //    total[j] += total_temp[j];
    //}
   // }
    return total;
}

// [[Rcpp::export(".coeff_var")]]
NumericVector coeff_var(const NumericMatrix& degree,  
                        const NumericVector& f,
                        const NumericVector& normalized_const,
                        const NumericVector& m_t,
                        const NumericMatrix& offset, 
                        const int            length_theta) {
    int nrow = degree.nrow();
    int ncol = degree.ncol();
    NumericMatrix temp(nrow,length_theta);

    NumericVector total(length_theta);

    for (int j = 0; j < ncol; j++) {
        for (int t = 0; t < nrow; t++)
           if (degree(t,j) >= 0) {
                temp(t,degree(t,j)) += f(j);
            }
    }
    #pragma omp parallel for
    for (int k = 0; k < length_theta; k++) {
        for (int t = 0; t < nrow; t++)
            if (normalized_const[t] != 0)
                total(k) += pow(temp(t,k) + offset(t,k),2)* m_t(t) / pow(normalized_const(t),2);
    }
    return total;
}

// [[Rcpp::export(".cal_var_f")]]
int cal_var_f(      NumericVector& cov_f, 
              const NumericVector& non_zero_f,
              const NumericMatrix& degree,
              const NumericVector& theta,
              const NumericVector& f,
              const NumericVector& z_j,
              const NumericVector& normalized_const,
              const NumericVector& m_t, 
              const double         shape) {
    int T    = degree.nrow();
    int N    = non_zero_f.size();
    #pragma omp parallel for
    for (int j = 0; j < N; j++) {
        double total = 0;
        for (int i = 0; i < T; i++)
            if ((degree(i,non_zero_f(j) - 1) >= 0) && (normalized_const(i) != 0)) {
                total += m_t(i) / pow(normalized_const(i),2) * pow(theta(degree(i,non_zero_f(j) - 1)),2);
            }
          cov_f(j) = 1/(z_j(non_zero_f(j) - 1)/pow(f(non_zero_f(j) - 1),2) + //- total +
                     (shape - 1)*pow(f(non_zero_f(j) - 1),2));
    }
    return 0;
}
