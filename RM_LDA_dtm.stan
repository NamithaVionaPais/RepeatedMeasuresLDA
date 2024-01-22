data {
  int<lower=2> J;                           // num topics
  int<lower=2> V;                           // Word levels
  int<lower=1> M;                           // num documents
  int<lower=1> n_d;                         // number of sub-documents in each document
  int w[M,n_d,V];          // word n or product n
  vector<lower=0>[J] mu;                 // mu
  int num_sub[M];
  
}
parameters {
  simplex[V] beta[J];       // word dist for topic j
  corr_matrix[J] sigma[M];  // doc specific var-cov matrix
  vector[J] lalpha[M,n_d];
  simplex[J] theta[M,n_d];  
}


model {
  for(d in 1:M){
    for(r in 1:num_sub[d]){
      lalpha[d,r] ~ multi_normal(mu,sigma[d]);
      theta[d, r] ~ dirichlet(exp(lalpha[d,r]));
      for (v in 1:V) {
        real gamma[J];
        for (j in 1:J)
        {
          gamma[j] =log(theta[d,r,j]) + w[d,r,v]*log(beta[j, v]);
        }
        target += log_sum_exp(gamma);  
      }
    }
  }
}


