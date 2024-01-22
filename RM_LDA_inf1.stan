data {
  int<lower=2> J;                           // num topics
  int<lower=2> V;                           // Word levels
  int<lower=1> M;                           // num documents
  int<lower=1> N;                           // total word instances
  int<lower=1> n_d;                         // number of sub-documents in each document
  int<lower=1,upper=V> w[M,n_d,N];          // word n or product n
  vector<lower=0>[J] mu;                 // mu
  simplex[V] beta[J];       // word dist for topic j
  corr_matrix[J] sigma[M];  // doc specific var-cov matrix
  
}
parameters {
  simplex[J] theta[M, n_d];        // theta is a simplex vector of size J
  vector[J] lalpha[M,n_d];
}


model {
  for(d in 1:M){
    for(r in 1:n_d){
      lalpha[d,r] ~ multi_normal(mu,sigma[d]);
      theta[d, r] ~ dirichlet(exp(lalpha[d,r]));
      for (n in 1:N) {
        real gamma[J];
        for (j in 1:J)
        {
          gamma[j] =log(theta[d,r,j]) + log(beta[j, w[d,r,n]]);
        }
        target += log_sum_exp(gamma);  
      }
    }
  }
}


