
data{
  int<lower = 1> M; // number of years
  int<lower = 1> Mgen; // number of years of posterior data to generate
  int<lower = 1> Nt; // number of days per year
  int<lower = 1> S; // number of sites
  real<lower=0> y[S,M,Nt]; // all rainfall values
  int<lower=0> N[S,M];  // number of events per year with PRCP>0 for site S
  matrix[S,3] X; //(1,lat,long)


  // iper-parametri
  real<lower = 0> sg0prior[2]; // a, b of inverse gamma for gamma
  real<lower = 0> sd0prior[2]; // a, b of inverse gamma for delta
  vector[3] ma0prior; // mean of the normal for gamma
  vector<lower = 0>[3] sa0prior; // sd of the normal for gamma 
  vector[3] mn0prior; //mean of the normal for delta
  vector<lower = 0>[3] sn0prior; // sd of the normal for delta
  vector[3] mb0prior; // mean of the normal for lambda
  vector<lower = 0>[3] sb0prior; // sd of the normal for lambda
}

parameters{
  vector<lower=0>[M] gamma[S]; //gamma_j(s) (shape)
  vector<lower=0>[M] delta[S]; //delta_j(s) (scale)
  
  real<lower=0> sg; // sigma_gamma
  real<lower=0> sdl; // sigma_delta

  vector[3] ag; // regression parameters for gamma
  vector[3] nd; // regression parameters for delta
  vector[3] bl; // regression parameters for lambda
}

transformed parameters{
  vector[S] mg; // mu_gamma
  vector[S] md; // mu_delta
  vector<lower=0, upper = 1>[S] lambda; //lamba(s)
  lambda = inv_logit(X*bl);
  mg = X*ag;
  md = X*nd;
}

model{
  //prior
  bl[1] ~ normal(mb0prior[1], sb0prior[1]);
  bl[2] ~ normal(mb0prior[2], sb0prior[2]);
  bl[3] ~ normal(mb0prior[3], sb0prior[3]);
  ag[1] ~ normal(ma0prior[1], sa0prior[1]);
  ag[2] ~ normal(ma0prior[2], sa0prior[2]);
  ag[3] ~ normal(ma0prior[3], sa0prior[3]);
  nd[1] ~ normal(mn0prior[1], sn0prior[1]);
  nd[2] ~ normal(mn0prior[2], sn0prior[2]);
  nd[3] ~ normal(mn0prior[3], sn0prior[3]);
  sg ~ inv_gamma(sg0prior[1], sg0prior[2]);
  sdl ~ inv_gamma(sd0prior[1], sd0prior[2]);

  //model per theta_j(s)
  for(s in 1:S){ 
   gamma[s] ~ gumbel(mg[s], sg);
   delta[s] ~ gumbel(md[s], sdl);
  }  

  //model for the number of events/year
  for (s in 1:S) {
    for (m in 1:M) {
    	target += binomial_lpmf(N[s,m] | Nt, lambda[s]);
      		for(j in 1:Nt) {
        		if (y[s,m,j] >1e-6) {
             	  y[s,m,j] ~ weibull( gamma[s,m], delta[s,m]);
      			}
    		}
  		}
  	}

}

generated quantities {
  int<lower=0> Ngen[S,Mgen]; // all observed values
  real<lower=0> ggen[S,Mgen]; // generated yearly parameters
  real<lower=0> dgen[S,Mgen]; // generated yearly parameters
  int count_d;
  real dtmp;
  int count_g;
  real gtmp;

  for (s in 1:S) {
   count_d = 0;
    while (count_d < Mgen) {
      dtmp =  gumbel_rng(md[s], sdl);
      if (dtmp > 0){
        count_d = count_d + 1;
        dgen[s,count_d] = dtmp;
      }
    }

    count_g = 0;
    while (count_g < Mgen){
      gtmp =  gumbel_rng(mg[s], sg);
      if (gtmp > 0){
        count_g = count_g + 1;
        ggen[s,count_g] = gtmp;
      }
    }
  }

for (s in 1:S) {
  for (m in 1:Mgen) {
    Ngen[s,m] = binomial_rng(Nt, lambda[s]);
   }
  }    
    
}
