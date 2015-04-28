//Integrated Model for Experiment 1 Data
data {
  int<lower=0> Trials;
  int<lower=0> PicConds;
  int<lower=0> IntConds;
  real<lower=1> NumPic[PicConds];
  int<lower=0> Int[IntConds];
  int<lower=0> P_s_s[PicConds*IntConds]; //Num Participants for Same
  int<lower=0> P_s_o[PicConds*IntConds]; //Num Participants for Switch
  int<lower=0> P_n_s[PicConds*IntConds]; //Num Participants for Same
  int<lower=0> P_n_o[PicConds*IntConds]; //Num Participants for Switch
  int<lower=0> B_s_s[PicConds,IntConds]; //Condition start indices for Same
  int<lower=0> B_s_o[PicConds,IntConds]; //Condition start indices for Switch
  int<lower=0> B_n_s[PicConds,IntConds]; //Condition start indices for Same
  int<lower=0> B_n_o[PicConds,IntConds]; //Condition start indices for Switch
  int<lower=0> E_s_s[PicConds,IntConds]; //Condition end indices for Same
  int<lower=0> E_s_o[PicConds,IntConds]; //Condition end indices for Switch
  int<lower=0> E_n_s[PicConds,IntConds]; //Condition end indices for Same
  int<lower=0> E_n_o[PicConds,IntConds]; //Condition end indices for Switch
  int<lower=0,upper=Trials> S_S[sum(P_s_s)]; //prop. correct for Same
  int<lower=0,upper=Trials> S_O[sum(P_s_o)]; //prop. correct for Switch
  int<lower=0,upper=Trials> N_S[sum(P_n_s)]; //prop. correct for Same
  int<lower=0,upper=Trials> N_O[sum(P_n_o)]; //prop. correct for Switch
}

parameters {
  real<lower=0,upper=1> sigma_s; //Social Attention allocation parameter
  real<lower=0,upper=1> sigma_n; //No-social Attention allocation parameter
  real<lower=0> lambda; //Memory decay parameter
  real<lower=0> gamma_s; //Initial encoding strength parameter
  real<lower=0> gamma_n; //Initial encoding strength parameter
}

transformed parameters {
  //Compute strength of encoding in each NumRef/Int condition 
  real<lower=0>c_s_s[PicConds,IntConds];
  real<lower=0>c_s_o[PicConds,IntConds];
  real<lower=0>c_n_s[PicConds,IntConds];
  real<lower=0>c_n_o[PicConds,IntConds];
  
  for (n in 1:PicConds)
    for (i in 1:IntConds) {
      c_s_s[n,i] <- sigma_s*pow(Int[i],-lambda)*gamma_s; //Same
      c_s_o[n,i] <- ((1-sigma_s)/(NumPic[n]-1))*pow(Int[i],-lambda)*gamma_s; //Switch
      
      c_n_s[n,i] <- sigma_n*pow(Int[i],-lambda)*gamma_n; //Same
      c_n_o[n,i] <- ((1-sigma_n)/(NumPic[n]-1))*pow(Int[i],-lambda)*gamma_n; //Switch
      
      //Prop. Correct must be <= 1
      if(c_s_s[n,i] > 1)
        c_s_s[n,i] <- 1;
      if(c_s_o[n,i] > 1)
        c_s_o[n,i] <- 1;
      if(c_n_s[n,i] > 1)
        c_n_s[n,i] <- 1;
      if(c_n_o[n,i] > 1)
        c_n_o[n,i] <- 1;
    }
}

model {
  
  for (n in 1:PicConds) {
    for (i in 1:IntConds) {
      //Compute binomial probabilities of success given encoding strengths
      for(s in B_s_s[n,i]:E_s_s[n,i]){
        S_S[s] ~ binomial(Trials, c_s_s[n,i] + (1-c_s_s[n,i])/NumPic[n]);
      }
      for(o in B_s_o[n,i]:E_s_o[n,i]){
        S_O[o] ~ binomial(Trials, c_s_o[n,i] + (1-c_s_o[n,i])/NumPic[n]);
      }
      for(s in B_n_s[n,i]:E_n_s[n,i]){
        N_S[s] ~ binomial(Trials, c_n_s[n,i] + (1-c_n_s[n,i])/NumPic[n]);
      }
      for(o in B_n_o[n,i]:E_n_o[n,i]){
        N_O[o] ~ binomial(Trials, c_n_o[n,i] + (1-c_n_o[n,i])/NumPic[n]);
      }
    }
  }
}