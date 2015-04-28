//Single Referent Tracking Model for Both Experiments
data {
  int<lower=0> Trials;
  int<lower=0> PicConds;
  int<lower=0> IntConds;
  real<lower=1> NumPic[PicConds];
  int<lower=0> Int[IntConds];
  int<lower=0> P_s[PicConds*IntConds]; //Num Participants for Same
  int<lower=0> P_o[PicConds*IntConds]; //Num Participants for Switch
  int<lower=0> P_w[PicConds*IntConds]; //Num Participants for New Label
  int<lower=0> B_s[PicConds,IntConds]; //Condition start indices for Same
  int<lower=0> B_o[PicConds,IntConds]; //Condition start indices for Switch
  int<lower=0> B_w[PicConds,IntConds]; //Condition start indices for New Label
  int<lower=0> E_s[PicConds,IntConds]; //Condition end indices for Same
  int<lower=0> E_o[PicConds,IntConds]; //Condition end indices for Switch
  int<lower=0> E_w[PicConds,IntConds]; //Condition end indices for New Label
  int<lower=0,upper=Trials> S[sum(P_s)]; //prop. correct for Same
  int<lower=0,upper=Trials> O[sum(P_o)]; //prop. correct for Switch
  int<lower=0,upper=Trials> W[sum(P_w)]; //prop. correct for New Label
}

parameters {
  real<lower=0> lambda; //Memory decay parameter
  real<lower=0> gamma; //Initial encoding strength parameter
}

transformed parameters {
  //Compute strength of encoding in each NumRef/Int condition 
  real<lower=0,upper=1>c_s[PicConds,IntConds];
  real<lower=0,upper=1>c_o[PicConds,IntConds];
  real<lower=0,upper=1>c_w[PicConds,IntConds];
  
  for (n in 1:PicConds)
    for (i in 1:IntConds) {
      c_s[n,i] <- pow(Int[i],-lambda)*gamma; //Same
      c_o[n,i] <- 0; //Switch: No memory for non-hypothesized referents
      c_w[n,i] <- 0; //New Label: No memory for non-hypothesized referents
    }
}

model {
  
  for (n in 1:PicConds) {
    for (i in 1:IntConds) {
      //Compute binomial probabilities of success given encoding strengths
      for(s in B_s[n,i]:E_s[n,i]){
        S[s] ~ binomial(Trials, c_s[n,i] + (1-c_s[n,i])/NumPic[n]);
      }
      for(o in B_o[n,i]:E_o[n,i]){
        O[o] ~ binomial(Trials, c_o[n,i] + (1-c_o[n,i])/NumPic[n]);
      }
      for(w in B_w[n,i]:E_w[n,i]){
        W[w] ~ binomial(Trials,(1-c_w[n,i])/NumPic[n] - c_w[n,i]);
      }
    }
  }
}