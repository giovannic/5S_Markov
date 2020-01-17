# # functions

# a function that draws all the model parameters, IF those 
# are not explicitly provided, and returns them as a list
draw_params = function(death_prob=NULL,
                      initial_distribution=NULL,
                      cons = NULL,
                      gamma = NULL,
                      transplant_rate = NULL,
                      tb1 =NULL,
                      tb2 = NULL,
                      tb3 = NULL,
                      tb4 = NULL,
                      ti1 =NULL,
                      ti2 = NULL,
                      ti3 = NULL,
                      ti4 = NULL,
                      base_iv_days = NULL,
                      int_iv_days = NULL,
                      s12_dis=NULL,
                      s22_dis = NULL,
                      s1_util=NULL,
                      s4_utl=NULL,
                      iv_excer_disutil=NULL,
                      s2_util = NULL,
                      s3_util = NULL,
                      state_utils = NULL,
                      int_costs_once =NULL,
                      int_costs_yearly=NULL,
                      band_costs=NULL,
                      transplant_costs=NULL,
                      b1=NULL,
                      b2=NULL,
                      b3=NULL,
                      drug_costs=NULL,
                      prop_drug_use=NULL,
                      drug_use_costs=NULL,
                      iv_drug_cost=NULL,
                      iv_day_costs=NULL,
                      iv_prop_hosp_days = NULL
                      ){
  
  # age-dependent survival function / death probabilities
  if(is.null(cons)){
    cons = rnorm(1,0.00446820,0.00044682)
  }
  if(is.null(gamma)){
    gamma = rnorm(1,0.06082920,0.00608292)
  }
  if(is.null(death_prob)){
    age_range = 0:100
    surv_probs = exp(-cons*(gamma^-1)*(exp(gamma*age_range)-1))
    cond_prob = rep(NA,100)
    cond_prob[1] = 1
    for(s in 2:length(surv_probs)){
      cond_prob[s] = surv_probs[s]/surv_probs[s-1]
    }
    death_prob = data.frame(age=age_range ,value=1-cond_prob) 
  }
  
  #  cohort starting distributon
  if(is.null(initial_distribution)){initial_distribution = as.numeric(rdirichlet(1,c(2892,1966,860,0)))  }
  
  # draw baseline state transistion probabilities
  # draw transplant probability  (s3 -> s4)
  if(is.null(transplant_rate)){
    transplant_rate = -1
    while(transplant_rate<0){
      transplant_rate = rnorm(n=1,mean=0.00432827833974037,sd=0.00432827833974037*0.1)
    }
  }
  # baseline transistion from s1
  if(is.null(tb1)){tb1 = as.numeric(rdirichlet(1,c(6057.33,921.69,22.98,0)))}
  # baseline transistion from s2
  if(is.null(tb2)){tb2 = as.numeric(rdirichlet(1,c(503.08,2861.77,381.15,0)))}
  # baseline transistion from s3 with transplant rate
  # will be problemativ to fix...
  if(is.null(tb3)){
    temp = c(36.37,186.16,1130.47)
    temp = as.numeric(rdirichlet(1,temp))
    temp = c(temp,transplant_rate*temp[3])
    temp[3] = (1-transplant_rate)*temp[3]
    tb3 = temp
  }
  # baseline transistions from s4
  if(is.null(tb4)){tb4 = c(0,0,0,1) }
  # Draw intervention transisiton probs
  # int transistion from s1
  if(is.null(ti1)){ti1 = as.numeric(rdirichlet(1,c(6101.83,878.38,21.79,0)))}
  # int transistion from s2
  if(is.null(ti2)){ti2 = as.numeric(rdirichlet(1,c(546.51,2850.04,349.45,0))) }
  # int transistion from s3 with transplant rate
  # will be problemativ to fix...
  if(is.null(ti3)){
    temp = c(35.65,183.04,1134.30)
    temp = c(temp,transplant_rate*temp[3])
    temp[3] = (1-transplant_rate)*temp[3]
    ti3 = as.numeric(rdirichlet(1,temp))
  }
  # int transistions from s4
  if(is.null(ti4)){ti4 = c(0,0,0,1)}
  # IV days         
  if(is.null(base_iv_days)){
    temp1 = rbeta(1,1802.22634780,42817.93952816)
    temp2 = rbeta(1,3134.65873055,32608.88484427)
    temp3 = rbeta(1,1184.12863848,6181.86183954)
    base_iv_days = c(temp1,temp2,temp3,0,0)
  }
  if(is.null(int_iv_days)){
    temp.rr = exp(rnorm(1,log(0.44658867)-(0.5*0.08928571^2),0.08928571) )
    # temp.rr = rnorm(1,1,1)
    int_iv_days = base_iv_days * temp.rr
  }   
  
  # UTILITIES        
  if(is.null(s12_dis)){s12_dis = rbeta(1,4.11028060769,72.00602694210)}
  if(is.null(s22_dis)){s22_dis = rbeta(1,12.59854497777,61.51054312676)}
  if(is.null(s1_util)){s1_util = rbeta(1,108.52029880066,17.08189888529)}
  if(is.null(s4_utl)){s4_utl = rbeta(1,319.30588235,65.40000000)}
  if(is.null(iv_excer_disutil)){iv_excer_disutil = -rbeta(1,3.4821053999,16.5299945993)}
  s2_util = s1_util - s12_dis
  s3_util = s2_util - s22_dis
  state_utils = c(s1_util,s2_util,s3_util,s4_utl,0)
  
  # DRUG USE
  # iv proportion hospital days
  if(is.null(iv_prop_hosp_days)){iv_prop_hosp_days = rbeta(1,93455.00,78452.00)}
  
  # band distributions
  if(is.null(b1)){b1 = as.numeric(rdirichlet(1,c(886,83,1063,1497,799,75,15))) }
  if(is.null(b2)){b2 = as.numeric(rdirichlet(1,c(115,28,250,807,787,247,56)))}
  if(is.null(b3)){b3 = as.numeric(rdirichlet(1,c(20,6,60,234,310,225,92)))}
  if(is.null(prop_drug_use)){ 
    prop_drug_use = c(rbeta(1,3949.00,2660.00),
                      rbeta(1,779.00,5830.00),
                      rbeta(1,235.00,6374.00),
                      rbeta(1,0.00,6609.00),
                      rbeta(1,2555.00,4054.00),
                      rbeta(1,1602.00,5007.00),
                      rbeta(1,0.00,6609.00),
                      rbeta(1,779.00,5830.00))}
  # COSTS
  # adherence intervention costs
  if(is.null(int_costs_once)){int_costs_once = c(583.44)}
  if(is.null(int_costs_yearly)){int_costs_yearly= c(583.44)}
  
  # band costs
  if(is.null(band_costs)){band_costs = c(5033.00,7447.00,7447.00,12036.00,18422.00,33224.00,40054.00)}
  # transplant costs
  if(is.null(transplant_costs)){transplant_costs = rnorm(1,40000.00,4000.00)}
  if(is.null(drug_costs)){drug_costs = c(6044.04,6016.25,14228.64,0.00,2366.82,8181.60,0.00,11674.96)} # fixed
  # state band costs
  s1_band_costs = sum(b1 * band_costs)
  s2_band_costs = sum(b2 * band_costs)
  s3_band_costs = sum(b3 * band_costs)
  # drug costs
  drug_use_costs = sum(prop_drug_use*drug_costs)
  if(is.null(iv_drug_cost)){iv_drug_cost = c(71.99)} # fixed
  if(is.null(iv_day_costs)){iv_day_costs = rnorm(1,361.68,77.48)}
  
  iv_costs =  iv_drug_cost + (iv_prop_hosp_days * (iv_day_costs + iv_drug_cost))
  
  
  ### full state costs
  base_state_costs = c(s1_base_costs = s1_band_costs +drug_use_costs +(base_iv_days[1]   *365.25 * iv_costs),
                       s2_base_costs = s2_band_costs +drug_use_costs +(base_iv_days[2]   *365.25 * iv_costs),
                       s3_base_costs = s3_band_costs + drug_use_costs+ transplant_rate * transplant_costs + (base_iv_days[3]   *365.25 * iv_costs),
                       s4_costs = 0,
                       s5_costs = 0)
  
  int_state_costs = c(s1_int_costs = s1_band_costs +drug_use_costs +(int_iv_days[1]   *365.25 * iv_costs) + int_costs_yearly,
                      s2_int_costs = s2_band_costs +drug_use_costs +(int_iv_days[2]   *365.25 * iv_costs) + int_costs_yearly,
                      s3_int_costs = s3_band_costs + drug_use_costs+ transplant_rate * transplant_costs + (int_iv_days[3]   *365.25 * iv_costs) + int_costs_yearly,
                      s4_costs = 0+ int_costs_yearly,
                      s5_costs = 0)
  
  # combine and return results
  res_list = list(death_prob = death_prob,
                  age = death_prob$age,
                  value = death_prob$value,
                  initial_distribution = initial_distribution,
                  transplant_rate = transplant_rate,
                  tb1 = tb1,
                  cons = cons,
                  gamma = gamma,
                  tb2 = tb2,
                  tb3 = tb3,
                  tb4 = tb4,
                  ti1 = ti1,
                  ti2 = ti2,
                  ti3 = ti3,
                  ti4 = ti4,
                  base_iv_days = base_iv_days,
                  int_iv_days = int_iv_days,
                  s12_dis = s12_dis,
                  s22_dis = s22_dis,
                  s1_util = s1_util,
                  s4_utl = s4_utl,
                  iv_excer_disutil = iv_excer_disutil,
                  s2_util = s2_util,
                  s3_util = s3_util,
                  state_utils = state_utils,
                  int_costs_once = int_costs_once,
                  int_costs_yearly = int_costs_yearly,
                  band_costs = band_costs,
                  transplant_costs = transplant_costs,
                  b1 = b1,
                  b2 = b2,
                  b3 = b3,
                  drug_costs = drug_costs,
                  prop_drug_use = prop_drug_use,
                  drug_use_costs = drug_use_costs,
                  iv_drug_cost = iv_drug_cost,
                  iv_day_costs = iv_day_costs,
                  iv_prop_hosp_days = iv_prop_hosp_days,
                  iv_costs = iv_costs,
                  s1_band_costs = s1_band_costs,
                  s2_band_costs = s2_band_costs,
                  s3_band_costs = s3_band_costs,
                  base_state_costs = base_state_costs,
                  int_state_costs = int_state_costs
                  
  )
  return(res_list)
}

# draw_params()

# adjust transition mat for age-death rate
death_age_fmat = function(mat,x = 16,age = m$death_prob$age,value = m$death_prob$value){
  d_rate = value[age==x]
  mat = apply(t(mat),1,function(x) x * (1-d_rate))
  mat = cbind(mat,d_rate)
  mat = rbind(mat,c(0,0,0,0,1))
  colnames(mat)[5] = "death"
  rownames(mat)[5] = "death"
  return(mat)
}


# function to create a named transistion matrix
trans_mat = function(r1,r2,r3,r4,names = NULL){
  if(is.null(names)){names=c(">70% FEV","69-40% FEV","<40% FEV","Post trans")}
  mat = matrix(ncol=4,data=c(r1,r2,r3,r4),dimnames = list(names,names))
  mat = t(mat)
  return(mat)
}


# extract discounted cost+qalys from markov trace
get_CE = function(
  trace = markov_trace_int, 
  utils = state_utils,
  costs= int_state_costs,
  extra_costs = 0, #
  iv_rates=int_iv_days,
  iv_disutil = excer_disutil,
  DRQ = DR_QALY, DRC= DR_COSTS,
  return.sum = T){
  
  # half cycle correction
  trace_hcc = matrix(nrow=dim(trace)[1]-1,ncol=dim(trace)[2],data=NA)
  for(j in 2:dim(trace )[1]){
    trace_hcc[j-1,] = (trace[j-1,]+trace[j,])/2
  }
  
  disc_year = c(1:dim(trace_hcc)[1])
  
  state_costs_undisc = rowSums(t(t(trace_hcc) * costs))
  state_costs_disc = state_costs_undisc / (1+DRC)^disc_year
  
  iv_days  = rowSums(t(iv_rates*t(trace_hcc)))
  iv_QALY_undisc = iv_days *  iv_disutil 
  state_QALY_undisc  = rowSums(t(utils*t(trace_hcc)))
  QALY_undisc = state_QALY_undisc +iv_QALY_undisc
  QALY_disc = QALY_undisc / (1+DRQ)^disc_year
  
  
  if(return.sum){
    res = list(QALY_disc=sum(QALY_disc),
               state_costs_disc=sum(state_costs_disc)+extra_costs)
  } else {
    res = list(QALY_disc=(QALY_disc),
               state_costs_disc=(state_costs_disc)) 
  }
  return(res)
}


# compute Cost-effectiveness Acceptability Curve
get_ceac = function(iq = incr_Q, ic =incr_C,lambda = seq(0,100000,1000)){
  len = length(lambda)
  nb = rep(NA,len)
  for(l in 1:len){
    nb[l] = sum((incr_Q * lambda[l] -incr_C)>0) / length(incr_C)  
  }
  res = data.frame(prob_ce = nb,wtp = lambda)
  return(res)
}


