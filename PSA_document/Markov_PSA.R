#---------------------------------------------------------------------------#
# A sample markov chain model
# 
# Based on a paper by Paul Tappenden, Susannah Sadler, Martin Wildman: 
# 'An Early Health Economic Analysis of the Potential Cost Effectiveness 
# of an Adherence Intervention to Improve Outcomes for Patients with Cystic Fibrosis'
# https://doi.org/10.1007/s40273-017-0500-x
#
# Excel model translated into R by 
# Paul Schneider, p.schneider@sheffield.ac.uk
# http://bitowaqr.github.io/
#
#---------------------------------------------------------------------------#

# load packages 
# install.packages("MCMCpack")
library(MCMCpack)
library(ggplot2)
# no scientific notations
options(scipen = 99)

# load functions
  source("./support_functions.R")

plot_plane <- function(int_costs_yearly) {
### Setup
  cycle_length = 84          # observational period - 84 years
  start_age = 16             # age at which agents enter the model - 16
  DR_QALY = DR_COSTS = 0.035 # discount rates for costs and QALYs - 3.5%
  PSA_length = 100           # how many runs for the probabilistic sensitivity analysis - 1000
  
  # create empty matrix to store results from PSA 
  # each row stores the main results of one iteration
  PSA_res = matrix(data=NA, ncol=5, nrow=PSA_length, 
                   dimnames = list(1:PSA_length,c("QALY base", "Costs base","QALY int","Costs int","Icer")))
  
## Draw parameters outside of loop and store as a list (of lists).
  # # ORIGNIAL RESULTS:
  #psa.input <- replicate(n=PSA_length,expr = draw_params(),simplify = F)                     
  
  # # OR, SIMULATE BORDERLINE COST-EFFECTIVENESS
  # # + increase the annual costs to the intervention
   psa.input <- replicate(n=PSA_length,expr = draw_params(int_costs_yearly = int_costs_yearly), simplify = F)                                

  
## RUN PSA LOOP 
  for(r in 1:PSA_length){
    # transistion matrices
    m       = psa.input[[r]] # take the PSA input for the r'th iteration 
    mat_b   = trans_mat(m$tb1,m$tb2,m$tb3,m$tb4) # baseline transition matrix 
    mat_int = trans_mat(m$ti1,m$ti2,m$ti3,m$ti4) # intervention transition matrix
    
    # create empty matrix traces
    markov_trace_base = matrix(data=NA,ncol=dim(mat_b)[1]+1,nrow=1+cycle_length)
    markov_trace_base[1,] = c(m$initial_distribution,0)
    markov_trace_int = markov_trace_base
    aging = 0
    
    ### INNER MARKOV LOOP
      for(i in 1:cycle_length){
        # adjust transisiton probs for age-dependent death risk
        step_matrix_b = death_age_fmat(mat=mat_b,x=start_age+aging,age = m$death_prob$age,value = m$death_prob$value)
        step_matrix_int = death_age_fmat(mat=mat_int,x=start_age+aging,age = m$death_prob$age,value = m$death_prob$value)
        # Move one step in the markov chain
        markov_trace_base[i+1,]= markov_trace_base[i,] %*% step_matrix_b
        markov_trace_int[i+1,]= markov_trace_int[i,] %*% step_matrix_int
        # increase age
        aging = aging + 1
        }
    
    ### EXTRACT QALY+COST RESULTS FROM MARKOV LOOP
    # baseline
    res_base  = get_CE(trace=markov_trace_base,
                       utils = m$state_utils,
                       costs = m$base_state_costs,
                       iv_rates = m$base_iv_days,
                       iv_disutil = m$iv_excer_disutil,
                       DRQ = DR_QALY,
                       DRC = DR_COSTS
                       )
    # intervention 
    res_int  = get_CE(trace=markov_trace_int,
                       utils = m$state_utils,
                       costs = m$int_state_costs,
                       extra_costs = m$int_costs_once,
                       iv_rates = m$int_iv_days,
                       iv_disutil = m$iv_excer_disutil,
                       DRQ = DR_QALY,
                       DRC = DR_COSTS
                       )
    # calculate ICER
    Icer = (res_base$state_costs_disc-res_int$state_costs_disc)/ (res_base$QALY_disc - res_int$QALY_disc)
    
    ## store results in PSA trace
    PSA_res[r,] = c(res_base$QALY_disc,
                    res_base$state_costs_disc,
                    res_int$QALY_disc,
                    res_int$state_costs_disc,
                    Icer)
  }

####################
### PSA RESULTS
###########
    
  # CE-PLANE PLOT
  incr_Q = PSA_res[,"QALY int"] - PSA_res[,"QALY base"] 
  incr_C = PSA_res[,"Costs int"] - PSA_res[,"Costs base"] 
    
    ce_plane_plot = ggplot() +
      # plot iteration results and mean estimate
      geom_point(aes(x=incr_Q,y=incr_C),col="cadetblue4",alpha=0.8,size=0.8) +
      geom_point(aes(x=mean(incr_Q),y=mean(incr_C)),size=2,col="blue") +
      # coordinate systems
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      geom_line(aes(x=c(-1,1),y=c(-20000,20000))) +
      # symmetric axes
      coord_cartesian(xlim=c(-max(abs(incr_Q)),max(abs(incr_Q))), 
                      ylim=c(-max(abs(incr_C)),max(abs(incr_C))))  +
      # labels
      xlab("Incremental QALY") +
      ylab("Incremental Costs") +
      theme_minimal()  + 
      NULL
    
    ce_plane_plot
}
  
  ## CEAC PLOT
    #ceac_data = get_ceac()
    #ceac_plot = ggplot(ceac_data) +
      #geom_line(aes(x=wtp,y=prob_ce,col="Intervention")) +
      #geom_line(aes(x=wtp,y=1-prob_ce,col="Baseline")) +
      #ylim(c(0,1)) +
      #ylab("Prob. cost effective") +
      #xlab("Willingness to pay threshold") +
      #theme_minimal()+
      #ggtitle("Cost acceptibility curve") +
      #theme(legend.position = "bottom") +
      #scale_color_manual(name="",values = c("blue","magenta"))
    ## ceac_plot
    
  ### PSA results summary
    #res_summary = data.frame(QALY_base = mean(PSA_res[,"QALY base"] ),
                             #QALY_int = mean(PSA_res[,"QALY int"]),
                             #incr_Qaly = mean(PSA_res[,"QALY int"] - PSA_res[,"QALY base"] ),
                             #cost_base = mean(PSA_res[,"Costs base"]),
                             #cost_int = mean(PSA_res[,"Costs int"]),
                             #incr_costs = mean(PSA_res[,"Costs int"] - PSA_res[,"Costs base"]),
                             #netbenefit = mean(incr_Q*20000 - incr_C),
                             #row.names = "Point estimate")
    #t(res_summary)

## Done.
