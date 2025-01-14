
library(ggplot2)
library(gridExtra)
library(psych)
library(MetBrewer)
library(reshape2)
library(plyr)
library(patchwork)
library(RColorBrewer)


############## ############## ############################ ############## ############################ ############## ##############
#################### ############## ############## I. Aggregations ############################ ### ############## ##############
############## ############## ############################ ############## ############################ ############## ##############
n_Sex = 2
n_Age = 8
n_Time = 10
n_subsize = 100

# Case I (not depending in Y)
DEPEND_ON_Y = FALSE

# Case II (depending in Y)
DEPEND_ON_Y = TRUE

agg_overall_BL = agg_overall_Fixed = agg_overall_Indep = agg_overall_RW = rep(NA, 50)
agg_time_BL = agg_time_Fixed = agg_time_Indep = agg_time_RW = matrix(NA, nrow = 50, ncol = n_Time)
agg_time_RW_noAge = agg_time_Indep_noAge = agg_time_Fixed_noAge = agg_time_BL_noAge = matrix(NA, nrow = 50, ncol = n_Time)
agg_overall_PopLevel = rep(NA, 50)


for(i in 1:50){
  name1 <- paste0("../simulation_fitted_model_itr_", i, "_age_",  n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y, ".rda")
  name2 <- paste0("../simulation_fitted_model_noAge_itr_", i, "_age_",  n_Age,"_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y,".rda")
  name3 <- paste0("../simulation_fitted_model_PopLevel_itr_", i, "_age_",  n_Age,"_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y,".rda")
 if(file.exists(name1)){
    print(i)
    load(name1) 
    agg_overall_param  = mean(models_fit$sim_data$param$pi)
    agg_overall_BL[i] = mean(apply(models_fit$BL$poster_sample$sample_pis, 1, mean))-agg_overall_param
    agg_overall_Fixed[i] = mean(apply(models_fit$Fixed$poster_sample$sample_pis, 1, mean))-agg_overall_param
    agg_overall_Indep[i] = mean(apply(models_fit$Indep$poster_sample$sample_pis, 1, mean))-agg_overall_param
    agg_overall_RW[i] = mean(apply(models_fit$RW$poster_sample$sample_pis, 1, mean))-agg_overall_param
    
    pi_dim = array(NA, dim = c(n_Sex, n_Time, n_Age))
    pi_BL_dim = pi_Indep_dim = pi_Fixed_dim = pi_RW_dim = array(NA, dim = c(n_Sex, n_Time, n_Age))
    for (nSex in c(0, 1)){
      for(nTime in 0:(n_Time-1)){
        for(nAge in 1:n_Age){
          pi_dim[nSex+1, nTime+1, nAge] = models_fit$sim_data$param$pi[nSex * n_Time * n_Age + nTime * n_Age + nAge]
          pi_BL_dim[nSex+1, nTime+1, nAge] = apply(models_fit$BL$poster_sample$sample_pis, 2, mean)[nSex * n_Time * n_Age + nTime * n_Age + nAge]
          pi_Indep_dim[nSex+1, nTime+1, nAge] = apply(models_fit$Indep$poster_sample$sample_pis, 2, mean)[nSex * n_Time * n_Age + nTime * n_Age + nAge]
          pi_Fixed_dim[nSex+1, nTime+1, nAge] = apply(models_fit$Fixed$poster_sample$sample_pis, 2, mean)[nSex * n_Time * n_Age + nTime * n_Age + nAge]
          pi_RW_dim[nSex+1, nTime+1, nAge] = apply(models_fit$RW$poster_sample$sample_pis, 2, mean)[nSex * n_Time * n_Age + nTime * n_Age + nAge]
        }
      }
    }
    agg_time_param  = apply(pi_dim, 2, mean)
    agg_time_BL[i, ] = apply(pi_BL_dim, 2, mean)-agg_time_param
    agg_time_Fixed[i, ] =  apply(pi_Fixed_dim, 2, mean)-agg_time_param
    agg_time_Indep[i, ] =  apply(pi_Indep_dim, 2, mean)-agg_time_param
    agg_time_RW[i, ] =  apply(pi_RW_dim, 2, mean)-agg_time_param
  }
  
  if(file.exists(name2)){
    # print(i)
    load(name2)
    agg_time_RW_noAge[i, ] = apply(models_fit_noAge$RW$poster_sample$sample_pis, 2, mean) - agg_time_param
    agg_time_Indep_noAge[i, ] = apply(models_fit_noAge$Indep$poster_sample$sample_pis, 2, mean) - agg_time_param
    agg_time_Fixed_noAge[i, ] = apply(models_fit_noAge$Fixed$poster_sample$sample_pis, 2, mean) - agg_time_param
    agg_time_BL_noAge[i, ] = apply(models_fit_noAge$BL$poster_sample$sample_pis, 2, mean) - agg_time_param
    
  }
  if(file.exists(name3)){
    # print(i)
    load(name3)
    agg_overall_PopLevel[i] = mean(model_fit_PopLevel$poster_sample$sample_pi) - agg_overall_param
  }
}

# save to local
agg_time_bias = list()
agg_time_bias$RW = agg_time_RW
agg_time_bias$Indep = agg_time_Indep
agg_time_bias$Fixed = agg_time_Fixed
agg_time_bias$BL = agg_time_BL
agg_time_bias$agg_time_param = agg_time_param
save(agg_time_bias, file = paste0("../simulation_fitted_model_agg_time_bias","_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y, ".rda"))

models_fit_noAge_bias = list()
models_fit_noAge_bias$RW = agg_time_RW_noAge
models_fit_noAge_bias$Indep = agg_time_Indep_noAge
models_fit_noAge_bias$Fixed = agg_time_Fixed_noAge
models_fit_noAge_bias$BL = agg_time_BL_noAge

models_fit_PopLevel_bias = agg_overall_PopLevel

agg_overall_bias = list()
agg_overall_bias$BL = agg_overall_BL
agg_overall_bias$Indep = agg_overall_Indep
agg_overall_bias$Fixed = agg_overall_Fixed
agg_overall_bias$RW = agg_overall_RW

save(models_fit_noAge_bias, file = paste0("../simulation_fitted_model_noAge_bias", "_age_", n_Age, "_time_", n_Time,  "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y, ".rda"))
save(models_fit_PopLevel_bias, file = paste0("../simulation_fitted_model_PopLevel_bias", "_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y, ".rda"))    
save(agg_overall_bias, file = paste0("../simulation_fitted_model_overall_bias", "_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y, ".rda"))  


out_agg_time = NULL
out_agg_overall = NULL
for(DY in c(TRUE, FALSE)){

  load(paste0("../simulation_fitted_model_agg_time_bias","_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DY, ".rda"))
  load(paste0("../simulation_fitted_model_noAge_bias","_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DY, ".rda"))
  load(paste0("../simulation_fitted_model_PopLevel_bias","_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DY, ".rda"))
  load(paste0("../simulation_fitted_model_overall_bias","_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DY, ".rda"))
  
  
  agg_time_RW_noAge_df =data.frame(agg_bias = as.vector(models_fit_noAge_bias$RW),
                                   agg_prev = as.vector(sweep(models_fit_noAge_bias$RW, 2, agg_time_bias$agg_time_param, `+`)),
                                   time = rep(1:n_Time, each = 50),
                                   itr = rep(seq(50), n_Time),
                                   upper_sampling = US,
                                   type = "RW1 RE Model",
                                   group = "Time-only Model",
                                   spec = "RW1 RE Model")
  agg_time_Indep_noAge_df =data.frame(agg_bias = as.vector(models_fit_noAge_bias$Indep),
                                      agg_prev =  as.vector(sweep(models_fit_noAge_bias$Indep, 2, agg_time_bias$agg_time_param, `+`)),
                                      time = rep(1:n_Time, each = 50),
                                      itr = rep(seq(50), n_Time),
                                      upper_sampling = US,
                                      type = "Indep RE Model",
                                      group = "Time-only Model", 
                                      spec = "Indep RE Model")
  agg_time_Fixed_noAge_df =data.frame(agg_bias = as.vector(models_fit_noAge_bias$Fixed),
                                      agg_prev =  as.vector(sweep(models_fit_noAge_bias$Fixed, 2, agg_time_bias$agg_time_param, `+`)),
                                      time = rep(1:n_Time, each = 50),
                                      itr = rep(seq(50), n_Time),
                                      upper_sampling = US,
                                      type = "Fixed Effect Model", 
                                      group = "Time-only Model",
                                      spec = "Fixed Effect Model")
  agg_time_BL_noAge_df =data.frame(agg_bias = as.vector(models_fit_noAge_bias$BL),                                   
                                   agg_prev =  as.vector(sweep(models_fit_noAge_bias$BL, 2, agg_time_bias$agg_time_param, `+`)),
                                   time = rep(1:n_Time, each = 50),
                                   itr = rep(seq(50), n_Time),
                                   upper_sampling = US,
                                   type = "Unstructured Baseline",
                                   group = "Time-only Model",
                                   spec = "Unstructured Baseline")
  
  
  agg_time_BL_df = data.frame(agg_bias = as.vector(agg_time_bias$BL),
                              agg_prev = as.vector(sweep(agg_time_bias$BL, 2, agg_time_bias$agg_time_param, `+`)),
                              time = rep(1:n_Time, each = 50),
                              itr = rep(seq(50), n_Time),
                              upper_sampling = US,
                              type = "Unstructured Baseline", 
                              group = "Age-Sex-Time Model",
                              spec = "Unstructured Baseline")
  agg_time_Fixed_df = data.frame(agg_bias = as.vector(agg_time_bias$Fixed), 
                                 agg_prev = as.vector(sweep(agg_time_bias$Fixed, 2, agg_time_bias$agg_time_param, `+`)),
                                 time = rep(1:n_Time, each = 50),
                                 itr = rep(seq(50), n_Time),
                                 upper_sampling = US,
                                 type = "Fixed Effect Model",
                                 group = "Age-Sex-Time Model",
                                 spec = "Fixed Effect Model")
  agg_time_Indep_df = data.frame(agg_bias = as.vector(agg_time_bias$Indep), 
                                 agg_prev = as.vector(sweep(agg_time_bias$Indep, 2, agg_time_bias$agg_time_param, `+`)),
                                 time = rep(1:n_Time, each = 50),
                                 itr = rep(seq(50), n_Time),
                                 upper_sampling = US,
                                 type = "Indep RE Model",
                                 group = "Age-Sex-Time Model",
                                 spec = "Indep RE Model")
  agg_time_RW_df = data.frame(agg_bias = as.vector(agg_time_bias$RW), 
                              agg_prev = as.vector(sweep(agg_time_bias$RW, 2, agg_time_bias$agg_time_param, `+`)),
                              time = rep(1:n_Time, each = 50),
                              itr = rep(seq(50), n_Time),
                              upper_sampling = US,
                              type = "RW1 RE Model", 
                              group = "Age-Sex-Time Model",
                              spec = "RW1 RE Model")
  
  
  
  out_agg_time <- rbind(out_agg_time, agg_time_BL_noAge_df, agg_time_Fixed_noAge_df, agg_time_Indep_noAge_df, agg_time_RW_noAge_df, agg_time_BL_df, agg_time_Fixed_df, agg_time_Indep_df, agg_time_RW_df)
  
  
  agg_overall_BL_df = data.frame(agg_prev = agg_overall_bias$BL, upper_sampling = US, type = "Unstructured\nBaseline")
  agg_overall_Fixed_df = data.frame(agg_prev = agg_overall_bias$Fixed, upper_sampling = US, type = "Fixed Effect\nModel")
  agg_overall_Indep_df = data.frame(agg_prev = agg_overall_bias$Indep, upper_sampling = US, type = "Indep RE\nModel")
  agg_overall_RW_df = data.frame(agg_prev = agg_overall_bias$RW, upper_sampling = US, type = "RW1 RE\nModel")
  
  agg_overall_PopLevel_df = data.frame(agg_prev = models_fit_PopLevel_bias, upper_sampling = US, type = "Unstratified\nBaseline")
  
  out_agg_overall <- rbind(out_agg_overall, agg_overall_PopLevel_df, agg_overall_BL_df, agg_overall_Fixed_df, agg_overall_Indep_df,agg_overall_RW_df)
}


out_agg_time <- subset(out_agg_time, !is.na(agg_prev))
# Check all equal to 50*10
table(out_agg_time$type, out_agg_time$group)
out_agg_time$type <- factor(out_agg_time$spec, levels = c("Unstructured Baseline", "Fixed Effect Model", "Indep RE Model", "RW1 RE Model"))


time_index = month.abb[1:10]
out_agg_time$time <- time_index[out_agg_time$time]
out_agg_time$time <- factor(out_agg_time$time, levels = time_index)
out_agg_time$itr = factor(out_agg_time$itr, levels = 1:500)
out_agg_time$upper_sampling = factor(out_agg_time$upper_sampling, levels = c(12, 14), labels = c("Verification Case I", "Verification Case II"))

out_agg_time_sub = subset(out_agg_time, spec %in% c("RW1 RE Model", "Unstructured Baseline"))
out_agg_time_sub$spec = factor(out_agg_time_sub$spec , levels = c("Unstructured Baseline", "RW1 RE Model"))


out_agg_time_sub_CaseI <- out_agg_time_sub[which(out_agg_time_sub$upper_sampling == "Verification Case I"), ]
out_agg_time_sub_CaseII <- out_agg_time_sub[which(out_agg_time_sub$upper_sampling == "Verification Case II"), ]

g_CaseI_aggregate_time <- ggplot(out_agg_time_sub_CaseI) +  aes(x = time, y = agg_bias) + 
  geom_jitter(aes(color = spec), width = 0.1, alpha = 0.4, size = 1.8) + 
  geom_boxplot(alpha = 0.01, outlier.size = 0.8, linewidth = 0.8, color = "grey30", position = position_dodge(.9)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.9, color = "red", linewidth = 1) + 
  scale_color_manual("", values = c("Unstructured Baseline" = "#d95f02","RW1 RE Model"= "#1b9e77")) + 
  facet_grid(group~type) +
  xlab("Time periods") + ylab("Bias") + 
  theme_bw() + 
  ylim(-0.2, 0.12)+
  theme(legend.position = 'none', strip.text = element_text(size = 10, face = "bold"))

g_CaseI_aggregate_time



g_CaseII_aggregate_time <- ggplot(out_agg_time_sub_CaseII) +  aes(x = time, y = agg_bias) + 
  geom_jitter(aes(color = spec), width = 0.1, alpha = 0.4, size = 1.8) + 
  geom_boxplot(alpha = 0.01, outlier.size = 0.8, linewidth = 0.8, color = "grey30", position = position_dodge(.9)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.9, color = "red", linewidth = 1) + 
  scale_color_manual("", values = c("Unstructured Baseline" = "#d95f02","RW1 RE Model"= "#1b9e77")) + 
  facet_grid(group~type) +
  xlab("Time periods") + ylab("Bias") + 
  theme_bw() + 
  ylim(-0.2, 0.12)+
  theme(legend.position = 'none', strip.text = element_text(size = 10, face = "bold"))

g_CaseII_aggregate_time


###### Aggregate to all boxplot #########
out_agg_overall$type <- factor(out_agg_overall$type, levels = c("Unstratified\nBaseline", 
                                                                "Unstructured\nBaseline", "Fixed Effect\nModel", "Indep RE\nModel", "RW1 RE\nModel"))
out_agg_overall <- subset(out_agg_overall, !is.na(agg_prev))
# check all equal to 50
table(out_agg_overall$type)
out_agg_overall$upper_sampling = factor(out_agg_overall$upper_sampling, levels = c(12, 14), labels = c("Verification Case I", "Verification Case II"))

out_agg_overall_CaseI <- out_agg_overall[which(out_agg_overall$upper_sampling == "Verification Case I"), ]
out_agg_overall_CaseII <- out_agg_overall[which(out_agg_overall$upper_sampling == "Verification Case II"), ]


cols <- met.brewer("Juarez", n = 6, type = "discrete")[-4]

g_CaseI_aggregate_all <-  ggplot(out_agg_overall_CaseI, aes(x = type, y = agg_prev)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 1, color = "red", linewidth = 0.7) + 
  geom_violin(aes(fill = type), alpha = 0.8, linewidth = 0.7, draw_quantiles = c(0.1, 0.5, 0.9), color = "gray20") + 
  scale_fill_manual("", values = c("Unstratified\nBaseline" = cols[5], "Unstructured\nBaseline" = cols[1],"Fixed Effect\nModel" = cols[2],"Indep RE\nModel" = cols[3],"RW1 RE\nModel"= cols[4])) + 
  theme_bw() + 
  ylim(-0.15, 0.05)+
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1),) + xlab("") + ylab(("Bias"))

g_CaseII_aggregate_all <-  ggplot(out_agg_overall_CaseII, aes(x = type, y = agg_prev)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 1, color = "red", linewidth = 0.7) + 
  geom_violin(aes(fill = type), alpha = 0.8, linewidth = 0.7, draw_quantiles = c(0.1, 0.5, 0.9), color = "gray20") + 
  scale_fill_manual("", values = c("Unstratified\nBaseline" = cols[5], "Unstructured\nBaseline" = cols[1],"Fixed Effect\nModel" = cols[2],"Indep RE\nModel" = cols[3],"RW1 RE\nModel"= cols[4])) + 
  theme_bw() + 
  ylim(-0.15, 0.05)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 25, hjust = 1),) + xlab("") + ylab(("Bias"))




g_sim_caseI <- (g_CaseI_aggregate_all + g_CaseI_aggregate_time + plot_layout(widths = c(1.5, 3)) + plot_annotation(title = 'Verification Mechanism Case I')) 
g_sim_caseII <- (g_CaseII_aggregate_all + g_CaseII_aggregate_time + plot_layout(widths = c(1.5, 3))+ plot_annotation(title = 'Verification Mechanism Case I')) 
  
 # plot_annotation(tag_levels = c('I', 'II'), tag_prefix = 'Case ') 
g_sim_cases <- g_sim_caseI / g_sim_caseII
g_sim_cases

g_sim_cases <- (g_CaseI_aggregate_all + g_CaseI_aggregate_time + plot_layout(widths = c(1.5, 3))) / 
  (g_CaseII_aggregate_all + g_CaseII_aggregate_time + plot_layout(widths = c(1.5, 3))) + 
  plot_annotation(tag_levels = 'A', tag_prefix = '(', tag_suffix = ')') 
g_sim_cases
ggsave(g_sim_cases, file = "../figures/final/simulation_aggregated_time_and_all_boxplots.pdf", width = 12, height = 9)



############## ############## ########################################## ############## ################################
######################## II. Simulation Prevalence Posterior mean and CI ########################################## ############## ##############
############## ############## ########################################## ############## ################################

n_Age = 8
n_Time = 10
n_Sex = 2

n_subsize = 100

for(case in 1:2){
  if(case == 1){
    ### Case I ###
    DEPEND_ON_Y = FALSE 
    load(paste0("../simulation_fitted_model_itr_", i, "_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y, ".rda"))
    filetitle <- "simulation_1"
  }

  if(case == 2){
    DEPEND_ON_Y = TRUE 
    load(paste0("../simulation_fitted_model_itr_", i, "_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y, ".rda"))
    filetitle <- "simulation_2"
  }
  ####### plot #####
  fitted_model_BL = models_fit$BL
  fitted_model_Fixed = models_fit$Fixed 
  fitted_model_Indep = models_fit$Indep
  fitted_model_RW = models_fit$RW
  sim_data = models_fit$sim_data
 
  est0 <- data.frame(sex = c(rep("F", n_Time * n_Age), rep("M", n_Time * n_Age)),
                     time = rep(rep(1:n_Time, each = n_Age), 2), 
                     age = rep(rep(1:n_Age, n_Time), 2), 
                     prev = fitted_model_BL$bias+sim_data$param$pi, 
                     lower = apply(fitted_model_BL$poster_sample$sample_pis,2, function(x) quantile(x, 0.025)),
                     upper = apply(fitted_model_BL$poster_sample$sample_pis,2, function(x) quantile(x, 0.975)),
                     type = "BL")
  est1 <- data.frame(sex = c(rep("F", n_Time * n_Age), rep("M", n_Time * n_Age)),
                     time = rep(rep(1:n_Time, each = n_Age), 2), 
                     age = rep(rep(1:n_Age, n_Time), 2), 
                     prev = fitted_model_Fixed$bias+sim_data$param$pi, 
                     lower = apply(fitted_model_Fixed$poster_sample$sample_pis,2, function(x) quantile(x, 0.025)),
                     upper = apply(fitted_model_Fixed$poster_sample$sample_pis,2, function(x) quantile(x, 0.975)),
                     type = "Fixed")
  est2 <-  data.frame(sex = c(rep("F", n_Time * n_Age), rep("M", n_Time * n_Age)),
                      time = rep(rep(1:n_Time, each = n_Age), 2), 
                      age = rep(rep(1:n_Age, n_Time), 2), 
                      prev = fitted_model_Indep$bias+sim_data$param$pi, 
                      lower = apply(fitted_model_Indep$poster_sample$sample_pis,2, function(x) quantile(x, 0.025)),
                      upper = apply(fitted_model_Indep$poster_sample$sample_pis,2, function(x) quantile(x, 0.975)),
                      type = "Indep")
  est3 <-  data.frame(sex = c(rep("F", n_Time * n_Age), rep("M", n_Time * n_Age)),
                      time = rep(rep(1:n_Time, each = n_Age), 2), 
                      age = rep(rep(1:n_Age, n_Time), 2), 
                      prev = fitted_model_RW$bias+sim_data$param$pi, 
                      lower = apply(fitted_model_RW$poster_sample$sample_pis,2, function(x) quantile(x, 0.025)),
                      upper = apply(fitted_model_RW$poster_sample$sample_pis,2, function(x) quantile(x, 0.975)),
                      type = "RW")


  param <- data.frame(sex = c(rep("F", n_Time * n_Age), rep("M", n_Time * n_Age)),
                      time = rep(rep(1:n_Time, each = n_Age), 2), 
                      age = rep(rep(1:n_Age, n_Time), 2), 
                      prev = sim_data$param$pi, 
                      lower = NA,
                      upper = NA)

  emp.obs <- aggregate(sim_data$data$Y.t ~sim_data$data$A + sim_data$data$Ti + sim_data$data$S  , FUN = mean, na.rm = TRUE)

  observed <- data.frame(sex = c(rep("F", n_Time * n_Age), rep("M", n_Time * n_Age)),
                         time = rep(rep(1:n_Time, each = n_Age), 2), 
                         age = rep(rep(1:n_Age, n_Time), 2), 
                         prev = emp.obs$`sim_data$data$Y.t`, 
                         lower = NA,
                         upper = NA)


  out <- rbind(est0, est1, est2, est3)

  time_index = month.abb[1:10]
  out$time <- time_index[out$time]
  out$time <- factor(out$time, levels = time_index)
  param$time <- time_index[param$time]
  param$time <- factor(param$time, levels = time_index)

  observed$time <- factor(observed$time, levels = 1:n_Time)
  out$age <- factor(out$age, levels = 1:n_Age)
  param$age <- factor(param$age, levels = 1:n_Age)
  observed$age <- factor(observed$age, levels = 1:n_Age)
  out$sex <- factor(out$sex, levels = c("F", "M"))
  param$sex <- factor(param$sex, levels = c("F", "M"))
  observed$sex <- factor(observed$sex,  levels = c("F", "M"))
  out$type <- factor(out$type, levels = c("BL", "Fixed", "Indep","RW"))


  out$type <- revalue(out$type, c("BL" = "Unstructured", 
                                  "Fixed" = "Fixed Effect", 
                                  "Indep" = "Indep RE", 
                                  "RW" = "RW1 RE"))
  out$sex <- revalue(out$sex, c("M" = "Male", 
                                "F" = "Female"))
  param$sex <- revalue(param$sex, c("M" = "Male", 
                                "F" = "Female"))

  cols <- met.brewer("Juarez", n = 6, type = "discrete")[-4]
  g_posterior <- ggplot(out, aes(x = age, y = prev))  + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = type, group = interaction(time, type)), alpha = 0.5) + 
    geom_line(aes(color = type, group = interaction(time, type)), linewidth = 1.05) + 
    geom_point(data = param, aes(x = age, y = prev), col = "black", size = 1.2) + 
    facet_grid(sex * type  ~ time) +
    scale_fill_manual("", values = c("Unstructured" = cols[1],"Fixed Effect" = cols[2],"Indep RE" = cols[3], "RW1 RE" = cols[4]))+  
    scale_color_manual("", values = c("Unstructured" = cols[1],"Fixed Effect" = cols[2],"Indep RE" = cols[3], "RW1 RE" = cols[4])) + 
    theme_bw() +
    xlab("Age Group") + ylab("Prevalence") +
    theme(axis.text.x = element_text(size = 10), # Make X axis labels (age) bold
          strip.text = element_text(face = "bold",size = 10),  # Make facet labels (time) bold
          legend.text = element_text(size = 10))+ # Make legend text bold and increase font size

    theme(legend.position = "none", strip.text = element_text(size = 12))

  ggsave(g_posterior, file = paste0("../figures/final/", filetitle, "_posterior.pdf"), width = 12, height = 10)


  g_posterior2 <- ggplot(subset(out, type %in% c("Unstructured", "RW1 RE")), aes(x = age, y = prev))  + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = type, group = interaction(time, type)), alpha = 0.4) + 
    geom_line(aes(color = type, group = interaction(time, type)), linewidth = 1.05) + 
    geom_point(data = param, aes(x = age, y = prev), col = "black", size = 1.2) + 
    facet_grid(sex ~ time) +
    scale_color_manual("Model", values = c("#d95f02", "#1b9e77"))+ 
    scale_fill_manual("Model", values = c("#d95f02", "#1b9e77"))+ 
    theme_bw() +
    xlab("Age Group") + ylab("Prevalence") +
    theme(axis.text.x = element_text(size = 10), # Make X axis labels (age) bold
          strip.text = element_text(face = "bold",size = 10),  # Make facet labels (time) bold
          legend.text = element_text(size = 10))+ 
    theme(legend.position = "bottom", strip.text = element_text(size = 12))

  ggsave(g_posterior2, file = paste0("../figures/final/", filetitle, "_posterior2.pdf"), width = 12, height = 6)

}




#################################### ############## ####################################################################
##################################### III. Simulation true prevalance  ################################
################################### ############## ####################################################################

n_Age = 8
n_Time = 10
n_Sex = 2

DEPEND_ON_Y = FALSE # or TRUE, corresponding to case i and ii of whether sampling with dependency on Y
n_subsize = 50

S_ind = c(rep("Female", n_Time * n_Age), rep("Male", n_Time*n_Age))
T_ind =  rep(rep(1:n_Time, each = n_Age), n_Sex)
A_ind =  rep(rep(1:n_Age,  n_Time), n_Sex)

i = 1
load(paste0("../simulation_fitted_model_itr_", i, "_age_", n_Age, "_time_", n_Time, "_n_Subsize_", n_subsize, "_dependOnY_", DEPEND_ON_Y, ".rda"))


prev_simulated_df = data.frame(pi = models_fit$sim_data$param$pi,
                               sex = S_ind,
                               time = T_ind,
                               age = A_ind)

time_index = month.abb[1:10]
prev_simulated_df$time <- time_index[prev_simulated_df$time]
prev_simulated_df$time <- factor(prev_simulated_df$time, levels = time_index)

prev_simulated_df$sex = as.factor(prev_simulated_df$sex)
prev_simulated_df$age = as.factor(prev_simulated_df$age)

cols <- met.brewer("Juarez", n = 6, type = "discrete")[-2]
g_simulated_prev <- ggplot(prev_simulated_df, aes(x = age, y = pi)) +
  geom_line(aes(color = sex, group = interaction(sex, time)), linewidth = 0.85) + 
  geom_point(aes(color = sex, alpha = 0.7))+
  facet_wrap(~time, nrow = 1, ncol = 10) +
  labs(title = " ",
       x = "Age Group",
       y = "Prevalence",
       color = "Sex") +
  scale_color_manual(values = c("Male" = cols[1], "Female" = cols[5]))+
  theme(axis.text.x = element_text(size = 10), # Make X axis labels (age) bold
        axis.text.y = element_text(face = "bold",size = 10), # Make Y axis labels bold
        strip.text = element_text(face = "bold",size = 10),  # Make facet labels (time) bold
        legend.text = element_text(size = 10))+ # Make legend text bold and increase font size
  guides(alpha = "none") # Remove the alpha legend

g_simulated_prev

