###
#
# VAKI_diagnosis_functions.R
#
###

library(openxlsx)
library(plyr)
library(tidyverse)
library(randomForestSRC)
library(stringr)
library(ggplot2)
library(cobalt)
library(Metrics)
library(MLmetrics)
library(WeightIt)
library(predtools)
library(rms)
library(parallel)
library(pbapply)
library(plotrix)
library(CalibratR)
library(BART)

#### function to load RData and return it ####
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#### function for parallel cluster bootstrap ####

cluster_bootstrap = function(train, R, cluster, statistics, test, test_CA, variables, output_folder, IC_names, seed){
  
  set.seed(seed)
  
  # print status
  print("Resampling the data...")
  
  # sample clusters with replacement, store the samples in a list
  clusters = as.numeric(unique(train[,cluster]))
  n_clusters = length(clusters)
  samples = list()
  for(i in 1:R){samples[[i]] = sample(clusters, n_clusters, replace = TRUE)}
  
  # obtain the data of the resampled clusters and store in the same list
  for(i in 1:R){samples[[i]] = do.call(rbind, lapply(samples[[i]], FUN = function(x, train, cluster){train[train[[cluster]] == x,]}, train, cluster))}
  
  # save the bootstrap data
  save(samples, file = paste0(output_folder, "/RData/bootstrap_data.RData"))
  
  # make cluster (CPU 'cluster', n.b.) for parallel processing
  cl = makeCluster(detectCores() - 1)
  
  # load libraries in cluster
  clusterEvalQ(cl, library(openxlsx))
  clusterEvalQ(cl, library(plyr))
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(randomForestSRC))
  clusterEvalQ(cl, library(stringr))
  clusterEvalQ(cl, library(ggplot2))
  clusterEvalQ(cl, library(cobalt))
  clusterEvalQ(cl, library(Metrics))
  clusterEvalQ(cl, library(MLmetrics))
  clusterEvalQ(cl, library(WeightIt))
  clusterEvalQ(cl, library(predtools))
  clusterEvalQ(cl, library(rms))
  clusterEvalQ(cl, library(parallel))
  clusterEvalQ(cl, library(pbapply))
  clusterEvalQ(cl, library(plotrix))
  clusterEvalQ(cl, library(BART))
  clusterEvalQ(cl, library(CalibratR))
  clusterExport(cl, list("calibration", "PPV", "seed"))
  clusterEvalQ(cl, set.seed(seed))
  
  # print status
  print("Modeling...")
  
  # for each sample: run the supplied statistics function parallelized over the CPUs
  #results = pblapply(samples, FUN = statistics, test = test, test_CA = test_CA, variables = variables, output_folder = output_folder, IC_names = IC_names, seed = seed, bootstrap = T, cl = cl)
  results = pblapply(samples, FUN = function(x, statistics, ...){tryCatch(statistics(x, ...), error = function(e){return(NA)})}, statistics, test = test, test_CA = test_CA, variables = variables, output_folder = output_folder, IC_names = IC_names, seed = seed, bootstrap = T, cl = cl)
  
  # stop the cluster
  stopCluster(cl)
  
  # return the results
  return(results)
  
}

#### function for calibration ####

calibration = function(predicted_train, label_train, predicted_test, label_test, seed, file, subset = NA, bootstrap = F, method = "GUESS"){
  
  library(stats)
  library(CalibratR)
  set.seed(seed)
  
  # placeholders for plots
  plot_before = NA
  plot_after_platt = NA
  plot_after_isotonic = NA
  plot_after_GUESS = NA
  
  # if subset, then use subset to fit and test
  if (length(subset) > 1){
    predicted_train_set = predicted_train[subset[[1]]]
    label_train_set = label_train[subset[[1]]]
    predicted_test_set = predicted_test[subset[[2]]]
    label_test_set = label_test[subset[[2]]]
  } else{
    predicted_train_set = predicted_train
    label_train_set = label_train
    predicted_test_set = predicted_test
    label_test_set = label_test
  }
  
  if (!bootstrap){
    
    # plot before calibration on test set
    png(filename = paste0(file, "before.png"), width = 800, height = 600)
    plot_before = val.prob(predicted_test_set, label_test_set, logistic.cal = F, cex = 1.15, legendloc = c(.65, .37))
    dev.off()
  }  
  
  # calibrate: platt
  model_platt = glm(y~x, data = data.frame(x = predicted_train_set, y = label_train_set), family = "binomial")
  predicted_train_new_platt = predict(object = model_platt, newdata = data.frame(x = predicted_train_set), type = "response")
  predicted_test_new_platt = predict(object = model_platt, newdata = data.frame(x = predicted_test_set), type = "response")
  
  # calibrate: isotonic
  model_isotonic = as.stepfun(isoreg(x = predicted_train_set, y = label_train_set))
  predicted_train_new_isotonic = model_isotonic(predicted_train_set)
  predicted_test_new_isotonic = model_isotonic(predicted_test_set)
  
  # calibrate: GUESS (https://doi.org/10.1093/bioinformatics/bty984)
  model_GUESS = CalibratR::calibrate(label_train_set, predicted_train_set, model_idx = 5, folds = 1, evaluate_no_CV_error = FALSE, evaluate_CV_error = FALSE, nCores = 1, n_seeds = 1)
  predicted_train_new_GUESS = predict_GUESS(model_GUESS$calibration_models$models_final$GUESS, new = predicted_train_set)[["predictions"]]
  predicted_test_new_GUESS = predict_GUESS(model_GUESS$calibration_models$models_final$GUESS, new = predicted_test_set)[["predictions"]]
  
  if (!bootstrap){
    
    # save the after plots for all methods
    png(filename = paste0(file, "platt_after.png"), width = 800, height = 600)
    plot_after_platt = val.prob(predicted_test_new_platt, label_test_set, logistic.cal = FALSE, cex = 1.15, legendloc = c(.65, .37))
    dev.off()
    
    png(filename = paste0(file, "isotonic_after.png"), width = 800, height = 600)
    plot_after_isotonic = val.prob(predicted_test_new_isotonic, label_test_set, logistic.cal = FALSE, cex = 1.15, legendloc = c(.65, .37))
    dev.off()
    
    png(filename = paste0(file, "GUESS_after.png"), width = 800, height = 600)
    plot_after_GUESS = val.prob(predicted_test_new_GUESS, label_test_set, logistic.cal = FALSE, cex = 1.15, legendloc = c(.65, .37))
    dev.off()
    
  }
  
  # make the full calibrated predictions for the complete set in case of a specified subset
  if (length(subset) > 1){
    predicted_train_new_platt = predict(object = model_platt, newdata = data.frame(x = predicted_train), type = "response")
    predicted_test_new_platt = predict(object = model_platt, newdata = data.frame(x = predicted_test), type = "response")
    predicted_train_new_isotonic = model_isotonic(predicted_train)
    predicted_test_new_isotonic = model_isotonic(predicted_test)
    predicted_train_new_GUESS = predict_GUESS(model_GUESS$calibration_models$models_final$GUESS, new = predicted_train)[["predictions"]]
    predicted_test_new_GUESS = predict_GUESS(model_GUESS$calibration_models$models_final$GUESS, new = predicted_test)[["predictions"]]
  }

  # return model, new predictions and plot based on choice of calibration method
  if (method == "platt"){ 
    
    predicted_train_new = predicted_train_new_platt
    predicted_test_new = predicted_test_new_platt
    plot_after = plot_after_platt
    model = model_platt
    
  } else if(method == "isotone"){
    
    predicted_train_new = predicted_train_new_isotonic
    predicted_test_new = predicted_test_new_isotonic
    plot_after = plot_after_isotonic
    model = model_isotonic
    
  } else if (method == "GUESS"){
    
    predicted_train_new = predicted_train_new_GUESS
    predicted_test_new = predicted_test_new_GUESS
    plot_after = plot_after_GUESS
    model = model_GUESS
    
  }
  
  return(list("plot_before" = plot_before, "model" = model, "predicted_train_new" = predicted_train_new, "predicted_test_new" = predicted_test_new, "plot_after" = plot_after))

}

#### function for the PPV ####
PPV = function(true, pred){
  
  if (length(true) != length(pred)) {
    stop("Unequal true and pred vector lengths")
  }
  
  tp = sum(true == 1 & pred == 1)
  fp = sum(true == 0 & pred == 1)
  
  # calculate ppv
  ppv = tp / (tp + fp)
  
  return(ppv)
  
}
#### function for modeling and plotting ####

modeling = function(train, test, test_CA, variables, output_folder, IC_names, seed, bootstrap = F, nodesize = 25){
  
  # RF hyperparameter nodesize to 25, for propensity scores: https://doi.org/10.1016%2Fj.cct.2015.12.012
  # ntree may also be set to 500: https://doi.org/10.1016%2Fj.cct.2015.12.012 
  
  set.seed(seed)
  start = Sys.time()
  
  #### scale and center continuous vars ####
  
  # define continuous vars
  cont = c("age", "nice_ap4_prob", "creatinine_baseline_first", "eGFR", "MAP_min", "UO_rate_min", "serum_albumin_min", "leukocytes_max", "temperature_max", "SOFA_max", "serum_creatinine_max", "time_to_init")
  
  # scale and center
  for(var in cont){
    train[[paste0(var, "_old")]] = train[[var]]
    train[[var]] = scale(train[[var]])[,1]
    test[[paste0(var, "_old")]] = test[[var]]
    test[[var]] = scale(test[[var]])[,1]
  }
  
  #### variable selection ####
  
  # AKI risk factors, SCr baseline, APACHE IV mortality probability, eGFR
  # for sepsis, only use the summarized version
  acute_AKI_risk_factors = c("hypoalbuminemia_max", "acute_heart_failure", "graft_and_transplant_surgery", "hypotension_max", "hypovolemia", "major_surgery", "mechanical_ventilation", "sepsis", "trauma", "burns")
  chronic_AKI_risk_factors = c("age", "alcohol_abuse", "cardiovascular_disease", "chronic_kidney_disease", "chronic_pulmonary_disease", "diabetes_mellitus", "gender", "liver_disease", "malignancy", "obesity")
  additional_AKI_risk_factors = c("nice_ap4_prob", "creatinine_baseline_first", "eGFR")
  
  # other variables
  other_variables = c(paste0(c("MAP", "UO_rate", "serum_albumin"),"_min"),
                      paste0(c("leukocytes", "temperature", "SOFA"),"_max"),
                      "time_to_init")
  
  # consensus nephrotoxins according to Gray et al; https://doi.org/10.1007/s40264-022-01173-4
  nephrotoxins = variables$variables[grepl("_max", variables$variables)]
  nephrotoxins = nephrotoxins[!nephrotoxins %in% c(other_variables, acute_AKI_risk_factors, "AKI_max", "serum_creatinine_max", "sepsis_max")]
  
  # define selected variables and remember all
  selected_variables = c(acute_AKI_risk_factors, chronic_AKI_risk_factors, additional_AKI_risk_factors, other_variables, nephrotoxins)
  selected_variables_all = selected_variables
  
  # drop variables with < 2% prevalence; Patrick et al; https://doi.org/10.1002/pds.2098, they used 5%
  selected_variables = selected_variables[!(colMeans(train[,selected_variables]) < 0.02 & !selected_variables %in% cont)]
  
  # select variables with p <= 0.2 in univariate association with AKI for other variables and nephrotoxins; Boorkhart et al; https://doi.org/10.1093/aje/kwj149
  p = unlist(lapply(selected_variables, function(x){cor.test(train[[x]], train[["AKI_max_num"]], method = "spearman")[["p.value"]]}))
  selected_variables = selected_variables[p <= 0.2 | selected_variables %in% c(acute_AKI_risk_factors, chronic_AKI_risk_factors, additional_AKI_risk_factors)]
  
  # if J01CR05_max (piperacillin+tazobactam) dropped, include it again given past evidence of interaction
  selected_variables = unique(c(selected_variables, "J01CR05_max"))
  
  # write the variables
  vars = data.frame("Variable" = selected_variables_all,
                    "Type" = c(rep("Acute AKI risk factor", length(acute_AKI_risk_factors)), rep("Chronic AKI risk factor", length(chronic_AKI_risk_factors)), rep("Additional AKI risk factor", length(additional_AKI_risk_factors)), rep("Additional variable", length(other_variables)), rep("Nephrotoxin", length(nephrotoxins))),
                    "Status" = ifelse(selected_variables_all %in% selected_variables, "Kept", "Dropped"))
  
  # only save if not bootstrap run
  if(!bootstrap){
    write.xlsx(vars, paste0(output_folder, "/tables/selected_variables.xlsx"))
  }
    
  #### propensity models ####
  
  # GLM
  model_prop_GLM = glm(as.formula(paste("arm~", paste(selected_variables, collapse="+",sep=""), sep = "")),
                       data  = train,
                       family = binomial())
  
  # RF for imbalanced data: random forests quantile-classifer (RFQ) approach of O'Brien and Ishwaran (2017): https://doi.org/10.1016/j.patcog.2019.01.036
  model_prop_RF = imbalanced(as.formula(paste("arm~", paste(selected_variables, collapse="+",sep=""), sep = "")),
                        data = train,
                        nodesize = nodesize,
                        seed = seed)
  
  # BART
  model_prop_BART = gbart(x.train = train[, selected_variables], y.train = train$arm_num, type = "pbart", printevery = 100000, seed = seed)
  
  #### predict propensity scores ####
  
  # predict the propensity scores on train
  propensity_GLM = predict(object = model_prop_GLM, newdata = train, type = "response")
  propensity_RF = predict(object = model_prop_RF, newdata = train, type = "response")
  propensity_BART = predict(object = model_prop_BART, newdata = train[, names(model_prop_BART$treedraws$cutpoints)])[["prob.test.mean"]]
  results = data.frame(propensity_GLM = propensity_GLM, propensity_RF = propensity_RF$predicted[,2], propensity_BART = propensity_BART)
  results$arm = train$arm
  results$arm_num = train$arm_num
  
  # predict the propensities on test
  test_propensity_RF = predict(object = model_prop_RF, newdata = test)
  test_propensity_GLM = predict(object = model_prop_GLM, newdata = test, type = "response")
  test_propensity_BART = predict(object = model_prop_BART, newdata = test[, names(model_prop_BART$treedraws$cutpoints)])[["prob.test.mean"]]
  test$propensity_RF = test_propensity_RF$predicted[,2]
  test$propensity_GLM = test_propensity_GLM
  test$propensity_BART = test_propensity_BART
  
  #### calibrate propensity score models ####
  
  # calibrate (and save plots)
  propensity_RF_calibration_results = calibration(results$propensity_RF, results$arm_num, test$propensity_RF, test$arm_num, seed, paste0(output_folder, "/plots/calibration/RF_propensity_calibration_"), bootstrap = bootstrap)
  propensity_GLM_calibration_results = calibration(results$propensity_GLM, results$arm_num, test$propensity_GLM, test$arm_num, seed, paste0(output_folder, "/plots/calibration/GLM_propensity_calibration_"), bootstrap = bootstrap)
  propensity_BART_calibration_results = calibration(results$propensity_BART, results$arm_num, test$propensity_BART, test$arm_num, seed, paste0(output_folder, "/plots/calibration/BART_propensity_calibration_"), bootstrap = bootstrap)
 
  # assign new predicted probabilities to train and test
  results$propensity_RF_cal = propensity_RF_calibration_results$predicted_train_new
  test$propensity_RF_cal = propensity_RF_calibration_results$predicted_test_new
  
  results$propensity_GLM_cal = propensity_GLM_calibration_results$predicted_train_new
  test$propensity_GLM_cal = propensity_GLM_calibration_results$predicted_test_new
  
  results$propensity_BART_cal = propensity_BART_calibration_results$predicted_train_new
  test$propensity_BART_cal = propensity_BART_calibration_results$predicted_test_new
  
  #### logistic regression assumptions plots for propensity model ####
  
  # only if not bootstrap run
  if(!bootstrap){
    
    # select continuous vars, join logit of the predicted probs and plot
    train %>% mutate(logit = log(results$propensity_GLM_cal/(1 - results$propensity_GLM_cal))) %>% select(c(selected_variables[selected_variables %in% cont], "logit")) %>% gather(predictors, value, -logit) %>%
      ggplot(aes(logit, value))+
      geom_point(size = 0.5, alpha = 0.5) +
      geom_smooth(method = "loess") + 
      theme_bw() + 
      facet_wrap(~predictors, scales = "free_y")
    
    # save the plot
    ggsave(paste0(output_folder, "/plots/LR_assumptions/LR_linearity_propensity.png"))
    
  }
  
  #### propensity plots ####
  
  # only if not bootstrap run
  if(!bootstrap){
    propensity_GLM_plot = ggplot(results, aes(x = propensity_GLM_cal, fill = arm)) +
                                 geom_density(alpha=.25) + xlab("Propensity") + 
                                 ylab("Density") + guides(fill = guide_legend(title = "Arm")) + 
                                 ggtitle("Propensities: logistic regression") + 
                                 xlim(c(0,1)) +
                                 scale_fill_discrete(labels = c('Alternative', 'Vancomycin'))
    
    ggsave(paste0(output_folder, "/plots/balance/propensities_GLM.png"))
    
    
    propensity_RF_plot = ggplot(results, aes(x = propensity_RF_cal, fill = arm)) + 
                                geom_density(alpha=.25) + xlab("Propensity") + 
                                ylab("Density") + guides(fill = guide_legend(title = "Arm")) +
                                ggtitle("Propensities: random forests") + 
                                xlim(c(0,1)) + 
                                scale_fill_discrete(labels = c('Alternative', 'Vancomycin'))
    
    ggsave(paste0(output_folder, "/plots/balance/propensities_RF.png"))
    
    
    propensity_BART_plot = ggplot(results, aes(x = propensity_BART_cal, fill = arm)) + 
                                geom_density(alpha=.25) + xlab("Propensity") + 
                                ylab("Density") + guides(fill = guide_legend(title = "Arm")) +
                                ggtitle("Propensities: BART") + 
                                xlim(c(0,1)) + 
                                scale_fill_discrete(labels = c('Alternative', 'Vancomycin'))
    
    ggsave(paste0(output_folder, "/plots/balance/propensities_BART.png"))
    
    
    propensity_all_plot = results %>% mutate(BART = propensity_BART_cal, RF = propensity_RF_cal, LR = propensity_GLM_cal) %>% select(arm, BART, RF, LR) %>% gather(model, propensity, -arm) %>% mutate(model = factor(model, levels = c("BART", "RF", "LR"))) %>%
                            ggplot(aes(x = propensity, fill = arm)) + 
                            geom_density(alpha=.25) + xlab("Propensity") + 
                            ylab("Density") + guides(fill = guide_legend(title = "Arm")) +
                            xlim(c(0,1)) + 
                            scale_fill_discrete(labels = c('Alternative', 'Vancomycin')) +
                            facet_wrap(~model) +
                            theme(legend.position = c(0.07, 0.87), legend.background = element_blank())
    
    ggsave(paste0(output_folder, "/plots/balance/propensities_all.png"), width = 12, height = 4)
    
    
    propensity_all_plot_2 = results %>% mutate(BART = propensity_BART_cal, RF = propensity_RF_cal, LR = propensity_GLM_cal) %>% select(arm, BART, RF, LR) %>% gather(model, propensity, -arm) %>% mutate(model = factor(model, levels = c("BART", "LR", "RF"))) %>%
      ggplot(aes(x = propensity, fill = arm)) + 
      geom_density(alpha=.25) + xlab("Propensity") + 
      ylab("Density") + guides(fill = guide_legend(title = "Arm")) +
      scale_x_continuous(breaks = seq(0.00, 1.00, 0.2)) +
      scale_fill_discrete(labels = c('Alternative', 'Vancomycin')) +
      facet_wrap(~model) + theme(panel.spacing.x = unit(1, "lines")) +
      theme(legend.position = c(0.07, 0.87), legend.background = element_blank(), text = element_text(size = 16))
    
    ggsave(paste0(output_folder, "/plots/balance/propensities_all_order.png"), width = 12, height = 4)
  }
    
  #### check overlap and subset ####
  
  # find the area of overlap
  min = max(min(results[results$arm == 0,]$propensity_BART_cal), min(results[results$arm == 1,]$propensity_BART_cal))
  max = min(max(results[results$arm == 0,]$propensity_BART_cal), max(results[results$arm == 1,]$propensity_BART_cal))
  
  # subset of data within overlap
  subset = train[results$propensity_BART_cal < max & results$propensity_BART_cal > min,]
  
  # obtain the results with the propensity scores for the subset with overlap
  results = results[results$propensity_BART_cal < max & results$propensity_BART_cal > min,]
  
  # save overlap report
  # only if not bootstrap run
  if(!bootstrap){
    
    overlap = data.frame("admissions_before" = nrow(train),
                         "exposed cases before" = nrow(train[train$arm == 1 & train$AKI_max == 1,]),
                         "admissions_after" = nrow(subset),
                         "exposed cases after" = nrow(subset[subset$arm == 1 & subset$AKI_max == 1,]),
                         "p_min" = min,
                         "p_max" = max)
    
    write.xlsx(overlap, file = paste0(output_folder, "/tables/overlap_report.xlsx"))
  }
  
  #### outcome models ####
  
  # for treated: RF for imbalanced data: random forests quantile-classifer (RFQ) approach of O'Brien and Ishwaran (2017): https://doi.org/10.1016/j.patcog.2019.01.036
  model_treat_RF = imbalanced(as.formula(paste("AKI_max~", paste(selected_variables, collapse="+",sep=""), sep = "")),
                            data = subset[subset$arm == 1,],
                            nodesize = nodesize,
                            seed = seed)
  
  model_treat_GLM = glm(as.formula(paste("AKI_max_num~", paste(selected_variables, collapse="+",sep=""), sep = "")),
                        data  = subset[subset$arm == 1,],
                        family = binomial())
  
  model_treat_BART = gbart(x.train = subset[subset$arm == 1, selected_variables], y.train = subset[subset$arm == 1,]$AKI_max_num, type = "pbart", printevery = 100000, seed = seed)
  
  # for control: RF for imbalanced data: random forests quantile-classifer (RFQ) approach of O'Brien and Ishwaran (2017): https://doi.org/10.1016/j.patcog.2019.01.036
  model_control_RF = imbalanced(as.formula(paste("AKI_max~", paste(selected_variables, collapse="+",sep=""), sep = "")),
                              data  = subset[subset$arm == 0,],
                              nodesize = nodesize,
                              seed = seed)
  
  model_control_GLM = glm(as.formula(paste("AKI_max_num~", paste(selected_variables, collapse="+",sep=""), sep = "")),
                        data  = subset[subset$arm == 0,],
                        family = binomial())
  
  model_control_BART = gbart(x.train = subset[subset$arm == 0, selected_variables], y.train = subset[subset$arm == 0,]$AKI_max_num, type = "pbart", printevery = 100000, seed = seed)
  
  #### predict outcomes ####
  
  # train set (subset)
  p_treat_RF = predict(object = model_treat_RF, newdata = subset)
  p_control_RF = predict(object = model_control_RF, newdata = subset)
  p_treat_GLM = predict(object = model_treat_GLM, newdata = subset, type = "response")
  p_control_GLM = predict(object = model_control_GLM, newdata = subset, type = "response")
  p_treat_BART = predict(object = model_treat_BART, newdata = subset[, names(model_treat_BART$treedraws$cutpoints)])[["prob.test.mean"]]
  p_control_BART = predict(object = model_control_BART, newdata = subset[, names(model_control_BART$treedraws$cutpoints)])[["prob.test.mean"]]
  
  results$treat_RF = p_treat_RF$predicted[,2]
  results$control_RF = p_control_RF$predicted[,2]
  results$treat_GLM = p_treat_GLM
  results$control_GLM = p_control_GLM
  results$treat_BART = p_treat_BART
  results$control_BART = p_control_BART
  
  results$AKI_max = subset$AKI_max
  results$AKI_max_num = subset$AKI_max_num
  
  # test set
  p_treat_RF = predict(object = model_treat_RF, newdata = test)
  p_control_RF = predict(object = model_control_RF, newdata = test)
  p_treat_GLM = predict(object = model_treat_GLM, newdata = test, type = "response")
  p_control_GLM = predict(object = model_control_GLM, newdata = test, type = "response")
  p_treat_BART = predict(object = model_treat_BART, newdata = test[, names(model_treat_BART$treedraws$cutpoints)])[["prob.test.mean"]]
  p_control_BART = predict(object = model_control_BART, newdata = test[, names(model_control_BART$treedraws$cutpoints)])[["prob.test.mean"]]
  
  test$treat_RF = p_treat_RF$predicted[,2]
  test$control_RF = p_control_RF$predicted[,2]
  test$treat_GLM = p_treat_GLM
  test$control_GLM = p_control_GLM
  test$treat_BART = p_treat_BART
  test$control_BART = p_control_BART
  
  #### calibrate outcome models ####
  
  # calibrate
  treat_RF_calibration_results = calibration(results$treat_RF, results$AKI_max_num, test$treat_RF, test$AKI_max_num, seed, file = paste0(output_folder, "/plots/calibration/RF_treat_calibration_"), subset = list(results$arm == 1, test$arm == 1), bootstrap = bootstrap)
  control_RF_calibration_results = calibration(results$control_RF, results$AKI_max_num, test$control_RF, test$AKI_max_num, seed, file = paste0(output_folder, "/plots/calibration/RF_control_calibration_"), subset = list(results$arm == 0, test$arm == 0), bootstrap = bootstrap)
  treat_GLM_calibration_results = calibration(results$treat_GLM, results$AKI_max_num, test$treat_GLM, test$AKI_max_num, seed, file = paste0(output_folder, "/plots/calibration/GLM_treat_calibration_"), subset = list(results$arm == 1, test$arm == 1), bootstrap = bootstrap)
  control_GLM_calibration_results = calibration(results$control_GLM, results$AKI_max_num, test$control_GLM, test$AKI_max_num, seed, file = paste0(output_folder, "/plots/calibration/GLM_control_calibration_"), subset = list(results$arm == 0, test$arm == 0), bootstrap = bootstrap)
  treat_BART_calibration_results = calibration(results$treat_BART, results$AKI_max_num, test$treat_BART, test$AKI_max_num, seed, file = paste0(output_folder, "/plots/calibration/BART_treat_calibration_"), subset = list(results$arm == 1, test$arm == 1), bootstrap = bootstrap)
  control_BART_calibration_results = calibration(results$control_BART, results$AKI_max_num, test$control_BART, test$AKI_max_num, seed, file = paste0(output_folder, "/plots/calibration/BART_control_calibration_"), subset = list(results$arm == 0, test$arm == 0), bootstrap = bootstrap)
  
  # assign new predicted probabilities to train and test: RF
  results$treat_RF_cal = treat_RF_calibration_results$predicted_train_new
  results$control_RF_cal = control_RF_calibration_results$predicted_train_new
  test$treat_RF_cal = treat_RF_calibration_results$predicted_test_new
  test$control_RF_cal = control_RF_calibration_results$predicted_test_new
  
  # assign new predicted probabilities to train and test: GLM
  results$treat_GLM_cal = treat_GLM_calibration_results$predicted_train_new
  results$control_GLM_cal = control_GLM_calibration_results$predicted_train_new
  test$treat_GLM_cal = treat_GLM_calibration_results$predicted_test_new
  test$control_GLM_cal = control_GLM_calibration_results$predicted_test_new
  
  # assign new predicted probabilities to train and test: BART
  results$treat_BART_cal = treat_BART_calibration_results$predicted_train_new
  results$control_BART_cal = control_BART_calibration_results$predicted_train_new
  test$treat_BART_cal = treat_BART_calibration_results$predicted_test_new
  test$control_BART_cal = control_BART_calibration_results$predicted_test_new
  
  #### logistic regression assumptions plots for outcome models ####
  
  # only if not bootstrap run
  if(!bootstrap){
    
    # for treat: select continuous vars, join logit of the predicted probs and plot
    subset %>% mutate(logit = log(results$treat_GLM_cal/(1 - results$treat_GLM_cal))) %>% select(c(selected_variables[selected_variables %in% cont], "logit")) %>% gather(predictors, value, -logit) %>%
      ggplot(aes(logit, value))+
      geom_point(size = 0.5, alpha = 0.5) +
      geom_smooth(method = "loess") + 
      theme_bw() + 
      facet_wrap(~predictors, scales = "free_y")
    
    # save the plot
    ggsave(paste0(output_folder, "/plots/LR_assumptions/LR_linearity_treat.png"))
    
    # for control: select continuous vars, join logit of the predicted probs and plot
    subset %>% mutate(logit = log(results$control_GLM_cal/(1 - results$control_GLM_cal))) %>% select(c(selected_variables[selected_variables %in% cont], "logit")) %>% gather(predictors, value, -logit) %>%
      ggplot(aes(logit, value))+
      geom_point(size = 0.5, alpha = 0.5) +
      geom_smooth(method = "loess") + 
      theme_bw() + 
      facet_wrap(~predictors, scales = "free_y")
    
    # save the plot
    ggsave(paste0(output_folder, "/plots/LR_assumptions/LR_linearity_control.png"))
    
  }
  #### effect estimands ####
  
  # ARD and ERR
  results$ARD_RF = results$treat_RF_cal - results$control_RF_cal
  results$ARD_GLM = results$treat_GLM_cal - results$control_GLM_cal
  results$ARD_BART = results$treat_BART_cal - results$control_BART_cal
  results$ERR_RF = unlist(lapply(1 - (1/(results$treat_RF_cal / results$control_RF_cal)), FUN = function(x){max(0, x)}))
  results$ERR_GLM = unlist(lapply(1 - (1/(results$treat_GLM_cal / results$control_GLM_cal)), FUN = function(x){max(0, x)}))
  results$ERR_BART = unlist(lapply(1 - (1/(results$treat_BART_cal / results$control_BART_cal)), FUN = function(x){max(0, x)}))
  
  # IPTW to compare: ATE
  results$weight_GLM = ifelse(results$arm == 1, 1/results$propensity_GLM_cal, 1/(1-results$propensity_GLM_cal))
  results$weight_RF = ifelse(results$arm == 1, 1/results$propensity_RF_cal, 1/(1-results$propensity_RF_cal))
  results$weight_BART = ifelse(results$arm == 1, 1/results$propensity_BART_cal, 1/(1-results$propensity_BART_cal))
  iptw_GLM = results %>% group_by(arm) %>% mutate(r = AKI_max_num * weight_GLM) %>% summarise(r = sum(r), n = sum(weight_GLM)) %>% mutate(risk = r/n) %>% select(arm, risk)
  iptw_RF = results %>% group_by(arm) %>% mutate(r = AKI_max_num * weight_RF) %>% summarise(r = sum(r), n = sum(weight_RF)) %>% mutate(risk = r/n) %>% select(arm, risk)
  iptw_BART = results %>% group_by(arm) %>% mutate(r = AKI_max_num * weight_BART) %>% summarise(r = sum(r), n = sum(weight_BART)) %>% mutate(risk = r/n) %>% select(arm, risk)
  
  # IPTW to compare: ATT
  results$weight_GLM_ATT = ifelse(results$arm == 1, 1, results$propensity_GLM_cal/(1-results$propensity_GLM_cal))
  results$weight_RF_ATT = ifelse(results$arm == 1, 1, results$propensity_RF_cal/(1-results$propensity_RF_cal))
  results$weight_BART_ATT = ifelse(results$arm == 1, 1, results$propensity_BART_cal/(1-results$propensity_BART_cal))
  iptw_GLM_ATT = results %>% group_by(arm) %>% mutate(r = AKI_max_num * weight_GLM_ATT) %>% summarise(r = sum(r), n = sum(weight_GLM_ATT)) %>% mutate(risk = r/n) %>% select(arm, risk)
  iptw_RF_ATT = results %>% group_by(arm) %>% mutate(r = AKI_max_num * weight_RF_ATT) %>% summarise(r = sum(r), n = sum(weight_RF_ATT)) %>% mutate(risk = r/n) %>% select(arm, risk)
  iptw_BART_ATT = results %>% group_by(arm) %>% mutate(r = AKI_max_num * weight_BART_ATT) %>% summarise(r = sum(r), n = sum(weight_BART_ATT)) %>% mutate(risk = r/n) %>% select(arm, risk)  
  
  # add the admission data to the results, don't add cols already there
  results = cbind(results, subset[,!colnames(subset) %in% colnames(results)])
  
  # summarise in a table: ATE
  effects_ATE = data.frame("Method" = c("Unadjusted", "T-learner - BART", "T-learner - RF", "T-learner - LR", "IPTW - BART", "IPTW - RF", "IPTW - LR"),
                       "AKI risk - alternative" = c(mean(results[results$arm == 0,]$AKI_max_num), mean(results$control_BART_cal), mean(results$control_RF_cal), mean(results$control_GLM_cal), iptw_BART$risk[1], iptw_RF$risk[1], iptw_GLM$risk[1]),
                       "AKI risk - vancomycin" = c(mean(results[results$arm == 1,]$AKI_max_num), mean(results$treat_BART_cal), mean(results$treat_RF_cal), mean(results$treat_GLM_cal), iptw_BART$risk[2], iptw_RF$risk[2], iptw_GLM$risk[2]),
                       "Absolute risk difference" = c(mean(results[results$arm == 1,]$AKI_max_num), mean(results$treat_BART_cal), mean(results$treat_RF_cal), mean(results$treat_GLM_cal), iptw_BART$risk[2], iptw_RF$risk[2], iptw_GLM$risk[2]) - c(mean(results[results$arm == 0,]$AKI_max_num), mean(results$control_BART_cal), mean(results$control_RF_cal), mean(results$control_GLM_cal), iptw_BART$risk[1], iptw_RF$risk[1], iptw_GLM$risk[1]),
                       "Risk ratio" = c(mean(results[results$arm == 1,]$AKI_max_num), mean(results$treat_BART_cal), mean(results$treat_RF_cal), mean(results$treat_GLM_cal), iptw_BART$risk[2], iptw_RF$risk[2], iptw_GLM$risk[2]) / c(mean(results[results$arm == 0,]$AKI_max_num), mean(results$control_BART_cal), mean(results$control_RF_cal), mean(results$control_GLM_cal), iptw_BART$risk[1], iptw_RF$risk[1], iptw_GLM$risk[1]),
                       "Excess risk ratio" = 1 - 1/(c(mean(results[results$arm == 1,]$AKI_max_num), mean(results$treat_BART_cal), mean(results$treat_RF_cal), mean(results$treat_GLM_cal), iptw_BART$risk[2], iptw_RF$risk[2], iptw_GLM$risk[2]) / c(mean(results[results$arm == 0,]$AKI_max_num), mean(results$control_BART_cal), mean(results$control_RF_cal), mean(results$control_GLM_cal), iptw_BART$risk[1], iptw_RF$risk[1], iptw_GLM$risk[1])),
                       check.names = F)
  
  # summarise in a table: ATT
  effects_ATT = data.frame("Method" = c("Unadjusted", "T-learner - BART", "T-learner - RF", "T-learner - LR", "IPTW - BART", "IPTW - RF", "IPTW - LR"),
                           "AKI risk - alternative" = c(mean(results[results$arm == 0,]$AKI_max_num), mean(results[results$arm == 1,]$control_BART_cal), mean(results[results$arm == 1,]$control_RF_cal), mean(results[results$arm == 1,]$control_GLM_cal), iptw_BART_ATT$risk[1], iptw_RF_ATT$risk[1], iptw_GLM_ATT$risk[1]),
                           "AKI risk - vancomycin" = c(mean(results[results$arm == 1,]$AKI_max_num), mean(results[results$arm == 1,]$treat_BART_cal), mean(results[results$arm == 1,]$treat_RF_cal), mean(results[results$arm == 1,]$treat_GLM_cal), iptw_BART_ATT$risk[2], iptw_RF_ATT$risk[2], iptw_GLM_ATT$risk[2]),
                           "Absolute risk difference" = c(mean(results[results$arm == 1,]$AKI_max_num), mean(results[results$arm == 1,]$treat_BART_cal), mean(results[results$arm == 1,]$treat_RF_cal), mean(results[results$arm == 1,]$treat_GLM_cal), iptw_BART_ATT$risk[2], iptw_RF_ATT$risk[2], iptw_GLM_ATT$risk[2]) - c(mean(results[results$arm == 0,]$AKI_max_num), mean(results[results$arm == 1,]$control_BART_cal), mean(results[results$arm == 1,]$control_RF_cal), mean(results[results$arm == 1,]$control_GLM_cal), iptw_BART_ATT$risk[1], iptw_RF_ATT$risk[1], iptw_GLM_ATT$risk[1]),
                           "Risk ratio" = c(mean(results[results$arm == 1,]$AKI_max_num), mean(results[results$arm == 1,]$treat_BART_cal), mean(results[results$arm == 1,]$treat_RF_cal), mean(results[results$arm == 1,]$treat_GLM_cal), iptw_BART_ATT$risk[2], iptw_RF_ATT$risk[2], iptw_GLM_ATT$risk[2]) / c(mean(results[results$arm == 0,]$AKI_max_num), mean(results[results$arm == 1,]$control_BART_cal), mean(results[results$arm == 1,]$control_RF_cal), mean(results[results$arm == 1,]$control_GLM_cal), iptw_BART_ATT$risk[1], iptw_RF_ATT$risk[1], iptw_GLM_ATT$risk[1]),
                           "Excess risk ratio" = 1 - 1/(c(mean(results[results$arm == 1,]$AKI_max_num), mean(results[results$arm == 1,]$treat_BART_cal), mean(results[results$arm == 1,]$treat_RF_cal), mean(results[results$arm == 1,]$treat_GLM_cal), iptw_BART_ATT$risk[2], iptw_RF_ATT$risk[2], iptw_GLM_ATT$risk[2]) / c(mean(results[results$arm == 0,]$AKI_max_num), mean(results[results$arm == 1,]$control_BART_cal), mean(results[results$arm == 1,]$control_RF_cal), mean(results[results$arm == 1,]$control_GLM_cal), iptw_BART_ATT$risk[1], iptw_RF_ATT$risk[1], iptw_GLM_ATT$risk[1])),
                           check.names = F)
  
  # summarise in a table: ATTY
  effects_ATTY = data.frame("Method" = c("T-learner - BART", "T-learner - RF", "T-learner - LR"),
                           "AKI risk - alternative" = c(mean(results[results$arm == 1 & results$AKI_max_num == 1,]$control_BART_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$control_RF_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$control_GLM_cal)),
                           "AKI risk - vancomycin" = c(mean(results[results$arm == 1 & results$AKI_max_num == 1,]$treat_BART_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$treat_RF_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$treat_GLM_cal)),
                           "Absolute risk difference" = c(mean(results[results$arm == 1 & results$AKI_max_num == 1,]$treat_BART_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$treat_RF_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$treat_GLM_cal)) - c(mean(results[results$arm == 1 & results$AKI_max_num == 1,]$control_BART_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$control_RF_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$control_GLM_cal)),
                           "Risk ratio" = c(mean(results[results$arm == 1 & results$AKI_max_num == 1,]$treat_BART_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$treat_RF_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$treat_GLM_cal)) / c(mean(results[results$arm == 1 & results$AKI_max_num == 1,]$control_BART_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$control_RF_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$control_GLM_cal)),
                           "Excess risk ratio" = 1 - 1/(c(mean(results[results$arm == 1 & results$AKI_max_num == 1,]$treat_BART_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$treat_RF_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$treat_GLM_cal)) / c(mean(results[results$arm == 1 & results$AKI_max_num == 1,]$control_BART_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$control_RF_cal), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$control_GLM_cal))),
                           check.names = F)
  
  # summarise in a table: PC
  effects_PC = data.frame("Method" = c("T-learner - BART", "T-learner - RF", "T-learner - LR"),
                          "ATE" = c(mean(results$ERR_BART), mean(results$ERR_RF), mean(results$ERR_GLM)),
                          "ATT" = c(mean(results[results$arm == 1,]$ERR_BART), mean(results[results$arm == 1,]$ERR_RF), mean(results[results$arm == 1,]$ERR_GLM)),
                          "ATTY" = c(mean(results[results$arm == 1 & results$AKI_max_num == 1,]$ERR_BART), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$ERR_RF), mean(results[results$arm == 1 & results$AKI_max_num == 1,]$ERR_GLM)))
  
  #### AUC on train and test, save test vanco AKIs ####
  
  # AUCs on train, test and validation
  performance = data.frame("Model" = c("Propensity - BART", "Propensity - RF", "Propensity - LR", "Vancomycin - BART", "Vancomycin - RF", "Vancomycin - LR", "Control - BART", "Control - RF", "Control - LR"),
                           "Train AUC (OOB for RF)" = c(auc(train$arm_num, model_prop_BART$prob.train.mean),
                                                        auc(model_prop_RF$yvar, model_prop_RF$predicted.oob[,2]),
                                                        auc(model_prop_GLM$y, model_prop_GLM$fitted.values),
                                                        auc(subset[subset$arm == 1,]$AKI_max_num, model_treat_BART$prob.train.mean),
                                                        auc(model_treat_RF$yvar, model_treat_RF$predicted.oob[,2]),
                                                        auc(model_treat_GLM$y, model_treat_GLM$fitted.values),
                                                        auc(subset[subset$arm == 0,]$AKI_max_num, model_control_BART$prob.train.mean),
                                                        auc(model_control_RF$yvar, model_control_RF$predicted.oob[,2]),
                                                        auc(model_control_GLM$y, model_control_GLM$fitted.values)),
                           "Test AUC" = c(auc(test$arm, test$propensity_BART_cal),
                                          auc(test$arm, test$propensity_RF_cal),
                                          auc(test$arm, test$propensity_GLM_cal),
                                          auc(test[test$arm == 1,]$AKI_max, test[test$arm == 1,]$treat_BART_cal),
                                          auc(test[test$arm == 1,]$AKI_max, test[test$arm == 1,]$treat_RF_cal),
                                          auc(test[test$arm == 1,]$AKI_max, test[test$arm == 1,]$treat_GLM_cal),
                                          auc(test[test$arm == 0,]$AKI_max, test[test$arm == 0,]$control_BART_cal),
                                          auc(test[test$arm == 0,]$AKI_max, test[test$arm == 0,]$control_RF_cal),
                                          auc(test[test$arm == 0,]$AKI_max, test[test$arm == 0,]$control_GLM_cal)),
                           "Test AUC-PR" = c(PRAUC(test$propensity_BART_cal, test$arm),
                                          PRAUC(test$propensity_RF_cal, test$arm),
                                          PRAUC(test$propensity_GLM_cal, test$arm),
                                          PRAUC(test[test$arm == 1,]$treat_BART_cal, test[test$arm == 1,]$AKI_max),
                                          PRAUC(test[test$arm == 1,]$treat_RF_cal, test[test$arm == 1,]$AKI_max),
                                          PRAUC(test[test$arm == 1,]$treat_GLM_cal, test[test$arm == 1,]$AKI_max),
                                          PRAUC(test[test$arm == 0,]$control_BART_cal, test[test$arm == 0,]$AKI_max),
                                          PRAUC(test[test$arm == 0,]$control_RF_cal, test[test$arm == 0,]$AKI_max),
                                          PRAUC(test[test$arm == 0,]$control_GLM_cal, test[test$arm == 0,]$AKI_max)),
                           check.names = FALSE)
  
  # set rownames as first col
  rownames(performance) = performance[,1]
  performance = performance[,-1]
  
  # calculate causal effect estimands for test set
  test$ARD_RF = test$treat_RF_cal - test$control_RF_cal
  test$ARD_GLM = test$treat_GLM_cal - test$control_GLM_cal
  test$ARD_BART = test$treat_BART_cal - test$control_BART_cal
  test$ERR_RF = unlist(lapply(1 - (1/(test$treat_RF_cal / test$control_RF_cal)), FUN = function(x){max(0, x)}))
  test$ERR_GLM = unlist(lapply(1 - (1/(test$treat_GLM_cal / test$control_GLM_cal)), FUN = function(x){max(0, x)}))
  test$ERR_BART = unlist(lapply(1 - (1/(test$treat_BART_cal / test$control_BART_cal)), FUN = function(x){max(0, x)}))
  
  # gather info about the admissions with vanco and AKI in the test set
  test_vanco_AKIs = data.frame("ID" = test[test$arm == 1 & test$AKI_max == 1,]$ID,
                               "IC" = do.call(rbind,strsplit(test[test$arm == 1 & test$AKI_max == 1,]$ID, "_"))[,1],
                               "IC_name" = NA,
                               "admno" = do.call(rbind,strsplit(test[test$arm == 1 & test$AKI_max == 1,]$ID, "_"))[,2],
                               "nice_adm_icu" = test[test$arm == 1 & test$AKI_max == 1,]$nice_adm_icu,
                               "nice_dis_icu" = test[test$arm == 1 & test$AKI_max == 1,]$nice_dis_icu,
                               "vancomycin_init" = test[test$arm == 1 & test$AKI_max == 1,]$nice_adm_icu + (test[test$arm == 1 & test$AKI_max == 1,]$time_to_init_old * 24 * 60 * 60),
                               "AKI_time" = test[test$arm == 1 & test$AKI_max == 1,]$AKI_SCr_static,
                               "ERR_RF" = test[test$arm == 1 & test$AKI_max == 1,]$ERR_RF,
                               "ERR_GLM" = test[test$arm == 1 & test$AKI_max == 1,]$ERR_GLM,
                               "ERR_BART" = test[test$arm == 1 & test$AKI_max == 1,]$ERR_BART)
  
  test_vanco_AKIs$IC_name = join(test_vanco_AKIs, IC_names, by = "IC")[["IC.name"]]
  
  # save test set performance and test set ICs / admno's with vanco and AKI (for causality assessments)
  # only if not bootstrap run
  if(!bootstrap){
    save(performance, file = paste0(output_folder, "/RData/performance_point.RData"))
    write.xlsx(test_vanco_AKIs[,1:(ncol(test_vanco_AKIs) - 3)], file = paste0(output_folder, "/tables/test_vanco_AKIs.xlsx"))
    write.xlsx(test_vanco_AKIs, file = paste0(output_folder, "/tables/test_vanco_AKIs_PC.xlsx"))
  }
  
  #### AUC, PPV and MSE of PC predictions vs test set's causality assessments ####
  
  # initiate a dataframe to store the metrics
  performance_CA = data.frame()
  
  # define a threshold for the PPV
  threshold = 0.5
  
  # for all models, calculate the AUC, PPV and MSE
  for (model in c("ERR_BART", "ERR_RF", "ERR_GLM")){
    
    # for the PPV, set the predictions to 1 if > 0.5, else 0
    ppv_naranjo = round(PPV(test_CA[["Naranjo"]], ifelse(test[test$arm == 1 & test$AKI_max == 1,][[model]] > threshold, 1, 0)), 2)
    ppv_who = round(PPV(test_CA[["WHO"]], ifelse(test[test$arm == 1 & test$AKI_max == 1,][[model]] > threshold, 1, 0)), 2)
    ppv_awdishu = round(PPV(test_CA[["Awdishu"]], ifelse(test[test$arm == 1 & test$AKI_max == 1,][[model]] > threshold, 1, 0)), 2)
    
    performance_CA = rbind(performance_CA, data.frame("Model" = str_replace(model, "ERR_", ""),
                                                      "AUC - Naranjo" = round(auc(test_CA[["Naranjo"]], test[test$arm == 1 & test$AKI_max == 1,][[model]]), 2),
                                                      "AUC - WHO" = round(auc(test_CA[["WHO"]], test[test$arm == 1 & test$AKI_max == 1,][[model]]), 2),
                                                      "AUC - Awdishu" = round(auc(test_CA[["Awdishu"]], test[test$arm == 1 & test$AKI_max == 1,][[model]]), 2),
                                                      "PPV - Naranjo" = ppv_naranjo,
                                                      "PPV - WHO" = ppv_who,
                                                      "PPV - Awdishu" = ppv_awdishu,
                                                      "MSE - Awdishu" = round(mean((test[test$arm == 1 & test$AKI_max == 1,][[model]] - test_CA[["Awdishu_prob"]])^2), 2),
                                                      check.names = F))
  }
  
  # set rownames as first col
  rownames(performance_CA) = performance_CA[,1]
  performance_CA = performance_CA[,-1]
  
  # save performance vs causality assessments
  # only if not bootstrap run
  if(!bootstrap){
    save(performance_CA, file = paste0(output_folder, "/RData/performance_CA_point.RData"))
  }
  
  #### correlations between PC predictions ####
  
  # dataframe with the correlations
  corrs = as.data.frame(cor(test[test$arm == 1 & test$AKI_max == 1, c("ERR_BART", "ERR_RF", "ERR_GLM")], method = "pearson"))
  
  # edit the colnames and rownames
  colnames(corrs) = c("BART", "RF", "LR")
  rownames(corrs) = c("BART", "RF", "LR")
  
  # save
  # only if not bootstrap run
  if(!bootstrap){
    save(corrs, file = paste0(output_folder, "/RData/PC_corrs_point.RData"))
  }
  
  #### save outcome model results ####
  
  # set first col as rowname
  rownames(effects_ATE) = effects_ATE[,1]
  effects_ATE = effects_ATE[,-1]
  rownames(effects_ATT) = effects_ATT[,1]
  effects_ATT = effects_ATT[,-1]
  rownames(effects_ATTY) = effects_ATTY[,1]
  effects_ATTY = effects_ATTY[,-1]
  rownames(effects_PC) = effects_PC[,1]
  effects_PC = effects_PC[,-1]
  
  # only if not bootstrap run
  if(!bootstrap){
    
    # save
    save(effects_ATE, file = paste0(output_folder, "/RData/effects_ATE_point.RData"))
    save(effects_ATT, file = paste0(output_folder, "/RData/effects_ATT_point.RData"))
    save(effects_ATTY, file = paste0(output_folder, "/RData/effects_ATTY_point.RData"))
    save(effects_PC, file = paste0(output_folder, "/RData/effects_PC_point.RData"))
  }
    
  #### balance assessment by SMDs and variance ratios ####
    
  # only if not bootstrap run
  if(!bootstrap){
    # cutoffs: SMD < 0.25, var < 2; https://doi.org/10.1016%2Fj.jclinepi.2013.01.013
    
    # make the weight objects
    weights_GLM = weightit(formula = as.formula(paste("arm~", paste(selected_variables, collapse="+",sep=""), sep = "")), data = results, ps = results$propensity_GLM_cal)
    weights_RF = weightit(formula = as.formula(paste("arm~", paste(selected_variables, collapse="+",sep=""), sep = "")), data = results, ps = results$propensity_RF_cal)
    weights_BART = weightit(formula = as.formula(paste("arm~", paste(selected_variables, collapse="+",sep=""), sep = "")), data = results, ps = results$propensity_BART_cal)
    
    # create the balance tables
    bal_GLM = bal.tab(weights_GLM, stats = c("mean.diffs", "variance.ratios"), thresholds = c(m = 0.25, v = 2), binary = "std", un = TRUE)
    bal_RF = bal.tab(weights_RF, stats = c("mean.diffs", "variance.ratios"), thresholds = c(m = 0.25, v = 2), binary = "std", un = TRUE)
    bal_BART = bal.tab(weights_BART, stats = c("mean.diffs", "variance.ratios"), thresholds = c(m = 0.25, v = 2), binary = "std", un = TRUE)
    
    # create and save the plots
    png(paste0(output_folder, "/plots/balance/balance_GLM.png"), height = 800, width = 600)
    love.plot(weights_GLM, stats = c("mean.diffs", "variance.ratios"), thresholds = c(m = .25, v = 2), binary = "std")
    dev.off()
    
    png(paste0(output_folder, "/plots/balance/balance_RF.png"), height = 800, width = 600)
    love.plot(weights_RF, stats = c("mean.diffs", "variance.ratios"), thresholds = c(m = .25, v = 2), binary = "std")
    dev.off()
    
    png(paste0(output_folder, "/plots/balance/balance_BART.png"), height = 800, width = 600)
    love.plot(weights_BART, stats = c("mean.diffs", "variance.ratios"), thresholds = c(m = .25, v = 2), binary = "std")
    dev.off()
  }
  
  #### causal effect plots: ARDs ####
  
  # only if not bootstrap run
  if(!bootstrap){
  
    # set ARD to BART estimate
    results$ARD = results$ARD_BART
    
    # mean ARD
    results_mean = results %>% group_by(arm) %>% summarise(ARD_mean = mean(ARD))
    
    # ARD per arm
    ARD_plot = ggplot(results, aes(x = ARD, fill = arm)) +
      geom_density(alpha=.25) + xlab("ARD") + 
      ylab("Density") + guides(fill = guide_legend(title = "Arm")) + 
      scale_fill_discrete(labels = c(paste0("Alternative: ", round(results_mean$ARD_mean[1],2)), paste0("Vancomycin: ", round(results_mean$ARD_mean[2],2)))) +
      ggtitle("Absolute risk difference, stratified by treatment arm") +
      geom_vline(data = results_mean, aes(xintercept = ARD_mean, color = arm), show.legend = F)
    
    ggsave(paste0(output_folder, "/plots/effects/ARD.png"))
    
    # ARD across AKI risk factors: only for categorical ones
    binary = selected_variables[apply(results[,selected_variables], 2, FUN = function(x){length(unique(x))}) == 2]
    binary_ARDs = results %>% select(ARD, binary) %>% gather(binary, value, -ARD) %>% group_by(binary, value) %>% summarise(ARD = mean(ARD)) %>% group_by(binary) %>% mutate(diff = diff(ARD))
    
    # create risk factor levels
    binary_ARDs$binary = factor(binary_ARDs$binary, levels = arrange(binary_ARDs, diff) %>% pull(binary) %>% unique())
    binary_ARDs$value = as.factor(binary_ARDs$value)
    levels(binary_ARDs$value) = c("0", "1")
    
    binary_plot = ggplot(binary_ARDs, aes(x = ARD, y = binary)) +
      geom_point(size = 2.5, aes(color = value)) +
      ggtitle("Absolute risk difference, stratified by AKI risk factors") +
      xlab("Absolute risk difference") +
      ylab("AKI risk factor") +
      scale_colour_discrete("Value")
    
    ggsave(paste0(output_folder, "/plots/effects/ARD_binary_AKI_risk_factors.png"))
    
    # ARD across age
    age_plot = ggplot(results, aes(x = age_old, y = ARD)) +
      geom_point() +
      geom_smooth(method = "lm") + 
      ggtitle("Absolute risk difference across age") +
      xlab("Age") +
      ylab("Absolute risk difference")
    
    ggsave(paste0(output_folder, "/plots/effects/ARD_age.png"))
    
    # ARD across age groups
    age_group_plot = results %>% mutate(age_groups = cut(age_old, breaks = c(17.9, 40, 64, 80, 100))) %>% 
        ggplot(aes(x = age_groups, y = ARD, fill = "red")) +
        stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.5), width = 0.2) +
        geom_boxplot(width = 0.5, position = position_dodge(width = 0.5)) +
        theme(legend.position = "none") +
        xlab("Age group") + ylab("Absolute risk difference") +
        ggtitle("Absolute risk difference across age groups") +
        scale_x_discrete(labels = c("18-40", "41-64", "65-80", "81-100"))
    
    ggsave(paste0(output_folder, "/plots/effects/ARD_age_groups.png"))
    
    # ARD across age groups - ATT
    age_group_plot = results %>% subset(arm == 1) %>% mutate(age_groups = cut(age_old, breaks = c(17.9, 40, 64, 80, 100))) %>% 
      ggplot(aes(x = age_groups, y = ARD, fill = "red")) +
      stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.5), width = 0.2) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 0.5)) +
      theme(legend.position = "none") +
      xlab("Age group") + ylab("Absolute risk difference") +
      ggtitle("Absolute risk difference across age groups") +
      scale_x_discrete(labels = c("18-40", "41-64", "65-80", "81-100"))
    
    ggsave(paste0(output_folder, "/plots/effects/ARD_age_groups_ATT.png"))
    
    # box plot for sex
    sex_plot = results %>% mutate(gender = ifelse(gender == 1, "Male", "Female")) %>% 
      ggplot(aes(x = gender, y = ARD, fill = "red")) +
      stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.5), width = 0.2) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 0.5)) +
      theme(legend.position = "none") +
      xlab("Sex") + ylab("Absolute risk difference") +
      ggtitle("Absolute risk difference stratified by sex")
      
    ggsave(paste0(output_folder, "/plots/effects/ARD_sex.png"))
    
    # box plot for sex - ATT
    sex_plot = results %>% subset(arm == 1) %>% mutate(gender = ifelse(gender == 1, "Male", "Female")) %>% 
      ggplot(aes(x = gender, y = ARD, fill = "red")) +
      stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.5), width = 0.2) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 0.5)) +
      theme(legend.position = "none") +
      xlab("Sex") + ylab("Absolute risk difference") +
      ggtitle("Absolute risk difference stratified by sex")
    
    ggsave(paste0(output_folder, "/plots/effects/ARD_sex_ATT.png"))
    
    # box plot for CKD
    CKD_plot = results %>% mutate(chronic_kidney_disease = ifelse(chronic_kidney_disease == 1, "CKD", "No CKD")) %>% 
      ggplot(aes(x = chronic_kidney_disease, y = ARD, fill = "red")) +
      stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.5), width = 0.2) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 0.5)) +
      theme(legend.position = "none") +
      xlab("CKD") + ylab("Absolute risk difference") +
      ggtitle("Absolute risk difference stratified by CKD")
    
    ggsave(paste0(output_folder, "/plots/effects/ARD_CKD.png"))
    
    # box plot for CKD - ATT
    CKD_plot = results %>% subset(arm == 1) %>% mutate(chronic_kidney_disease = ifelse(chronic_kidney_disease == 1, "CKD", "No CKD")) %>% 
      ggplot(aes(x = chronic_kidney_disease, y = ARD, fill = "red")) +
      stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.5), width = 0.2) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 0.5)) +
      theme(legend.position = "none") +
      xlab("CKD") + ylab("Absolute risk difference") +
      ggtitle("Absolute risk difference stratified by CKD")
    
    ggsave(paste0(output_folder, "/plots/effects/ARD_CKD_ATT.png"))
  }
    
  #### causal effect plots: PCs for vanco admissions with AKI ####
  
  # only if not bootstrap run
  if(!bootstrap){
  
    # set ERR to BART estimate
    results$ERR = results$ERR_BART
    
    # ERR for exposed admissions with AKI
    ERR = results %>% subset(arm == 1 & AKI_max == 1)
    
    # ERR: density
    ERR_plot = ggplot(ERR, aes(x = ERR)) +
      geom_density(alpha = .25, color = "blue", fill = "lightblue") + xlab("Probability of causation") + 
      ylab("Density") + 
      ggtitle("Probability of causation\nFor admissions exposed to vancomycin with subsequent AKI") +
      xlim(c(0,1))
    
    ggsave(paste0(output_folder, "/plots/effects/PC.png"))
    
    # ERR: histogram
    ERR_plot = ggplot(ERR, aes(x = ERR)) +
      geom_histogram(alpha = .25, color = "blue", fill = "lightblue", binwidth = 0.1, center = 0.05) + xlab("Probability of causation") + 
      ylab("Frequency") + 
      scale_x_continuous(expand = c(0.01,0.03), breaks = seq(0.00, 1.00, 0.1)) +
      ylim(c(0,40))
    
    ggsave(paste0(output_folder, "/plots/effects/PC_hist.png"), height = 4, width = 4)
    
    # ERR: histograms for all models
    ERR_plot = results %>% subset(arm == 1 & AKI_max == 1) %>% mutate(BART = ERR_BART, RF = ERR_RF, LR = ERR_GLM) %>% select(BART, RF, LR) %>% gather(model, ERR) %>% mutate(model = factor(model, levels = c("BART", "RF", "LR"))) %>%
                      ggplot(aes(x = ERR)) +
                      geom_histogram(alpha = .25, color = "blue", fill = "lightblue", binwidth = 0.1, center = 0.05) + xlab("Probability of causation") + 
                      ylab("Frequency") + 
                      scale_x_continuous(expand = c(0.01,0.03), breaks = seq(0.00, 1.00, 0.1)) +
                      ylim(c(0,40)) +
                      facet_wrap(~model)
      
    ggsave(paste0(output_folder, "/plots/effects/PC_hist_all.png"), height = 4, width = 12)
    
    # ERR: histograms for all models with alternative order
    ERR_plot = results %>% subset(arm == 1 & AKI_max == 1) %>% mutate(BART = ERR_BART, RF = ERR_RF, LR = ERR_GLM) %>% select(BART, RF, LR) %>% gather(model, ERR) %>% mutate(model = factor(model, levels = c("BART", "LR", "RF"))) %>%
      ggplot(aes(x = ERR)) +
      geom_histogram(alpha = .25, color = "blue", fill = "lightblue", binwidth = 0.1, center = 0.05) + xlab(bquote(PC[low])) + 
      ylab("Frequency") + 
      scale_x_continuous(expand = c(0.01,0.03), breaks = seq(0.0, 1.0, 0.2)) +
      ylim(c(0,40)) +
      facet_wrap(~model) + theme(text = element_text(size = 16))
    
    ggsave(paste0(output_folder, "/plots/effects/PC_hist_all_order.png"), height = 4, width = 12)
    
    # ERR across AKI risk factors: only for categorical ones
    binary = selected_variables[apply(results[,selected_variables], 2, FUN = function(x){length(unique(x))}) == 2]
    binary_ERRs = ERR %>% select(ERR, binary) %>% gather(binary, value, -ERR) %>% group_by(binary, value) %>% summarise(ERR = mean(ERR)) %>% group_by(binary) %>% mutate(diff = ifelse(length(value) == 2, diff(ERR),0))
    
    # create risk factor levels
    binary_ERRs$binary = factor(binary_ERRs$binary, levels = arrange(binary_ERRs, diff) %>% pull(binary) %>% unique())
    binary_ERRs$value = as.factor(binary_ERRs$value)
    levels(binary_ERRs$value) = c("0", "1")
    
    binary_plot = ggplot(binary_ERRs, aes(x = ERR, y = binary)) +
      geom_point(size = 2.5, aes(color = value)) +
      ggtitle("Probability of causation stratified by AKI risk factors\nFor admissions exposed to vancomycin with subsequent AKI") +
      theme(plot.title = element_text(size = 11)) +
      xlab("Probability of causation") +
      ylab("AKI risk factor") +
      scale_colour_discrete("Value")
    
    ggsave(paste0(output_folder, "/plots/effects/PC_binary_AKI_risk_factors.png"))
    
    # ERR across age
    age_plot = ggplot(ERR, aes(x = age_old, y = ERR)) +
      geom_point() +
      geom_smooth(method = "lm") + 
      ggtitle("Probability of causation across age\nFor admissions exposed to vancomycin with subsequent AKI") +
      xlab("Age") +
      ylab("Probability of causation")
    
    ggsave(paste0(output_folder, "/plots/effects/PC_age.png"))
    
    # ERR across age groups
    age_group_plot = ERR %>% mutate(age_groups = cut(age_old, breaks = c(17.9, 40, 64, 80, 100))) %>% 
        ggplot(aes(x = age_groups, y = ERR, fill = "red")) +
        stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.5), width = 0.2) +
        geom_boxplot(width = 0.5, position = position_dodge(width = 0.5)) +
        theme(legend.position = "none") +
        xlab("Age group") + ylab("Probability of causation") +
        ggtitle("Probability of causation across age groups\nFor admissions exposed to vancomycin with subsequent AKI") +
        scale_x_discrete(labels = c("18-40", "41-64", "65-80", "81-100"))
    
    ggsave(paste0(output_folder, "/plots/effects/PC_age_groups.png"))
    
    # box plot for sex
    sex_plot = ERR %>% mutate(gender = ifelse(gender == 1, "Male", "Female")) %>% 
      ggplot(aes(x = gender, y = ERR, fill = "red")) +
      stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.5), width = 0.2) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 0.5)) +
      theme(legend.position = "none") +
      xlab("Sex") + ylab("Probability of causation") +
      ggtitle("Probability of causation stratified by sex\nFor admissions exposed to vancomycin with subsequent AKI")
    
    ggsave(paste0(output_folder, "/plots/effects/PC_sex.png"))
    
    # box plot for CKD
    CKD_plot = ERR %>% mutate(chronic_kidney_disease = ifelse(chronic_kidney_disease == 1, "CKD", "No CKD")) %>% 
      ggplot(aes(x = chronic_kidney_disease, y = ERR, fill = "red")) +
      stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.5), width = 0.2) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 0.5)) +
      theme(legend.position = "none") +
      xlab("CKD") + ylab("Probability of causation") +
      ggtitle("Probability of causation stratified by CKD\nFor admissions exposed to vancomycin with subsequent AKI")
    
    ggsave(paste0(output_folder, "/plots/effects/PC_CKD.png"))
  }
  
  #### timing ####
  
  end = Sys.time()
  elapsed = end - start
  
  #### return models and data if not bootstrap ####
  
  # return the env
  if(!bootstrap){
   
    env = list2env(mget(ls()))
    return(env)
     
  }
  
  #### return bootstrap statistics of causal effects and model performance ####
  
  if(bootstrap){
    
    # convert to matrix
    effects_ATE = as.matrix(effects_ATE)
    effects_ATT = as.matrix(effects_ATT)
    effects_ATTY = as.matrix(effects_ATTY)
    effects_PC = as.matrix(effects_PC)
    performance = as.matrix(performance)
    performance_CA = as.matrix(performance_CA)
    corrs = as.matrix(corrs)
    
    output = list(effects_ATE, effects_ATT, effects_ATTY, effects_PC, performance, performance_CA, corrs)
    return(output)
  }
  

  
}

#### function for plotting effects of different models ####
effects_plot = function(data, estimand, center, output){
  
  library(ggplot2)
  
  estimates = data.frame(Estimator = names(data), point = NA, lower = NA, upper = NA)
  
  # split the estimates and CIs
  for (estimate in 1:length(data)){
    estimates[estimate, c("point", "lower", "upper")] = do.call(rbind, str_split(str_replace(str_replace(str_replace_all(str_replace(str_replace(data[estimate], " - ", "#"), " ", ""), "\\*", ""), "\\)", "#"), "\\(", "#"), "#"))[,1:3]  
  }
  
  # convert to numeric
  for (col in c("point", "lower", "upper")){
    estimates[[col]] = as.numeric(estimates[[col]])
  }
  
  # estimator to factor
  estimates$Estimator = factor(estimates$Estimator, levels = estimates$Estimator)
  
  # plot
  plot = ggplot(estimates, aes(x = Estimator, y = point, ymin = lower, ymax = upper)) + 
          geom_hline(aes(yintercept = center), lty = 1) +
          geom_point(size = 3, shape = 21, fill = "#008fd5", colour = "black", stroke = 1) +
          scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
          coord_flip() +
          xlab("Estimator") +
          ylab(estimand) +
          geom_errorbar(aes(ymin = lower, ymax = upper), width = .2, position = position_dodge(.9)) 
        
  ggsave(output)
  
}
