###
#
# VAKI_diagnosis.R
#
###

#### set the input, code and output folders ####

input_folder_data = ""
input_folder_misc = ""
code_folder = ""
output_folder = ""

#### load the libraries and functions ####

library(openxlsx)
library(plyr)
library(tidyverse)
library(coxed)
library(stringr)
source(paste0(code_folder, "VAKI_diagnosis_functions.R"))

#### load the data ####

# load the IC names
IC_names = read.xlsx(paste0(input_folder_misc, "ICs.xlsx"))

# load the variable names
variables = read.xlsx(paste0(input_folder_misc, "variables.xlsx"))

# load the train and test data
train = loadRData(paste0(input_folder_data, "data_imputed_train.RData"))
test = loadRData(paste0(input_folder_data, "data_imputed_test.RData"))

# load the causality assessments
test_CA = read.xlsx(paste0(output_folder, "/tables/test_vanco_AKIs_gold_prob.xlsx"))
test_CA = arrange(test_CA, ID)

# set causality assessment outcome categories 'unassessible', 'possible', 'probable' and 'nearly certain' to 1 (coded as "yes - vancomycin (may have) caused the AKI" in the excel sheet)
for (col in 5:7){test_CA[[col]] = ifelse(test_CA[[col]] == "yes - vancomycin (may have) caused the AKI", 1, 0)}

#### data preprocessing ####

# select the variables for vancomycin
variables = variables[!is.na(variables$name) & variables$name == "Vancomycin",]

# preprocess all sets
for (set in c("train", "test")){
  
  # print status
  print(paste0("Working on ", set , " set..."))
  
  # get the data
  tmp = get(set)
  
  # only keep relevant variables
  tmp = tmp %>% subset(variable %in% c(variables$variables, "AKI_SCr_static"))
  
  # convert dosages to binary
  tmp[grepl("_DDDs_sum", tmp$variable) & tmp$value > 0,]$value = 1
  
  # rename vars with _DDDs_sum to _max
  tmp$variable = str_replace(tmp$variable, "_DDDs_sum", "_max")
  
  # keep and summarize the baseline characteristics and outcome (periods up to 0 for baseline vars)
  tmp = tmp %>% subset(period < 1 | variable == "AKI_max") %>% group_by(ID, nice_adm_icu, nice_dis_icu, variable) %>% summarise(min = min(value), max = max(value))
  
  # for variables with _min: keep min, else keep max
  tmp = tmp %>% ungroup() %>% mutate(value = ifelse(grepl("_min", variable), min, max)) %>% select(ID, nice_adm_icu, nice_dis_icu, variable, value)
  
  # spread the variables 
  tmp = tmp %>% spread(variable, value)
  
  # edit the AKI SCr timestamps: make POSIX
  tmp[["AKI_SCr_static"]] = as.POSIXct(tmp[["AKI_SCr_static"]], origin = "1970-01-01")
  
  # summarise sepsis and sepsis_max by taking max
  tmp$sepsis = as.numeric(rowSums(tmp[,c("sepsis", "sepsis_max")]) > 0)
  
  # store numeric AKI and arm
  tmp$arm_num = tmp$arm
  tmp$AKI_max_num = tmp$AKI_max
  
  # convert arm, AKI and ICU to factors
  tmp$arm = as.factor(tmp$arm)
  tmp$AKI_max = as.factor(tmp$AKI_max)
  tmp$ICU = as.factor(tmp$ICU)
  
  # convert to dataframe
  tmp = as.data.frame(tmp)
  
  # assign to set
  assign(set, tmp)
  
}

# rename vars with _DDDs_sum to _max in variables index
variables$variables = str_replace(variables$variables, "_DDDs_sum", "_max")

#### make the output folders ####

if(!dir.exists(paste0(output_folder, "/tables/"))){dir.create(paste0(output_folder, "/tables/"))}
if(!dir.exists(paste0(output_folder, "/plots/"))){dir.create(paste0(output_folder, "/plots/"))}
if(!dir.exists(paste0(output_folder, "/RData/"))){dir.create(paste0(output_folder, "/RData/"))}
if(!dir.exists(paste0(output_folder, "/plots/calibration/"))){dir.create(paste0(output_folder, "/plots/calibration/"))}
if(!dir.exists(paste0(output_folder, "/plots/effects/"))){dir.create(paste0(output_folder, "/plots/effects/"))}
if(!dir.exists(paste0(output_folder, "/plots/balance/"))){dir.create(paste0(output_folder, "/plots/balance/"))}
if(!dir.exists(paste0(output_folder, "/plots/LR_assumptions/"))){dir.create(paste0(output_folder, "/plots/LR_assumptions/"))}

#### run the modeling, plotting and bootstrap ####

# define seed
seed = 123

# modeling
output = modeling(train, test, test_CA, variables, output_folder, IC_names, seed)

# bootstrap
R = 500
results = cluster_bootstrap(train, R, "ICU", modeling, test, test_CA, variables, output_folder, IC_names, seed)

# save the results
save(results, file = paste0(output_folder, "/RData/bootstrap_results.RData"))

#### process and save the estimates and performance ####

# obtain BCA CI; for effects-ATE
effects_ATE = lapply(results, first)
effects_ATE = array(unlist(effects_ATE), dim = c(7, 5, R))
effects_ATE_CI = apply(effects_ATE, c(1,2), function(x){paste0(round(bca(x), 2), collapse = "$")})

# obtain BCA CI; for effects-ATT
effects_ATT = lapply(results, nth, n = 2)
effects_ATT = array(unlist(effects_ATT), dim = c(7, 5, R))
effects_ATT_CI = apply(effects_ATT, c(1,2), function(x){paste0(round(bca(x), 2), collapse = "$")})

# obtain BCA CI; for effects-ATTY
effects_ATTY = lapply(results, nth, n = 3)
effects_ATTY = array(unlist(effects_ATTY), dim = c(3, 5, R))
effects_ATTY_CI = apply(effects_ATTY, c(1,2), function(x){paste0(round(bca(x), 2), collapse = "$")})

# obtain BCA CI; for effects-PC
effects_PC = lapply(results, nth, n = 4)
effects_PC = array(unlist(effects_PC), dim = c(3, 3, R))
effects_PC_CI = apply(effects_PC, c(1,2), function(x){paste0(round(bca(x), 2), collapse = "$")})

# obtain BCA CI; for performances
performance = lapply(results, nth, n = 5)
performance = array(unlist(performance), dim = c(9, 3, R))
performance_CI = apply(performance, c(1,2), function(x){paste0(round(bca(x), 2), collapse = "$")})

# obtain BCA CI; for performances for causality assessments
performance_CA = lapply(results, nth, n = 6)
performance_CA = array(unlist(performance_CA), dim = c(3, 7, R))
performance_CA_CI = apply(performance_CA, c(1,2), function(x){paste0(round(bca(x), 2), collapse = "$")})

# obtain BCA CI; for corrs between predicted PCs
corrs = lapply(results, last)
corrs = array(unlist(corrs), dim = c(3, 3, R))
corrs_CI = apply(corrs, c(1,2), function(x){paste0(round(bca(x), 2), collapse = "$")})

# get the point/mean effects and performance
effects_ATE_mean = as.matrix(loadRData(paste0(output_folder, "/RData/effects_ATE_point.RData")))
effects_ATT_mean = as.matrix(loadRData(paste0(output_folder, "/RData/effects_ATT_point.RData")))
effects_ATTY_mean = as.matrix(loadRData(paste0(output_folder, "/RData/effects_ATTY_point.RData")))
effects_PC_mean = as.matrix(loadRData(paste0(output_folder, "/RData/effects_PC_point.RData")))
performance_mean = apply(performance, c(1,2), function(x){mean(x)})
performance_CA_mean = apply(performance_CA, c(1,2), function(x){mean(x)})
corrs_mean = apply(corrs, c(1,2), function(x){mean(x)})

# round everything to two decimals and add CIs
effects_ATE_final = paste0(format(round(effects_ATE_mean, 2), nsmall = 2), " (", matrix(str_replace(format(as.numeric(lapply((strsplit(effects_ATE_CI,"\\$")), first)), nsmall = 2), " ", ""), nrow = 7, ncol = 5), " - ", matrix(str_replace(format(as.numeric(lapply((strsplit(effects_ATE_CI,"\\$")), last)), nsmall = 2), " ", ""), nrow = 7, ncol = 5), ")") %>% matrix(nrow = nrow(effects_ATE_mean))
effects_ATT_final = paste0(format(round(effects_ATT_mean, 2), nsmall = 2), " (", matrix(str_replace(format(as.numeric(lapply((strsplit(effects_ATT_CI,"\\$")), first)), nsmall = 2), " ", ""), nrow = 7, ncol = 5), " - ", matrix(str_replace(format(as.numeric(lapply((strsplit(effects_ATT_CI,"\\$")), last)), nsmall = 2), " ", ""), nrow = 7, ncol = 5), ")") %>% matrix(nrow = nrow(effects_ATT_mean))
effects_ATTY_final = paste0(format(round(effects_ATTY_mean, 2), nsmall = 2), " (", matrix(str_replace(format(as.numeric(lapply((strsplit(effects_ATTY_CI,"\\$")), first)), nsmall = 2), " ", ""), nrow = 3, ncol = 5), " - ", matrix(str_replace(format(as.numeric(lapply((strsplit(effects_ATTY_CI,"\\$")), last)), nsmall = 2), " ", ""), nrow = 3, ncol = 5), ")") %>% matrix(nrow = nrow(effects_ATTY_mean))
effects_PC_final = paste0(format(round(effects_PC_mean, 2), nsmall = 2), " (", matrix(str_replace(format(as.numeric(lapply((strsplit(effects_PC_CI,"\\$")), first)), nsmall = 2), " ", ""), nrow = 3, ncol = 3), " - ", matrix(str_replace(format(as.numeric(lapply((strsplit(effects_PC_CI,"\\$")), last)), nsmall = 2), " ", ""), nrow = 3, ncol = 3), ")") %>% matrix(nrow = nrow(effects_PC_mean))
performance_final = paste0(format(round(performance_mean, 2), nsmall = 2), " (", matrix(str_replace(format(as.numeric(lapply((strsplit(performance_CI,"\\$")), first)), nsmall = 2), " ", ""), nrow = 9, ncol = 3), " - ", matrix(str_replace(format(as.numeric(lapply((strsplit(performance_CI,"\\$")), last)), nsmall = 2), " ", ""), nrow = 9, ncol = 3), ")") %>% matrix(nrow = nrow(performance_mean))
performance_CA_final = paste0(format(round(performance_CA_mean, 2), nsmall = 2), " (", matrix(str_replace(format(as.numeric(lapply((strsplit(performance_CA_CI,"\\$")), first)), nsmall = 2), " ", ""), nrow = 3, ncol = 7), " - ", matrix(str_replace(format(as.numeric(lapply((strsplit(performance_CA_CI,"\\$")), last)), nsmall = 2), " ", ""), nrow = 3, ncol = 7), ")") %>% matrix(nrow = nrow(performance_CA_mean))
corrs_final = paste0(format(round(corrs_mean, 2), nsmall = 2), " (", matrix(str_replace(format(as.numeric(lapply((strsplit(corrs_CI,"\\$")), first)), nsmall = 2), " ", ""), nrow = 3, ncol = 3), " - ", matrix(str_replace(format(as.numeric(lapply((strsplit(corrs_CI,"\\$")), last)), nsmall = 2), " ", ""), nrow = 3, ncol = 3), ")") %>% matrix(nrow = nrow(corrs_mean))

# add colnames and rownames
colnames(effects_ATE_final) = colnames(effects_ATE_mean)
rownames(effects_ATE_final) = rownames(effects_ATE_mean)
colnames(effects_ATT_final) = colnames(effects_ATT_mean)
rownames(effects_ATT_final) = rownames(effects_ATT_mean)
colnames(effects_ATTY_final) = colnames(effects_ATTY_mean)
rownames(effects_ATTY_final) = rownames(effects_ATTY_mean)
colnames(effects_PC_final) = colnames(effects_PC_mean)
rownames(effects_PC_final) = rownames(effects_PC_mean)
colnames(performance_final) = colnames(results[[1]][[5]])
rownames(performance_final) = rownames(results[[1]][[5]])
colnames(performance_CA_final) = colnames(results[[1]][[6]])
rownames(performance_CA_final) = rownames(results[[1]][[6]])
colnames(corrs_final) = colnames(results[[1]][[7]])
rownames(corrs_final) = rownames(results[[1]][[7]])

# plots with bars
effects_plot(effects_ATE_final[,3], "ATE - Absolute risk difference", 0, paste0(output_folder, "/plots/effects/ATE_ARD.png"))
effects_plot(effects_ATE_final[,4], "ATE - Risk ratio", 1, paste0(output_folder, "/plots/effects/ATE_RR.png"))
effects_plot(effects_ATE_final[,5], "ATE - Excess risk ratio", 0, paste0(output_folder, "/plots/effects/ATE_ERR.png"))
effects_plot(effects_ATT_final[,3], "ATT - Absolute risk difference", 0, paste0(output_folder, "/plots/effects/ATT_ARD.png"))
effects_plot(effects_ATT_final[,4], "ATT - Risk ratio", 1, paste0(output_folder, "/plots/effects/ATT_RR.png"))
effects_plot(effects_ATT_final[,5], "ATT - Excess risk ratio", 0, paste0(output_folder, "/plots/effects/ATT_ERR.png"))
effects_plot(effects_PC_final[,1], "Probability of causation - ATE", 0, paste0(output_folder, "/plots/effects/PC_ATE.png"))
effects_plot(effects_PC_final[,2], "Probability of causation - ATT", 0, paste0(output_folder, "/plots/effects/PC_ATT.png"))
effects_plot(effects_PC_final[,3], "Probability of causation - ATTY", 0, paste0(output_folder, "/plots/effects/PC_ATTY.png"))

# save effects and performance with CIs
write.xlsx(effects_ATE_final, paste0(output_folder, "/tables/effects_ATE.xlsx"), row.names = TRUE)
write.xlsx(effects_ATT_final, paste0(output_folder, "/tables/effects_ATT.xlsx"), row.names = TRUE)
write.xlsx(effects_PC_final, paste0(output_folder, "/tables/effects_PC.xlsx"), row.names = TRUE)
write.xlsx(performance_final, paste0(output_folder, "/tables/performance.xlsx"), row.names = TRUE)
write.xlsx(performance_CA_final, paste0(output_folder, "/tables/performance_CA.xlsx"), row.names = TRUE)
write.xlsx(corrs_final, paste0(output_folder, "/tables/PC_corrs.xlsx"), row.names = TRUE)

#### save modeling output and session ####

# modeling environment
save(list = ls(output), file = paste0(output_folder, "/RData/session_modeling.RData"), envir = output)

# session
save.image(paste0(output_folder, "/RData/session.RData"))

