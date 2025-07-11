## -- interaction analysis presented in Figure 6A to evaluate the interaction term between night shift and the Polygenic Risk scores

# -- data -- this is the polygenic risk scores, phenotype file from DeepPheWAS package and a file to subset the phenotypes
phe=read.table("PHE_SCORE_noACCEL_inclCOV_inclSCORE.txt",header=T)
dat1=read.table("pheno_data_ready/SCORE_POPULATION_FINAL_PHENO.txt", header=T, sep="\t")
sum.tab1=read.table("summary_table_for_SCORE_set_gt20.txt",sep="\t",header=T)

##### -------------------------------------------------  creating R-objects restriced to specfic sample sets:
day=phe.dat[phe.dat$night_shift %in% c(0),]
some=phe.dat[phe.dat$night_shift %in% c(1),]
night=phe.dat[phe.dat$night_shift %in% c(2,3),]
ever=phe.dat[phe.dat$night_shift %in% c(1,2,3), ]
phe.dat$ever=rep(NA,nrow(phe.dat))
phe.dat$ever[phe.dat$night_shift %in% c(0)]=c(0)
phe.dat$ever[phe.dat$night_shift %in% c(1,2,3)]=c(1)
table(phe.dat$ever)
table(is.na(phe.dat$ever))
phe.dat$night=rep(NA,nrow(phe.dat))
phe.dat$night[phe.dat$night_shift %in% c(0)]=c(0)
phe.dat$night[phe.dat$night_shift %in% c(1)]=c(1)
phe.dat$night[phe.dat$night_shift %in% c(2)]=c(2)
phe.dat$night[phe.dat$night_shift %in% c(3)]=c(3)

####    -------------------------------------   sumstat wrangling function
glm_summary_stats <- function(formula, data, family, id, set_name, type) {
  
  # Ensure family is a valid GLM family object
  if (!inherits(family, "family")) {
    stop("Error: 'family' must be a valid GLM family function (e.g., binomial(), gaussian()).")
  }
  
  # Run the GLM model
  glm_model <- glm(formula, data = data, family = family)
  coeff <- summary(glm_model)$coefficients[2,]
  
  # Extract values
  estimate <- coeff[1]
  std_error <- coeff[2]
  z_value <- coeff[3]
  p_value <- coeff[4]
  
  # Handle binomial and gaussian cases
  if (family$family == "binomial") {
    or <- exp(estimate)  # Odds Ratio
    ci_low <- exp(estimate - 1.96 * std_error)  # Lower bound of 95% CI
    ci_up <- exp(estimate + 1.96 * std_error)  # Upper bound of 95% CI
    
    # Extract outcome variable name from formula
    outcome_var <- all.vars(formula)[1]  
    
    # Compute case and control counts
    cases <- sum(data[[outcome_var]] == 1, na.rm = TRUE)
    controls <- sum(data[[outcome_var]] == 0, na.rm = TRUE)
    
  } else if (family$family == "gaussian") {
    or <- NA
    ci_low <- NA
    ci_up <- NA
    cases <- NA
    controls <- NA
  }
  
  # Prepare final named result vector with correct column names
  result <- c(
    ID = id,
    outcome = all.vars(formula)[1],  # Extract outcome name dynamically
    type = type,
    set = set_name,
    Estimate = estimate,
    `Std. Error` = std_error,
    `z value` = z_value,
    `Pr(>|z|)` = p_value,
    `Odds Ratio` = or,
    `CI Lower Bound` = ci_low,
    `CI Upper Bound` = ci_up,
    Cases = cases,
    Controls = controls
  )
  
  return(result)
}

###  --------------  function for interactioh analysis
glm_interaction_stats <- function(formula, data, family, id, set_name, type, interaction_term) {
  
  # Ensure family is a valid GLM family object
  if (!inherits(family, "family")) {
    stop("Error: 'family' must be a valid GLM family function (e.g., binomial(), gaussian()).")
  }
  
  # Run GLM model
  glm_model <- glm(formula, data = data, family = family)
  interaction_coeffs <- summary(glm_model)$coefficients
  
  # Check if interaction term exists in the model output
  interaction_name <- paste0(id, ":", interaction_term)
  
  if (interaction_name %in% rownames(interaction_coeffs)) {
    coeff <- interaction_coeffs[interaction_name,]
    
    estimate <- coeff[1]
    std_error <- coeff[2]
    z_value <- coeff[3]
    p_value <- coeff[4]
  } else {
    # If the interaction term is not found, set values to NA
    estimate <- NA
    std_error <- NA
    z_value <- NA
    p_value <- NA
  }
  
  # Since OR, CI, and case/control numbers are not applicable, set them to NA
  or <- NA
  ci_low <- NA
  ci_up <- NA
  cases <- NA
  controls <- NA
  
  # Prepare output vector with correct column names
  result <- c(
    ID = id,
    outcome = all.vars(formula)[1],  # Extract outcome variable dynamically
    type = type,
    set = set_name,
    Estimate = estimate,
    `Std. Error` = std_error,
    `z value` = z_value,
    `Pr(>|z|)` = p_value,
    `Odds Ratio` = or,
    `CI Lower Bound` = ci_low,
    `CI Upper Bound` = ci_up,
    Cases = cases,
    Controls = controls
  )
  
  return(result)
}
###   --------------------------- 
k=c("sc_F2_score")
phenotypes = as.character(pheno.bin.sig$ID[pheno.bin.sig$sc_F2_score %in% c(1)])
rm(pgs.out)
for (t in phenotypes) {
  print (paste ("++++++++ outcome is ", t, "now +++++++++"))
  
  formula1 <- as.formula(paste(t, "~", k, "+ age + sex + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  
  day.add2 =  glm_summary_stats(formula = formula1, data = day, family = binomial(), 
                                id = k, set_name = "day", type = "sens")
  ## ----- ever night
  ever.add2 =  glm_summary_stats(formula = formula1, data = ever, family = binomial(), 
                                 id = k, set_name = "ever nightshift", type = "sens")
  ## ----    night full
  
  some.add2=glm_summary_stats(formula = formula1, data = some, family = binomial(), 
                              id = k, set_name = "some nightshift", type = "sens")
  ##--- night
  night.add2=glm_summary_stats(formula = formula1, data = night, family = binomial(), 
                               id = k, set_name = "usually_always nigthshift", type = "sens")
  ###   --------------------------------------------------------------------------------------
  ### --------- interaction analysis: ever
  i.formula <- as.formula(paste0( t,"  ~ " ,k,"*ever"  , " + age + sex + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
  ever.interact2 <- glm_interaction_stats(formula = i.formula, data = phe.dat, family = binomial(),
                                          id = k, set_name = "ever nightshift", type = "interaction", interaction_term = "ever")
  ###   -----------------------------------
  ### --------- interaction analysis: night
  i.formula2 <- as.formula(paste0( t,"  ~ " ,k,"*night"  , " + age + sex  + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
  
  night.interact2 <- glm_interaction_stats(formula = i.formula2, data = phe.dat, family = binomial(),
                                           id = k, set_name = "usually_always nightshift", type = "interaction", interaction_term = "night")
  
  ###  --- make output
  if (!exists("pgs.out")){
    pgs.out=rbind(day.add2,ever.add2,night.add2,some.add2,ever.interact2,night.interact2)
  }else {
    pgs.out=as.data.frame(rbind(pgs.out,day.add2,ever.add2,night.add2,some.add2,ever.interact2,night.interact2))
  }
} # pgs loop

fin.inter.out2=merge(pgs.out,anno,by.x="outcome",by.y="PheWAS_ID",all.x=T )

options(scipen = 999)
write.table(fin.inter.out2,file=paste0("Interaction_",k,"_sigOnly_dichotomous_only.txt"),sep="\t",col.names=T,row.names=F,quote=F)
