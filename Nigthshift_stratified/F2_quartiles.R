## this script is supposed to calculate quartiles form the F2* score and evaluate the association of these quartiles to clinically relevant outcomes

# -- data -- this is the polygenic risk scores, phenotype file from DeepPheWAS package and a file to subset the phenotypes
phe=read.table("PHE_SCORE_noACCEL_inclCOV_inclSCORE.txt",header=T)
dat1=read.table("pheno_data_ready/SCORE_POPULATION_FINAL_PHENO.txt", header=T, sep="\t")
sum.tab1=read.table("summary_table_for_SCORE_set_gt20.txt",sep="\t",header=T)



# Create a new variable with quartile labels (1 = lowest, 4 = highest)
library(dplyr)
phe.dat <- phe.dat %>%
  mutate(sc_F2inclPSYCH_quartile = ntile(sc_F2inclPSYCH_score, 4))

#### ----------------  stratify night shift

###   -------------------------------------------  ---- day ---  some -- night 
phe.dat$night[phe.dat$night_shift %in% c(0)]=c("day")
phe.dat$night[phe.dat$night_shift %in% c(1)]=c("some")
phe.dat$night[phe.dat$night_shift %in% c(2)]=c("night")
phe.dat$night[phe.dat$night_shift %in% c(3)]=c("night")
### ----------------------------------------------------------  ever
phe.dat$night[phe.dat$night_shift %in% c(0)]=c("day")
phe.dat$night[phe.dat$night_shift %in% c(1)]=c("ever")
phe.dat$night[phe.dat$night_shift %in% c(2)]=c("ever")
phe.dat$night[phe.dat$night_shift %in% c(3)]=c("ever")
### -----------------------------------------------------------------------------
###  ----------------------- have a qquuick look
hist(phe.dat$sc_F2inclPSYCH_score)
table(phe.dat$sc_F2inclPSYCH_quartile)
table(phe.dat$night, phe.dat$P2035.2, phe.dat$sc_F2inclPSYCH_quartile)

## -----------------------  analysis
library(dplyr)
library(broom)
library(purrr)
library(tibble)
library(stringr)

##    --------- prep
k=c("sc_F2inclPSYCH")
###
phe.dat <- phe.dat %>%
  mutate(
    PRS_quartile = factor(sc_F2inclPSYCH_quartile, levels = 1:4),
    Shift_work = factor(night, levels = c("day","ever")),
    strata = interaction(PRS_quartile, Shift_work, sep = "_")
  )
###
phe.dat$strata <- relevel(phe.dat$strata, ref = "1_day")
# Create output container
results_all <- list()

# Loop over traits 
for (i in seq_len(nrow(trait_table))) {
  trait_code <- trait_table$plot.trait[i]
  trait_label <- trait_table$row_labels[i]
  
  # Skip if variable doesn't exist in dataset
  if (!trait_code %in% names(phe.dat)) next
  
  # Subset to relevant rows and remove missing values
  subdat <- phe.dat %>%
    select(
      FID, strata, !!sym(trait_code),
      age, sex, centre,
      PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10
    ) %>%
    filter(!is.na(!!sym(trait_code)), !is.na(strata)) %>%
    mutate(strata = droplevels(factor(strata)))  # ensure correct factor
  
  
  # Fit logistic model
  model <- glm(as.formula(paste(trait_code, "~ strata + age + sex  + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data = subdat, family = binomial)
  
  
  # Get log-odds (β), SE, z-values
  model_log <- tidy(model, conf.int = FALSE, exponentiate = FALSE) %>%
    filter(term != "(Intercept)") %>%
    mutate(strata = gsub("strata", "", term)) %>%
    select(strata, beta = estimate, SE = std.error, z = statistic, p.value)
  
  # Get odds ratios and confidence intervals
  model_or <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(strata = gsub("strata", "", term)) %>%
    select(strata, OR = estimate, CI_lower = conf.low, CI_upper = conf.high)
  
  # Merge log-odds and OR output
  model_tidy <- left_join(model_log, model_or, by = "strata")
  
  # Add the reference row (1_day)
  ref_row <- tibble(
    strata = "1_day",
    beta = NA_real_,
    SE = NA_real_,
    z = NA_real_,
    p.value = NA_real_,
    OR = 1,
    CI_lower = NA_real_,
    CI_upper = NA_real_
  )
  
  model_tidy <- bind_rows(ref_row, model_tidy)
  
  # Case/control counts
  counts <- subdat %>%
    group_by(strata) %>%
    summarise(
      cases = sum(!!sym(trait_code) == 1),
      controls = sum(!!sym(trait_code) == 0),
      odds = round(cases / controls, 4),
      .groups = "drop"
    )
  
  # Merge and finalize
  merged <- full_join(model_tidy, counts, by = "strata") %>%
    mutate(
      outcome = trait_label,
      ID = trait_code
    ) %>%
    select(ID, outcome, strata, beta, SE, z, OR, CI_lower, CI_upper, p.value, cases, controls, odds)
  
  # Store result
  results_all[[trait_code]] <- merged
}


# Combine all results into a single dataframe
final_results <- bind_rows(results_all)

options(scipen = 999)
write.table(final_results,file=paste0("Strat_Forest_",k,"_BASIS_1_EVER.txt"),sep="\t",col.names=T,row.names=F,quote=F)











