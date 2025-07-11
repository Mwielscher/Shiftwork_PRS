## this script is supposed to calculate quartiles form the F2* score and evaluate the association of these quartiles to clinically relevant outcomes in an obesity stratified analysis

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

## -------------- run analysis
##-- prepare
k=c("sc_F2inclPSYCH")
###
phe.dat <- phe.dat %>%
  mutate(sc_F2inclPSYCH_quartile = ntile(sc_F2inclPSYCH_score, 4))

phe.dat <- phe.dat %>%
  mutate(obese = factor(if_else(BMI >= 30, "obese", "non-obese")))

###
phe.dat <- phe.dat %>%
  mutate(
    PRS_quartile = factor(sc_F2inclPSYCH_quartile, levels = 1:4),
    Shift_work = factor(night, levels = c("day","ever")),
    strata = interaction(PRS_quartile, Shift_work, sep = "_")
  )

phe.dat <- phe.dat %>%
  mutate(
    PRS_quartile = factor(sc_F2inclPSYCH_quartile, levels = 1:4),
    Shift_work = factor(night, levels = c("day", "ever")),
    strata = interaction(PRS_quartile, Shift_work, obese, sep = "_")
  )

table(phe.dat$strata)
phe.dat$strata <- relevel(phe.dat$strata, ref = "1_day_non-obese")
### 
# Create output container
results_all <- list()

# Loop over each trait
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
      PC1:PC10
    ) %>%
    filter(!is.na(!!sym(trait_code)), !is.na(strata)) %>%
    mutate(strata = droplevels(factor(strata)))  # ensure clean factor levels
  
  # Fit logistic model
  model <- glm(
    as.formula(paste(trait_code, "~ strata + age + sex + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")),
    data = subdat,
    family = binomial
  )
  
  # Log-odds output
  model_log <- tidy(model, conf.int = FALSE, exponentiate = FALSE) %>%
    filter(term != "(Intercept)") %>%
    mutate(strata = gsub("strata", "", term)) %>%
    select(strata, beta = estimate, SE = std.error, z = statistic, p.value)
  
  # Odds ratios
  model_or <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(strata = gsub("strata", "", term)) %>%
    select(strata, OR = estimate, CI_lower = conf.low, CI_upper = conf.high)
  
  # Merge results
  model_tidy <- left_join(model_log, model_or, by = "strata")
  
  # Add the reference row (1_day_non-obese)
  ref_row <- tibble(
    strata = "1_day_non-obese",
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
      cases = sum(!!sym(trait_code) == 1, na.rm = TRUE),
      controls = sum(!!sym(trait_code) == 0, na.rm = TRUE),
      odds = round(cases / controls, 4),
      .groups = "drop"
    )
  
  # Merge final result
  merged <- full_join(model_tidy, counts, by = "strata") %>%
    mutate(
      outcome = trait_label,
      ID = trait_code
    ) %>%
    select(ID, outcome, strata, beta, SE, z, OR, CI_lower, CI_upper, p.value, cases, controls, odds)
  
  # Store
  results_all[[trait_code]] <- merged
}


# Combine all results into a single dataframe
final_results <- bind_rows(results_all)

final_res1=merge(final_results,trait_table,all.x=T, by.x="ID",by.y="plot.trait")



