# Phewas phenotype file was produced using DeepPheWAS R package
# https://richard-packer.github.io/DeepPheWAS_site/

## -----------------
phe=read.table("PHE_SCORE_noACCEL_inclCOV_inclSCORE.txt",header=T)
dat=read.table("pheno_data_ready/all_pop_FINAL_PHEWAS_PHENO",header=T, sep="\t")
phe.map=read.csv("pheno_data_ready/FINAL_PHE_stats_relate_remove.csv")

# Function to summarize the dataset
summarize_data <- function(dat) {
  summary_dt <- data.table(
    Variable = colnames(dat),
    NA_Count = sapply(dat, function(x) sum(is.na(x))),
    Zero_Count = sapply(dat, function(x) if (is.numeric(x)) sum(x == 0, na.rm = TRUE) else NA),
    One_Count = sapply(dat, function(x) if (is.numeric(x)) sum(x == 1, na.rm = TRUE) else NA),
    Non_NA_Count = sapply(dat, function(x) sum(!is.na(x))),
    Type = sapply(dat, function(x) {
      if (is.numeric(x)) {
        unique_vals <- unique(na.omit(x))
        if (all(unique_vals %in% c(0, 1))) {
          return("Binary")
        } else {
          return("Continuous")
        }
      } else {
        return("Other")
      }
    })
  )
  
  setnames(summary_dt, c("Variable", "NA_Count", "Zero_Count", "One_Count", "Non_NA_Count", "Type"))
  return(summary_dt)
}

summary_table <- summarize_data(dat1)
sum.tab1=summary_table[summary_table$One_Count > 20 |summary_table$Type %in% c("Continuous") , ]
write.table(summary_table,file="summary_table_for_SCORE_set.txt",sep="\t",col.names=T,row.names = F, quote =F)

## ------- use scaled scores only
k=c("sc_F2_score")
k=c("sc_F4_score")
k=c("sc_F2inclPSYCH_score")
k=c("sc_F4inclPSYCH_score")
####  ----- run biary traits
res=as.data.frame(matrix(NA,nrow=length(phenotypes.bin),ncol=11))
colnames(res)=c("score","trait", "Estimate" ,"SE","z_value","P","OR","CI_up","CI_low","n_case","n_cont" )
u=0
for (t in phenotypes.bin) {
  print (paste ("++++++++ outcome is ", t, "now +++++++++"))
  u= u+1
  formula <- as.formula(paste( t,"  ~",  k, "+ age + sex + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) ##  + BMI
  coeff=summary(glm(formula ,data = phe.dat, family = "binomial" ))$coefficients[2,]
  estimate <- coeff[1]
  std_error <- coeff[2]
  or <- exp(estimate)  # Odds Ratio
  ci_low <- exp(estimate - 1.96 * std_error)  # Lower bound of 95% CI
  ci_up <- exp(estimate + 1.96 * std_error)  # Upper bound of 95% CI
  # Calculate number of cases and controls
  cases <- sum(phe.dat[[t]] == 1, na.rm = TRUE)  # Assuming 't' is the target binary outcome
  controls <- sum(phe.dat[[t]] == 0, na.rm = TRUE)
  res[u,c(1:2)]=c(k,t)
  res[u,c(3:6)]= coeff
  res[u,c(7:11)] = c(or,ci_up,ci_low,cases, controls)
  
}

####  ---------------------------- quantitiative traits
phenotypes.cont=as.character(sum.tab2$Variable[sum.tab2$Type %in% c("Continuous")])
res_2=as.data.frame(matrix(NA,nrow=length(phenotypes.cont),ncol=11))
colnames(res_2)=c("score","trait", "Estimate" ,"SE","z_value","P","OR","CI_up","CI_low","n_case","n_cont" )
u=0
for (t in phenotypes.cont) {  
  print(paste("++++++++ Outcome is", t, "now +++++++++"))
  u = u + 1
  
  # Define regression formula
  formula <- as.formula(paste(t, "~", k, "+ age + sex + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")) 
  
  # Fit linear model
  coeff <- summary(lm(formula, data = phe.dat))$coefficients[2,]
  
  estimate <- coeff[1]  # Regression coefficient
  std_error <- coeff[2]  # Standard error
  t_value <- coeff[3]  # T-statistic
  p_value <- coeff[4]  # P-value
  
  # Calculate additional summary statistics
  mean_val <- mean(phe.dat[[t]], na.rm = TRUE)
  sd_val <- sd(phe.dat[[t]], na.rm = TRUE)
  n_samples <- sum(!is.na(phe.dat[[t]]))  # Number of non-missing samples
  
  # Store results in res_2
  res_2[u, c(1:2)] <- c(k, t)
  res_2[u, c(3:6)] <- coeff  # Estimate, Std. Error, T-value, P-value
  res_2[u, c(7:11)] <- c(NA,NA,NA, n_samples,"continous trait")  # Store mean, SD, and case count
}

res.fin=as.data.frame(rbind(res,res_2))
#ann.dat=read.csv("pheno_data_ready/FINAL_PHE_stats_relate_remove.csv")
#anno=ann.dat[,c("PheWAS_ID","phenotype") ]
res.fin2=merge(res.fin, anno,by.x="trait",by.y="PheWAS_ID",all.x=T )
res.fin2=res.fin2[,c("score","trait","phenotype",colnames(res.fin2)[3:11])]

options(scipen = 999)
write.table(res.fin2,file="PHEWAS_sc_F4_score_FULL_Set_noSCALE.txt",sep="\t",col.names=T,row.names=F,quote=F )
