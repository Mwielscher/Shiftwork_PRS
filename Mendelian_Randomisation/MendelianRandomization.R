library(TwoSampleMR)
library(MRPRESSO)
library(RadialMR)
library(data.table)
library(httr)
library(jsonlite)
library(dplyr)
library(ggplot2)

###   -----------------  find proxies function 

find_proxies_ensembl <- function(snp, population = "1000GENOMES:phase_3:EUR", r2_threshold = 0.8) {
  url <- paste0("https://rest.ensembl.org/ld/human/", snp, "/1000GENOMES:phase_3:EUR?content-type=application/json")
  response <- httr::GET(url)
  
  if (response$status_code != 200) {
    message("Error fetching data for SNP: ", snp, " (HTTP status code: ", response$status_code, ")")
    return(NULL)
  }
  
  content <- httr::content(response, as = "text")
  # Parse JSON content
  parsed_content <- jsonlite::fromJSON(content)
  # Check if 'data' field is present
  if (!is.data.frame(parsed_content)) {
    message("No valid data found for SNP: ", snp)
    return(NULL)
  }
  proxies <- as.data.frame(parsed_content)
  # Check if proxies data contains expected structure
  if (!"variation2" %in% colnames(proxies) || !"r2" %in% colnames(proxies)) {
    message("Expected columns not found for SNP: ", snp)
    return(NULL)
  }
  # Filter proxies by r2 threshold
  proxies <- proxies %>% filter(r2 >= r2_threshold)
  
  # Debugging: Print the first few rows of proxies data
  message("Proxies for SNP: ", snp)
  return(proxies)
}


####  ----------- safe extract function
safe_extract <- function(expr) {
  tryCatch({
    val <- eval(expr)
    if (is.null(val)) {
      NA
    } else {
      val
    }
  }, error = function(e) {
    NA
  })
}


sumstat.outcome=read.table("sumstats/outcomes/all_outcome_sumstat.txt",header=T)


##   ------------------------------------------------------------------------------------
## loop though cardiovascular traits
## produce one file per factor

### ------------------- outcome will stay the same 
exp=c("factor1")
exp_file=paste0("independent_loci/factors/",exp,"_GCTA_sentinels_clean_clean.txt")
exp.dataf=fread(exp_file)

exp.sum.full=fread(paste0("sumstats/factors/",exp,"_multiGWAS_final.txt"))

for (i in 1:nrow(exp.dataf)) {
  snp=exp.dataf$RSID[i]
  allele=as.character(exp.sum.full$A2[exp.sum.full$rsid %in% snp])
  if (length(allele) <1 ){
    next
  }
  exp.dataf$A2[i]=allele
}

#exp.dataf=exp.dataf[1:31,]

rm(final.output)


for (out in sumstat.outcome$trait){
  #for (exp in sumstat.outcome$trait[11:length(sumstat.outcome$trait)]){
  print(c(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  "))
  print(paste("this is exposure",exp, "vs outcome",out))
  print(c(" +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  "))
  
  ### ---------------------------------------   find proxies
  trait.file=sumstat.outcome$file[sumstat.outcome$trait %in% out]
  out_file=paste0("sumstats/outcomes/",trait.file )
  out.dataf=fread(out_file)
  
  exp.proxies=exp.dataf$RSID[!exp.dataf$RSID %in% out.dataf$rsid ]
  outcome_snps = unique(out.dataf$rsid)
  
  # Initialize a data frame to store results
  proxy_results <- data.frame(exposure_snp = character(),
                              proxy_snp = character(),
                              r2 = numeric(),
                              stringsAsFactors = FALSE)
  
  # Find proxies for each SNP in the exposure dataset
  for (snp in exp.proxies) {
    message("Finding proxies for SNP: ", snp)
    proxies <- find_proxies_ensembl(snp)
    if (!is.null(proxies)) {
      # Filter proxies that are present in the outcome dataset
      proxies_present <- proxies %>% filter(variation2 %in% outcome_snps)
      if (nrow(proxies_present) > 0) {
        proxy_results <- rbind(proxy_results, 
                               data.frame(exposure_snp = snp, 
                                          proxy_snp = proxies_present$variation2, 
                                          r2 = proxies_present$r2))
      }
    }
  }
  ####      --------------------------------------------------------------
  ####  ------------------------   update exposure dataframe
  exp.match=exp.dataf[! exp.dataf$RSID %in% exp.proxies,]
  proxi.found=unique(proxy_results$exposure_snp)
  not.possible=exp.dataf$RSID[! exp.dataf$RSID %in% c(exp.match$RSID,proxi.found)]
  
  print(paste0("exact match exposure outcome for: ", length(exp.match$RSID) ))
  print(paste0("exposure proxy matched to outcome for: ", length(proxi.found) ))
  print(paste0("exposure SNP could NOT be matched to outcome for: ", length(not.possible) ))
  
  
  inter.exp=fread(paste0("sumstats/factors/",exp,"_multiGWAS_final.txt"))
  
  ## get proxy with lowest P-value
  rm(exp.new)
  for (k in proxi.found) {
    extr=proxy_results[proxy_results$exposure_snp %in% k,]
    candidates=inter.exp[inter.exp$rsid %in% extr$proxy_snp,]
    if (nrow(candidates) <1) {
      print(paste("proxies not found in full exposure dataset for SNP", k))
      not.possible=c(not.possible,k)
    } else {
      add.on=candidates[candidates$P %in% min(candidates$P),] 
      add.on1=add.on[1,]
      if (!exists("exp.new")) {
        exp.new=as.data.frame(add.on1)
      }else {
        exp.new=as.data.frame(rbind(exp.new,add.on1) )
      }
    }
  }
  
  
  if (!exists("exp.new")){
    final.exp = exp.match
    exp.new=data.frame(Name = character(0), Age = integer(0), Gender = character(0))
  } else {
    
    
    print(paste0("exact match exposure outcome for: ", length(exp.match$RSID) ))
    print(paste0("exposure proxy matched to outcome FINAL: ", dim(exp.new)[1] ))
    print(paste0("exposure SNP could NOT be matched to outcome FINAL: ", length(not.possible) ))
    
    exp.new$CHR=rep(NA,nrow(exp.new))
    exp.new$POS=rep(NA,nrow(exp.new))
    colnames(exp.new)=as.character(colnames(exp.match))
    
    ## ---------------
    final.exp=as.data.frame(rbind(exp.match,exp.new))
  }    ## else loop for new proxies
 
  
  write.table(final.exp, file=paste0("independent_loci/factors/",exp,"_incl_proxy_for_",out,".txt"),sep="\t",
              col.names=TRUE,row.names=FALSE,quote=FALSE)
  ### ------------------------------------  
  ##### get data intp R-package
  expo_file=paste0("independent_loci/factors/",exp,"_incl_proxy_for_",out,".txt")
  
  ex.dat <- read_exposure_data(
    filename = expo_file,
    sep = "\t",
    snp_col = "RSID",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A2",
    other_allele_col = "A1",
    eaf_col = "Freq",
    pval_col = "P_value",
    samplesize_col = "N")
  #####   ---------------------------  subset and fix BP SE
  out.dataf2=out.dataf[out.dataf$rsid %in% final.exp$RSID,]
  
  ####--- get SE for blood pressure traits
  if (out %in% c("DBP","SBP")) {
    for(k in 1:nrow(out.dataf2)) {
      # Degrees of freedom (assuming simple linear regression with one predictor)
      df <- out.dataf2$N[k] - 2
      # Calculate the t-statistic from the p-value
      t_stat <- qt(out.dataf2$P[k] / 2, df, lower.tail = FALSE)
      # Calculate the standard error
      SE_beta <- out.dataf2$beta[k] / t_stat
      out.dataf2$se[k]=abs(SE_beta)
    }
  }
  
  write.table(out.dataf2, file=paste0("sumstats/outcomes/matched_sumstat_",exp,"and",out,".txt"),sep="\t",
              col.names=TRUE,row.names=FALSE,quote=FALSE) 
  
  out_file2=paste0("sumstats/outcomes/matched_sumstat_",exp,"and",out,".txt")
 
   out.dat <- read_outcome_data(
    filename = out_file2,
    sep = "\t",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "A2",
    other_allele_col = 'A1',
    eaf_col = "Freq",
    pval_col = "P",
    samplesize_col = "N")
  ####    ---------------------------------------------------
  #### -----------------------------   harmonize
  dat.harm = harmonise_data(ex.dat,out.dat)
  
  ##  ------------------------------------------------------------------
  ## ------------------------------------        assess instruments:
  F = dat.harm$beta.exposure[dat.harm$mr_keep %in% c("TRUE")]^2/dat.harm$se.exposure[dat.harm$mr_keep %in% c("TRUE")]^2
  mF = mean(F)
  mF
  
  if (nrow(dat.harm)<2){
    steiger=directionality_test(dat.harm) # The Steiger test compares the variance explained by the genetic instruments in the exposure versus the outcome. The hypothesis being tested is whether the genetic instruments explain more variance in the exposure than in the outcome.
    # Significant Steiger P Value indicates that the genetic instruments are more strongly associated with the exposure than with the outcome. This supports the hypothesis that the exposure causes the outcome.
    sin = mr_singlesnp(dat.harm)
    ### now fix het and plt
    het=as.data.frame(matrix(NA, nrow = 2, ncol = 8))
    plt=as.data.frame(matrix(NA, nrow = 2, ncol = 8))
    
  }else {
    steiger=directionality_test(dat.harm) # The Steiger test compares the variance explained by the genetic instruments in the exposure versus the outcome. The hypothesis being tested is whether the genetic instruments explain more variance in the exposure than in the outcome.
    # Significant Steiger P Value indicates that the genetic instruments are more strongly associated with the exposure than with the outcome. This supports the hypothesis that the exposure causes the outcome.
    
    het = mr_heterogeneity(dat.harm)  # Null Hypothesis: All instruments estimate the same causal effect. Alternative Hypothesis: There is variability in the causal estimates from different instruments.
    #Indicates that the estimates from different instruments are not consistent. Suggesting possible horizontal pleiotropy or other violations of MR assumptions.
    
    plt = mr_pleiotropy_test(dat.harm) ## this is a MR-Egger Intercept Test
    # If the intercept term is significantly different from zero, it suggests the presence of directional pleiotropy.
    sin = mr_singlesnp(dat.harm)
  }
  #### ------------------------------------------------------------
  ## ---------------------------------------   run analysis
  all_res_empty=as.data.frame(matrix(NA,nrow=8,ncol=14))
  colnames(all_res_empty)=c("id.exposure","id.outcome","outcome","exposure","method","nsnp","b","se","pval","lo_ci","up_ci","or","or_lci95","or_uci95")
  
  tryCatch({
    all_res <- mr(dat.harm, method_list = c("mr_ivw", "mr_simple_median", "mr_weighted_median", "mr_weighted_mode", 
                                            "mr_ivw_radial", "mr_egger_regression", "mr_raps", "mr_sign"))
    ## get odds ratios and CI
    all_res_OR = generate_odds_ratios(all_res)
  }, error = function(e) {
    # If an error occurs, fill the data frame with NA values
    all_res_OR <- all_res_empty
    message("An error occurred: ", e$message)
  })
  
  
  # MR-RAPS is designed to handle situations where some of the genetic instruments are weak, meaning they have a small effect on the exposure.
  # IVW radial down-weights the influence of genetic variants that are outliers, which can potentially reduce bias in the causal estimate. This is achieved by re-weighting the contribution of each variant based on its residual. It is more robust to heterogeneity among genetic instruments compared to traditional IVW, as it adjusts the weights of the variants.
  if (nrow(dat.harm)<5) {
    presso[1:2,2:6] = NA
    mr_pre[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue =NA
    mr_pre[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue =NA
    mr_pre[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices` = NULL
    
    
  }else {
    mr_pre = run_mr_presso(dat.harm, NbDistribution = 1000, SignifThreshold = 0.05)
  }
  ## MR-Presso:
  #MR-PRESSO performs instrumental variable regression to estimate the causal effect of the exposure on the outcome. It calculates the residuals for each genetic variant.
  ## outlier TEST: MR-PRESSO identifies outliers by comparing the residuals to a distribution of residuals generated under the null hypothesis of no horizontal pleiotropy. 
  # Genetic variants with residuals that are significantly larger than expected are flagged as outliers.
  ## Distortion TEST: The outlier test provides a list of outliers, and these outliers are removed or down-weighted in the analysis to correct the causal estimate. The distortion test compares the causal estimates before and after the removal of outliers. A significant distortion test indicates that the outliers were substantially influencing the causal estimate, and the corrected estimate should be considered more reliable.
  ## Global Test for Horizontal Pleiotropy: It evaluates whether the observed residual sum of squares (RSS) is significantly higher than what would be expected by chance. It performs a large number of simulations to generate a null distribution of RSS values. The observed RSS is then compared to this distribution to calculate the global P-value.
  global_pvalue <- safe_extract(quote(mr_pre[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue))
  distortion_pvalue <- safe_extract(quote(mr_pre[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue))
  num_outliers <- safe_extract(quote(length(mr_pre[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)))
  
  
  ####  -------------------------------------------------------------------------------------
  ####   --------              organize output
  
  output=as.data.frame(all_res_OR)
  output=output[,3:ncol(output)]
  output$outcome=rep(out,nrow(output))
  output$exposure=rep(exp,nrow(output))
  presso.head=c("outcome","exposure","MR_Analysis","Causal_Estimate",'Sd',"T-stat","P_value", "Presso_global_test_P","Presso_Distortion Test","#Presso_outlier","NA","NA")
  output=as.data.frame(rbind(output,presso.head))
  presso=as.data.frame(mr_pre[[1]][1])
  output[10,1:2]=c(out,exp)
  output[10,3:7]=presso[1,2:6]
  output[11,1:2]=c(out,exp)
  output[11,3:7]=presso[2,2:6]
  output[10,8:10]=c(global_pvalue,distortion_pvalue,num_outliers )
  output[11,8:10]=c(global_pvalue,distortion_pvalue,num_outliers )
  
  snp.head=c("outcome","exposure","totalSNPs","#matched",'#proxy',"#not_possible","#harmonisation_excl", "matched","poxy","not_possible","harmonisation_excl","inAnalysis")
  output=as.data.frame(rbind(output,snp.head))
  output[13,1:2]=c(NA,exp)
  output[13,3:7]=c(length(exp.dataf$RSID),length(exp.match$RSID),dim(exp.new)[1],length(not.possible),length(dat.harm$SNP[dat.harm$mr_keep %in% c("FALSE")]))
  output[13,8]=paste(exp.match$RSID, collapse = "|")
  output[13,9]=paste(exp.new$RSID, collapse = "|") 
  output[13,10]=paste(not.possible, collapse = "|") 
  output[13,11]=paste(dat.harm$SNP[dat.harm$mr_keep %in% c("FALSE")], collapse = "|") 
  output[13,12]=paste(dat.harm$SNP[dat.harm$mr_keep %in% c("TRUE")], collapse = "|") 
  snp2.head=c("outcome","exposure","het_Q","het_Q_df",'het_Q_pval',"plt_Egger_intercept", "plt_Egger_se", "plt_Egger_inter_Pval", "steiger_exposure_R2","steiger_outcome_R2","steiger_causal_direction","steiger_pval")
  output=as.data.frame(rbind(output,snp2.head))
  output[15,1:2]=c(out,exp)
  output[15,3:5]=het[2,6:8]
  output[15,6:8]=plt[1,5:7]
  output[15,9:12]=steiger[1,5:8]
  
  
  ##### -------------------------------------------------------------
  ### --------------------  collect and combine output
  if (!exists("final.output")) {
    final.output=as.data.frame(output)
  }else {
    final.output=as.data.frame(rbind(final.output,output) )
  }
  
  ###  ---------------------------------------------------------------------------------
  ####      plotting !!!!
  
  tryCatch({
    plot_res = mr(dat.harm, method_list=c("mr_ivw" , "mr_ivw_radial",
                                          "mr_egger_regression",  "mr_raps"))
    scatter_plots <- mr_scatter_plot(plot_res, dat.harm)
    
    # Access the first plot if it's a list
    if (is.list(scatter_plots)) {
      scatter_plot <- scatter_plots[[1]]
    } else {
      scatter_plot <- scatter_plots
    }
    
    # Customize the scatter plot
    scatter_plot <- scatter_plot +
      ggtitle("Mendelian Randomization Scatter Plot") +
      xlab(paste("Effect size exposure:", exp)) +
      ylab(paste("Effect size on outcome:", out)) +
      theme_minimal() +
      theme(legend.position = "bottom") +
      theme(legend.title = element_blank())
    scatter_plot
    ggsave(paste0("results/plots/MRscatter_exp_",exp,"_out_",out,".pdf"), plot = scatter_plot, device = "pdf", width = 5, height = 5, units = "in")
    
  }, error = function(e) {
    # If an error occurs, fill the data frame with NA values
    
    message("An plotting error occurred: ", e$message)
  })    
  
  
  funnel_plots=mr_funnel_plot(sin)
  # Access the first plot if it's a list
  if (is.list(funnel_plots)) {
    funnel_plot <- funnel_plots[[1]]
  } else {
    funnel_plot <- funnel_plots
  }
  # Customize the funnel plot to match the scatter plot
  funnel_plot <- funnel_plot +
    ggtitle(paste("MR Funnel Plot for exposure:",exp) ) +
    xlab(paste("Standard Error of SNP Effect on outcome:",out)) +
    ylab(paste("Effect Size of SNP on outcome:",out)) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    theme(legend.title = element_blank())
  # Print the customized funnel plot
  ggsave(paste0("results/plots/MRfunnel_exp_",exp,"_out_",out,".pdf"), plot = funnel_plot, device = "pdf", width = 5, height = 4, units = "in")
  
  
}  ## colse exposure for loop


write.table(final.output, file=paste0("results/MR",exp,"_vs_all_exposures.txt"),sep="\t",
            col.names=TRUE,row.names=FALSE,quote=FALSE)
