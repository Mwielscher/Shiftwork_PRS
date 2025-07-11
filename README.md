## Table of contents
1. [About this Repository](#About-this-Repository)
2. [LDSC](#LDSC)
3. [MultiVar_GWAS](#MultiVar_GWAS)
4. [Mendelian_Randomisation](#Mendelian_Randomisation)
5. [PRS_CS](#PRS_CS)
4. [UKB_PheWAS](#UKB_PheWAS)
5. [All_of_Us_PheWAS](#All_of_Us_PheWAS)
6. [Nigthshift_stratified](#Nigthshift_stratified)

## About this Repository
This repository accompanies the manuscript [__"Modeling Sleep Genetics Uncovers Distinct Disease Risks in Night vs. Day Workers"__](https://www.nature.) 
<br/><br/> <br/>

<p align="center">
<img src="/img/figure_1.png" alt="Overview Figure" width="600"/>
<br/><br/>



> **_Abstract:_**  Shift work has adverse health consequences, but the genetic basis of this vulnerability remains underexplored. Using Genomic structural equation modelling, we derived six latent sleep factors reflecting core RU-SATED domains from large-scale GWAS data. Polygenic scores for sleep regularity and daytime alertness showed strong genetic correlations with metabolic and behavioral traits and were causally linked to depression, well-being, and insulin regulation in Mendelian randomization. A phenome-wide association study in UK Biobank subset of night shift workers and controls (N = 46,211) revealed associations with T2D, hypertension, and COPD, and other outcomes, which were independently replicated in the All of Us cohort (N = 131,729). Incorporating behavioral traits into extended PRS improved predictive power. Stratified analyses showed that genetic risk for poor sleep regularity had the strongest differential effect, increasing T2D risk by over twofold in night shift workers compared to day workers. These findings highlight how genetically encoded sleep traits interact with environmental exposures to shape disease risk. 
<p>
<br/>


The scripts to perform the analyses presented in the manuscript are deposited in this repository. These are custom scripts, which were shared with collaborators for the manuscript or run on my local destop computer. Most shell scripts are designed to run on a [SGE cluster](http://gridscheduler.sourceforge.net/htmlman/manuals.html).  

## LDSC
LDSC generates scores reflecting whether the GWAS test statistic of a biologically relevant variant correlates with nearby variants in high linkage disequilibrium. The z statistic for the genetic association of each variant with trait 1 are multiplied with the z statistic for the genetic association with trait 2, followed by regression of this product of statistics against the LD scores. The slope (coefficient) represents genetic correlation. When large, the same genetic variants impact both the traits.

## MultiVar_GWAS  
We incorporated 21 sets of summary statistics, totaling 3.5 million SNPs, into our structural multivariable regression model. All models were based on samples of European ancestry or trans-ethnic meta-analyses, with European linkage disequilibrium maps as references.
To mitigate type 1 error inflation, we used the built-in genomic control function of the Genomic SEM package. Additionally, we employed the "smooth_check" option to exclude SNPs requiring excessive smoothing, enhancing the reliability of the GWAS results. The "fix_measurement" option was also applied to ensure stable factor loadings, preventing implausible estimates and Heywood cases.
__GCTA:__ We used the UK Biobank genotyping array data as our reference sample (see UKB genetic data QC section). A genome-wide significance threshold of 5x10⁻⁸ was applied. We analyzed a 10 MB window around each SNP, assuming that SNPs outside this window are in complete linkage equilibrium. Additionally, we applied a collinearity threshold of 0.9 to exclude highly correlated SNPs.


## Mendelian_Randomisation
Scripts include: 
__Proxy search__ - we used proxies in high linkage disequilibrium (LD) with the original instruments, applying an LD threshold of 0.8 in the European sample of the 1000 Genomes Project Phase 3. Proxy SNPs were retrieved using the Ensembl REST API. We then selected the SNP with the lowest P-value from our summary statistics.
__Pleiotropy evaluation;__ we also calculate heterogeneity indicators specific to each outcome, including Cochran's Q test results, which assess variability in causal estimates and potential horizontal pleiotropy. We also evaluated the Egger intercept to flag directional pleiotropy in the instruments. Scripts also return Steiger P values for directionality  
__Sensitivity analysis:__ radial IVW method, designed to handle outliers within the instrument set by re-weighting the contribution of each variant based on its residual, making it more robust to heterogeneity among genetic instruments. Additionally, we included MR-RAPS to handle cases where some genetic instruments are weak, mitigating weak instrument bias.
Finally, we reported results from the MR-PRESSO method, including the MR-PRESSO Global P-value, which tests for horizontal pleiotropy by simulating a null distribution of residual sum of squares and comparing it to the observed values. This identifies whether the observed residual sum of squares falls within the expected variance or suggests horizontal pleiotropy. We also provided the MR-PRESSO Distortion P-value, which compares the causal estimate before and after outlier correction, indicating the influence of outliers on the causal estimate.


## PRS_CS
We calculated scores using the PRS-CS approach, which employs a Bayesian framework that incorporates linkage disequilibrium (LD) information from large reference panels. This method adaptively shrinks SNP effect sizes, allowing for the flexible handling of both large and small genetic effects.  
__MTAG__ We ran MTAG using its default settings, thus obtained the same number of summary statistics as the input traits, with each set containing updated meta-analyzed effect estimates based on the genetic correlations among the included traits. We then focused on the MTAG results where effect estimates were driven by the sleep factor, rather than the psychiatric or behavioral traits. These adjusted factors were subsequently used in PheWAS analysis and polygenic score (PGS) analysis, conducted separately in night shift and day shift workers.  

## UKB_PheWAS
The following data fields were used: 20002 (Self-reported non-cancer diagnosis), 20008 (Year of self-reported diagnosis), 20004 (Self-reported operation code), 40006 (Cancer diagnosis from cancer registry), and 40005 (Date of cancer diagnosis). Additionally, we retrieved record-level data from 42040 (GP clinical event records), 42039 (GP prescription records), 41259 (Hospital Episode Statistics [HES] inpatient main dataset), 41234 (HES inpatient diagnosis), 41149 (HES inpatient operations), 40023 (Cause of death), and 40023 (Date of death). Using DeepPheWAS, we mapped these data to phecodes, created quantitative phenotypes, and generated combined phenotypes. This resulted in a dataset containing 1,234 variables, each with at least 20 entries in the analysed UKB subpopulation 


## [All_of_Us_PheWAS]( https://github.com/LeoVincenzi/AllOfUs_to_UKBB)  
The replication was conducted in the All of Us Research Program, using survey and EHR data mapped to UK Biobank trait definitions. After careful harmonisation, we were able to replicate 933 of the original 1,233 phenotypes tested in the UK Biobank. The entire process and all code to harmonise between cohorts is stored in a [dedicated repository]( https://github.com/LeoVincenzi/AllOfUs_to_UKBB) 


## Nigthshift_stratified

