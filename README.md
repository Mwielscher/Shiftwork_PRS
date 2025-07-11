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


## Mendelian_Randomisation


## PRS_CS

## UKB_PheWAS


## All_of_Us_PheWAS


## Nigthshift_stratified

