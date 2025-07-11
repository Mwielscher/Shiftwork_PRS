wdir=c("/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/genomicSEM/multiVAR_GWAS/")
wld<-c("/gpfs/data/fs71407/epiwielsch/ukb_2/application/initial_4_traits/genomic_SEM/eur_w_ld_chr/")


library(GenomicSEM,lib.loc = "/gpfs/data/fs71407/epiwielsch/ukb_2/application/initial_4_traits/genomic_SEM")

traits_dat=read.table(paste0(wdir,"sleep_traits_factor_analysis.txt"),header=T,sep=" ")
traits=as.character(traits_dat$file)
trait_names=as.character(traits_dat$name)


#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev<-rep(.5,length(traits))
#vector of population prevalences
population.prev<-rep(NA,length(traits))

#the folder of LD weights [typically the same as folder of LD scores]
ld=as.character(wld)

#run LDSC
LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait_names)

save(LDSCoutput,file=paste0(wdir,"LDSCoutput_for_GWAS_20240523.RData"))

