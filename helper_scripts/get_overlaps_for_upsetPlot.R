####  -----------------------------  factor 1 --- threshold is BON
rm(list=ls())
dat=read.table("sleep_factor1_GCTA_sentinels_clean.txt",header=T)

indep_1=read.table("sleep_trait_independent/SR_chronotype_significant_indepenedent.txt",header=T)
indep_2=read.table("sleep_trait_independent/SR_morning_person_significant_indepenedent.txt",header=T)
indep_3=read.table("sleep_trait_independent/ACC_L5_TIME_RAW_significant_indepenedent.txt",header=T)
indep_4=read.table("sleep_trait_independent/ACC_M10_TIME_RAW_significant_indepenedent.txt",header=T)
indep_5=read.table("sleep_trait_independent/ACC_SLEEP_MIDP_RAW_significant_indepenedent.txt",header=T)

### for Plink extraction
plink.dat=rbind(indep_1[,c(1:2)],indep_2[,c(1:2)],indep_3[,c(1:2)],indep_4[,c(1:2)],indep_5[,c(1:2)])
write.table(plink.dat,file="../FACTOR_1_sumstat_forPLINK.txt",sep="\t",col.names=T,row.names=F,quote=F)

base.dat=list(indep_1,indep_2,indep_3,indep_4,indep_5)


sub.traits=c("SR_chronotype","SR_morning_person","ACC_L5_TIME","ACC_M10_TIME","ACC_SLEEP_MIDP")
names(base.dat) = sub.traits
dat1=dat
for (k in 1:5){
  dat1=as.data.frame(cbind(dat1,rep(0,nrow(dat1))))
  colnames(dat1)[ncol(dat1)]=as.character(sub.traits[k])
  inter=as.data.frame(base.dat[[k]])
  print (paste("this is sumstat" ,sub.traits[k], "now"))
 ## now check overlap
 for (t in 1:nrow(dat1)) {
    senti=dat1[t,]
    rsid=senti$SNP
    chr=senti$Chr
    upper=senti$bp + 500000
    lower = senti$bp - 500000
     for (w in 1:nrow(inter)) {
       n.senti=inter[w,]
       if (n.senti$CHROM %in% chr & n.senti$POS > lower & n.senti$POS < upper) {
         dat1[t,sub.traits[k]] = n.senti$RSID
         print (paste("there is overlap to", sub.traits[k], "sumstats; new ID is ",rsid, "old ID was", inter[w,"RSID"]))
         inter[w,"RSID"]=as.character(rsid)
         ## update base.dat
         base.dat[[k]][w,"RSID"]=as.character(rsid)
     }
  }#close  lookup
 }#  close dat loop    
} # close master trait loop
###  --------------------------------------------------------------------------
write.table(dat1,file="factor1_sentinels_overlaps_to_input.txt",sep="\t",col.names=T,row.names=F,quote=F)
##    ----      now rest  just run up p and k manually !!
p=1
for (k in 2:5) {
  print (paste("!!!!!  reference is ",sub.traits[p],"now!!!!!!"))
  new.ref=as.data.frame(base.dat[[p]])
  inter=as.data.frame(base.dat[[k]])
  print (paste("this is sumstat" ,sub.traits[k], "now"))
  ## now check overlap
  for (t in 1:nrow(new.ref)) {
    senti=new.ref[t,]
    rsid=senti$RSID
    chr=senti$CHROM
    upper=senti$POS + 500000
    lower = senti$POS - 500000
    for (w in 1:nrow(inter)) {
      n.senti=inter[w,]
      if (n.senti$CHROM %in% chr & n.senti$POS > lower & n.senti$POS < upper) {
        dat1[t,sub.traits[k]] = n.senti$RSID
        print (paste("there is overlap to", sub.traits[k], "sumstats; new ID is ",rsid, "old ID was", inter[w,"RSID"]))
        inter[w,"RSID"]=as.character(rsid)
        ## update base.dat
        base.dat[[k]][w,"RSID"]=as.character(rsid)
      }
    }
  }#close outer lookup
} # close master trait loop

factor1.dat=base.dat
save(factor1.dat,file="factor1_for_upset.RData")