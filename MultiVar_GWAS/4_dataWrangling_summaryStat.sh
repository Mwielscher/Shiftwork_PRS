#!/bin/bash

wdir='/gpfs/data/fs71407/epiwielsch/ukb_2/application/exploratory_analysis/genomicSEM/multiVAR_GWAS/'
#mkdir -p ${wdir}sum_result


while read -r fac; do
echo "------ this is ${fac} ----"

## --- for now we take the header line from the new directory change diretory names
sed -n 1p ${wdir}/result/${fac}_chr21_20240523.txt > ${wdir}sum_result/${fac}_multiGWAS_incl19.txt
tail -n +2 ${wdir}/result/${fac}_chr*.txt | sed -e '1d' >> ${wdir}sum_result/${fac}_multiGWAS_incl19.txt

cat <<EOF >plot_and_sort_${fac}.R

dat=read.table("${wdir}sum_result/${fac}_multiGWAS_incl19.txt",sep="\t",header=T,fill=T)
dat=dat[!is.na(dat\$Pval_Estimate),]
dat=dat[!is.na(dat\$est),]
dat1=dat[,!colnames(dat) %in% c("warning","lhs","op", "rhs", "free", "label")]
dat2=dat1[order(dat1\$Pval_Estimate),]
write.table(dat2,file="${wdir}sum_result/${fac}_multiGWAS_sorted_incl19.txt",sep="\t",col.names=T,row.names=F,quote=F)
dat1=dat2
##  make manhattan and qq-plot
library(ggplot2)
library(dplyr)
ci=0.95


colnames(dat1)[10] =c("p")


res <- dat1 %>% filter(!is.na(p)) %>%
    arrange(p) %>%
    mutate(r=rank(p, ties.method = "random"),
           pexp=r/length(p),
           clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = length(p)-r)),
           cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = length(p)-r)))

p1 <- res %>%
    ggplot(aes(x=-log10(pexp),y=-log10(p))) +
    geom_ribbon(mapping = aes(x = -log10(pexp), ymin = clower, ymax = cupper),
      alpha = 0.1,color="darkgray") +
    geom_point() +
    geom_abline(slope=1,intercept=0) +
   ggtitle("${fac} GWAS") +
   theme(plot.title = element_text(hjust = 0.5))+
    xlab(expression(Expected -log[10](p))) +
    ylab(expression(Observed -log[10](p))) + 
    theme_bw()



ggsave("${wdir}/sum_result/qqplot_${fac}_incl19.jpeg", plot = p1, width = 6, height = 5, dpi = 300)

#### now manhattan plot

dat1=dat1[!is.na(dat1\$BP),]
dat1\$BP=as.numeric(as.character(dat1\$BP ))
dat1\$CHR=factor(dat1\$CHR, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22" ))
colnames(dat1)[10]= c("P") 

don <- dat1 %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dat1, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BP1=BP+tot)

don\$CHR=factor(don\$CHR, levels=c("1","2","3","4","5","6","7","8","9","10",
                                       "11","12","13","14","15","16","17","18","19","20","21","22" ))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BP1) + min(BP1) ) / 2 )
axisdf\$CHR=factor(axisdf\$CHR, levels=c("1","2","3","4","5","6","7","8","9","10",
"11","12","13","14","15","16","17","18","19","20","21","22" ))



don\$CHR=factor(don\$CHR, levels=c("1","2","3","4","5","6","7","8","9","10",
"11","12","13","14","15","16","17","18","19","20","21","22" ))
print(length(don\$P))

don\$P=as.numeric(as.character(don\$P))
man <- ggplot(don, aes(x=BP1, y=-log10(P))) +
    
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    ylim(0,c(max(-log10(don\$P))+1) ) +
    # custom X axis:
    scale_x_continuous( label = axisdf\$CHR, breaks= axisdf\$center ) +
    #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    ggtitle("${fac} multivariate GWAS") +
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

ggsave("${wdir}/sum_result/manhattan_${fac}_incl19.jpeg", plot = man, width = 10, height = 5, dpi = 300)


EOF
chmod +x plot_and_sort_${fac}.R

cat <<EOFSUBMIT >submit_qc_plots_${fac}.sh
#!/bin/bash
#!/bin/bash
#SBATCH -J score1
#SBATCH --ntasks=15
#SBATCH --mem=30GB
#SBATCH -p skylake_0096
#SBATCH --qos=skylake_0096

module load singularity
singularity exec --bind /gpfs/data/fs71407/epiwielsch/ukb_2/application /gpfs/data/fs71407/epiwielsch/shortcake.sif Rscript plot_and_sort_${fac}.R


EOFSUBMIT

echo " --------------  submit plotting job for ${fac} now !!" 

chmod +x submit_qc_plots_${fac}.sh
sbatch submit_qc_plots_${fac}.sh
 

done<FAC_LIST
