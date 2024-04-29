####7. Mediation analyses#############################################

library(RMediation)
##"finngen_R9_I9_HYPTENS.gz" was downloaded from https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_I9_HYPTENS.gz

#####7.1. Identification of potential mediators###############################

##############################online datasets

lipcardio.gwasid=c("ieu-a-300","ieu-a-299","ieu-a-302","ieu-b-107","ieu-b-108",
                   "ieu-a-835","ebi-a-GCST006867","ieu-b-35","ieu-a-7")
iv.lipcardio.pre<-c()
for (i in 1:length(lipcardio.gwasid)) {
  gwas.lipcardio<-extract_instruments(lipcardio.gwasid[i],p1 = 5e-08,r2 = 0.001,access_token = NULL)
  iv.lipcardio.pre<-rbind(iv.lipcardio.pre,gwas.lipcardio)
}
iv.lipcardio.pre$F.statistic=(iv.lipcardio.pre$beta.exposure/iv.lipcardio.pre$se.exposure)^2
iv.lipcardio.pre<-subset(iv.lipcardio.pre,F.statistic>10)
table(iv.lipcardio.pre$id.exposure)
rm(gwas.lipcardio)

##############################local dataset: hypertension

hp.gwas<-fread("finngen_R9_I9_HYPTENS.gz",header = T,sep = "\t",stringsAsFactors = F)
head(hp.gwas)
iv.hp.pre<-subset(hp.gwas,pval<5e-08) %>%
  format_data(type = "exposure",snp_col = "rsids",effect_allele_col = "alt",other_allele_col = "ref",
              eaf_col = "af_alt",beta_col = "beta",se_col = "sebeta",pval_col = "pval")
clump.data<-ld_clump(dplyr::tibble(rsid=iv.hp.pre$SNP, pval=iv.hp.pre$pval.exposure, id=iv.hp.pre$id.exposure),
                     clump_kb = 10000,
                     clump_r2 = 0.001,
                     plink_bin = get_plink_exe(),
                     bfile = "./g1000_eur")
iv.hp.pre<-subset(iv.hp.pre,SNP %in% clump.data$rsid)
iv.hp.pre$F.statistic=(iv.hp.pre$beta.exposure/iv.hp.pre$se.exposure)^2
iv.hp.pre<-subset(iv.hp.pre,F.statistic>10)

################################MR analyses for AF before removing pleiotropic SNPs

MRresults.lipcardio.afnie<-c()
for (i in 1:length(lipcardio.gwasid)) {
  iv.lipcardio<-subset(iv.lipcardio.pre,id.exposure==lipcardio.gwasid[i])
  lipcardio.afnie<-extract_outcome_data(outcomes = "ebi-a-GCST006414",snps = iv.lipcardio$SNP,proxies = F,access_token = NULL) %>%
    harmonise_data(iv.lipcardio,.,action = 2) %>%
    subset(mr_keep==T & pval.exposure<pval.outcome)
  MRresults.lipcardio.afnie<-rbind(MRresults.lipcardio.afnie,
                                   cbind(traits=lipcardio.gwasid[i],generate_odds_ratios(mr(lipcardio.afnie,method_list=c("mr_wald_ratio","mr_ivw")))))
}

rm(lipcardio.afnie)
write.table(MRresults.lipcardio.afnie,"MRresults_lipcardio_afnie.txt",quote = F,sep="\t",row.names = F)

hp.afnie<-extract_outcome_data(outcomes = "ebi-a-GCST006414",snps = iv.hp.pre$SNP,proxies = F,access_token = NULL) %>%
  harmonise_data(iv.hp.pre,.,action = 2) %>%
  subset(mr_keep==T & pval.exposure<pval.outcome)
generate_odds_ratios(mr(hp.afnie,method_list=c("mr_wald_ratio","mr_ivw")))
rm(hp.afnie)

#####7.2. MR analyses for AF after removing pleiotropic SNPs###############################

##potential pleiotropic SNPs were removed using PhenoscannerV2 (BMI, Hypertension, and CHD)

removed.snps<-read.table("pheno_removed_snps.txt",header = T,sep = "\t",stringsAsFactors = F)
iv.lipcardio.final<-list(iv.bmi.pre,iv.hp.pre,iv.chd.pre)
names(iv.lipcardio.final)=c("BMI","Hypertension","CHD")
for (i in 1:length(iv.lipcardio.final)) {
  rm.id=match(subset(removed.snps,Trait==names(iv.lipcardio.final)[i])$SNP,iv.lipcardio.final[[i]]$SNP)
  iv.lipcardio.final[[i]]<-iv.lipcardio.final[[i]][-rm.id,]
}

MRresults.lipcardiofinal.afnie<-c()
for (i in 1:length(iv.lipcardio.final)) {
  iv.lipcardio<-iv.lipcardio.final[[i]]
  lipcardio.afnie<-extract_outcome_data(outcomes = "ebi-a-GCST006414",snps = iv.lipcardio$SNP,proxies = F,access_token = NULL) %>%
    harmonise_data(iv.lipcardio,.,action = 2) %>%
    subset(mr_keep==T & pval.exposure<pval.outcome)
  MRresults.lipcardiofinal.afnie<-rbind(MRresults.lipcardiofinal.afnie,
                                        cbind(traits=names(iv.lipcardio.final)[i],pleiotropy=mr_pleiotropy_test(lipcardio.afnie)[7],
                                              generate_odds_ratios(mr(lipcardio.afnie,method_list=c("mr_wald_ratio","mr_ivw")))))
}
rm(iv.lipcardio,lipcardio.afnie)
write.table(MRresults.lipcardiofinal.afnie,"MRresults_lipcardiofinal_afnie.txt",quote = F,sep="\t",row.names = F)

#####7.3. Two-step mediation analyses for BMI, hypertension, and CHD###################

mr_mediation<-function(x,y){
  MR.A<-generate_odds_ratios(mr(x,method_list=c("mr_ivw")))
  MR.B<-generate_odds_ratios(mv_multiple(y)$result)[2,]
  mediation<-medci(mu.x = MR.A$b,mu.y = MR.B$b,type = "asymp",
                   se.x = MR.A$se,se.y = MR.B$se)
  Beta=mediation$Estimate  ###indirect effect
  Se=mediation$SE
  Proportion=Beta/MR.C$b   ###mediation proportion
  Zvalue=Beta/Se
  Pvalue=2*pnorm(abs(Zvalue),lower.tail = F)
  return(c(Beta=Beta,Se=Se,Prop=Proportion,Zval=Zvalue,Pval=Pvalue))
}
MR.C=subset(MRresults.afnie,genes=="LPL" & method=="Inverse variance weighted")   ###total effect

############################BMI
lpl.bmi<-extract_outcome_data(outcomes = "ieu-a-835",snps = iv.lpl$SNP,proxies = F,access_token = NULL) %>%
  harmonise_data(iv.lpl,.,action = 2) %>% 
  subset(mr_keep==T & pval.exposure<pval.outcome)
generate_odds_ratios(mr(lpl.bmi,method_list = c("mr_ivw")))
####the above result shows that BMI is not a mediator

############################Hypertension
lpl.hp<-subset(hp.gwas,rsids %in% iv.lpl$SNP) %>%
  format_data(type = "outcome",snps = iv.lpl$SNP,snp_col = "rsids",effect_allele_col = "alt",other_allele_col = "ref",
              eaf_col = "af_alt",beta_col = "beta",se_col = "sebeta",pval_col = "pval") %>%
  harmonise_data(iv.lpl,.,action = 2) %>% 
  subset(mr_keep==T & pval.exposure<pval.outcome)
generate_odds_ratios(mr(lpl.hp,method_list = c("mr_ivw"))) ###β1

lplhp.snp=unique(c(iv.lpl$SNP,iv.lipcardio.final[["Hypertension"]]$SNP))
iv.tsmv1<-extract_outcome_data("ieu-a-302",snps = lplhp.snp,proxies = F,access_token = NULL) %>%
  convert_outcome_to_exposure()
iv.tsmv1$beta.exposure=-iv.tsmv1$beta.exposure
iv.tsmv2<-subset(hp.gwas,rsids %in% lplhp.snp) %>%
  format_data(type = "exposure",snp_col = "rsids",effect_allele_col = "alt",other_allele_col = "ref",
              eaf_col = "af_alt",beta_col = "beta",se_col = "sebeta",pval_col = "pval",chr_col = "#chrom",pos_col = "pos")
iv.tsmv2<-iv.tsmv2[,names(iv.tsmv1)]
iv.tsmv<-rbind(iv.tsmv1,iv.tsmv2)
lplhp.afnie<-extract_outcome_data(outcomes = "ebi-a-GCST006414",snps = iv.tsmv$SNP,proxies = F,access_token = NULL) %>%
  mv_harmonise_data(iv.tsmv,.)
generate_odds_ratios(mv_multiple(lplhp.afnie)$result) ###β2

hp.mediation<-mr_mediation(lpl.hp,lplhp.afnie) ###get β1*β2
hp.mediation
exp(hp.mediation[1])
exp(hp.mediation[1]-1.96*hp.mediation[2])
exp(hp.mediation[1]+1.96*hp.mediation[2])

############################CHD
lpl.chd<-extract_outcome_data(outcomes = "ieu-a-7",snps = iv.lpl$SNP,proxies = F,access_token = NULL) %>%
  harmonise_data(iv.lpl,.,action = 2) %>% 
  subset(mr_keep==T & pval.exposure<pval.outcome)
generate_odds_ratios(mr(lpl.chd,method_list = c("mr_ivw"))) ###β1

lplchd.snp=unique(c(iv.lpl$SNP,iv.lipcardio.final[["CHD"]]$SNP))
iv.tsmv1<-extract_outcome_data("ieu-a-302",snps = lplchd.snp,proxies = F,access_token = NULL) %>%
  convert_outcome_to_exposure()
iv.tsmv1$beta.exposure=-iv.tsmv1$beta.exposure
iv.tsmv2<-extract_outcome_data("ieu-a-7",snps = lplchd.snp,proxies = F,access_token = NULL) %>%
  convert_outcome_to_exposure()
identical(names(iv.tsmv1),names(iv.tsmv2))
iv.tsmv<-rbind(iv.tsmv1,iv.tsmv2)
lplchd.afnie<-extract_outcome_data(outcomes = "ebi-a-GCST006414",snps = iv.tsmv$SNP,proxies = F,access_token = NULL) %>%
  mv_harmonise_data(iv.tsmv,.)
generate_odds_ratios(mv_multiple(lplchd.afnie)$result) ###β2

chd.mediation<-mr_mediation(lpl.chd,lplchd.afnie) ###get β1*β2
chd.mediation
exp(chd.mediation[1])
exp(chd.mediation[1]-1.96*chd.mediation[2])
exp(chd.mediation[1]+1.96*chd.mediation[2])

############################Hypertension+CHD
lplhpchd.snp=unique(c(iv.lpl$SNP,iv.lipcardio.final[["Hypertension"]]$SNP,iv.lipcardio.final[["CHD"]]$SNP))
iv.tsmv1<-extract_outcome_data("ieu-a-302",snps = lplhpchd.snp,proxies = F,access_token = NULL) %>%
  convert_outcome_to_exposure()
iv.tsmv1$beta.exposure=-iv.tsmv1$beta.exposure
iv.tsmv2<-subset(hp.gwas,rsids %in% lplhpchd.snp) %>%
  format_data(type = "exposure",snp_col = "rsids",effect_allele_col = "alt",other_allele_col = "ref",
              eaf_col = "af_alt",beta_col = "beta",se_col = "sebeta",pval_col = "pval",chr_col = "#chrom",pos_col = "pos")
iv.tsmv2<-iv.tsmv2[,names(iv.tsmv1)]
iv.tsmv3<-extract_outcome_data("ieu-a-7",snps = lplhpchd.snp,proxies = F,access_token = NULL) %>%
  convert_outcome_to_exposure()
identical(names(iv.tsmv1),names(iv.tsmv3))
iv.tsmv<-rbind(iv.tsmv1,iv.tsmv2,iv.tsmv3)
lplhpchd.afnie<-extract_outcome_data(outcomes = "ebi-a-GCST006414",snps = iv.tsmv$SNP,proxies = F,access_token = NULL) %>%
  mv_harmonise_data(iv.tsmv,.)

b.direct=generate_odds_ratios(mv_multiple(lplhpchd.afnie)$result)[2,]$b
se.direct=generate_odds_ratios(mv_multiple(lplhpchd.afnie)$result)[2,]$se
b.total=MR.C$b
se.total=MR.C$se
b.indirect=b.total-b.direct

prop.combined=b.indirect/b.total
se.indirect=sqrt(se.total^4*MR.C$nsnp+se.direct^4*generate_odds_ratios(mv_multiple(lplhpchd.afnie)$result)[2,]$nsnp)
p.indirect=2*pnorm(abs(b.indirect/se.indirect),lower.tail = F)
exp(b.indirect)
exp(b.indirect-1.96*se.indirect)
exp(b.indirect+1.96*se.indirect)


