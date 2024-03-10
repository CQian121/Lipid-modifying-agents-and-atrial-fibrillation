####1. Environment and big data###################################################

library(TwoSampleMR)
library(data.table)
library(plinkbinr)
library(ieugwasr)
##"g1000_eur" and "g1000_eas" were downloaded from https://ctg.cncr.nl/software/magma
##"finngen_R9_I9_AF.gz" was downloaded from https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_I9_AF.gz
##"logTG_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz" was downloaded from https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/
##"nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl.gz" was downloaded from https://csg.sph.umich.edu/willer/public/afib2018/

####2. Selecting IVs for drug targets############################################

regional.snps<-read.table("100kb_snps_EUR.csv",header = T,sep=',',stringsAsFactors = F)
regional.snps$genes<-ifelse(regional.snps$genes %in% c("ABCG5","ABCG8"),"ABCG5/ABCG8",regional.snps$genes)
ldlc.targets=c("LDLR","HMGCR","PCSK9","NPC1L1","APOB","ABCG5/ABCG8")
hdlc.targets="CETP"
tg.targets=c("LPL","ANGPTL3","APOC3")
lipid.genes=c(ldlc.targets,hdlc.targets,tg.targets)

iv.targets<-c()
for (i in 1:length(lipid.genes)) {
  targets.lipid<-subset(regional.snps,genes==lipid.genes[i]) %>% {
    if(lipid.genes[i] %in% ldlc.targets) {
      extract_outcome_data(outcomes="ieu-a-300",snps =.$snps,proxies = F,access_token = NULL)
    } else if(lipid.genes[i] %in% tg.targets) {
      extract_outcome_data(outcomes="ieu-a-302",snps =.$snps,proxies = F,access_token = NULL)
    } else {
      extract_outcome_data(outcomes="ieu-a-299",snps =.$snps,proxies = F,access_token = NULL)
    }
  } %>%
    convert_outcome_to_exposure() %>% subset(pval.exposure<5e-08)   #####significant SNPs with p<5e-08
  clump.data<-ld_clump(dplyr::tibble(rsid=targets.lipid$SNP, pval=targets.lipid$pval.exposure, id=targets.lipid$id.exposure),
                       clump_kb = 10000,
                       clump_r2 = 0.3,  ##clumped using r2<0.3 within ¡À10Mb
                       plink_bin = get_plink_exe(),
                       bfile = "./g1000_eur")  ##Reference data for Europeans from from Phase 3 of 1000 Genomes
  targets.lipid<-subset(targets.lipid,SNP %in% clump.data$rsid)
  if(lipid.genes[i] %in% c(ldlc.targets,tg.targets)){targets.lipid$beta.exposure=-targets.lipid$beta.exposure}
  iv.targets<-rbind(iv.targets,cbind(targets.lipid,genes=lipid.genes[i]))
}
rm(clump.data,targets.lipid)
iv.targets$F.statistic=(iv.targets$beta.exposure/iv.targets$se.exposure)^2   ##calculate F statistic
write.csv(iv.targets,"iv_targets.csv",quote = F,row.names = F)

####3. Positive control analyses############################################

MRresults.chd<-c()
for (i in 1:length(lipid.genes)) {
  targets.lipid<-subset(iv.targets,genes==lipid.genes[i])
  targets.chd<-extract_outcome_data(outcomes = "ieu-a-7",snps = targets.lipid$SNP,proxies = F,access_token = NULL) %>%
    harmonise_data(targets.lipid,.,action = 2) %>%
    subset(mr_keep==T & pval.exposure<pval.outcome)
  MRresults.chd<-rbind(MRresults.chd,
                       cbind(genes=lipid.genes[i],generate_odds_ratios(mr(targets.chd,method_list=c("mr_wald_ratio","mr_ivw")))))
}
rm(targets.lipid,targets.chd)
write.table(MRresults.chd,"MRresults_chd.txt",quote = F,row.names = F,sep="\t")

####4. MR analyses for AF in Nielsen et al (discovery cohort)##########################################

#####4.1. Main analyses#######################

MRresults.afnie<-c()
MRsensi.afnie<-c()
for (i in 1:length(lipid.genes)) {
  targets.lipid<-subset(iv.targets,genes==lipid.genes[i])
  targets.afnie<-extract_outcome_data(outcomes = "ebi-a-GCST006414",snps = targets.lipid$SNP,proxies = F,access_token = NULL) %>%
    harmonise_data(targets.lipid,.,action = 2) %>%
    subset(mr_keep==T & pval.exposure<pval.outcome)
  MRresults.afnie<-rbind(MRresults.afnie,
                         cbind(genes=lipid.genes[i],
                               generate_odds_ratios(mr(targets.afnie,method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_weighted_median")))))
  MRsensi.afnie<-rbind(MRsensi.afnie,
                       cbind(genes=lipid.genes[i],
                             hetero=mr_heterogeneity(targets.afnie)[nrow(mr_heterogeneity(targets.afnie)),8], ##heterogeneity
                             pleio=mr_pleiotropy_test(targets.afnie)[,7]))  ##pleiotropy
}
rm(targets.lipid,targets.afnie)
write.table(MRresults.afnie,"MRresults_afnie.txt",quote = F,row.names = F,sep = "\t")
write.table(MRsensi.afnie,"MRsensi_afnie.txt",quote = F,row.names = F,sep = "\t")

#####4.2. MR_PRESSO and leave-one-out analyses#################

iv.lpl<-subset(iv.targets,genes=="LPL")
lpl.afnie<-extract_outcome_data(outcomes = "ebi-a-GCST006414",snps = iv.lpl$SNP,proxies = F,access_token = NULL) %>%
  harmonise_data(iv.lpl,.,action = 2) %>%
  subset(mr_keep==T & pval.exposure<pval.outcome)
run_mr_presso(lpl.afnie)
leaveoneout.lpl.afnie<-generate_odds_ratios(mr_leaveoneout(lpl.afnie))
write.table(leaveoneout.lpl.afnie,"Leaveoneout_lpl_afnie.txt",quote = F,row.names = F,sep = "\t")

####5. MR analyses for AF in FinnGen (validation cohort)#####################################

#####5.1. Main analyses#######################

af.finngen<-fread("finngen_R9_I9_AF.gz",header = T,sep = "\t",stringsAsFactors = F)
af.fin<-subset(af.finngen,rsids %in% iv.lpl$SNP)

lpl.affin<-format_data(af.fin,type="outcome",snps = iv.lpl$SNP,snp_col = "rsids",effect_allele_col = "alt",other_allele_col = "ref",
                       eaf_col = "af_alt",beta_col = "beta",se_col = "sebeta",pval_col = "pval") %>%
  harmonise_data(iv.lpl,.,action = 2) %>%
  subset(mr_keep==T & pval.exposure<pval.outcome)
MRresults.affin<-generate_odds_ratios(mr(lpl.affin,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median")))

MRsensi.affin<-data.frame(hetero=mr_heterogeneity(lpl.affin)[nrow(mr_heterogeneity(lpl.affin)),8],
                          pleio=mr_pleiotropy_test(lpl.affin)[,7])
write.table(MRresults.affin,"MRresults_affin.txt",quote = F,row.names = F,sep="\t")
write.table(MRsensi.affin,"MRsensi_affin.txt",quote = F,row.names = F,sep = "\t")

#####5.2. MR_PRESSO and leave-one-out analyses################

run_mr_presso(lpl.affin)
leaveoneout.lpl.affin<-generate_odds_ratios(mr_leaveoneout(lpl.affin))
write.table(leaveoneout.lpl.affin,"Leaveoneout_lpl_affin.txt",quote = F,row.names = F,sep = "\t")

####6. Colocalization##################################

#####6.1 Nielsen's dataset############################

ciseQTL.lpl<-read.table("ciseqtl_LPL.txt",header = T,sep = "\t",stringsAsFactors = F)  ##Cis-eQTL data for LPL in subcutaneous adipose tissues from GETxV8
lpl.pos=c(8,19759228,19824769)
eqtl.data<-subset(ciseQTL.lpl,BP<=lpl.pos[3]+250e03 & BP>=lpl.pos[2]-250e03)
af.nie<-fread("nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl.gz",sep = "\t",header = T,check.names = F)
head(af.nie)
share.id=intersect(eqtl.data$SNP,af.nie$rs_dbSNP147)
gwas.data<-af.nie[match(share.id,af.nie$rs_dbSNP147),]
gwas.data$MAF=ifelse(gwas.data$Freq_A2<0.5,gwas.data$Freq_A2,1-gwas.data$Freq_A2)
gwas.data$MAF=ifelse(gwas.data$MAF==0,gwas.data$MAF+1e-08,gwas.data$MAF)
eqtl.data<-eqtl.data[eqtl.data$SNP %in% share.id,]
dataset1=list(pvalues=as.numeric(gwas.data$Pvalue),type="cc", s=60620/1030836, N=1030836)
dataset2=list(pvalues=eqtl.data$p, type="quant", N=581)
coloc.abf(dataset1,dataset2, MAF=gwas.data$MAF)

#####6.2 FinnGen dataset############################

eqtl.data<-subset(ciseQTL.lpl,BP<=lpl.pos[3]+250e03 & BP>=lpl.pos[2]-250e03)
share.id=intersect(eqtl.data$SNP,af.finngen$rsids)
gwas.data<-af.finngen[match(share.id,af.finngen$rsids),]
gwas.data$MAF=ifelse(gwas.data$af_alt<0.5,gwas.data$af_alt,1-gwas.data$af_alt)
gwas.data$MAF=ifelse(gwas.data$MAF==0,gwas.data$MAF+1e-08,gwas.data$MAF)
eqtl.data<-eqtl.data[eqtl.data$SNP %in% share.id,]
dataset1=list(pvalues=as.numeric(gwas.data$pval),type="cc", s=45766/237690, N=237690)
dataset2=list(pvalues=eqtl.data$p, type="quant", N=581)
coloc.abf(dataset1,dataset2, MAF=gwas.data$MAF)

####7. Validation in East Asian population##################################

regional.snps<-read.table("100kb_snps_EAS.csv",header = T,sep=',',stringsAsFactors = F)
ldlc.asian<-fread("logTG_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz",header = T,sep = "\t",stringsAsFactors = F)
head(ldlc.asian)
ldlc.asian<-subset(ldlc.asian,rsID %in% regional.snps$snps)
ldlc.asian<-format_data(ldlc.asian,type="exposure",snp_col = "rsID",effect_allele_col = "ALT",other_allele_col = "REF",
                        eaf_col = "POOLED_ALT_AF",beta_col = "EFFECT_SIZE",se_col = "SE",pval_col = "pvalue")

iv.lpl.asian<-subset(ldlc.asian,SNP %in% subset(regional.snps,genes=="LPL")$snps) %>%
  subset(pval.exposure<5e-08)
clump.data<-ld_clump(dplyr::tibble(rsid=iv.lpl.asian$SNP, pval=iv.lpl.asian$pval.exposure, id=iv.lpl.asian$id.exposure),
                     clump_kb = 10000,
                     clump_r2 = 0.3,
                     plink_bin = get_plink_exe(),
                     bfile = "./g1000_eas")  ##Reference data for East Asians from from Phase 3 of 1000 Genomes
iv.lpl.asian<-subset(iv.lpl.asian,SNP %in% clump.data$rsid)
iv.lpl.asian$beta.exposure=-iv.lpl.asian$beta.exposure
iv.lpl.asian$F.statistic=(iv.lpl.asian$beta.exposure/iv.lpl.asian$se.exposure)^2 

lpl.af.asian<-extract_outcome_data(outcomes = "bbj-a-71",snps = iv.lpl.asian$SNP,proxies = F,access_token = NULL) %>%
  harmonise_data(iv.lpl.asian,.,action = 2) %>%
  subset(mr_keep==T & pval.exposure<pval.outcome)
MRresults.afasian<-data.frame(generate_odds_ratios(mr(lpl.af.asian,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))))
mr_heterogeneity(lpl.af.asian)
mr_pleiotropy_test(lpl.af.asian)
write.table(MRresults.afasian,"MRresults_lplasian_af.txt",quote = F,row.names = F,sep="\t")
