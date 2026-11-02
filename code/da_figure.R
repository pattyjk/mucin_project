#!/usr/bin/env Rscript
# ANCOM-BC2 mucin differential-abundance figure.
# Usage: Rscript da_figure.R results taxonomy.tsv   (results = directory written by analysis_20C_mucin.R)
suppressMessages({library(tidyverse);library(patchwork)})
args <- commandArgs(trailingOnly=TRUE); R <- args[1]; tx <- read.delim(args[2],row.names=1)
X <- read.delim(file.path(R,"asv_table_20C_carbon_antibiotic.txt"),row.names=1,check.names=FALSE)
tax <- tx[rownames(X),"Taxon"]
# same lowest-named-rank label as the pipeline, plus the rank it came from
lab <- t(sapply(strsplit(tax,"; "),function(p){p<-p[!grepl("__$|uncultured",p)]; x<-tail(p,1); c(sub("^[a-z]__","",x),substr(x,1,1))}))
rk <- tibble(taxon=lab[,1],rank=lab[,2]) %>% count(taxon,rank) %>% group_by(taxon) %>% slice_max(n,n=1,with_ties=FALSE) %>% ungroup() %>% select(-n)
pretty <- function(t,r){ t2 <- gsub("_"," ",t); ifelse(r %in% c("g","s"), t2, ifelse(t=="Bacteria","Unclassified Bacteria",paste0("Unclassified ",t2)))}
mk <- function(f,keep=NULL){ d <- read.csv(file.path(R,f))
  if(!is.null(keep)) d <- d %>% filter(taxon %in% keep)
  d %>% filter(diff,pss) %>% left_join(rk,by="taxon") %>%
    mutate(label=pretty(taxon,rank), italic=rank %in% c("g","s") & grepl("^[A-Z][a-z]+( [a-z]+)?$",label),  # placeholder clades (WD260, BSV26) stay upright
           dir=ifelse(lfc>0,"Higher with mucin","Lower with mucin"),
           raw_agrees=factor(ifelse(sign(meanRA_Yes_pct-meanRA_No_pct)==sign(lfc),"agrees","opposite"),levels=c("agrees","opposite")),
           lo=lfc-1.96*se, hi=lfc+1.96*se)}
sens <- read.csv(file.path(R,"ancombc2_antibiotic_PABD4_sensitivity.csv"))
dA <- mk("ancombc2_genus_mucin_antibiotic.csv",keep=sens$taxon[sens$robust])
dC <- mk("ancombc2_genus_mucin_carbon.csv")
cols <- c("Higher with mucin"="#114477","Lower with mucin"="#AA4455")
panel <- function(d,title){
  d <- d %>% arrange(lfc) %>% mutate(lab_expr=ifelse(italic,paste0("italic('",label,"')"),paste0("'",label,"'")),
                                     y=factor(lab_expr,levels=lab_expr))
  ggplot(d,aes(lfc,y,colour=dir))+geom_vline(xintercept=0,colour="grey60",linewidth=0.3)+
    geom_errorbar(aes(xmin=lo,xmax=hi),width=0,linewidth=0.5,orientation="y")+
    geom_point(aes(shape=raw_agrees),size=2,fill="white",stroke=0.7)+
    scale_shape_manual(values=c(agrees=16,opposite=21),limits=c("agrees","opposite"),name="Direction of raw mean\nabundance change",drop=FALSE)+
    scale_colour_manual(values=cols,limits=names(cols),name=NULL)+scale_y_discrete(labels=function(x) parse(text=x))+
    labs(x="Log fold change, mucin vs. no mucin\n(ANCOM-BC2, bias-corrected, 95% CI)",y=NULL,title=title)+
    theme_classic(base_size=9)+theme(plot.title=element_text(size=9,face="bold"),axis.text.y=element_text(size=7.5),legend.position="bottom",legend.box="vertical",legend.margin=margin(0,0,0,0))}
pA <- panel(dA,sprintf("Antibiotics (Synthetic Amphibian Skin Media)\n%d genera robust to PABD4 exclusion",nrow(dA)))
pC <- panel(dC,sprintf("Carbon sources\n%d genera, adjusted for carbon source",nrow(dC)))
P <- pA | (pC+theme(legend.position="none"))  # single legend under the antibiotic panel (it carries every level)
ggsave(file.path(R,"differential_abundance_mucin.png"),plot=P,width=9,height=6.2,dpi=300)
write.csv(bind_rows(dA %>% mutate(Experiment="Antibiotics"),dC %>% mutate(Experiment="CarbonTemp")) %>%
  select(Experiment,taxon,label,rank,lfc,se,lo,hi,q,meanRA_No_pct,meanRA_Yes_pct,raw_agrees),file.path(R,"differential_abundance_mucin_plotted.csv"),row.names=FALSE)
cat("antibiotic",nrow(dA),"(opposite raw:",sum(dA$raw_agrees=="opposite"),") carbon",nrow(dC),"(opposite raw:",sum(dC$raw_agrees=="opposite"),")\n")
