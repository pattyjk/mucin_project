#!/usr/bin/env Rscript
# 20C mucin x antibiotic / carbon-source 16S analysis
# Usage: Rscript analysis_20C_mucin.R asv_table.txt synthcom_map.txt taxonomy.tsv [outdir]
# Requires: vegan, tidyverse, phyloseq, ANCOMBC, ggh4x, patchwork
suppressMessages({library(vegan);library(tidyverse);library(phyloseq);library(ANCOMBC);library(ggh4x);library(patchwork)})
args <- commandArgs(trailingOnly=TRUE)
f_asv <- args[1]; f_map <- args[2]; f_tax <- args[3]; out <- ifelse(length(args)>3,args[4],"results")
dir.create(out,showWarnings=FALSE,recursive=TRUE)
RARE <- 1000; set.seed(42)

## ---- 1. Sample selection: 20C, CarbonTemp + Antibiotics plates, no BlankWell, no Diversity plate
sm <- read.delim(f_map,check.names=FALSE)
md <- sm %>% filter(Temp==20, Experiment %in% c("CarbonTemp","Antibiotics"), is.na(Antibiotic) | Antibiotic!="BlankWell") %>%
  mutate(Treatment=ifelse(Experiment=="Antibiotics",Antibiotic,CarbonSource), Mucin=factor(Mucin,c("No","Yes")))
A <- read.delim(f_asv,row.names=1,check.names=FALSE); A <- A[,!grepl("^Unnamed",colnames(A)),drop=FALSE]
tx <- read.delim(f_tax,row.names=1)
X <- A[,md$SampleID]; X <- X[rowSums(X)>0,]
## ---- 2. Taxonomic filter (decontamination was done upstream)
tax <- tx[rownames(X),"Taxon"]
bad <- is.na(tax)|grepl("Chloroplast|Mitochondria|Eukaryota|Unassigned|Archaea",tax)|!grepl("d__Bacteria",tax)
X <- X[!bad,]; tax <- tax[!bad]
md$depth <- colSums(X)[md$SampleID]
cat("Samples:",nrow(md)," ASVs kept:",nrow(X)," removed:",sum(bad),"\n")
write.table(cbind(`OTU ID`=rownames(X),X),file.path(out,"asv_table_20C_carbon_antibiotic.txt"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(md %>% select(-depth),file.path(out,"map_20C_carbon_antibiotic.txt"),sep="\t",quote=FALSE,row.names=FALSE,na="")

## ---- 3. Rarefy + alpha diversity
keep <- md$SampleID[md$depth>=RARE]; cat("Dropped <",RARE," reads:",setdiff(md$SampleID,keep),"\n")
XR <- t(rrarefy(t(X[,keep]),RARE))
al <- tibble(SampleID=keep,Observed=colSums(XR>0),Shannon=diversity(t(XR))) %>% left_join(md,by="SampleID")
write.csv(al %>% select(SampleID,Experiment,Treatment,Mucin,depth,Observed,Shannon),file.path(out,"alpha_rarefied1000.csv"),row.names=FALSE)
st <- al %>% filter(Treatment!="Mucin")   # mucin-as-sole-carbon wells have no -mucin pair

res_alpha <- list(); res_perm <- list(); disp <- list(); dtests <- list(); sep <- list()
for(ex in c("Antibiotics","CarbonTemp")){
  d <- st %>% filter(Experiment==ex)
  for(v in c("Observed","Shannon")){ an <- anova(lm(as.formula(paste0("rank(",v,")~Mucin*Treatment")),d))
    res_alpha[[length(res_alpha)+1]] <- tibble(Experiment=ex,metric=v,term=rownames(an)[1:3],F=an$`F value`[1:3],p=an$`Pr(>F)`[1:3])}
  for(dist in c("bray","jaccard")){
    D <- vegdist(t(XR[,d$SampleID]),dist,binary=(dist=="jaccard"))
    p <- adonis2(D~Mucin*Treatment,data=d,permutations=9999,by="terms")
    res_perm[[length(res_perm)+1]] <- tibble(Experiment=ex,distance=dist,term=rownames(p)[1:3],R2=p$R2[1:3],F=p$F[1:3],p=p$`Pr(>F)`[1:3])
    ## ---- beta dispersion: (a) mucin-level centroid, (b) treatment x mucin centroid
    b1 <- betadisper(D,d$Mucin,type="centroid"); t1 <- permutest(b1,permutations=9999)$tab
    b2 <- betadisper(D,interaction(d$Mucin,d$Treatment,drop=TRUE),type="centroid")
    dd <- tibble(SampleID=d$SampleID,Experiment=ex,distance=dist,Treatment=d$Treatment,Mucin=d$Mucin,
                 dist_mucin_centroid=b1$distances,dist_group_centroid=b2$distances)
    an2 <- anova(lm(rank(dist_group_centroid)~Mucin+Treatment,dd)); disp[[length(disp)+1]] <- dd
    dtests[[length(dtests)+1]] <- tibble(Experiment=ex,distance=dist,
      mean_No=mean(b1$distances[d$Mucin=="No"]),mean_Yes=mean(b1$distances[d$Mucin=="Yes"]),F_mucinCentroid=t1$F[1],p_mucinCentroid=t1$`Pr(>F)`[1],
      grp_mean_No=mean(dd$dist_group_centroid[dd$Mucin=="No"]),grp_mean_Yes=mean(dd$dist_group_centroid[dd$Mucin=="Yes"]),
      F_groupCentroid_mucin=an2$`F value`[1],p_groupCentroid_mucin=an2$`Pr(>F)`[1])}
  ## ---- treatment separation within each mucin level
  for(mu in c("No","Yes")){ dm <- d %>% filter(Mucin==mu); D <- vegdist(t(XR[,dm$SampleID]),"bray")
    p <- adonis2(D~Treatment,data=dm,permutations=9999); M <- as.matrix(D); same <- outer(dm$Treatment,dm$Treatment,"=="); ut <- upper.tri(M)
    sep[[length(sep)+1]] <- tibble(Experiment=ex,Mucin=mu,R2_treatment=p$R2[1],p=p$`Pr(>F)`[1],
      mean_BC_between_treatments=mean(M[ut&!same]),mean_BC_within_treatment=mean(M[ut&same]))}
  ## ---- PCoA coordinates for plotting
  dA <- al %>% filter(Experiment==ex); pc <- cmdscale(vegdist(t(XR[,dA$SampleID]),"bray"),k=2,eig=TRUE); pe <- pc$eig[pc$eig>0]
  write.csv(data.frame(SampleID=dA$SampleID,PC1=pc$points[,1],PC2=pc$points[,2],Treatment=dA$Treatment,Mucin=dA$Mucin,
                       ve1=100*pc$eig[1]/sum(pe),ve2=100*pc$eig[2]/sum(pe)),file.path(out,paste0("pcoa_",ex,".csv")),row.names=FALSE)
}
write.csv(bind_rows(res_alpha),file.path(out,"alpha_rank_anova.csv"),row.names=FALSE)
write.csv(bind_rows(res_perm),file.path(out,"permanova_by_experiment.csv"),row.names=FALSE)
write.csv(bind_rows(disp),file.path(out,"betadisper_distances.csv"),row.names=FALSE)
write.csv(bind_rows(dtests),file.path(out,"betadisper_tests.csv"),row.names=FALSE)
write.csv(bind_rows(sep),file.path(out,"betadisper_treatment_separation_by_mucin.csv"),row.names=FALSE)

## combined model, treatment nested in experiment
D <- vegdist(t(XR[,st$SampleID]),"bray")
pc_all <- adonis2(D~Experiment+Mucin+Experiment:Mucin+Experiment:Treatment+Experiment:Treatment:Mucin,data=st,permutations=9999,by="terms")
write.csv(data.frame(term=rownames(pc_all),pc_all,check.names=FALSE),file.path(out,"permanova_combined.csv"),row.names=FALSE)

## per-treatment mucin effect
pt <- bind_rows(lapply(split(st,list(st$Experiment,st$Treatment),drop=TRUE),function(d){
  r <- adonis2(vegdist(t(XR[,d$SampleID]),"bray")~Mucin,data=d,permutations=9999)
  tibble(Experiment=d$Experiment[1],Treatment=d$Treatment[1],n_No=sum(d$Mucin=="No"),n_Yes=sum(d$Mucin=="Yes"),
         Obs_No=median(d$Observed[d$Mucin=="No"]),Obs_Yes=median(d$Observed[d$Mucin=="Yes"]),
         p_Obs=wilcox.test(Observed~Mucin,d,exact=FALSE)$p.value,R2_bray=r$R2[1],p_bray=r$`Pr(>F)`[1])})) %>%
  group_by(Experiment) %>% mutate(q_Obs=p.adjust(p_Obs,"BH"),q_bray=p.adjust(p_bray,"BH")) %>% ungroup()
write.csv(pt,file.path(out,"mucin_effect_per_treatment.csv"),row.names=FALSE)

## ---- 4. Taxon labels: nearest named rank
rank_lab <- function(tx, lvl){ sapply(strsplit(tx,"; "),function(p){ p <- sub("^[a-z]__","",p); p <- c(p,rep("",7))[1:7]
  ok <- function(x) x!="" & !grepl("^uncultured|^unidentified|^metagenome",x,ignore.case=TRUE)
  if(ok(p[lvl])) return(p[lvl]); k <- max(c(0,which(sapply(p[1:(lvl-1)],ok)))); if(k==0) "Unassigned" else paste0("Unclassified ",p[k])})}
glab <- sapply(strsplit(tax,"; "),function(p){p<-p[!grepl("__$|uncultured",p)]; sub("^[a-z]__","",tail(p,1))})

## ---- 5. ANCOM-BC2 genus-level mucin effect, per experiment (covariate = treatment)
for(ex in c("Antibiotics","CarbonTemp")){
  d <- md %>% filter(Experiment==ex,Treatment!="Mucin",depth>=RARE)
  G <- rowsum(as.matrix(X[,d$SampleID]),glab); sd <- d %>% column_to_rownames("SampleID")
  ps <- phyloseq(otu_table(G,taxa_are_rows=TRUE),sample_data(sd)); set.seed(1)
  r <- ancombc2(data=ps,fix_formula="Mucin + Treatment",group=NULL,p_adj_method="BH",prv_cut=0.10,lib_cut=0,
                struc_zero=FALSE,alpha=0.05,n_cl=1,verbose=FALSE)$res %>%
    select(taxon,lfc=lfc_MucinYes,se=se_MucinYes,q=q_MucinYes,diff=diff_MucinYes,pss=passed_ss_MucinYes) %>% arrange(q)
  Gr <- sweep(G,2,colSums(G),"/")
  r$meanRA_No_pct <- round(rowMeans(Gr[r$taxon,sd$Mucin=="No"])*100,3); r$meanRA_Yes_pct <- round(rowMeans(Gr[r$taxon,sd$Mucin=="Yes"])*100,3)
  write.csv(as.data.frame(r),file.path(out,paste0("ancombc2_genus_mucin_",ifelse(ex=="Antibiotics","antibiotic","carbon"),".csv")),row.names=FALSE)
  cat(ex,"ANCOM-BC2: tested",nrow(r),"significant & passed sensitivity",sum(r$diff&r$pss),"\n")}

## ---- 6. Stacked bars: independent top-N per experiment
pal <- c("#771155","#CC99BB","#114477","#4477AA","#117777","#44AAAA","#77CCCC","#117744","#44AA77","#88CCAA","#777711","#AAAA44","#DDDD77","#774411","#AA7744","#DDAA77","#771122","#AA4455","#DD7788","#41AB5D","#252525","#525252","#737373","#969696")
shortn <- c(NoAB="No AB",Ampicillin="Amp",Chloramphenicol="Chlor",Gentamycin="Gent",Neomycin="Neo",PSN="PSN",Streptomycin="Strep",Sulfamethazine="Sulfa",Trimethoprim="Trim",
  `N-acetyl glucosamine`="GlcNAc",`N-acetylneuraminic acid`="Neu5Ac",VitaminMix="Vit. mix",YE="Yeast ext.",Mucin="Mucin only")
ab <- c("NoAB","Ampicillin","Chloramphenicol","Gentamycin","Neomycin","PSN","Streptomycin","Sulfamethazine","Trimethoprim")
cs <- c("Fucose","Galactose","Glucose","Glycerol","Mannitol","Methanol","N-acetyl glucosamine","N-acetylneuraminic acid","Tryptone","VitaminMix","YE","Mucin")
lab2 <- function(x) ifelse(x %in% names(shortn), shortn[x], x)
lv <- unname(lab2(c(ab,cs)))
mdk <- md %>% filter(depth>=RARE); rel <- sweep(as.matrix(X[,mdk$SampleID]),2,colSums(X[,mdk$SampleID]),"/")
mkpanel <- function(R,ex,k,nm,title){
  ids <- mdk$SampleID[mdk$Experiment==ex]; top <- names(sort(rowMeans(R[,ids]),decreasing=TRUE))[1:k]
  R2 <- rbind(R[top,ids],Other=colSums(R[!rownames(R)%in%top,ids,drop=FALSE]))
  df <- as.data.frame(t(R2)) %>% rownames_to_column("SampleID") %>% pivot_longer(-SampleID,names_to="Taxon",values_to="RA") %>% left_join(mdk,by="SampleID") %>%
    group_by(Treatment,Mucin,Taxon) %>% summarise(RA=mean(RA),n=n(),.groups="drop") %>%
    mutate(Treatment=factor(unname(lab2(Treatment)),levels=lv),Mucin=factor(ifelse(Mucin=="Yes","+","-"),levels=c("-","+")),
           Taxon=factor(Taxon,levels=c(top,"Other")),Experiment=title)
  cols <- setNames(c(pal[1:k],"#D9D9D9"),c(top,"Other"))
  p <- ggplot(df,aes(Mucin,RA,fill=Taxon))+geom_col(width=0.9)+facet_grid(.~Treatment,scales="free_x",space="free_x",switch="x")+
    scale_fill_manual(values=cols,name=paste0("Top ",k," ",ifelse(nm=="Genus","genera","families")))+scale_y_continuous(expand=c(0,0),labels=scales::percent)+
    labs(x="Mucin (- / +)",y="Mean relative abundance",title=title)+theme_classic(base_size=9)+
    theme(strip.placement="outside",strip.background=element_blank(),strip.text.x=element_text(size=7,angle=90,hjust=1,vjust=0.5),panel.spacing.x=unit(1.5,"pt"),
          plot.title=element_text(size=9,face="bold"),legend.key.size=unit(8,"pt"),legend.text=element_text(size=6.5,face=ifelse(nm=="Genus","italic","plain")))
  list(p=p,df=df)}
for(L in c(Family=5,Genus=6)){ nm <- names(which(c(Family=5,Genus=6)==L)); R <- rowsum(rel,rank_lab(tax,L))
  a1 <- mkpanel(R,"Antibiotics",20,nm,"Antibiotics (Synthetic Amphibian Skin Media)"); c1 <- mkpanel(R,"CarbonTemp",15,nm,"Carbon sources")
  ggsave(file.path(out,paste0("stackedbar_independent_",tolower(nm),".png")),plot=a1$p+c1$p+plot_layout(widths=c(18,23)),width=12,height=4.6,dpi=300)
  write.csv(bind_rows(a1$df,c1$df),file.path(out,paste0("stackedbar_independent_",tolower(nm),"_means.csv")),row.names=FALSE)}
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo_R.txt"))
cat("Done ->",out,"\n")
