#!/usr/bin/env Rscript
# Compare R2A control wells (no cortisol) with Synthetic Amphibian Skin Media (SASM) no-antibiotic wells, -/+ mucin.
# Usage: Rscript r2a_vs_sasm.R asv_table_synth_com.txt cort_synth_com_map.txt taxonomy_r2a.tsv \
#                              results/asv_table_20C_carbon_antibiotic.txt results/map_20C_carbon_antibiotic.txt taxonomy.tsv outdir
suppressMessages({library(vegan);library(tidyverse);library(patchwork)})
a <- commandArgs(trailingOnly=TRUE); out <- a[7]; dir.create(out,showWarnings=FALSE,recursive=TRUE); set.seed(42); RARE <- 1000
rm2 <- read.delim(a[2],check.names=FALSE)
Ar <- read.delim(a[1],row.names=1,check.names=FALSE); Ar <- Ar[,!grepl("^Unnamed",colnames(Ar)),drop=FALSE]
ctl <- intersect(rm2$SampleID[rm2$Type=="Control"],colnames(Ar))          # exclude cortisol wells
Ar <- Ar[,ctl]; Ar <- Ar[rowSums(Ar)>0,]
Xs <- read.delim(a[4],row.names=1,check.names=FALSE); ms <- read.delim(a[5],check.names=FALSE)
sid <- ms$SampleID[ms$Experiment=="Antibiotics" & ms$Treatment=="NoAB"]  # SASM, no antibiotic
Xs <- Xs[,sid]; Xs <- Xs[rowSums(Xs)>0,]
txr <- read.delim(a[3],row.names=1); txs <- read.delim(a[6],row.names=1)
ids <- union(rownames(Ar),rownames(Xs))                                    # ASV IDs are sequence hashes -> directly comparable
X <- matrix(0,length(ids),length(ctl)+length(sid),dimnames=list(ids,c(ctl,sid)))
X[rownames(Ar),ctl] <- as.matrix(Ar); X[rownames(Xs),sid] <- as.matrix(Xs)
tax <- ifelse(ids %in% rownames(txr), txr[ids,"Taxon"], txs[ids,"Taxon"])
bad <- is.na(tax)|grepl("Chloroplast|Mitochondria|Eukaryota|Unassigned",tax)|!grepl("d__Bacteria|d__Archaea",tax)
X <- X[!bad,,drop=FALSE]; tax <- tax[!bad]
md <- tibble(SampleID=colnames(X), Group=c(rep("R2A",length(ctl)), ifelse(ms$Mucin[match(sid,ms$SampleID)]=="Yes","SASM + mucin","SASM"))) %>%
  mutate(Group=factor(Group,c("R2A","SASM","SASM + mucin")), depth=colSums(X)[SampleID])
inR <- rowSums(X[,ctl,drop=FALSE])>0; inS <- rowSums(X[,sid,drop=FALSE])>0
ov <- tibble(R2A_only=sum(inR&!inS),SASM_only=sum(!inR&inS),shared=sum(inR&inS),
             R2A_reads_in_shared=sum(X[inR&inS,ctl])/sum(X[,ctl]),SASM_reads_in_shared=sum(X[inR&inS,sid])/sum(X[,sid]))
write.csv(ov,file.path(out,"r2a_sasm_asv_overlap.csv"),row.names=FALSE)
## alpha
XR <- t(rrarefy(t(X[,md$depth>=RARE]),RARE)); md <- md %>% filter(depth>=RARE)
al <- md %>% mutate(Observed=colSums(XR>0)[SampleID],Shannon=diversity(t(XR))[SampleID])
write.csv(al,file.path(out,"r2a_sasm_alpha.csv"),row.names=FALSE)
kw <- bind_rows(lapply(c("Observed","Shannon"),function(v){k<-kruskal.test(al[[v]],al$Group)
  pw <- pairwise.wilcox.test(al[[v]],al$Group,p.adjust.method="BH",exact=FALSE)$p.value
  tibble(metric=v,KW_chi2=k$statistic,KW_p=k$p.value,
         p_R2A_vs_SASM=pw["SASM","R2A"],p_R2A_vs_SASMmucin=pw["SASM + mucin","R2A"],p_SASM_vs_SASMmucin=pw["SASM + mucin","SASM"])}))
write.csv(kw,file.path(out,"r2a_sasm_alpha_tests.csv"),row.names=FALSE)
## composition: ASV-level Bray/Jaccard and genus-level Bray (robust to ASV-calling differences between runs)
glab <- sapply(strsplit(tax,"; "),function(p){p<-p[!grepl("__$|uncultured",p)]; sub("^[a-z]__","",tail(p,1))})
GR <- rowsum(XR,glab[match(rownames(XR),rownames(X))])
Ds <- list(`Bray-Curtis (ASV)`=vegdist(t(XR),"bray"),`Jaccard (ASV)`=vegdist(t(XR>0),"jaccard"),`Bray-Curtis (genus)`=vegdist(t(GR),"bray"))
pairs <- list(c("R2A","SASM"),c("R2A","SASM + mucin"),c("SASM","SASM + mucin"))
perm <- bind_rows(lapply(names(Ds),function(n){D<-Ds[[n]]; g<-adonis2(D~Group,data=al,permutations=9999)
  bd <- permutest(betadisper(D,al$Group),permutations=9999)$tab
  pw <- bind_rows(lapply(pairs,function(p){k<-al$Group %in% p; s<-al[k,]
    r<-adonis2(as.dist(as.matrix(D)[k,k])~Group,data=s,permutations=9999); tibble(comparison=paste(p,collapse=" vs "),R2=r$R2[1],p=r$`Pr(>F)`[1])})) %>% mutate(p_BH=p.adjust(p,"BH"))
  bind_rows(tibble(distance=n,comparison="overall (3 groups)",R2=g$R2[1],p=g$`Pr(>F)`[1],p_BH=NA,betadisper_p=bd$`Pr(>F)`[1]),pw %>% mutate(distance=n))}))
write.csv(perm,file.path(out,"r2a_sasm_permanova.csv"),row.names=FALSE)
disp <- bind_rows(lapply(names(Ds),function(n){b<-betadisper(Ds[[n]],al$Group); tibble(distance=n,SampleID=al$SampleID,Group=al$Group,dist_to_centroid=b$distances)}))
write.csv(disp,file.path(out,"r2a_sasm_dispersion.csv"),row.names=FALSE)
## figures
cols <- c(R2A="#771155",SASM="#B0B0B0",`SASM + mucin`="#114477")
theme_set(theme_classic(base_size=9)+theme(plot.title=element_text(size=9,face="bold"),legend.position="none"))
dot <- function(v,yl) ggplot(al,aes(Group,.data[[v]],fill=Group))+geom_jitter(width=0.12,height=0,shape=21,size=2,stroke=0.3)+
  stat_summary(fun=median,geom="crossbar",width=0.45,linewidth=0.3)+scale_fill_manual(values=cols)+labs(x=NULL,y=yl)
pc <- cmdscale(Ds[[1]],k=2,eig=TRUE); ve <- round(100*pc$eig[1:2]/sum(pc$eig[pc$eig>0]))
pcd <- al %>% mutate(PC1=pc$points[SampleID,1],PC2=pc$points[SampleID,2])
pP <- ggplot(pcd,aes(PC1,PC2,fill=Group))+geom_point(shape=21,size=2.3,stroke=0.3)+scale_fill_manual(values=cols,name=NULL)+
  labs(x=paste0("PCoA 1 (",ve[1],"%)"),y=paste0("PCoA 2 (",ve[2],"%)"),title="Bray-Curtis (ASV)")+theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.position="right")
pD <- ggplot(disp %>% filter(distance=="Bray-Curtis (ASV)"),aes(Group,dist_to_centroid,fill=Group))+geom_jitter(width=0.12,height=0,shape=21,size=2,stroke=0.3)+
  stat_summary(fun=median,geom="crossbar",width=0.45,linewidth=0.3)+scale_fill_manual(values=cols)+labs(x=NULL,y="Distance to group centroid",title="Dispersion")
fig1 <- (dot("Observed","Observed ASVs (1,000 reads)")+ggtitle("Richness") | dot("Shannon","Shannon index")+ggtitle("Shannon")) / (pP | pD) + plot_annotation(tag_levels="a")
ggsave(file.path(out,"r2a_sasm_diversity_composition.png"),fig1,width=7.5,height=6,dpi=300)
## stacked bars: top 15 genera across the three groups, per well and group mean
pal <- c("#771155","#CC99BB","#114477","#4477AA","#117777","#44AAAA","#77CCCC","#117744","#44AA77","#88CCAA","#777711","#AAAA44","#DDDD77","#774411","#AA7744","#DDAA77","#771122","#AA4455","#DD7788","#41AB5D")
rank_lab <- function(tx,lvl) sapply(strsplit(tx,"; "),function(p){ p<-sub("^[a-z]__","",p); p<-c(p,rep("",7))[1:7]
  ok<-function(x) x!="" & !grepl("^uncultured|^unidentified|^metagenome",x,ignore.case=TRUE)
  if(ok(p[lvl])) return(p[lvl]); k<-max(c(0,which(sapply(p[1:(lvl-1)],ok)))); if(k==0) "Unassigned" else paste0("Unclassified ",p[k])})
rel <- sweep(X[,md$SampleID],2,colSums(X[,md$SampleID]),"/"); R <- rowsum(rel,rank_lab(tax,6))
gm <- sapply(levels(md$Group),function(g) rowMeans(R[,md$SampleID[md$Group==g],drop=FALSE]))
top <- names(sort(rowMeans(gm),decreasing=TRUE))[1:15]                     # rank by mean of group means (each medium weighted equally)
R2 <- rbind(R[top,],Other=colSums(R[!rownames(R)%in%top,,drop=FALSE]))
df <- as.data.frame(t(R2)) %>% rownames_to_column("SampleID") %>% pivot_longer(-SampleID,names_to="Taxon",values_to="RA") %>% left_join(md,by="SampleID") %>%
  mutate(Taxon=factor(Taxon,levels=c(top,"Other")))
ord <- df %>% filter(Taxon==top[1]) %>% arrange(Group,desc(RA)) %>% pull(SampleID); df$SampleID <- factor(df$SampleID,levels=ord)
fc <- setNames(c(pal[1:15],"#D9D9D9"),c(top,"Other"))
sb <- function(d,x,ttl) ggplot(d,aes({{x}},RA,fill=Taxon))+geom_col(width=0.9)+scale_fill_manual(values=fc,name="Top 15 genera")+
  scale_y_continuous(expand=c(0,0),labels=scales::percent)+labs(x=NULL,y="Relative abundance",title=ttl)+
  theme(legend.position="right",legend.key.size=unit(8,"pt"),legend.text=element_text(size=6.5,face="italic"))
pW <- sb(df,SampleID,"Individual wells")+facet_grid(.~Group,scales="free_x",space="free_x")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),strip.background=element_blank(),strip.text=element_text(size=8,face="bold"))
dm <- df %>% group_by(Group,Taxon) %>% summarise(RA=mean(RA),n=n(),.groups="drop")
pM <- sb(dm,Group,"Group means")+theme(legend.position="none",axis.text.x=element_text(angle=30,hjust=1))
ggsave(file.path(out,"r2a_sasm_stackedbar_genus.png"),(pM|pW)+plot_layout(widths=c(1,4)),width=10,height=4.5,dpi=300)
write.csv(dm,file.path(out,"r2a_sasm_stackedbar_genus_means.csv"),row.names=FALSE)
cat("Groups:",paste(names(table(md$Group)),table(md$Group),collapse="; "),"\n")
