# Influences of plant genotype, chemotype and environment on the leaf bacterial community

**A. Malacrinò, R. Jakobs, S. Xu, C. Müller**

## Abstract

# Disclaimer

This repository contains the main components used to process the raw data and to analyze it. Raw data is available at NCBI SRA under the BioProject number `PRJNA1080585`.

Our pipeline included:
* nf-core/ampliseq v2.1.7 [https://github.com/nf-core/ampliseq/](https://github.com/nf-core/ampliseq/)
* MAFFT [https://academic.oup.com/nar/article/30/14/3059/2904316](https://academic.oup.com/nar/article/30/14/3059/2904316)
* FastTree [https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490)
* R v4.3.2 [https://www.R-project.org/](https://www.R-project.org/)

# Data processing

Run ampliseq

```bash
nextflow run nf-core/ampliseq -r 2.7.1 -profile singularity \
--input $INDIR \
--FW_primer "AACMGGATTAGATACCCKG" \
--RV_primer "ACGTCATCCCCACCTTCC" \
--outdir $OUTDIR \
--extension "/*_{1,2}.fastq.gz" \
--trunclenf 0 \
--trunclenr 0 \
--skip_qiime \
--skip_barplot \
--skip_abundance_tables \
--skip_alpha_rarefaction \
--skip_diversity_indices \
--skip_ancom \
--max_cpus 16 \
--max_memory '128.GB'

mafft --thread $NTHREADS ASV_seqs.fasta > asv_aligned.fasta

FastTree -gtr -nt < asv_aligned.fasta > tree.tre
```

Build a phyloseq object

```r
library("phyloseq")
library("Biostrings")

creatPSobj <- function(mdata, asvtab, tax, tree, refseqs){
  metadata <- read.table(mdata, sep = '\t', header = T, row.names = 1)
  rownames(metadata) <- gsub('[-]', '.', rownames(metadata))
  data <- read.table(asvtab, sep = '\t', header = T, row.names = 1)
  tax <- read.table(tax, sep = '\t', header = T, row.names = 1)[,c(1:7)]
  ref_tree <- ape::read.tree(tree) 
  refseq <- read.table(refseqs, sep = '\t', header = T)[,c(1,10)]
  seqs.vec <- refseq$sequence 
  names(seqs.vec) <- refseq$ASV_ID
  refseq <- DNAStringSet(seqs.vec)
  ps <- phyloseq(sample_data(metadata),
                 otu_table(data, taxa_are_rows = T),
                 tax_table(as.matrix(tax)),
                 phy_tree(ref_tree),
                 refseq)
  return(ps)
}

ps <- creatPSobj(mdata = "metadata.txt",
                     asvtab = "ASV_table.tsv",
                     tax = "ASV_tax.silva_138.tsv",
                     tree = 'tree.tre',
                     refseqs = "ASV_tax.silva_138.tsv")

rm(list=setdiff(ls(), c("ps")))

save.image(file = 'data.rds')
```

## Data analysis

### Load libraries

```r
library("tidyverse")
library("phyloseq")
library("DESeq2")
library("ggrepel")
library("emmeans")
library("car")
library("lme4")
library("microbiome")
library("ggvenn")
library("decontam")
library("Wrench")
library("RVAideMemoire")
library("picante")
library("decontam")
library("reshape2")
library("ggvenn")
library("MOFA2")
library("RColorBrewer")
library("data.table")
library("psych")
```

### Load and clean microbiome data

```r
load(file = 'data/data.rds')

ps <- subset_taxa(ps, Class !="Chloroplast")
ps <- subset_taxa(ps, Order !="Chloroplast")
ps <- subset_taxa(ps, Family !="Mitochondria")

remove.cont <- function(ps){
  sample_data(ps)$is.neg <- sample_data(ps)$location == "neg_ctrl"
  contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold = 0.05)
  cont.remove <- subset(contamdf.prev, contaminant == "TRUE")
  cont.remove <- row.names(cont.remove)
  allTaxa <- taxa_names(ps)
  allTaxa <- allTaxa[!(allTaxa %in% cont.remove)]
  ps <- prune_taxa(allTaxa, ps)
  return(ps)
}

ps <- remove.cont(ps)
ps <- phyloseq::subset_samples(ps, location != "neg_ctrl")
ps <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE)

mdata <- sample_data(ps) %>% as.matrix() %>% as.data.frame()
mdata <- mdata[order(row.names(mdata)),]

wrench.norm.ps <- function(ps){
  count_tab <- as.matrix(data.frame(otu_table(ps)))
  group <- paste0(sample_data(ps)$location)
  W <- wrench(count_tab, condition=group)
  norm_factors <- W$nf
  norm_counts <- sweep(count_tab, 2, norm_factors, FUN = '/')
  norm_counts_trim <- data.frame(t(data.frame(norm_counts)))                                                  
  norm_counts_trim[] <- lapply(norm_counts_trim, function(x) DescTools::Winsorize(x, probs = c(0, 0.97), type = 1))
  norm_counts_trim <- data.frame(t(norm_counts_trim))
  norm_counts_trim[norm_counts_trim == 0] <- 1
  norm_counts_trim[norm_counts_trim < 0] <- 1
  norm_counts_trim[is.na(norm_counts_trim)] <- 1
  norm_counts_trim <- log2(norm_counts_trim)
  colnames(norm_counts_trim) <- gsub("\\.", "-", colnames(norm_counts_trim))
  ps_norm <- ps
  otu_table(ps_norm) <- otu_table(norm_counts_trim, taxa_are_rows =  TRUE)
  return(ps_norm)
}

ps_n <- wrench.norm.ps(ps)
```

### Load terpenoid data

```r
data.metabolome <- read.table("data/terpenoid.txt", header = T, sep = "\t", row.names = 1)
data.metabolome <- as.data.frame(t(data.metabolome))
data.metabolome <- data.metabolome[order(row.names(data.metabolome)),]
```

### Analyze terpenoid data

```r
test.dist <- vegdist(data.metabolome, method='bray', na.rm = T)
set.seed(100)
test.div <- adonis2(test.dist ~ chemotype * maternal_ID * location, data=mdata, permutations = 999)
test.div
```

### NMDS

```r
vare.mds <- metaMDS(data.metabolome) 
data.scores <- as.data.frame(scores(vare.mds)$sites)
data.scores <- cbind(mdata, data.scores)
data.scores.ab <- data.scores

nmds <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2, color = location)) +
  theme_bw(base_size = 20) +
  stat_ellipse(mapping = aes(color = location),
               alpha = 0.4,
               type = "norm",
               show.legend=F) +
  geom_point(mapping = aes(color = location, shape = location), size = 5) +
  theme(legend.title= element_blank(), 
        legend.background = element_rect(color = NA),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = c(0.99, 0.01), legend.justification = c(0.99, 0.01),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "mm")) +
  scale_color_manual(name = "Legend", values=c("#1b9e77", "#d95f02"), labels = c("field", "greenhouse"), breaks = c("field", "greenhouse")) +
  scale_shape_manual(name = "Legend", values=c(19,17), labels = c("field", "greenhouse"), breaks = c("field", "greenhouse")) +
  xlim(-2, 2) + ylim(-2, 2)
nmds
```
![image](https://github.com/amalacrino/malacrino_et_al_chemodiversity/assets/21124426/3c50d380-8990-48b4-beb3-46b52ebc3495)

```r
nmds <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2, color = chemotype)) +
  theme_bw(base_size = 20) +
  stat_ellipse(mapping = aes(color = chemotype),
               alpha = 0.4,
               type = "norm",
               show.legend=F) +
  geom_point(mapping = aes(color = chemotype, shape = chemotype), size = 5) +
  theme(legend.title= element_blank(), 
        legend.background = element_rect(color = NA),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = c(0.99, 0.01), legend.justification = c(0.99, 0.01),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "mm")) +
  scale_color_manual(name = "Legend", values=c("#FFC466", "#F9623E", "#B481FF"), labels = c("Keto", "ABThu", "Myrox"), breaks = c("Keto", "ABThu", "Myrox")) +
  scale_shape_manual(name = "Legend", values=c(17, 19, 15), labels = c("Keto", "ABThu", "Myrox"), breaks = c("Keto", "ABThu", "Myrox")) +
  xlim(-2, 2) + ylim(-2, 2)
nmds
```
![image](https://github.com/amalacrino/malacrino_et_al_chemodiversity/assets/21124426/f3b16682-9b76-4516-8023-2eb8cc8291d8)


```r
nmds <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2, color = maternal_ID)) +
  theme_bw(base_size = 20) +
  stat_ellipse(mapping = aes(color = maternal_ID),
               alpha = 0.4,
               type = "norm",
               show.legend=F) +
  geom_point(mapping = aes(color = maternal_ID, shape = maternal_ID), size = 5) +
  theme(legend.title= element_blank(), 
        legend.background = element_rect(color = NA),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = c(0.99, 0.01), legend.justification = c(0.99, 0.01),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "mm")) +
  scale_color_manual(name = "Legend", values=c("#e41a1c", "#377eb8", "#4daf4a"), labels = c("M09", "M18", "M26"), breaks = c("M09", "M18", "M26")) +
  scale_shape_manual(name = "Legend", values=c(19, 17, 15), labels = c("M09", "M18", "M26"), breaks = c("M09", "M18", "M26")) +
  xlim(-2, 2) + ylim(-2, 2)
ggsave(nmds, filename = "figures/chem_nmds_mID.pdf", dpi = 600,  width = 6, height = 5, units = "in")
nmds
```
![image](https://github.com/amalacrino/malacrino_et_al_chemodiversity/assets/21124426/c384d6c2-b555-4d52-a63d-b92d9e71e184)

### Barplot

```r
df <- merge(mdata, data.metabolome, by=0, all=TRUE)

mycolors <- c("#a6cee3","#1f78b4","#b2df8a","#F9623E","#fb9a99","#FFC466","#33a02c","#ff7f00","#e31a1c","#6a3d9a","#ffff99","#b15928","#8dd3c7","#bebada","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f","#8c510a","#c51b7d","#e08214","#B481FF")

dat <- melt(df)
dat <- dat %>% group_by(location, chemotype, maternal_ID, variable) %>% dplyr::summarize(cs = mean(value)) %>% mutate(cs = cs/sum(cs)) 
dat$chemotype = factor(dat$chemotype, levels=c("Keto", "ABThu", "Myrox"))

taxa_plot <- ggplot(dat, aes(x = as.factor(maternal_ID), y = cs, fill = variable)) +
  facet_grid(location~chemotype) +
  theme_bw(base_size = 14) +
  geom_bar(stat="identity") +
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill="transparent"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.title=element_blank(),
        legend.text.align = 0) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = mycolors,
                    labels = c("1,8-cineole",
                               expression(paste(alpha, "-pinene", sep = "")),
                               expression(paste(alpha, "-thujene", sep = "")),
                               expression(paste(alpha, "-thujone", sep = "")),
                               "artemisia alcohol",
                               "artemisia ketone",
                               "artemisyl acetate",
                               expression(paste(beta, "-pinene", sep = "")),
                               expression(paste(beta, "-thujone", sep = "")),
                               "borneol",
                               "camphene",
                               "camphor",
                               "chrysanthenyl acetate",
                               "eugenol",
                               "lavandulol",
                               "limonene",
                               expression(paste(italic(p), "-cymene", sep = "")),
                               "sabinene",
                               "sabinene hydrate",
                               "santolina triene",
                               "trans-sabinol",
                               "umbellulone",
                               "unknown monoterpenoid 1",
                               "unknown monoterpenoid 4",
                               "unknown monoterpenoid 5",
                               "unknown monoterpenoid 9",
                               "unknown terpineol 1",
                               "yomogi alcohol",
                               expression(paste("(", italic(Z), ")-myroxide", sep = ""))
                    )) +
  labs(y = "relative abundance", x="")
taxa_plot
```
![image](https://github.com/amalacrino/malacrino_et_al_chemodiversity/assets/21124426/8dd3e766-47f3-4f46-9614-552c4f5242ae)

### Analyze microbiome data

### PERMANOVA

```r
sampledf <- data.frame(sample_data(ps_n))
dist.mat <- phyloseq::distance(ps, method = "wunifrac")
perm <- how(nperm = 999)
set.seed(100)
pmv <- adonis2(dist.mat ~ maternal_ID * chemotype * location, data = sampledf, permutations = perm)
pmv

pairwise.perm.manova(dist.mat, paste0(sampledf$maternal_ID), nperm = 999, progress = TRUE, p.method = "fdr", F = T, R2 = T)
```

### NMDS

```r
cap_ord <- ordinate(physeq = ps, method = "NMDS", distance = dist.mat, formula = ~ 1)
cap_plot <- plot_ordination(physeq = ps, ordination = cap_ord, axes = c(1,2)) +
  theme_bw(base_size = 20) +
  stat_ellipse(mapping = aes(color = location),
               alpha = 0.4,
               type = "norm",
               show.legend=F) +
  geom_point(mapping = aes(color = location, shape = location), size = 5) +
  theme(legend.title= element_blank(), 
        legend.background = element_rect(color = NA),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = c(0.99, 0.99), legend.justification = c(0.99, 0.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "mm")) +
  scale_color_manual(name = "Legend", values=c("#1b9e77", "#d95f02"), labels = c("field", "greenhouse"), breaks = c("field", "greenhouse")) +
  scale_shape_manual(name = "Legend", values=c(19,17), labels = c("field", "greenhouse"), breaks = c("field", "greenhouse")) +
  xlim(-0.4, 0.4) + ylim(-0.3, 0.3)
cap_plot
```
![image](https://github.com/amalacrino/malacrino_et_al_chemodiversity/assets/21124426/d37e7edd-a302-403d-9cd3-4601560564f6)

```r
cap_ord <- ordinate(physeq = ps, method = "NMDS", distance = dist.mat, formula = ~ 1)
cap_plot <- plot_ordination(physeq = ps, ordination = cap_ord, axes = c(1,2)) +
  theme_bw(base_size = 20) +
  stat_ellipse(mapping = aes(color = chemotype),
               alpha = 0.4,
               type = "norm",
               show.legend=F) +
  geom_point(mapping = aes(color = chemotype, shape = chemotype), size = 5) +
  theme(legend.title= element_blank(), 
        legend.background = element_rect(color = NA),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = c(0.99, 0.99), legend.justification = c(0.99, 0.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "mm")) +
  scale_color_manual(name = "Legend", values=c("#FFC466", "#F9623E", "#B481FF"), labels = c("Keto", "ABThu", "Myrox"), breaks = c("Keto", "ABThu", "Myrox")) +
  scale_shape_manual(name = "Legend", values=c(17, 19, 15), labels = c("Keto", "ABThu", "Myrox"), breaks = c("Keto", "ABThu", "Myrox")) +
  xlim(-0.4, 0.4) + ylim(-0.3, 0.3)
cap_plot
```
![image](https://github.com/amalacrino/malacrino_et_al_chemodiversity/assets/21124426/df0a4f8d-df1b-42fb-a1ce-00a52d32923b)

```r
cap_ord <- ordinate(physeq = ps, method = "NMDS", distance = dist.mat, formula = ~ 1)
cap_plot <- plot_ordination(physeq = ps, ordination = cap_ord, axes = c(1,2)) +
  theme_bw(base_size = 20) +
  stat_ellipse(mapping = aes(color = maternal_ID),
               alpha = 0.4,
               type = "norm",
               show.legend=F) +
  geom_point(mapping = aes(color = maternal_ID, shape = maternal_ID), size = 5) +
  theme(legend.title= element_blank(), 
        legend.background = element_rect(color = NA),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position = c(0.99, 0.99), legend.justification = c(0.99, 0.99),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "mm")) +
  scale_color_manual(name = "Legend", values=c("#e41a1c", "#377eb8", "#4daf4a"), labels = c("M09", "M18", "M26"), breaks = c("M09", "M18", "M26")) +
  scale_shape_manual(name = "Legend", values=c(19, 17, 15), labels = c("M09", "M18", "M26"), breaks = c("M09", "M18", "M26")) +
  xlim(-0.4, 0.4) + ylim(-0.3, 0.3)
cap_plot
```
![image](https://github.com/amalacrino/malacrino_et_al_chemodiversity/assets/21124426/d733965c-6e20-4e1d-a336-e74e10af56da)

### Taxaplot

```r
glom <- microbiome::aggregate_taxa(ps_n, "Genus")
glom <- microbiome::transform(glom, "compositional")
dat <- psmelt(glom)
filt.gen <- dat %>% group_by(Genus) %>% summarize(mean = mean(Abundance)) %>% filter(mean <= 0.01)
dat <- subset(dat, !(OTU %in% filt.gen$Genus))
dat <- dat %>% group_by(location, maternal_ID, Genus) %>% summarize(cs = mean(Abundance)) %>% mutate(cs = cs/sum(cs)) %>% mutate(Genus = replace(Genus, Genus == "", "unidentified taxa"))
nb.cols <- length(unique(dat$Genus))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

taxa_plot <- ggplot(dat, aes(x = as.factor(maternal_ID), y = cs, fill = Genus)) +
  facet_wrap(vars(location)) +
  theme_bw(base_size = 14) +
  geom_bar(stat="identity") +
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill="transparent"),
        legend.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank()) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = mycolors) +
  labs(y = "relative abundance", x="")
taxa_plot
```
![image](https://github.com/amalacrino/malacrino_et_al_chemodiversity/assets/21124426/47347c98-088f-41ca-a5df-a4bb68c785b9)

### Venn diagrams

```r
ps.1 <- subset_samples(ps, location == "field")
ps.1 <- filter_taxa(ps.1, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps.1 <- rownames(as.data.frame(otu_table(ps.1)))
ps.2 <- subset_samples(ps, location == "greenhouse")
ps.2 <- filter_taxa(ps.2, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps.2 <- rownames(as.data.frame(otu_table(ps.2)))

a <- list(Field = ps.1,
          Greenhouse = ps.2)
ggvenn(a) 

ps.1 <- subset_samples(ps, chemotype == "Myrox")
ps.1 <- filter_taxa(ps.1, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps.1 <- rownames(as.data.frame(otu_table(ps.1)))
ps.2 <- subset_samples(ps, chemotype == "Keto")
ps.2 <- filter_taxa(ps.2, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps.2 <- rownames(as.data.frame(otu_table(ps.2)))
ps.3 <- subset_samples(ps, chemotype == "ABThu")
ps.3 <- filter_taxa(ps.3, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps.3 <- rownames(as.data.frame(otu_table(ps.3)))

a <- list(Myrox = ps.1,
          Keto = ps.2,
          ABThu = ps.3)
ggvenn(a) 

ps.1 <- subset_samples(ps, maternal_ID == "M18")
ps.1 <- filter_taxa(ps.1, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps.1 <- rownames(as.data.frame(otu_table(ps.1)))
ps.2 <- subset_samples(ps, maternal_ID == "M26")
ps.2 <- filter_taxa(ps.2, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps.2 <- rownames(as.data.frame(otu_table(ps.2)))
ps.3 <- subset_samples(ps, maternal_ID == "M09")
ps.3 <- filter_taxa(ps.3, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps.3 <- rownames(as.data.frame(otu_table(ps.3)))

a <- list(`9` = ps.3,
          `18` = ps.1,
          `26` = ps.2)
ggvenn(a) 
```

### Volcano plots

```r
cal.diff.taxa <- function(object){
  diagdds <- phyloseq_to_deseq2(object, ~ 1)
  ts <- counts(diagdds)
  geoMeans = apply(ts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  diagdds = estimateSizeFactors(diagdds, geoMeans=geoMeans)
  diagdds = estimateDispersions(diagdds)
  diagdds$group <- factor(paste0(diagdds$location))
  design(diagdds) <- ~ group
  dds <-DESeq(diagdds, betaPrior=FALSE, parallel = T)
  c1 <- results(dds, contrast=c("group", "field", "greenhouse"), parallel = T)
  c1 <- as.data.frame(c1)
  c1 <- setDT(c1, keep.rownames = TRUE)[]
  c1 <- c1[,c("rn", "log2FoldChange", "padj")]
  tax.table <- as.data.frame(tax_table(object))
  tax.table <- setDT(tax.table, keep.rownames = TRUE)[]
  tx <- merge(c1, tax.table, by = "rn")
  tx$diffexpressed <- "no changes"
  tx$diffexpressed[tx$log2FoldChange > 0 & tx$padj < 0.05] <- "field"
  tx$diffexpressed[tx$log2FoldChange < 0 & tx$padj < 0.05] <- "greenhouse"
  return(tx)
}

df.diff <- cal.diff.taxa(ps)

plot <- ggplot(data=df.diff) +
  theme_bw(base_size = 14) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = diffexpressed), size = 3) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(diffexpressed != "no changes", Genus,"")), position = "dodge", size = 5, max.overlaps = Inf, force = 20, box.padding = 0.2, min.segment.length = 0, seed = 42) +
  geom_vline(xintercept=0, col="black", linetype = "longdash") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = "longdash") +
  theme(legend.justification=c(0.99,0.99), legend.position=c(0.99,0.99),
        panel.grid = element_blank(),
        legend.title=element_blank(), 
        legend.background = element_rect(color = NA),
        legend.key = element_rect(color = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  xlab(expression(paste(Log[2], " Fold Changes"))) +
  ylab(expression(paste(-Log[10], " P"))) +
  scale_color_manual(name = "Legend", values=c("#1b9e77", "#d95f02", "#000000"), labels = c("field", "greenhouse", "no changes"), breaks = c("field", "greenhouse", "no changes"))
plot
```
![image](https://github.com/amalacrino/malacrino_et_al_chemodiversity/assets/21124426/464ea58b-a339-402b-b48d-832db0062c88)


```r
cal.diff.taxa2 <- function(object, g1, g2){
  diagdds <- phyloseq_to_deseq2(object, ~ 1)
  ts <- counts(diagdds)
  geoMeans = apply(ts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  diagdds = estimateSizeFactors(diagdds, geoMeans=geoMeans)
  diagdds = estimateDispersions(diagdds)
  diagdds$group <- factor(paste0(diagdds$maternal_ID))
  design(diagdds) <- ~ group
  dds <-DESeq(diagdds, betaPrior=FALSE, parallel = T)
  c1 <- results(dds, contrast=c("group", g1, g2), parallel = T)
  c1 <- as.data.frame(c1)
  c1 <- setDT(c1, keep.rownames = TRUE)[]
  c1 <- c1[,c("rn", "log2FoldChange", "padj")]
  tax.table <- as.data.frame(tax_table(object))
  tax.table <- setDT(tax.table, keep.rownames = TRUE)[]
  tx <- merge(c1, tax.table, by = "rn")
  tx$diffexpressed <- "no changes"
  tx$diffexpressed[tx$log2FoldChange > 0 & tx$padj < 0.05] <- paste0(g1)
  tx$diffexpressed[tx$log2FoldChange < 0 & tx$padj < 0.05] <- paste0(g2)
  return(tx)
}

df.diff <- cal.diff.taxa2(ps, "M18", "M26")

plot <- ggplot(data=df.diff) +
  theme_bw(base_size = 14) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = diffexpressed), size = 3) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(diffexpressed != "no changes", Genus,"")), position = "dodge", size = 5, max.overlaps = Inf, force = 10, box.padding = 0.2, min.segment.length = 0, seed = 42) +
  geom_vline(xintercept=0, col="black", linetype = "longdash") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = "longdash") +
  theme(legend.justification=c(0.99,0.99), legend.position=c(0.99,0.99),
        panel.grid = element_blank(),
        legend.title=element_blank(), 
        legend.background = element_rect(color = NA),
        legend.key = element_rect(color = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  xlab(expression(paste(Log[2], " Fold Changes"))) +
  ylab(expression(paste(-Log[10], " P"))) +
  scale_color_manual(name = "Legend", values=c("#377eb8", "#4daf4a", "#000000"), labels = c("M18", "M26", "no changes"), breaks = c("M18", "M26", "no changes")) +
  ylim(0, 25) + xlim(-25,25)
plot
```
![image](https://github.com/amalacrino/malacrino_et_al_chemodiversity/assets/21124426/892be366-6ef8-4b78-9c17-a68fcc4f570b)

```r
df.diff <- cal.diff.taxa2(ps, "M09", "M18")

plot <- ggplot(data=df.diff) +
  theme_bw(base_size = 14) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = diffexpressed), size = 3) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(diffexpressed != "no changes", Genus,"")), position = "dodge", size = 5, max.overlaps = Inf, force = 10, box.padding = 0.2, min.segment.length = 0, seed = 42) +
  geom_vline(xintercept=0, col="black", linetype = "longdash") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = "longdash") +
  theme(legend.justification=c(0.99,0.99), legend.position=c(0.99,0.99),
        panel.grid = element_blank(),
        legend.title=element_blank(), 
        legend.background = element_rect(color = NA),
        legend.key = element_rect(color = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  xlab(expression(paste(Log[2], " Fold Changes"))) +
  ylab(expression(paste(-Log[10], " P"))) +
  scale_color_manual(name = "Legend", values=c("#e41a1c", "#377eb8", "#000000"), labels = c("M09", "M18", "no changes"), breaks = c("M09", "M18", "no changes")) +
  ylim(0, 25) + xlim(-25,25)
plot
```
![image](https://github.com/amalacrino/malacrino_et_al_chemodiversity/assets/21124426/82a0e08d-ea87-437e-9d1a-7ffc4ed6d934)

```r
df.diff <- cal.diff.taxa2(ps, "M09", "M26")

plot <- ggplot(data=df.diff) +
  theme_bw(base_size = 14) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = diffexpressed), size = 3) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(diffexpressed != "no changes", Genus,"")), position = "dodge", size = 5, max.overlaps = Inf, force = 10, box.padding = 0.2, min.segment.length = 0, seed = 42) +
  geom_vline(xintercept=0, col="black", linetype = "longdash") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = "longdash") +
  theme(legend.justification=c(0.99,0.99), legend.position=c(0.99,0.99),
        panel.grid = element_blank(),
        legend.title=element_blank(), 
        legend.background = element_rect(color = NA),
        legend.key = element_rect(color = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  xlab(expression(paste(Log[2], " Fold Changes"))) +
  ylab(expression(paste(-Log[10], " P"))) +
  scale_color_manual(name = "Legend", values=c("#e41a1c", "#4daf4a", "#000000"), labels = c("M09", "M26", "no changes"), breaks = c("M09", "M26", "no changes")) +
  ylim(0, 25) + xlim(-25,25)
ggsave(plot, filename = "figures/volcano_26vs9.pdf", dpi = 600,  width = 8, height = 6, units = "in")
plot
```
![image](https://github.com/amalacrino/malacrino_et_al_chemodiversity/assets/21124426/1d983741-7702-453e-938e-a6baf77be51a)

```r
cal.diff.taxa3 <- function(object, g1, g2){
  diagdds <- phyloseq_to_deseq2(object, ~ 1)
  ts <- counts(diagdds)
  geoMeans = apply(ts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  diagdds = estimateSizeFactors(diagdds, geoMeans=geoMeans)
  diagdds = estimateDispersions(diagdds)
  diagdds$group <- factor(paste0(diagdds$chemotype))
  design(diagdds) <- ~ group
  dds <-DESeq(diagdds, betaPrior=FALSE, parallel = T)
  c1 <- results(dds, contrast=c("group", g1, g2), parallel = T)
  c1 <- as.data.frame(c1)
  c1 <- setDT(c1, keep.rownames = TRUE)[]
  c1 <- c1[,c("rn", "log2FoldChange", "padj")]
  tax.table <- as.data.frame(tax_table(object))
  tax.table <- setDT(tax.table, keep.rownames = TRUE)[]
  tx <- merge(c1, tax.table, by = "rn")
  tx$diffexpressed <- "no changes"
  tx$diffexpressed[tx$log2FoldChange > 0 & tx$padj < 0.05] <- paste0(g1)
  tx$diffexpressed[tx$log2FoldChange < 0 & tx$padj < 0.05] <- paste0(g2)
  return(tx)
}

cal.diff.taxa3(ps, "Keto", "ABThu") %>% filter(padj < 0.05)
cal.diff.taxa3(ps, "Keto", "Myrox") %>% filter(padj < 0.05)
cal.diff.taxa3(ps, "ABThu", "Myrox") %>% filter(padj < 0.05)
```

### Mantel test

```r
micro.dist <- phyloseq::distance(ps, method = "wunifrac")
terp.dist <- vegdist(data.metabolome, method='bray', na.rm = T)
mantel(xdis = dist.mat, ydis = test.dist, method = "pearson", permutations = 999, na.rm = TRUE) 
```

### Correlate microbiome and terpenoid Shannon's diversity

```r
div <- microbiome::alpha(ps, index = "diversity_shannon")
div <- cbind(sample_data(ps), div)
H <- as.data.frame(diversity(data.metabolome, "shannon"))
colnames(H) <- c("shannon")
div$chem <- H$shannon
cor.test(div$diversity_shannon, div$chem)
```

MOFA2

```r
set.seed(100)
data.metabolome.mofa <- t(as.matrix(data.metabolome))
glom <- microbiome::transform(ps_n, "compositional")
otu.table.mofa <- as.matrix(otu_table(glom))
taxa.df <- tax_table(glom) %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var = "rn")
setdiff(colnames(data.metabolome.mofa), colnames(otu.table.mofa))
dataMOFA <- list(metabolome = data.metabolome.mofa, microbiota = otu.table.mofa)

MOFAobject <- create_mofa(dataMOFA)

data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10
train_opts <- get_default_training_options(MOFAobject)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

model <- run_mofa(MOFAobject, use_basilisk = TRUE)
model

factors <- get_factors(model, factors = "all")
weights <- get_weights(model, views = "all", factors = "all")
factors.l <- get_factors(model, as.data.frame = T)
weights.l <- get_weights(model, as.data.frame = T)
data.l <- get_data(model, as.data.frame = T)

print(plot_variance_explained(model, x="view", y="factor"))
plot_variance_explained(model, plot_total = T)[[2]]

model@cache$variance_explained$r2_per_factor[[1]]

corr.mat.fun <- function(x, perc1, perc2){
  set.seed(100)
  microb.w <-  weights.l %>%
    filter(factor == x) %>%
    filter(view == "microbiota") %>%
    slice_max(abs(value), prop = perc1)
  microb.w <-  merge(microb.w, taxa.df, by.x = "feature", by.y = "rn")
  metabol.w <-  weights.l %>%
    filter(factor == x) %>%
    filter(view == "metabolome") 
  mat.otu <- as.data.frame(t(otu.table.mofa))
  mat.met <- as.data.frame(t(data.metabolome.mofa))
  mat.otu <- mat.otu[which(names(mat.otu) %in% microb.w$feature)]
  mat.met <- mat.met[which(names(mat.met) %in% metabol.w$feature)]
  mat.otu <- mat.otu[which(row.names(mat.otu) %in% row.names(mdata)),]
  mat.met <- mat.met[which(row.names(mat.met) %in% row.names(mdata)),]  
  for(col in names(mat.otu)) {mat.otu[col] = mat.otu[col] / sum(mat.otu[col])}
  for(col in names(mat.met)) {mat.met[col] = mat.met[col] / sum(mat.met[col])}
  p.mat <- corr.test(mat.otu, mat.met)
  return(p.mat)
}

flattenCorrMatrix <- function(cormat, pmat) {
  set.seed(100)
  a <- melt(as.matrix(cormat))
  b <- melt(as.matrix(pmat))
  a$p <- b$value
  colnames(a) <- c("OTU", "Compound", "r", "p")
  flat.cor.mat <- a
  flat.cor.mat <- merge(flat.cor.mat, taxa.df, by.x = "OTU", by.y = "rn")
  return(flat.cor.mat)
}

p.mat <- corr.mat.fun("Factor1",  perc1 = 0.1, perc2 = 1)

flat.cor.mat <- flattenCorrMatrix(p.mat$r, p.mat$p) %>% filter(p < 0.05) %>% na.omit()

plot <- ggplot(flat.cor.mat, aes(OTU,Compound, fill=r)) + 
  geom_tile(height = 0.5, color = "black") + 
  scale_fill_viridis_c(direction=-1, option = "H") + 
  theme_bw(14) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_y_discrete(limits=rev,
                   labels=c("1_8_cineole"="1,8-cineole",
                            "a_pinene"= expression(paste(alpha, "-pinene", sep = "")),
                            "a_thujene"= expression(paste(alpha, "-thujene", sep = "")),
                            "a_thujone" = expression(paste(alpha, "-thujone", sep = "")),
                            "artemisia_alcohol" = "artemisia alcohol",
                            "artemisia_ketone" = "artemisia ketone",
                            "b_artemisia_acetate" = "artemisyl acetate",
                            "b_Pinene" = expression(paste(beta, "-pinene", sep = "")),
                            "b_thujone" = expression(paste(beta, "-thujone", sep = "")),
                            "borneol" = "borneol",
                            "camphene" = "camphene",
                            "camphor" = "camphor",
                            "chrysanthenyl_acetate" = "chrysanthenyl acetate",
                            "eugenol" = "eugenol",
                            "lavandulol" = "lavandulol",
                            "limonene" = "limonene",
                            "p-cymene" = expression(paste(italic(p), "-cymene", sep = "")),
                            "sabinene" = "sabinene",
                            "sabinene_hydrate" = "sabinene hydrate",
                            "santolina_triene" = "santolina triene",
                            "trans_Sabinol" = "trans-sabinol",
                            "umbellulone" = "umbellulone",
                            "unknown_monoterpene_1" = "unknown monoterpenoid 1",
                            "unknown_monoterpene_4" = "unknown monoterpenoid 4",
                            "unknown_monoterpene_5" = "unknown monoterpenoid 5",
                            "unknown_monoterpene_9" = "unknown monoterpenoid 9",
                            "unknown_terpineol_1" = "unknown terpineol 1",
                            "Yomogi_alcohol" = "yomogi alcohol",
                            "Z_myroxide" =  expression(paste("(", italic(Z), ")-myroxide", sep = "")))) 
plot
```
![image](https://github.com/amalacrino/malacrino_et_al_chemodiversity/assets/21124426/b61c3c68-359c-4f5b-b5c4-3048f1e64ac2)

### Output table

```r
write.table(flat.cor.mat, "tabs1.csv", sep = ",", row.names = F)
```
