# Taxonomic and Functional Analysis of 72 metagenomes from Andes Ag

#### 1. Setup #####
# Libraries
library(plyr)
library(tidyverse)
library(readxl)
library(writexl)
library(mctoolsr)
library(vegan)
library(microseq)
library(FSA)
library(car)
library(RColorBrewer)
library(scales)
library(ggrepel)
library(cowplot)
library(ggh4x)
library(pheatmap)
library(ANCOMBC)
library(phyloseq)
library(DESeq2)
library(KEGGREST)
library(corrplot)
library(VennDiagram)
library(ggvenn)

# Functions
`%notin%` <- Negate(`%in%`)
source("~/Documents/GitHub/SunflowerGxE/code/cliffplot_taxa_bars.R")
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]

# Effect size (from Jack Darcy, Science Advances 2018)
source("~/Documents/GitHub/Sunflower_AEM/code/effectSize.R")

# Working directory
setwd("~/Desktop/Fierer/AndesAg/")

d1 <- read.csv("SampleSheet_24471Fie_N25015_L003.csv") %>%
  separate(Sample, remove = F, sep = "_", into = c("Project", "sampleID")) %>%
  separate(sampleID, remove = F, sep = "-", into = c("Group", "Replicate", "Depth")) %>%
  mutate(sampleID = gsub("-", ".", sampleID)) %>%
  mutate(Sample = sampleID) %>%
  rename("PfClusterCount" = "Reads") %>%
  mutate(Treatment = ifelse(Group == "A", "+MP1",
                            ifelse(Group == "B", "+Basalt",
                                   ifelse(Group == "C", "+MP1+Basalt",
                                          ifelse(Group == "D", "Control", "Blank"))))) %>%
  mutate(MP1 = ifelse(Group == "A" | Group == "C", "MP1",
                      ifelse(Group == "B" | Group == "D", "Control", "Blank"))) %>%
  mutate(Basalt = ifelse(Group == "B" | Group == "C", "Basalt", 
                         ifelse(Group == "A" | Group == "D", "Control", "Blank"))) %>%
  filter(Treatment != "Blank") %>% # Go ahead and remove, few seqs
  dplyr::select(Sample, everything())
#write.table(d1, "metadata.txt", sep = "\t", row.names = F)
d1$sampleID
table(d1$Group)
table(d1$Replicate)
table(d1$Depth)
table(d1$Treatment)
table(d1$MP1)
table(d1$Basalt)

# Barcodes (thought were needed but didn't end up using)
bcf <- d1 %>%
  dplyr::select(Barcode1) %>%
  mutate(BC = "Barcode") %>%
  mutate(Num = seq(1:nrow(.))) %>%
  mutate(Header = paste0(BC, Num)) %>%
  dplyr::rename(Sequence = Barcode1) %>%
  dplyr::select(Header, Sequence)
bcr <- d1 %>%
  dplyr::select(Barcode2) %>%
  mutate(BC = "Barcode") %>%
  mutate(Num = seq(1:nrow(.))) %>%
  mutate(Header = paste0(BC, Num)) %>%
  dplyr::rename(Sequence = Barcode2) %>%
  dplyr::select(Header, Sequence) %>%
  add_case(Header = "Barcode76", Sequence = "GGGGGGGGGG")
#writeFasta(bcf, "barcodesF.fasta")
#writeFasta(bcr, "barcodesR.fasta")

# Check sequencing depth
ggplot(d1, aes(Treatment, Reads)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.25, pch = 16, size = 2) +
  theme_bw()

ggplot(d1, aes(Depth, Reads)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.25, pch = 16, size = 2) +
  theme_bw()

# Check sequencing depth of concatenated lanes, trimmed and filtered
# for n = 70 samples analyzed (dropped 2 samples and 3 blanks)
rc <- read.delim("read_counts.txt", header = F) %>%
  separate(V1, into = c("sampleID", "ReadsAll", "Junk"), sep = " ") %>%
  mutate(sampleID = gsub("_R1_combined.fastq.gz:", "", sampleID)) %>%
  mutate(sampleID = gsub("B1/", "", sampleID)) %>%
  mutate(sampleID = gsub("B2/", "", sampleID)) %>%
  mutate(sampleID = gsub("B3/", "", sampleID)) %>%
  mutate(sampleID = gsub("B4/", "", sampleID)) %>%
  mutate(sampleID = gsub("-", ".", sampleID)) %>%
  mutate(ReadsAll = as.integer(ReadsAll)) %>%
  dplyr::select(-Junk) %>%
  left_join(., d1, by = c("sampleID" = "Sample")) %>%
  mutate(SampleID = gsub("\\.", "_", sampleID))
mean(rc$ReadsAll)
se(rc$ReadsAll)
range(rc$ReadsAll)

ggplot(rc, aes(Treatment, ReadsAll)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.25, pch = 16, size = 2) +
  theme_bw()

ggplot(rc, aes(Depth, ReadsAll)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.25, pch = 16, size = 2) +
  theme_bw()
  


#### 2. Taxonomy ####
#### _mTAGs ####
# Import and make prokaryote and eukaryote mctroolsr datasets

# Prok
mt <- read.delim("merged_profileAll.otu.tsv") %>%
  dplyr::rename(taxonomy = X.taxpath) %>%
  filter(., !grepl("Eukaryota", taxonomy)) %>%
  mutate(taxonomy = gsub("root__Root;", "", taxonomy)) %>%
  mutate(taxonomy = gsub("domain__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("phylum__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("class__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("order__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("family__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("genus__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("otu__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("silva_138_complink_cons_", "", taxonomy)) %>%
  mutate(taxonomy = gsub("unknown", "NA", taxonomy)) %>%
  mutate(taxonomy = gsub("otu_", "NA;otu_", taxonomy)) %>%
  filter(taxonomy != "Unassigned") %>%
  filter(taxonomy != "Unaligned") %>%
  separate(taxonomy, into = c("a","s","d","f","g","h","j","otu"), remove = F, sep = ";") %>%
  dplyr::select(-c(a,s,d,f,g,h,j)) %>%
  dplyr::select(otu, everything()) %>%
  dplyr::select(-taxonomy, taxonomy)
n <- data.frame(name = names(mt)) %>%
  mutate(name = gsub(".bins", "", name))
mt <- mt %>%
  set_names(n$name)
out_fp <- "seqtab_wTax_mctoolsr_mtagProk_all.txt"
names(mt)[1] = "#OTU_ID"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(mt,
                             out_fp,
                             sep = "\t",
                             row.names = FALSE,
                             append = TRUE))
sum(n$name %in% d1$sampleID)

# Import mctoolsr
tax_table_fp <- "seqtab_wTax_mctoolsr_mtagProk_all.txt"
map_fp <- "metadata.txt"
input = load_taxa_table(tax_table_fp, map_fp) # 70 samples loaded

# Filter Chloroplast, Mitochondria, Domain NA
input_filt <- filter_taxa_from_input(input,
                                     taxa_to_remove = "Chloroplast") # 48 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "Mitochondria") # 41 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "NA",
                                     at_spec_level = 1) # 0 removed

# Remove singletons and doubletons
singdoub <- data.frame("count" = rowSums(input_filt$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.))

input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_IDs_to_remove = singdoub$ASV) # 9078 removed

# Removed C.R1.20 and B.R1.30, very few reads

# Rarefaction
rarecurve(t(input_filt$data_loaded), step = 500, label = F)
sort(colSums(input_filt$data_loaded))
mean(colSums(input_filt$data_loaded)) # 13603.51
se(colSums(input_filt$data_loaded)) # 176.8941
set.seed(530)
input_filt_rare <- single_rarefy(input_filt, 9299) # n = 70 still
sort(colSums(input_filt_rare$data_loaded))

# Add rarefied richness and Shannon
input_filt_rare$map_loaded$rich <- specnumber(input_filt_rare$data_loaded, MARGIN = 2)
input_filt_rare$map_loaded$shannon <- diversity(input_filt_rare$data_loaded, 
                                                index = "shannon", MARGIN = 2)

# Save
#saveRDS(input_filt_rare, "input_filt_rare_mtagsProk_all.rds")

# Euk
mt <- read.delim("merged_profileAll.otu.tsv") %>%
  rename(taxonomy = X.taxpath) %>%
  filter(., grepl("Eukaryota", taxonomy)) %>%
  mutate(taxonomy = gsub("root__Root;", "", taxonomy)) %>%
  mutate(taxonomy = gsub("domain__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("phylum__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("class__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("order__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("family__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("genus__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("otu__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("silva_138_complink_cons_", "", taxonomy)) %>%
  mutate(taxonomy = gsub("unknown", "NA", taxonomy)) %>%
  mutate(taxonomy = gsub("otu_", "NA;otu_", taxonomy)) %>%
  filter(taxonomy != "Unassigned") %>%
  filter(taxonomy != "Unaligned") %>%
  separate(taxonomy, into = c("a","s","d","f","g","h","j","otu"), remove = F, sep = ";") %>%
  dplyr::select(-c(a,s,d,f,g,h,j)) %>%
  dplyr::select(otu, everything()) %>%
  dplyr::select(-taxonomy, taxonomy)
n <- data.frame(name = names(mt)) %>%
  mutate(name = gsub(".bins", "", name))
mt <- mt %>%
  set_names(n$name)
out_fp <- "seqtab_wTax_mctoolsr_mtagEuk_all.txt"
names(mt)[1] = "#OTU_ID"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(mt,
                             out_fp,
                             sep = "\t",
                             row.names = FALSE,
                             append = TRUE))
sum(n$name %in% d1$sampleID)

# Import mctoolsr
tax_table_fp <- "seqtab_wTax_mctoolsr_mtagEuk_all.txt"
map_fp <- "metadata.txt"
input = load_taxa_table(tax_table_fp, map_fp) # 70 samples loaded

# Filter Prokaryotes and Domain NA
input_filt <- input
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = c("Archaea", "Bacteria")) # 0 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "NA",
                                     at_spec_level = 1) # 0 removed

# Remove singletons and doubletons
singdoub <- data.frame("count" = rowSums(input_filt$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.))

input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_IDs_to_remove = singdoub$ASV) # 795 removed

# Removed C.R1.20 and B.R1.30, very few reads

# Rarefaction
rarecurve(t(input_filt$data_loaded), step = 500, label = F)
sort(colSums(input_filt$data_loaded))
mean(colSums(input_filt$data_loaded)) # 267.5571
se(colSums(input_filt$data_loaded)) # 7.0515

# Don't rarefy. Very few sequences.

# Add  richness and Shannon
input_filt$map_loaded$rich <- specnumber(input_filt$data_loaded, MARGIN = 2)
input_filt$map_loaded$shannon <- diversity(input_filt$data_loaded, 
                                                index = "shannon", MARGIN = 2)

# Save
#saveRDS(input_filt, "input_filt_mtagsEuk_all.rds")



#### __Prok ####
input_filt_rare <- readRDS("input_filt_rare_mtagsProk_all.rds")
input_filt_rare$map_loaded <- input_filt_rare$map_loaded %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Control",
                                       "+Basalt",
                                       "+MP1",
                                       "+MP1+Basalt"))) %>%
  mutate(Treatment = case_match(Treatment,
                                "Control" ~ "UTC",
                                "+Basalt" ~ "UTC+B",
                                "+MP1" ~ "MP1",
                                "+MP1+Basalt" ~ "MP1+B")) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("UTC", "MP1", "UTC+B", "MP1+B"))) %>%
  mutate(Depth = as.factor(Depth)) %>%
  mutate(TrtDep = factor(paste(Treatment, Depth, sep = "_"),
                          levels = c("Control_10", "Control_20", "Control_30",
                                     "+Basalt_10", "+Basalt_20", "+Basalt_30",
                                     "+MP1_10", "+MP1_20", "+MP1_30",
                                     "+MP1+Basalt_10", "+MP1+Basalt_20",
                                     "+MP1+Basalt_30")))

# Got soil data later, import and merge that
# Note: Alk data is just 1 per column so repeated values
soil <- read.csv("andes_cs5_compiled_soil_data.csv") %>%
  mutate(Comb = paste(Treatment, Condition, sep = "_")) %>%
  mutate(Group = ifelse(Comb == "MP1_Basalt", "C",
                        ifelse(Comb == "MP1_NB", "A",
                               ifelse(Comb == "UTC_Basalt", "B", "D")))) %>%
  mutate(Comb = gsub("MP1_Basalt", "+MP1+Basalt", Comb)) %>%
  mutate(Comb = gsub("MP1_NB", "+MP1", Comb)) %>%
  mutate(Comb = gsub("UTC_Basalt", "+Basalt", Comb)) %>%
  mutate(Comb = gsub("UTC_NB", "Control", Comb)) %>%
  mutate(Treatment = factor(Comb,
                            levels = c("Control",
                                       "+Basalt",
                                       "+MP1",
                                       "+MP1+Basalt"))) %>%
  mutate(Treatment = case_match(Treatment,
                                "Control" ~ "UTC",
                                "+Basalt" ~ "UTC+B",
                                "+MP1" ~ "MP1",
                                "+MP1+Basalt" ~ "MP1+B")) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("UTC", "MP1", "UTC+B", "MP1+B"))) %>%
  mutate(sampleID = paste(Group, Column, Increment, sep = ".")) %>%
  dplyr::select(-Treatment, -Condition, -Increment, -Column, -Comb, -Group) %>%
  dplyr::select(sampleID, everything())

dim(input_filt_rare$map_loaded)
rownames(input_filt_rare$map_loaded)
input_filt_rare$map_loaded <- input_filt_rare$map_loaded %>%
  left_join(., soil, by = "sampleID")
dim(input_filt_rare$map_loaded)
rownames(input_filt_rare$map_loaded) <- input_filt_rare$map_loaded$sampleID

# Check soil correlations
soil <- input_filt_rare$map_loaded %>%
  dplyr::select(17:37)
m <- cor(soil)
corrplot(m, 
         method = "square",
         type = "lower",
         diag = FALSE,
         hclust.method = "ward.D2",
         tl.cex = 0.5)

# PCA of soil
pca <- rda(soil, scale = TRUE)  # `scale=TRUE` standardizes variables
summary(pca)

plot(pca, scaling = 2, type = "n")  # type = "n" plots empty frame
points(pca, display = "sites", scaling = 2, pch = 21, bg = "lightblue")
species_scores <- scores(pca, display = "species", scaling = 2)
arrows(0, 0, species_scores[, 1], species_scores[, 2], length = 0.1, col = "red")
text(species_scores[, 1], species_scores[, 2], labels = rownames(species_scores), 
     col = "red", cex = 0.8, pos = 3)

sites_scores <- scores(pca, display = "sites", scaling = 2)
species_scores <- scores(pca, display = "species", scaling = 2)
vec.df <- as.data.frame(species_scores) %>%
  mutate(shortnames = rownames(.))
pcaA1 <- paste("PC1: ", round((eigenvals(pca)/sum(eigenvals(pca)))[1]*100, 1), "%")
pcaA2 <- paste("PC2: ", round((eigenvals(pca)/sum(eigenvals(pca)))[2]*100, 1), "%")
input_filt_rare$map_loaded$Axis01 <- sites_scores[,1]
input_filt_rare$map_loaded$Axis02 <- sites_scores[,2]
micro.hulls <- ddply(input_filt_rare$map_loaded, c("Treatment"), find_hull)
pdf("SoilPCA.pdf", width = 7, height = 5)
ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, shape = 21, aes(fill = Treatment)) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text_repel(data = vec.df,
                  aes(x = PC1, y = PC2, label = shortnames),
                  size = 3, color = "red") +
  labs(x = pcaA1, 
       y = pcaA2,
       fill = "Treatment") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 23, size = 4))) +
  theme_bw() +  
  theme(legend.position = "right",
        legend.box.margin = margin(0,0,100,0, unit = "pt"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())
dev.off()


# Stats
range(input_filt_rare$map_loaded$rich)
range(input_filt_rare$map_loaded$shannon)
m1 <- lm(rich ~ Treatment * Depth, data = input_filt_rare$map_loaded)
shapiro.test(m1$residuals)
summary(m1) # R2 = 0.44, p < 0.001
Anova(m1, type = "III") # Treatment, Depth, no int.
m1 <- lm(rich ~ MP1 + Basalt + MP1:Basalt + Depth, data = input_filt_rare$map_loaded)
shapiro.test(m1$residuals)
summary(m1) # R2 = 0.36, p < 0.001
Anova(m1, type = "III") # Depth and interaction

m2 <- lm(shannon ~ Treatment, data = input_filt_rare$map_loaded)
shapiro.test(m2$residuals)
summary(m2) # R2 = 0.28, p < 0.001
m2 <- lm(shannon ~ MP1 + Basalt + MP1:Basalt + Depth, data = input_filt_rare$map_loaded)
shapiro.test(m2$residuals)
summary(m2) # R2 = 0.46, p < 0.001
Anova(m2, type = "III") # Depth and interaction

# Pies
eta_sq_m1 <- eta_sq_glm(m1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("MP1", "Basalt", "Depth***", "I*", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
set.seed(106)
g1 <- ggplot(eta_sq_m1, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA, alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m1,
                  aes(y = ypos, label = group), color = "black", size = 1.5,
                  box.padding = 0.01) +
  scale_fill_manual(values = c("#35B779FF", "grey30", "#FDE725FF", "#31688EFF", "grey70")) +
  theme_void() + 
  theme(legend.position = "none")
g1

eta_sq_m2 <- eta_sq_glm(m2) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("MP1", "Basalt", "Depth***", "I**", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
set.seed(100)
g2 <- ggplot(eta_sq_m2, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA, alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m2,
                  aes(y = ypos, label = group), color = "black", size = 1.5,
                  box.padding = 0.01) +
  scale_fill_manual(values = c("#35B779FF", "grey30", "#FDE725FF", "#31688EFF", "grey70")) +
  theme_void() + 
  theme(legend.position = "none")
g2

# Plots
alpha_long <- input_filt_rare$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
facet_names <- c("rich" = "A. OTU Richness",
                 "shannon" = "B. Shannon Diversity")
g3 <- ggplot(alpha_long, aes(Treatment, value, shape = Depth)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment)) +
  geom_point(size = 3, alpha = 1, position = position_jitterdodge(),
             aes(fill = Treatment)) +
  labs(x = NULL,
       y = NULL) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24, 22)) +
  facet_wrap(~ name, scales = "free_y", labeller = as_labeller(facet_names)) +
  guides(fill = "none") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid = element_blank())
g3

# Combine
plot.with.inset <-
  ggdraw() +
  draw_plot(g3) +
  draw_plot(g1, x = 0.21, y = 0.12, width = 0.21, height = 0.21) +
  draw_plot(g2, x = 0.65, y = 0.12, width = 0.21, height = 0.21)
plot.with.inset
pdf("FinalFigs/FigureS6.pdf", width = 7, height = 5)
plot.with.inset
dev.off()
png("FinalFigs/FigureS6.png", width = 7, height = 5, units = "in", res = 300)
plot.with.inset
dev.off()



# Beta diversity
bc <- calc_dm(input_filt_rare$data_loaded)
hist(bc)
range(bc)
mean(bc)

# PERMANOVA
set.seed(1150)
ado1 <- adonis2(bc ~ input_filt_rare$map_loaded$MP1 +
                  input_filt_rare$map_loaded$Basalt +
                  input_filt_rare$map_loaded$MP1:input_filt_rare$map_loaded$Basalt +
                  input_filt_rare$map_loaded$Depth,
                by = "terms")
ado1

# dbRDA
names(input_filt_rare$map_loaded)
env <- input_filt_rare$map_loaded %>%
  dplyr::select(17:37) # Careful - adjust if updates are made
mod0 <- dbrda(bc ~ 1, env)  # Model with intercept only
mod1 <- dbrda(bc ~ ., env)  # Model with all explanatory variables
set.seed(100)
dbmod <- ordistep(mod0, scope = formula(mod1))
dbmod$anova # Na, Fe, tot_Na_ppm

# Mantel
dist.env <- as.matrix(dist(env, method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.env) <- input_filt_rare$map_loaded$sampleID
colnames(dist.env) <- input_filt_rare$map_loaded$sampleID
set.seed(100)
mantel(bc, dist.env, permutations = 2000)
qplot(as.dist(dist.env), bc, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Environmental Distance",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))

# PCoA
pcoa <- cmdscale(bc, k = nrow(input_filt_rare$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input_filt_rare$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input_filt_rare$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]

set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier <- 0.2
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals <= 0.05) %>%
  mutate(shortnames = c("pH", "OM", "K", "Ca", "Mg", "Na", "Fe", "Al",
                        "Ni", "NO3", "NH4", "CEC", "tot_K", "tot_Ca", "tot_Mg",
                        "tot_Al", "tot_Na", "tot_Ni", "Alk"))

df_pcoa <- data.frame(SampleID = rownames(pcoa$points),
                      PCoA1 = pcoa$points[, 1],
                      PCoA2 = pcoa$points[, 2],
                      Treatment = input_filt_rare$map_loaded$Treatment)
centroids <- aggregate(cbind(PCoA1, PCoA2) ~ Treatment, data = df_pcoa, FUN = mean)

eta_sq_m3 <- eta_sq_adonis(ado1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("MP1", "Basalt***", "Depth***", "I*", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g4 <- ggplot(eta_sq_m3, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA, alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m3,
                  aes(y = ypos, label = group), color = "black", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#35B779FF", "grey30", "#FDE725FF", "#31688EFF", "grey70")) +
  theme_void() + 
  theme(legend.position = "none")
g4

micro.hulls <- ddply(input_filt_rare$map_loaded, c("Treatment"), find_hull)
g5 <- ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, shape = 21, aes(fill = Treatment)) +
  geom_point(data = centroids, aes(PCoA1, PCoA2, fill = Treatment),
             size = 5, shape = 23, stroke = 1) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       fill = "Treatment") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 23, size = 4))) +
  theme_bw() +  
  theme(legend.position = "right",
        legend.box.margin = margin(0,0,100,0, unit = "pt"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())
g5

# Combine
plot.with.inset <-
  ggdraw() +
  draw_plot(g5) +
  draw_plot(g4, x = 0.75, y = 0.1, width = 0.3, height = 0.3)
plot.with.inset
pdf("Beta.pdf", width = 7, height = 5)
plot.with.inset
dev.off()

# Depth
micro.hulls <- ddply(input_filt_rare$map_loaded, c("Depth"), find_hull)
g5 <- ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, aes(colour = Depth, fill = Depth),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, aes(shape = Depth, fill = Depth)) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) +
  geom_text_repel(data = vec.df,
                  aes(x = Dim1, y = Dim2, label = shortnames),
                  size = 3, color = "red") +
  labs(x = pcoaA1, 
       y = pcoaA2) +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_colour_manual(values = c("#DDAA33", "#BB5566", "#004488")) +
  scale_fill_manual(values = c("#DDAA33", "#BB5566", "#004488")) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(0, 40, 0, 0),
        panel.grid = element_blank())
g5

plot.with.inset <-
  ggdraw() +
  draw_plot(g5) +
  draw_plot(g4, x = 0.75, y = 0.1, width = 0.3, height = 0.3)
plot.with.inset
pdf("FinalFigs/Figure7.pdf", width = 7, height = 5)
plot.with.inset
dev.off()
png("FinalFigs/Figure7.png", width = 7, height = 5, units = "in", res = 300)
plot.with.inset
dev.off()

# TrtDep
micro.hulls <- ddply(input_filt_rare$map_loaded, c("TrtDep"), find_hull)
ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = TrtDep, fill = TrtDep),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.75) +
  labs(x = pcoaA1, 
       y = pcoaA2) +
  scale_colour_manual(values = c(rep("#440154FF", 3),
                                 rep("#31688EFF", 3),
                                 rep("#35B779FF", 3),
                                 rep("#FDE725FF", 3))) +
  scale_fill_manual(values = c(rep("#440154FF", 3),
                               rep("#31688EFF", 3),
                               rep("#35B779FF", 3),
                               rep("#FDE725FF", 3))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))

# Check taxa
cliffplot_taxa_bars(input_filt_rare, 2, variable = "sampleID")
cliffplot_taxa_bars(input_filt_rare, 3, variable = "sampleID")
cliffplot_taxa_bars(input_filt_rare, 4, variable = "sampleID")
cliffplot_taxa_bars(input_filt_rare, 5, variable = "sampleID")
cliffplot_taxa_bars(input_filt_rare, 6, variable = "sampleID")

# Phyla
tax_sum_phyla <- summarize_taxonomy(input = input_filt_rare, 
                                    level = 2, 
                                    report_higher_tax = F)
bars_phyla <- plot_taxa_bars(tax_sum_phyla,
                             input_filt_rare$map_loaded,
                             "sampleID",
                             num_taxa = 12,
                             data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare$map_loaded, by = c("group_by" = "sampleID"))
top_phyla <- bars_phyla %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_phyla <- bars_phyla %>%
  mutate(taxon = factor(taxon,
                        levels = c("Other", rev(top_phyla$taxon))))
pdf("Taxa_Phyla.pdf", width = 8, height = 6)
ggplot(bars_phyla, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_nested_wrap(~ Treatment + Depth, scales = "free_x", ncol = 12) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 12),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"))
dev.off()

tax_sum_phyla_t <- as.data.frame(t(tax_sum_phyla))
input_filt_rare$map_loaded$Cyanobacteria <- tax_sum_phyla_t$Cyanobacteria
leveneTest(input_filt_rare$map_loaded$Cyanobacteria ~ input_filt_rare$map_loaded$Treatment) # Homogeneous
leveneTest(input_filt_rare$map_loaded$Cyanobacteria ~ input_filt_rare$map_loaded$Depth) # Homogeneous
m <- aov(Cyanobacteria ~ MP1 + Basalt + MP1:Basalt + Depth, data = input_filt_rare$map_loaded)
Anova(m, type = "III", singular.ok = TRUE) # Depth
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
pdf("Cyanobacteria.pdf", width = 7, height = 5)
ggplot(input_filt_rare$map_loaded, aes(Treatment, Cyanobacteria)) +
  geom_boxplot(outlier.shape = NA, aes(colour = Treatment)) +
  geom_jitter(size = 3, alpha = 1, width = 0.2, aes(shape = Depth, fill = Treatment)) +
  labs(x = "Treatment", y = "Cyanobacteria relative abundance") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24, 22)) +
  guides(fill = "none",
         colour = "none",
         shape = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14))
dev.off()

# Test effect of depth on phyla
phy_kw <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_phyla, 
                                      metadata_map = input_filt_rare$map_loaded, 
                                      type_header = "Depth", 
                                      test_type = "KW")
phy_kw_sig <- phy_kw %>%
  filter(pvalsBon < 0.05) %>%
  dplyr::select(4:6)

# Genus
tax_sum_genus <- summarize_taxonomy(input = input_filt_rare, 
                                    level = 6, 
                                    report_higher_tax = F)
bars_genus <- plot_taxa_bars(tax_sum_genus,
                             input_filt_rare$map_loaded,
                             "sampleID",
                             num_taxa = 12,
                             data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare$map_loaded, by = c("group_by" = "sampleID"))
top_genus <- bars_genus %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_genus <- bars_genus %>%
  mutate(taxon = factor(taxon,
                        levels = c("Other", rev(top_genus$taxon))))
pdf("Taxa_Genus.pdf", width = 8, height = 6)
ggplot(bars_genus, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_nested_wrap(~ Treatment + Depth, scales = "free_x", ncol = 12) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 12),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"))
dev.off()

tax_sum_genus_t <- as.data.frame(t(tax_sum_genus))
input_filt_rare$map_loaded$Bacillus <- tax_sum_genus_t$Bacillus
leveneTest(input_filt_rare$map_loaded$Bacillus ~ input_filt_rare$map_loaded$Treatment) # Homogeneous
leveneTest(input_filt_rare$map_loaded$Bacillus ~ input_filt_rare$map_loaded$Depth) # Homogeneous
m <- aov(Bacillus ~ MP1 + Basalt + MP1:Basalt + Depth, data = input_filt_rare$map_loaded)
summary(m) # All
Anova(m, type = "III", singular.ok = TRUE) # Interaction
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
pdf("Bacillus.pdf", width = 7, height = 5)
ggplot(input_filt_rare$map_loaded, aes(Treatment, Bacillus)) +
  geom_boxplot(outlier.shape = NA, aes(colour = Treatment)) +
  geom_jitter(size = 3, alpha = 1, width = 0.2, aes(shape = Depth, fill = Treatment)) +
  labs(x = "Treatment", y = "Bacillus relative abundance") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24, 22)) +
  guides(fill = "none",
         colour = "none",
         shape = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14))
dev.off()

input_filt_rare$map_loaded$Pseudomonas <- tax_sum_genus_t$Pseudomonas
leveneTest(input_filt_rare$map_loaded$Pseudomonas ~ input_filt_rare$map_loaded$Treatment) # Homogeneous
leveneTest(input_filt_rare$map_loaded$Pseudomonas ~ input_filt_rare$map_loaded$Depth) # Homogeneous
m <- aov(Pseudomonas ~ MP1 + Basalt + MP1:Basalt + Depth, data = input_filt_rare$map_loaded)
summary(m) # All
Anova(m, type = "III", singular.ok = TRUE) # Interaction
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
ggplot(input_filt_rare$map_loaded, aes(Treatment, Pseudomonas)) +
  geom_boxplot(outlier.shape = NA, aes(colour = Treatment)) +
  geom_jitter(size = 3, alpha = 1, width = 0.2, height = 0, 
              aes(shape = Depth, fill = Treatment)) +
  labs(x = "Treatment", y = "Pseudomonas relative abundance") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24, 22)) +
  guides(fill = "none",
         colour = "none",
         shape = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14))

input_filt_rare$map_loaded$Burkholderia <- tax_sum_genus_t$`Burkholderia-Caballeronia-Paraburkholderia`
leveneTest(input_filt_rare$map_loaded$Burkholderia ~ input_filt_rare$map_loaded$Treatment) # Homogeneous
leveneTest(input_filt_rare$map_loaded$Burkholderia ~ input_filt_rare$map_loaded$Depth) # Homogeneous
m <- aov(Burkholderia ~ MP1 + Basalt + MP1:Basalt + Depth, data = input_filt_rare$map_loaded)
summary(m) # All
Anova(m, type = "III", singular.ok = TRUE) # Interaction
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
ggplot(input_filt_rare$map_loaded, aes(Treatment, Burkholderia)) +
  geom_boxplot(outlier.shape = NA, aes(colour = Treatment)) +
  geom_jitter(size = 3, alpha = 1, width = 0.2, height = 0, 
              aes(shape = Depth, fill = Treatment)) +
  labs(x = "Treatment", y = "Burkholderia relative abundance") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24, 22)) +
  guides(fill = "none",
         colour = "none",
         shape = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14))

# Known rock weathering genera
# Use list from literature
# Peribacillus (Zhang et al. 2025)
# Bacillus (Andes, Corbett et al. 2024)
# Gluconobacter (Marecos et al. 2024)
# Cupriavidus metallidurans (Corbett et al. 2024)
# Cyanobacteria (Jones et al. 2023, McCutcheon & Power 2023)
cyano <- filter_taxa_from_input(input_filt_rare,
                                taxa_to_keep = "Cyanobacteria") # Get genera
# Shewanella (Lunstrum et al. 2023)
# Geobacter, Geothrix (Wild et al. 2022)
# Pseudomonas, Streptomyces, Burkholderia (Wild)
# Acidithiobacillus, Sulfobacillus, Leptospirillum, Ferroplasma, Acidianus, Sulfolobus (Wild)
# Acidithiobacillus ferrooxidans (Sekerci & Balci 2022)
# Samuels review chapter mentions some of these
# Brady et al. Hawaii volcanoes
# Pseudogulbenkiania, Ralstonia, Rhodopseudomonas (Napieralski)
# Burkholderia fungorum (Stranghoener et al. 2018)
# Pseudomonas sp. HerB (Popa et al. 2012)
# Pseudomonas aureofaciens (Pokrovsky et al. 2009)

# Uroz et al. 2009 list:
# Agrobacterium, Aminobacter, Azospirillum, Labrys, Rhanella, Rhizobium, Sphingomonas
# Achromobacter, Burkholderia, Collimonas, Janthinobacterium
# Acinetobacter, Azotobacter, Geobacter
# Acidithiobacillus, Citrobacter, Dyella, Enterobacter, Frateuria, Pseudomonas, Serratia, Shewanella
# Arthrobacter, Bacillus, Mycobacterium, Paenebacillus, Staphylococcus, Streptomyces

weather_tax <- c("Bacillus", "Peribacillus", "Gluconobacter", "Cupriavidus",
                 "Shewanella", "Vampirovibrio", "PhormidiumSAG37.90", 
                 "LeptolyngbyaEcFYyyy-00", "LeptolyngbyaEs-Yyy1000", 
                 "TychonemaCCAP1459-11B", "Geobacter", "Geothrix", "Pseudomonas", 
                 "Streptomyces", "Bulkholderia", "Acidithiobacillus", "Sulfobacillus",
                 "Leptospirillum", "Ferroplasma", "Acidianus", "Sulfolobus",
                 "Pseudogulbenkiania", "Ralstonia", "Rhodopseudomonas", 
                 "Agrobacterium", "Aminobacter", "Azospirillum", "Labrys", "Rhanella", 
                 "Rhizobium", "Sphingomonas", "Achromobacter", "Collimonas", 
                 "Janthinobacterium", "Acinetobacter", "Azotobacter", "Citrobacter", 
                 "Dyella", "Enterobacter", "Frateuria", "Serratia", "Arthrobacter", 
                 "Mycobacterium", "Paenebacillus", "Staphylococcus")
length(weather_tax)

wt <- tax_sum_genus %>%
  filter(rownames(.) %in% weather_tax) # 22/45 of those genera present
bars_genus <- plot_taxa_bars(wt,
                             input_filt_rare$map_loaded,
                             "sampleID",
                             num_taxa = 22,
                             data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare$map_loaded, by = c("group_by" = "sampleID"))
top_genus <- bars_genus %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_genus <- bars_genus %>%
  mutate(taxon = factor(taxon,
                        levels = c(rev(top_genus$taxon))))
nb.cols <- 22
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
pdf("Taxa_Weathering.pdf", width = 8, height = 6)
ggplot(bars_genus, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(0.001, 0.001)) + 
  facet_nested_wrap(~ Treatment + Depth, scales = "free_x", ncol = 12) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 12),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"))
dev.off()

input_wt <- filter_taxa_from_input(input_filt_rare,
                                   taxa_to_keep = weather_tax,
                                   at_spec_level = 6)
nrow(input_wt$data_loaded)
input_wt$map_loaded$rich <- specnumber(input_wt$data_loaded, MARGIN = 2)
input_wt$map_loaded$shannon <- diversity(input_wt$data_loaded, index = "shannon", 
                                         MARGIN = 2)

# No difference in numbers of weathering genera OTUs
m1 <- lm(rich ~ MP1 * Basalt, data = input_wt$map_loaded)
shapiro.test(m1$residuals)
Anova(m1, type = "III") # NSD
m2 <- lm(shannon ~ MP1 * Basalt, data = input_wt$map_loaded)
shapiro.test(m2$residuals)
Anova(m2, type = "III") # NSD
alpha_long <- input_wt$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
facet_names <- c("rich" = "OTU Richness",
                 "shannon" = "Shannon Diversity")
g3 <- ggplot(alpha_long, aes(Treatment, value, shape = Depth)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment)) +
  geom_point(size = 3, alpha = 1, position = position_jitterdodge(),
             aes(fill = Treatment)) +
  labs(x = "Treatment",
       y = NULL) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24, 22)) +
  facet_wrap(~ name, scales = "free_y", labeller = as_labeller(facet_names)) +
  guides(fill = "none") +
  theme_bw() +
  theme(legend.position = "right")
g3

# No diff in weatherer commmunity comp either
bc <- calc_dm(input_wt$data_loaded)
set.seed(1150)
ado1 <- adonis2(bc ~ input_wt$map_loaded$MP1 +
                  input_wt$map_loaded$Basalt +
                  input_wt$map_loaded$MP1:input_wt$map_loaded$Basalt +
                  input_wt$map_loaded$Depth,
                by = "terms")
ado1 # NSD
pcoa <- cmdscale(bc, k = nrow(input_wt$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input_wt$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input_wt$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
micro.hulls <- ddply(input_wt$map_loaded, c("Treatment"), find_hull)
g5 <- ggplot(input_wt$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, shape = 21, aes(fill = Treatment)) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       fill = "Treatment") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 23, size = 4))) +
  theme_bw() +  
  theme(legend.position = "right",
        legend.box.margin = margin(0,0,100,0, unit = "pt"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())
g5



#### ___10 ####
# Rerun alpha and beta just on 10 cm samples
top <- filter_data(input_filt_rare,
                   filter_cat = "Depth",
                   keep_vals = 10) # 24 samples remaining
bac <- filter_taxa_from_input(top,
                              at_spec_level = 6,
                              taxa_to_keep = "Bacillus")
#View(bac$taxonomy_loaded)
#write.table(bac$taxonomy_loaded$taxonomy8, "my_otus.txt", row.names = F, col.names = F)

# Stats
m1 <- lm(rich ~ MP1 * Basalt, data = top$map_loaded)
shapiro.test(m1$residuals)
Anova(m1, type = "III") # Interaction

m2 <- lm(shannon ~ MP1 * Basalt, data = top$map_loaded)
shapiro.test(m2$residuals)
Anova(m2, type = "III") # Depth and interaction

# Pies
eta_sq_m1 <- eta_sq_glm(m1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("MP1", "Basalt", "I", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
set.seed(100)
g1 <- ggplot(eta_sq_m1, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA, alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m1,
                  aes(y = ypos, label = group), color = "black", size = 1.5,
                  box.padding = 0.01) +
  scale_fill_manual(values = c("#31688EFF", "#FDE725FF", "#35B779FF", "grey70")) +
  theme_void() + 
  theme(legend.position = "none")
g1

eta_sq_m2 <- eta_sq_glm(m2) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("MP1", "Basalt", "I", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
set.seed(100)
g2 <- ggplot(eta_sq_m2, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA, alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m2,
                  aes(y = ypos, label = group), color = "black", size = 1.5,
                  box.padding = 0.01) +
  scale_fill_manual(values = c("#31688EFF", "#FDE725FF", "#35B779FF", "grey70")) +
  theme_void() + 
  theme(legend.position = "none")
g2

# Plots
alpha_long <- top$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
facet_names <- c("rich" = "OTU Richness",
                 "shannon" = "Shannon Diversity")
g3 <- ggplot(alpha_long, aes(Treatment, value)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment)) +
  geom_point(size = 3, alpha = 1, pch = 21, aes(fill = Treatment)) +
  labs(x = "Treatment",
       y = NULL) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(~ name, scales = "free_y", labeller = as_labeller(facet_names)) +
  guides(colour = "none") +
  theme_bw() +
  theme(legend.position = "right")
g3

# Combine
plot.with.inset <-
  ggdraw() +
  draw_plot(g3) +
  draw_plot(g1, x = 0.03, y = 0.15, width = 0.2, height = 0.2) +
  draw_plot(g2, x = 0.47, y = 0.15, width = 0.2, height = 0.2)
plot.with.inset
pdf("Alpha_10.pdf", width = 7, height = 5)
plot.with.inset
dev.off()



# Beta diversity
bc <- calc_dm(top$data_loaded)
hist(bc)
range(bc)
mean(bc)

# PERMANOVA
set.seed(1150)
ado1 <- adonis2(bc ~ top$map_loaded$MP1 * top$map_loaded$Basalt, by = "terms")
ado1 # No interaction
set.seed(1150)
ado1 <- adonis2(bc ~ top$map_loaded$MP1 + top$map_loaded$Basalt, by = "terms")
ado1

# dbRDA
names(top$map_loaded)
env <- top$map_loaded %>%
  dplyr::select(17:37) # Careful - adjust if updates are made
mod0 <- dbrda(bc ~ 1, env)  # Model with intercept only
mod1 <- dbrda(bc ~ ., env)  # Model with all explanatory variables
set.seed(100)
dbmod <- ordistep(mod0, scope = formula(mod1))
dbmod$anova # tot_Na_ppm

# Mantel
dist.env <- as.matrix(dist(env, method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.env) <- top$map_loaded$sampleID
colnames(dist.env) <- top$map_loaded$sampleID
set.seed(100)
mantel(bc, dist.env, permutations = 2000)
qplot(as.dist(dist.env), bc, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Environmental Distance",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))

# PCoA
pcoa <- cmdscale(bc, k = nrow(top$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
top$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
top$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]

set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals <= 0.05) %>%
  mutate(shortnames = c("pH", "Mg", "Na", "Fe",
                        "NO3", "tot_Ca", "tot_Mg",
                        "tot_Na"))

df_pcoa <- data.frame(SampleID = rownames(pcoa$points),
                      PCoA1 = pcoa$points[, 1],
                      PCoA2 = pcoa$points[, 2],
                      Treatment = top$map_loaded$Treatment)
centroids <- aggregate(cbind(PCoA1, PCoA2) ~ Treatment, data = df_pcoa, FUN = mean)

eta_sq_m3 <- eta_sq_adonis(ado1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  filter(group != "top$map_loaded$MP1:top$map_loaded$Basalt") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("MP1", "Basalt***", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g4 <- ggplot(eta_sq_m3, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA, alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m3,
                  aes(y = ypos, label = group), color = "black", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#35B779FF", "#31688EFF", "grey70")) +
  theme_void() + 
  theme(legend.position = "none")
micro.hulls <- ddply(top$map_loaded, c("Treatment"), find_hull)
g5 <- ggplot(top$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, shape = 21, aes(fill = Treatment)) +
  geom_point(data = centroids, aes(PCoA1, PCoA2, fill = Treatment),
             size = 5, shape = 23, stroke = 1) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) +
  geom_text_repel(data = vec.df,
                  aes(x = Dim1, y = Dim2, label = shortnames),
                  size = 3, color = "red",
                  box.padding = 0.1,
                  min.segment.length = 1,
                  force = 0.5,
                  max.time = 2,
                  seed = 1) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       fill = "Treatment") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 23, size = 4))) +
  ggtitle("A. 0-10 cm") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())
g5
leg <- cowplot::get_legend(g5)
leg
g5 <- g5 +
  theme(legend.position = "none")
g5

# Combine
plot.with.inset.10 <-
  ggdraw() +
  draw_plot(g5) +
  draw_plot(g4, x = 0.22, y = 0.1, width = 0.3, height = 0.3)
plot.with.inset.10
pdf("Beta_10.pdf", width = 7, height = 5)
plot.with.inset.10
dev.off()
pdf("Beta_10_wVec.pdf", width = 7, height = 5)
plot.with.inset.10
dev.off()
png("Beta_10_wVec.png", width = 7, height = 5, units = "in", res = 300)
plot.with.inset.10
dev.off()



#### Phylum level
tax_sum_phy <- summarize_taxonomy(input = top,
                                  level = 2,
                                  relative = T,
                                  report_higher_tax = F)
cliffplot_taxa_bars(top, 2, "Treatment")
top_phy <- tax_sum_phy %>%
  mutate(MeanAbund = rowMeans(.)) %>%
  rownames_to_column(var = "Phylum") %>%
  dplyr::select(Phylum, MeanAbund) %>%
  arrange(desc(MeanAbund))

# Make nice graph of top 12 phyla
bars_phyla <- plot_taxa_bars(tax_table = tax_sum_phy, 
                             metadata_map = input_filt_rare$map_loaded, 
                             type_header = "sampleID",
                             num_taxa = 13,
                             data_only = T) %>%
  left_join(., input_filt_rare$map_loaded, by = c("group_by" = "sampleID"))
top_phyla <- bars_phyla %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_phyla <- bars_phyla %>%
  mutate(taxon = factor(taxon,
                        levels = c("Other", rev(top_phyla$taxon))))
phyplot <- ggplot(bars_phyla, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_grid(~ Treatment, space = "free", scales = "free_x") +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 7),
        strip.background = element_rect(linewidth = 0.5),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.4, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2), "lines"))
pdf("Phylum_Barplot_10.pdf", width = 8, height = 6)
phyplot
dev.off()
png("Phylum_Barplot_10.png", width = 8, height = 6, units = "in", res = 300)
phyplot
dev.off()

# Make another one focusing on phyla 14-45.
tax_sum_phy <- summarize_taxonomy(input = top,
                                  level = 2,
                                  relative = T,
                                  report_higher_tax = F)
cliffplot_taxa_bars(top, 2, "Treatment")
top_phy <- tax_sum_phy %>%
  mutate(MeanAbund = rowMeans(.)) %>%
  rownames_to_column(var = "Phylum") %>%
  dplyr::select(Phylum, MeanAbund) %>%
  arrange(desc(MeanAbund))

# Make nice graph of top 12 phyla
tax_sum_phy_other <- tax_sum_phy %>%
  filter(rownames(.) %in% top_phy$Phylum[14:45])
bars_phyla <- plot_taxa_bars(tax_table = tax_sum_phy_other, 
                             metadata_map = input_filt_rare$map_loaded, 
                             type_header = "sampleID",
                             num_taxa = 32,
                             data_only = T) %>%
  left_join(., input_filt_rare$map_loaded, by = c("group_by" = "sampleID"))
top_phyla <- bars_phyla %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  #filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_phyla <- bars_phyla %>%
  mutate(taxon = factor(taxon,
                        levels = c("NA", rev(top_phyla$taxon))))
colrs <- c(colorRampPalette(brewer.pal(12, "Paired"))(30), "black", "grey")
phyplot2 <- ggplot(bars_phyla, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = colrs) +
  scale_y_continuous(expand = c(0.001, 0.001)) + 
  facet_grid(~ Treatment, space = "free", scales = "free_x") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 7),
        strip.background = element_rect(linewidth = 0.5),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.4, "cm"),
        panel.spacing.x = unit(c(0.2, 0.2, 0.2), "lines"))
phyplot2
pdf("Phylum_Barplot_Other_10.pdf", width = 8, height = 6)
phyplot2
dev.off()
png("Phylum_Barplot_Other_10.png", width = 8, height = 6, units = "in", res = 300)
phyplot2
dev.off()

pdf("FinalFigs/Figure9.pdf", width = 8.5, height = 9)
plot_grid(phyplot, phyplot2, labels = "AUTO")
dev.off()
png("FinalFigs/Figure9.png", width = 8.5, height = 6, units = "in", res = 300)
plot_grid(phyplot, phyplot2, labels = "AUTO")
dev.off()


bc <- calc_dm(tax_sum_phy)
set.seed(1150)
ado1 <- adonis2(bc ~ top$map_loaded$MP1 * top$map_loaded$Basalt, by = "terms")
ado1 # Basalt. No int.
set.seed(1150)
ado1 <- adonis2(bc ~ top$map_loaded$MP1 + top$map_loaded$Basalt, by = "terms")
ado1
pcoa <- cmdscale(bc, k = nrow(top$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
top$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
top$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]

set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals <= 0.05) %>%
  mutate(shortnames = c("pH", "K", "Mg", "Al",
                        "Ni", "NO3", "tot_K",
                        "tot_Ca", "tot_Mg", "tot_Na", "Alk"))

df_pcoa <- data.frame(SampleID = rownames(pcoa$points),
                      PCoA1 = pcoa$points[, 1],
                      PCoA2 = pcoa$points[, 2],
                      Treatment = top$map_loaded$Treatment)
centroids <- aggregate(cbind(PCoA1, PCoA2) ~ Treatment, data = df_pcoa, FUN = mean)

eta_sq_m3 <- eta_sq_adonis(ado1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("MP1", "Basalt***", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g4 <- ggplot(eta_sq_m3, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA, alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m3,
                  aes(y = ypos, label = group), color = "black", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#35B779FF", "#31688EFF", "grey70")) +
  theme_void() + 
  theme(legend.position = "none")
g4

micro.hulls <- ddply(top$map_loaded, c("Treatment"), find_hull)
g5 <- ggplot(top$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, shape = 21, aes(fill = Treatment)) +
  geom_point(data = centroids, aes(PCoA1, PCoA2, fill = Treatment),
             size = 5, shape = 23, stroke = 1) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) +
  geom_text_repel(data = vec.df,
                  aes(x = Dim1, y = Dim2, label = shortnames),
                  size = 3, color = "red") +
  labs(x = pcoaA1, 
       y = pcoaA2,
       fill = "Treatment") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 23, size = 4))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())
g5

# Combine
plot.with.inset <-
  ggdraw() +
  draw_plot(g5) +
  draw_plot(g4, x = 0.1, y = 0.1, width = 0.3, height = 0.3)
plot.with.inset
pdf("FinalFigs/FigureS7.pdf", width = 7, height = 5)
plot.with.inset
dev.off()
png("FinalFigs/FigureS7.png", width = 7, height = 5, units = "in", res = 300)
plot.with.inset
dev.off()



#### ___20 ####
mid <- filter_data(input_filt_rare,
                   filter_cat = "Depth",
                   keep_vals = 20) # 23 samples remaining
bc <- calc_dm(mid$data_loaded)
hist(bc)
range(bc)
mean(bc)

# PERMANOVA
set.seed(1150)
ado1 <- adonis2(bc ~ mid$map_loaded$MP1 * mid$map_loaded$Basalt, by = "terms")
ado1 # No interaction
set.seed(1150)
ado1 <- adonis2(bc ~ mid$map_loaded$MP1 + mid$map_loaded$Basalt, by = "terms")
ado1

# dbRDA
names(mid$map_loaded)
env <- mid$map_loaded %>%
  dplyr::select(17:37) # Careful - adjust if updates are made
mod0 <- dbrda(bc ~ 1, env)  # Model with intercept only
mod1 <- dbrda(bc ~ ., env)  # Model with all explanatory variables
set.seed(100)
dbmod <- ordistep(mod0, scope = formula(mod1))
dbmod$anova # Fe_ppm

# Mantel
dist.env <- as.matrix(dist(env, method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.env) <- mid$map_loaded$sampleID
colnames(dist.env) <- mid$map_loaded$sampleID
set.seed(100)
mantel(bc, dist.env, permutations = 2000)
qplot(as.dist(dist.env), bc, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Environmental Distance",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))

# PCoA
pcoa <- cmdscale(bc, k = nrow(mid$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
mid$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
mid$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]

set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals <= 0.05) %>%
  mutate(shortnames = c("pH", "Ca", "Mg", "Fe",
                        "NO3", "NH4", "tot_Ca",
                        "tot_Na"))

df_pcoa <- data.frame(SampleID = rownames(pcoa$points),
                      PCoA1 = pcoa$points[, 1],
                      PCoA2 = pcoa$points[, 2],
                      Treatment = mid$map_loaded$Treatment)
centroids <- aggregate(cbind(PCoA1, PCoA2) ~ Treatment, data = df_pcoa, FUN = mean)
eta_sq_m3 <- eta_sq_adonis(ado1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("MP1*", "Basalt***", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g4 <- ggplot(eta_sq_m3, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA, alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m3,
                  aes(y = ypos, label = group), color = "black", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#35B779FF", "#31688EFF", "grey70")) +
  theme_void() + 
  theme(legend.position = "none")
micro.hulls <- ddply(mid$map_loaded, c("Treatment"), find_hull)
g5 <- ggplot(mid$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, shape = 24, aes(fill = Treatment)) +
  geom_point(data = centroids, aes(PCoA1, PCoA2, fill = Treatment),
             size = 5, shape = 23, stroke = 1) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) +
  geom_text_repel(data = vec.df,
                  aes(x = Dim1, y = Dim2, label = shortnames),
                  size = 3, color = "red",
                  box.padding = 0.1,
                  min.segment.length = 1,
                  force = 0.5,
                  max.time = 2,
                  seed = 1) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       fill = "Treatment") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 23, size = 4))) +
  ggtitle("B. 10-20 cm") +
  theme_bw() +  
  theme(legend.position = "none",
        legend.box.margin = margin(0,0,100,0, unit = "pt"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())
plot.with.inset.20 <-
  ggdraw() +
  draw_plot(g5) +
  draw_plot(g4, x = 0.65, y = 0.62, width = 0.3, height = 0.3)
plot.with.inset.20
pdf("Beta_20_wVec.pdf", width = 7, height = 5)
plot.with.inset.20
dev.off()
png("Beta_20_wVec.png", width = 7, height = 5, units = "in", res = 300)
plot.with.inset.20
dev.off()



#### ___30 ####
bot <- filter_data(input_filt_rare,
                   filter_cat = "Depth",
                   keep_vals = 30) # 23 samples remaining
bc <- calc_dm(bot$data_loaded)
hist(bc)
range(bc)
mean(bc)

# PERMANOVA
set.seed(1150)
ado1 <- adonis2(bc ~ bot$map_loaded$MP1 * bot$map_loaded$Basalt, by = "terms")
ado1 # Interaction

# dbRDA
names(bot$map_loaded)
env <- bot$map_loaded %>%
  dplyr::select(17:37) # Careful - adjust if updates are made
mod0 <- dbrda(bc ~ 1, env)  # Model with intercept only
mod1 <- dbrda(bc ~ ., env)  # Model with all explanatory variables
set.seed(100)
dbmod <- ordistep(mod0, scope = formula(mod1))
dbmod$anova # Fe_ppm

# Mantel
dist.env <- as.matrix(dist(env, method = "euclidean", diag = FALSE, upper = FALSE))
rownames(dist.env) <- bot$map_loaded$sampleID
colnames(dist.env) <- bot$map_loaded$sampleID
set.seed(100)
mantel(bc, dist.env, permutations = 2000)
qplot(as.dist(dist.env), bc, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Environmental Distance",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))

# PCoA
pcoa <- cmdscale(bc, k = nrow(bot$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
bot$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
bot$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]

set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
ordiplot(pcoa)
plot(ef, cex = 0.5, p.max = 0.05)
multiplier <- ordiArrowMul(ef)
multiplier
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * multiplier,
         Dim2 = Dim2 * multiplier) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals <= 0.05) %>%
  mutate(shortnames = c("pH", "K", "Ca", "Na", "Fe", "Ni",
                        "NO3", "CEC", "tot_K", "tot_Ca", "tot_Mg",
                        "tot_Na", "Alk"))

df_pcoa <- data.frame(SampleID = rownames(pcoa$points),
                      PCoA1 = pcoa$points[, 1],
                      PCoA2 = pcoa$points[, 2],
                      Treatment = bot$map_loaded$Treatment)
centroids <- aggregate(cbind(PCoA1, PCoA2) ~ Treatment, data = df_pcoa, FUN = mean)

eta_sq_m3 <- eta_sq_adonis(ado1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("MP1", "Basalt**", "I*", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g4 <- ggplot(eta_sq_m3, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA, alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m3,
                  aes(y = ypos, label = group), color = "black", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#35B779FF", "#FDE725FF", "#31688EFF", "grey70")) +
  theme_void() + 
  theme(legend.position = "none")
micro.hulls <- ddply(bot$map_loaded, c("Treatment"), find_hull)
g5 <- ggplot(bot$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, shape = 22, aes(fill = Treatment)) +
  geom_point(data = centroids, aes(PCoA1, PCoA2, fill = Treatment),
             size = 5, shape = 23, stroke = 1) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) +
  geom_text_repel(data = vec.df,
                  aes(x = Dim1, y = Dim2, label = shortnames),
                  size = 3, color = "red",
                  box.padding = 0.1,
                  min.segment.length = 1,
                  force = 0.5,
                  max.time = 2,
                  seed = 1) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       fill = "Treatment") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 23, size = 4))) +
  ggtitle("C. 20-30 cm") +
  theme_bw() +  
  theme(legend.position = "none",
        legend.box.margin = margin(0,0,100,0, unit = "pt"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())
g5
plot.with.inset.30 <-
  ggdraw() +
  draw_plot(g5) +
  draw_plot(g4, x = 0.22, y = 0.64, width = 0.3, height = 0.3)
plot.with.inset.30
pdf("Beta_30.pdf", width = 7, height = 5)
plot.with.inset.30
dev.off()
pdf("Beta_30_wVec.pdf", width = 7, height = 5)
plot.with.inset.30
dev.off()
png("Beta_30_wVec.png", width = 7, height = 5, units = "in", res = 300)
plot.with.inset.30
dev.off()



# Combine the 3 depths (Figure 8)
pdf("FinalFigs/Figure8.pdf", width = 10, height = 4)
plot_grid(plot.with.inset.10, plot.with.inset.20, plot.with.inset.30, leg,
          ncol = 4,
          rel_widths = c(0.3, 0.3, 0.3, 0.1))
dev.off()
png("FinalFigs/Figure8.png", width = 10, height = 4, units = "in", res = 300)
plot_grid(plot.with.inset.10, plot.with.inset.20, plot.with.inset.30, leg,
          ncol = 4,
          rel_widths = c(0.3, 0.3, 0.3, 0.1))
dev.off()



#### ___DA ####
# Need to do differential abundance analysis
# Many methods - SIMPER, Multipatt, ANCOM, DESeq2, Aldex2, LFse

# SIMPER (simple, give top drivers of variation)
sim <- simper(t(top$data_loaded), top$map_loaded$Treatment)
sim_sum <- summary(sim)
sim_df1 <- head(sim_sum$`+MP1_Control`, n = 10) %>%
  mutate(Comparison = "Control_+MP1",
         OTU = rownames(.)) %>%
  rename("MeanControl" = avb,
         "MeanTrt" = ava)
sim_df2 <- head(sim_sum$`+Basalt_Control`, n = 10) %>%
  mutate(Comparison = "Control_+Basalt",
         OTU = rownames(.)) %>%
  rename("MeanControl" = avb,
         "MeanTrt" = ava)
sim_df3 <- head(sim_sum$`+MP1+Basalt_Control`, n = 10) %>%
  mutate(Comparison = "Control_+MP1+Basalt",
         OTU = rownames(.)) %>%
  rename("MeanControl" = avb,
         "MeanTrt" = ava)
simper_results <- rbind(sim_df1, sim_df2, sim_df3) %>%
  left_join(., top$taxonomy_loaded, by = c("OTU" = "taxonomy8")) %>%
  mutate(Response = ifelse(MeanTrt > MeanControl, "Positive", "Negative")) %>%
  rename(c("Domain" = "taxonomy1",
           "Phylum" = "taxonomy2",
           "Class" = "taxonomy3",
           "Order" = "taxonomy4",
           "Family" = "taxonomy5",
           "Genus" = "taxonomy6",
           "Species" = "taxonomy7",
           "CumulativeContribution" = "cumsum")) %>%
  mutate(MeanTrt = round((MeanTrt/2315)*100, digits = 2),
         MeanControl = round((MeanControl/2315)*100, digits = 2),
         CumulativeContribution = round(CumulativeContribution*100, digits = 2)) %>%
  unite(Taxonomy, Phylum, Class, Order, Family, Genus, OTU,
        sep = "; ") %>%
  mutate(Taxonomy = make.unique(Taxonomy)) %>%
  dplyr::select(Comparison, Response, Domain, Taxonomy, MeanTrt,
                MeanControl, CumulativeContribution)
simper_meta <- simper_results %>%
  column_to_rownames(var = "Taxonomy") %>%
  dplyr::select(Response, Comparison)
simper_mat <- simper_results %>%
  column_to_rownames(var = "Taxonomy") %>%
  dplyr::select(MeanControl, MeanTrt, CumulativeContribution) %>%
  as.matrix()
ann_rows <- data.frame(row.names = rownames(simper_meta),
                       "Response" = simper_meta$Response,
                       "Comparison" = simper_meta$Comparison)
ann_colors <- list(Response = c(Positive = "red",
                                Negative = "blue"),
                   Comparison = c(`Control_+MP1` = "#31688EFF", 
                                  `Control_+Basalt` = "#35B779FF", 
                                  `Control_+MP1+Basalt` = "#FDE725FF"))
pheatmap(simper_mat,
         legend = T,
         legend_breaks = c(1, 2, 3, 4, 5, 6, max(na.omit(simper_mat))),
         legend_labels = c("1   ", "2   ", "3   ", "4   ", "5   ", "6   ", "%\n"),
         main = "",
         border_color = NA,
         scale = "none",
         angle_col = 315,
         fontsize = 6,
         fontsize_row = 6,
         fontsize_col = 8,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         gaps_row = c(10, 20),
         filename = "Simper_Prok.png",
         width = 5,
         height = 5)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# ANCOMBC (use input_filt, BC for small sample size)
input_filt$map_loaded$Treatment <- factor(input_filt$map_loaded$Treatment,
                                          levels = c("Control",
                                                     "+Basalt",
                                                     "+MP1",
                                                     "+MP1+Basalt"))
input_filt_top <- filter_data(input_filt,
                              filter_cat = "Depth",
                              keep_vals = 10) # 24 samples remaining
otu <- phyloseq::otu_table(input_filt_top$data_loaded, taxa_are_rows = T)
names(input_filt_top$taxonomy_loaded) <- c("Kingdom", "Phylum", "Class", "Order",
                                           "Family", "Genus", "Species", "OTU")
tax <- phyloseq::tax_table(as.matrix(input_filt_top$taxonomy_loaded))
map <- phyloseq::sample_data(input_filt_top$map_loaded)
input.phy <- phyloseq::phyloseq(otu, tax, map)
input.phy
levels(input.phy@sam_data$Treatment)

# OTU
set.seed(123)
out <- ancombc2(data = input.phy,
                assay_name = "counts",
                tax_level = NULL,
                fix_formula = "Treatment",
                rand_formula = NULL,
                p_adj_method = "holm",
                pseudo = 0,
                pseudo_sens = FALSE,
                prv_cut = 0.10,
                lib_cut = 0,
                s0_perc = 0.05,
                group = "Treatment", 
                struc_zero = TRUE,
                neg_lb = TRUE, 
                alpha = 0.05, 
                n_cl = 1, 
                verbose = TRUE)
saveRDS(out, file = "ancom_out_prok_all.rds")
out <- readRDS("ancom_out_prok_all.rds")
# Structural zeroes
sz <- out$zero_ind
res <- out$res %>%
  mutate(SigBasalt = ifelse(`p_Treatment+Basalt` < 0.05,
                            "Significant", "Not significant")) %>%
  mutate(SigMP1 = ifelse(`p_Treatment+MP1` < 0.05,
                         "Significant", "Not significant")) %>%
  mutate(SigComb = ifelse(`p_Treatment+MP1+Basalt` < 0.05,
                          "Significant", "Not significant"))
SigBasalt <- res %>%
  filter(`p_Treatment+Basalt` < 0.05)
SigMP1 <- res %>%
  filter(`p_Treatment+MP1` < 0.05)
SigComb <- res %>%
  filter(`p_Treatment+MP1+Basalt` < 0.05)
ggplot(res, aes(reorder(taxon, `lfc_Treatment+Basalt`), `lfc_Treatment+Basalt`,
                colour = SigBasalt)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1.5) +
  geom_point(pch = 16, alpha = 0.75) +
  geom_text_repel(data = SigBasalt, aes(taxon, `lfc_Treatment+Basalt`, label = taxon),
                  show.legend = F, color = "blue", size = 2) +
  labs(x = NULL, y = "LFC", colour = "Basalt") +
  scale_colour_manual(values = c("grey", "red")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank())
ggplot(res, aes(reorder(taxon, `lfc_Treatment+MP1`), `lfc_Treatment+MP1`,
                colour = SigMP1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1.5) +
  geom_point(pch = 16, alpha = 0.75) +
  geom_text_repel(data = SigMP1, aes(taxon, `lfc_Treatment+MP1`, label = taxon),
                  show.legend = F, color = "blue", size = 2) +
  labs(x = NULL, y = "LFC", colour = "MP1") +
  scale_colour_manual(values = c("grey", "red")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank())
ggplot(res, aes(reorder(taxon, `lfc_Treatment+MP1+Basalt`), `lfc_Treatment+MP1+Basalt`,
                colour = SigComb)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1.5) +
  geom_point(pch = 16, alpha = 0.75) +
  geom_text_repel(data = SigComb, aes(taxon, `lfc_Treatment+MP1+Basalt`, label = taxon),
                  show.legend = F, color = "blue", size = 2) +
  labs(x = NULL, y = "LFC", colour = "+MP1+Basalt") +
  scale_colour_manual(values = c("grey", "red")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank())
SigBasalt <- SigBasalt %>%
  mutate(Comparison = "Control_+Basalt") %>%
  mutate(LFC = `lfc_Treatment+Basalt`,
         se = `se_Treatment+Basalt`) %>%
  dplyr::select(taxon, Comparison, LFC, se)
SigMP1 <- SigMP1 %>%
  mutate(Comparison = "Control_+MP1") %>%
  mutate(LFC = `lfc_Treatment+MP1`,
         se = `se_Treatment+MP1`) %>%
  dplyr::select(taxon, Comparison, LFC, se)
SigComb <- SigComb %>%
  mutate(Comparison = "Control_+MP1+Basalt") %>%
  mutate(LFC = `lfc_Treatment+MP1+Basalt`,
         se = `se_Treatment+MP1+Basalt`) %>%
  dplyr::select(taxon, Comparison, LFC, se)

# Venn Diagram
venn_list <- list(
  "Control vs. +Basalt" = SigBasalt$taxon,
  "Control vs. +MP1" = SigMP1$taxon,
  "Control vs. +MP1+Basalt" = SigComb$taxon
)
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,   # NULL means draw to R plotting window
  fill = c("#31688EFF", "#35B779FF", "#FDE725FF"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-20, 20, 0)
)
grid::grid.draw(venn.plot)
png("VennDiagram.png", width = 7, height = 5, units = "in", res = 300)
ggvenn(venn_list,
       fill_color = c("#31688EFF", "#35B779FF", "#FDE725FF"),
       stroke_size = 0.5, set_name_size = 4, text_size = 4)
dev.off()

# Other plots
ancom_results <- rbind(SigBasalt, SigMP1, SigComb) %>%
  dplyr::left_join(., input_filt_top$taxonomy_loaded, by = c("taxon" = "OTU")) %>%
  mutate(Response = ifelse(LFC > 0, "Positive", "Negative"),
         OTU = taxon) %>%
  dplyr::rename(Domain = Kingdom) %>%
  unite(Taxonomy, Phylum, Class, Order, Family, Genus, OTU,
        sep = "; ") %>%
  mutate(Taxonomy = make.unique(Taxonomy)) %>%
  dplyr::select(Comparison, Response, Domain, Taxonomy, LFC)
ancom_meta <- ancom_results %>%
  column_to_rownames(var = "Taxonomy") %>%
  dplyr::select(Domain, Response, Comparison)
ancom_mat <- ancom_results %>%
  column_to_rownames(var = "Taxonomy") %>%
  dplyr::select(LFC) %>%
  as.matrix()
ann_rows <- data.frame(row.names = rownames(ancom_meta),
                       "Domain" = ancom_meta$Domain,
                       "Response" = ancom_meta$Response,
                       "Comparison" = ancom_meta$Comparison)
ann_colors <- list(Domain = c(Bacteria = "#440154FF",
                              Archaea = "#FCFDBFFF"),
                   Response = c(Positive = "red",
                                Negative = "blue"),
                   Comparison = c(`Control_+Basalt` = "#31688EFF", 
                                  `Control_+MP1` = "#35B779FF", 
                                  `Control_+MP1+Basalt` = "#FDE725FF"))
pheatmap(ancom_mat,
         legend = T,
         border_color = NA,
         scale = "none",
         angle_col = 315,
         fontsize = 6,
         fontsize_row = 3,
         fontsize_col = 8,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         gaps_row = c(71, 157),
         filename = "Ancom_Prok.png",
         width = 8,
         height = 12)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

ancom_results <- ancom_results %>%
  #separate(Taxonomy, into = c("Taxonomy", "Junk"), sep = "\\.") %>%
  #dplyr::select(-Junk) %>%
  mutate(Taxonomy = as.factor(Taxonomy)) %>%
  separate(Taxonomy, remove = F, into = c("Phylum", "Class", "Order", "Family",
                                          "Genus", "OTU"), sep = "; ") %>%
  mutate(Phylum = as.factor(Phylum))
  
pdf("AncomLFC_Prok.pdf", width = 8, height = 12)
ggplot(ancom_results, aes(reorder(Taxonomy, LFC), LFC)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "LFC") +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  facet_wrap(~ Comparison, ncol = 3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 4),
        strip.text = element_text(size = 5),
        panel.grid.minor = element_blank())
dev.off()
png("AncomLFC_Prok.png", width = 8, height = 12, units = "in", res = 300)
ggplot(ancom_results, aes(reorder(Taxonomy, LFC), LFC)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "LFC") +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  facet_wrap(~ Comparison, ncol = 3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 4),
        strip.text = element_text(size = 5),
        panel.grid.minor = element_blank())
dev.off()

pdf("AncomLFC_Prok_tax.pdf", width = 8, height = 12)
ggplot(ancom_results, aes(Taxonomy, LFC)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "LFC") +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  facet_grid(Phylum ~ Comparison, scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 4),
        strip.text.x = element_text(size = 5),
        strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        panel.grid.minor = element_blank())
dev.off()
png("AncomLFC_Prok_tax.png", width = 8, height = 12, units = "in", res = 300)
ggplot(ancom_results, aes(Taxonomy, LFC)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "LFC") +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  facet_grid(Phylum ~ Comparison, scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 4),
        strip.text.x = element_text(size = 5),
        strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        panel.grid.minor = element_blank())
dev.off()
write_xlsx(ancom_results, "ancom_results.xlsx", format_headers = F)

length(unique(ancom_results$Taxonomy))
# Need to plot these 33 to confirm the ANCOM results and see how abundant they are
ancom_results <- ancom_results %>%
  separate(Taxonomy, remove = F, into = c("Phylum", "Class", "Order",
                                          "Family", "Genus", "OTU"), sep = "; ")
ancom_results2 <- ancom_results %>%
  dplyr::select(Taxonomy, Comparison, Response) %>%
  pivot_wider(names_from = Comparison, values_from = Response) %>%
  mutate(Comb = paste(`Control_+Basalt`, `Control_+MP1`, `Control_+MP1+Basalt`, sep = "_"))
levels(as.factor(ancom_results2$Comb)) # 11 combinations, too many
tax_sum_otu <- summarize_taxonomy(top, level = 8, relative = T, report_higher_tax = F)
ancom_tax_sum <- tax_sum_otu %>%
  filter(rownames(.) %in% ancom_results$OTU) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Treatment = top$map_loaded$Treatment) %>%
  pivot_longer(cols = c(1:33))
ggplot(ancom_tax_sum, aes(Treatment, value*100)) +
  geom_boxplot(outlier.shape = NA, aes(colour = Treatment)) +
  geom_jitter(size = 1, pch = 16, alpha = 1, width = 0.2, height = 0, aes(colour = Treatment)) +
  labs(x = "Treatment", y = "% abundance") +
  scale_colour_viridis_d() +
  facet_wrap(~ name, ncol = 11) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 6))

# Too many, break into 3 groups
ancom_basalt <- ancom_results %>%
  filter(Comparison == "Control_+Basalt")
info <- ancom_basalt %>%
  dplyr::select(OTU, Response) %>%
  arrange(OTU) %>%
  mutate(Color = Response) %>%
  mutate(Color = gsub("Positive", "#619CFF", Color)) %>%
  mutate(Color = gsub("Negative", "#F8766D", Color))
strip_col <- strip_themed(background_x = elem_list_rect(fill = info$Color))
ancom_tax_sum <- tax_sum_otu %>%
  filter(rownames(.) %in% ancom_basalt$OTU) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Treatment = top$map_loaded$Treatment) %>%
  pivot_longer(cols = c(1:13)) %>%
  left_join(., info, by = c("name" = "OTU"))
pdf("SigTaxa_Abund_Basalt.pdf", width = 7, height = 5)
ggplot(ancom_tax_sum, aes(Treatment, value*100)) +
  geom_boxplot(outlier.shape = NA, aes(colour = Treatment), show.legend = F) +
  geom_jitter(size = 1, pch = 16, alpha = 1, width = 0.2, height = 0, 
              aes(colour = Treatment), show.legend = F) +
  labs(x = NULL, y = "% abundance") +
  scale_colour_viridis_d() +
  facet_wrap2(~ name, ncol = 7, strip = strip_col) +
  ggtitle("Significant Taxa: Control vs. +Basalt") +
  theme_bw() +
  theme(legend.position.inside = c(0.8, 0.1),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 7))
dev.off()

ancom_MP1 <- ancom_results %>%
  filter(Comparison == "Control_+MP1")
info <- ancom_MP1 %>%
  dplyr::select(OTU, Response) %>%
  arrange(OTU) %>%
  mutate(Color = Response) %>%
  mutate(Color = gsub("Positive", "#619CFF", Color)) %>%
  mutate(Color = gsub("Negative", "#F8766D", Color))
strip_col <- strip_themed(background_x = elem_list_rect(fill = info$Color))
ancom_tax_sum <- tax_sum_otu %>%
  filter(rownames(.) %in% ancom_MP1$OTU) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Treatment = top$map_loaded$Treatment) %>%
  pivot_longer(cols = c(1:13)) %>%
  left_join(., info, by = c("name" = "OTU"))
pdf("SigTaxa_Abund_MP1.pdf", width = 7, height = 5)
ggplot(ancom_tax_sum, aes(Treatment, value*100)) +
  geom_boxplot(outlier.shape = NA, aes(colour = Treatment), show.legend = F) +
  geom_jitter(size = 1, pch = 16, alpha = 1, width = 0.2, height = 0, 
              aes(colour = Treatment), show.legend = F) +
  labs(x = NULL, y = "% abundance") +
  scale_colour_viridis_d() +
  facet_wrap2(~ name, ncol = 7, strip = strip_col) +
  ggtitle("Significant Taxa: Control vs. +MP1") +
  theme_bw() +
  theme(legend.position.inside = c(0.8, 0.1),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 7))
dev.off()

ancom_MP1Basalt <- ancom_results %>%
  filter(Comparison == "Control_+MP1+Basalt")
info <- ancom_MP1Basalt %>%
  dplyr::select(OTU, Response) %>%
  arrange(OTU) %>%
  mutate(Color = Response) %>%
  mutate(Color = gsub("Positive", "#619CFF", Color)) %>%
  mutate(Color = gsub("Negative", "#F8766D", Color))
strip_col <- strip_themed(background_x = elem_list_rect(fill = info$Color))
ancom_tax_sum <- tax_sum_otu %>%
  filter(rownames(.) %in% ancom_MP1Basalt$OTU) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Treatment = top$map_loaded$Treatment) %>%
  pivot_longer(cols = c(1:13)) %>%
  left_join(., info, by = c("name" = "OTU"))
pdf("SigTaxa_Abund_MP1Basalt.pdf", width = 7, height = 5)
ggplot(ancom_tax_sum, aes(Treatment, value*100)) +
  geom_boxplot(outlier.shape = NA, aes(colour = Treatment), show.legend = F) +
  geom_jitter(size = 1, pch = 16, alpha = 1, width = 0.2, height = 0, 
              aes(colour = Treatment), show.legend = F) +
  labs(x = NULL, y = "% abundance") +
  scale_colour_viridis_d() +
  facet_wrap2(~ name, ncol = 7, strip = strip_col) +
  ggtitle("Significant Taxa: Control vs. +MP1+Basalt") +
  theme_bw() +
  theme(legend.position.inside = c(0.8, 0.1),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 7))
dev.off()

# Make heatmap of sig taxa from Control vs MP1
miss <- SigMP1 %>%
  filter(taxon %notin% ancom_MP1$OTU)
tax_sum_otu_DA <- tax_sum_otu %>%
  filter(rownames(.) %in% SigMP1$taxon) %>%
  dplyr::select(D.R1.10, D.R2.10, D.R3.10, D.R4.10, D.R5.10, D.R6.10,
                B.R1.10, B.R2.10, B.R3.10, B.R4.10, B.R5.10, B.R6.10,
                A.R1.10, A.R2.10, A.R3.10, A.R4.10, A.R5.10, A.R6.10,
                C.R1.10, C.R2.10, C.R3.10, C.R4.10, C.R5.10, C.R6.10)
meta10 <- top$map_loaded %>%
  arrange(SampleID = colnames(tax_sum_otu_DA)) %>%
  arrange(SampleID = colnames(tax_sum_otu_DA))
colSums(tax_sum_otu)
colSums(tax_sum_otu_DA)
colnames(tax_sum_otu_DA) == meta10$sampleID

ann_cols <- data.frame(row.names = colnames(tax_sum_otu_DA), 
                       Treatment = meta10$Treatment)
ann_colors <- list(Treatment = c("Control" = "#440154FF",
                                 "+Basalt" = "#31688EFF",
                                 "+MP1" = "#35B779FF",
                                 "+MP1+Basalt" = "#FDE725FF"))
max(tax_sum_otu_DA)
min(tax_sum_otu_DA[tax_sum_otu_DA > 0], na.rm = TRUE)
breaks <- c(0, seq(0.0001075384, 0.006774922, length.out = 100))
colors <- c("white", colorRampPalette(brewer.pal(9, "Reds"))(99))
pheatmap(tax_sum_otu_DA,
         color = colors,
         breaks = breaks,
         legend = T,
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         fontsize_row = 6,
         fontsize_col = 6,
         gaps_col = c(6, 12, 18),
         na_col = "white",
         border_color = "white",
         filename = "ANCOM_Heatmap_NoClust.png",
         width = 10,
         height = 8)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### __Euk ####
input_filt <- readRDS("input_filt_mtagsEuk_all.rds")
input_filt$map_loaded$Treatment <- factor(input_filt$map_loaded$Treatment,
                                          levels = c("Control",
                                                     "+Basalt",
                                                     "+MP1",
                                                     "+MP1+Basalt"))
input_filt$map_loaded$Depth <- as.factor(input_filt$map_loaded$Depth)
input_filt$map_loaded$TrtDep <- factor(paste(input_filt$map_loaded$Treatment,
                                             input_filt$map_loaded$Depth,
                                             sep = "_"),
                                       levels = c("Control_10", "Control_20", "Control_30",
                                                  "+Basalt_10", "+Basalt_20", "+Basalt_30",
                                                  "+MP1_10", "+MP1_20", "+MP1_30",
                                                  "+MP1+Basalt_10", "+MP1+Basalt_20",
                                                  "+MP1+Basalt_30"))

alpha_long <- input_filt$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
facet_names <- c("rich" = "OTU Richness",
                 "shannon" = "Shannon Diversity")
ggplot(alpha_long, aes(Treatment, value, shape = Depth)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment)) +
  geom_point(size = 3, alpha = 1, position = position_jitterdodge(),
             aes(fill = Treatment)) +
  labs(x = "Treatment",
       y = NULL) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24, 22)) +
  facet_wrap(~ name, scales = "free_y", labeller = as_labeller(facet_names)) +
  guides(fill = "none") +
  theme_bw() +
  theme(legend.position = "right")

bc <- calc_dm(input_filt$data_loaded)
pcoa <- cmdscale(bc, k = nrow(input_filt$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input_filt$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input_filt$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
df_pcoa <- data.frame(SampleID = rownames(pcoa$points),
                      PCoA1 = pcoa$points[, 1],
                      PCoA2 = pcoa$points[, 2],
                      Treatment = input_filt$map_loaded$Treatment)
centroids <- aggregate(cbind(PCoA1, PCoA2) ~ Treatment, data = df_pcoa, FUN = mean)
micro.hulls <- ddply(input_filt$map_loaded, c("Treatment"), find_hull)
ggplot(input_filt$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, shape = 21, aes(fill = Treatment)) +
  geom_point(data = centroids, aes(PCoA1, PCoA2, fill = Treatment),
             size = 5, shape = 23, stroke = 1) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       fill = "Treatment") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 23, size = 4))) +
  theme_bw() +  
  theme(legend.position = "right",
        legend.box.margin = margin(0,0,100,0, unit = "pt"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())

cliffplot_taxa_bars(input_filt, 2, variable = "sampleID")
cliffplot_taxa_bars(input_filt, 3, variable = "sampleID")
cliffplot_taxa_bars(input_filt, 4, variable = "sampleID")
cliffplot_taxa_bars(input_filt, 5, variable = "sampleID")
cliffplot_taxa_bars(input_filt, 6, variable = "sampleID")

tax_sum_phyla <- summarize_taxonomy(input = input_filt, 
                                    level = 2, 
                                    report_higher_tax = F)
bars_phyla <- plot_taxa_bars(tax_sum_phyla,
                             input_filt_rare$map_loaded,
                             "sampleID",
                             num_taxa = 12,
                             data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare$map_loaded, by = c("group_by" = "sampleID"))
top_phyla <- bars_phyla %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_phyla <- bars_phyla %>%
  mutate(taxon = factor(taxon,
                        levels = c("Other", rev(top_phyla$taxon))))
pdf("Taxa_Phyla_Euk.pdf", width = 8, height = 6)
ggplot(bars_phyla, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_nested_wrap(~ Treatment + Depth, scales = "free_x", ncol = 12) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 12),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"))
dev.off()



#### __Kaiju (L1) ####
# Don't analyze mTAGs 18S - not enough data
# Use Kaiju Fungal db classifications

# Taxa
# Knufia petricola and Aspergillus niger (Corbett et al. 2024, Gerrits et al. 2020)
# Trichoderma guizhouense (Wild et al. 2022)
# Paxillus involutus (Wild)
# Talaromyces flavus (Wild)

metadata <- read.delim("metadata.txt")

perc_fun <- read.delim("kaiju_fungi_genus.tsv") %>%
  filter(taxon_name == "unclassified") %>%
  separate(file, remove = T, sep = "_", into = c("sampleID", "Junk")) %>%
  separate(sampleID, remove = F, sep = "-", into = c("Group", "Replicate", "Depth")) %>%
  mutate(sampleID = gsub("-", ".", sampleID)) %>%
  mutate(perc_fun = 100 - percent) %>%
  dplyr::select(-taxon_id, -Junk) %>%
  mutate(Treatment = ifelse(Group == "A", "+MP1",
                            ifelse(Group == "B", "+Basalt",
                                   ifelse(Group == "C", "+MP1+Basalt",
                                          ifelse(Group == "D", "Control", "Blank"))))) %>%
  mutate(MP1 = ifelse(Group == "A" | Group == "C", "MP1",
                      ifelse(Group == "B" | Group == "D", "Control", "Blank"))) %>%
  mutate(Basalt = ifelse(Group == "B" | Group == "C", "Basalt", 
                         ifelse(Group == "A" | Group == "D", "Control", "Blank"))) %>%
  filter(Treatment != "Blank") %>%
  filter(perc_fun > 0)
ggplot(perc_fun, aes(Treatment, perc_fun)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment, shape = Depth)) +
  geom_point(size = 3, alpha = 1, position = position_dodge(width = 0.75),
             aes(fill = Treatment, shape = Depth)) +
  labs(x = "Treatment",
       y = "% Fungi") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24, 22)) +
  guides(fill = "none") +
  theme_bw() +
  theme(legend.position = "right")

fun_gen <- read.delim("kaiju_fungi_genus.tsv") %>%
  filter(taxon_name != "unclassified") %>%
  filter(taxon_name != "Viruses") %>%
  filter(taxon_name != "cannot be assigned to a (non-viral) genus") %>%
  separate(file, remove = T, sep = "_", into = c("sampleID", "Junk")) %>%
  mutate(sampleID = gsub("-", ".", sampleID)) %>%
  dplyr::select(-taxon_id, -Junk, -percent) %>%
  pivot_wider(names_from = sampleID, values_from = reads) %>%
  dplyr::rename(taxonomy = taxon_name) %>%
  dplyr::select(-Ext.Blank.1, -Ext.Blank.2, -Ext.Blank.3) %>%
  mutate(OTU_ID = paste0("OTU_", row_number())) %>%
  mutate(C.R1.20 = 0) %>%
  mutate(B.R1.30 = 0) %>%
  dplyr::select(OTU_ID, all_of(metadata$Sample), taxonomy) %>%
  as.data.frame()
sum(names(fun_gen) %in% metadata$Sample)
out_fp <- "seqtab_wTax_mctoolsr_kaijuGen.txt"
names(fun_gen)[1] = "#OTU_ID"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(fun_gen,
                             out_fp,
                             sep = "\t",
                             row.names = FALSE,
                             append = TRUE))

tax_table_fp <- "seqtab_wTax_mctoolsr_kaijuGen.txt"
map_fp <- "metadata.txt"
input = load_taxa_table(tax_table_fp, map_fp) # Error. Do manually
input = list()
input$map_loaded <- metadata
rownames(input$map_loaded) <- input$map_loaded$Sample
input$data_loaded <- as.data.frame(fun_gen[,2:73])
rownames(input$data_loaded) <- fun_gen$OTU_ID
input$taxonomy_loaded <- data.frame(taxonomy1 = "Fungi",
                                    taxonomy2 = "Phylum",
                                    taxonomy3 = "Class",
                                    taxonomy4 = "Order",
                                    taxonomy5 = "Family",
                                    taxonomy6 = fun_gen$taxonomy,
                                    taxonomy7 = fun_gen$OTU_ID)
rownames(input$taxonomy_loaded) <- fun_gen$OTU_ID

# Remove singletons and doubletons - none
singdoub <- data.frame("count" = rowSums(input$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.))

# Remove C.R1.20 and B.R1.30, very few reads, (already knew this, few raw reads too)
input_filt <- input

# Rarefaction
rarecurve(t(input_filt$data_loaded), step = 500, label = F)
sort(colSums(input_filt$data_loaded))
mean(colSums(input_filt$data_loaded)) # 263815.8
se(colSums(input_filt$data_loaded)) # 6343.258
set.seed(530)
input_filt_rare <- single_rarefy(input_filt, 184666) # n = 70 still
sort(colSums(input_filt_rare$data_loaded))

# Add rarefied richness and Shannon
input_filt_rare$map_loaded$rich <- specnumber(input_filt_rare$data_loaded, MARGIN = 2)
input_filt_rare$map_loaded$shannon <- diversity(input_filt_rare$data_loaded, 
                                                index = "shannon", MARGIN = 2)
input_filt_rare$map_loaded$rich # They all of 291
input_filt_rare$map_loaded$shannon

# Save
#saveRDS(input_filt_rare, "input_filt_rare_kaijuGen.rds")

input_filt_rare <- readRDS("input_filt_rare_kaijuGen.rds")
input_filt_rare$map_loaded$Treatment <- factor(input_filt_rare$map_loaded$Treatment,
                                               levels = c("Control",
                                                          "+Basalt",
                                                          "+MP1",
                                                          "+MP1+Basalt"))
input_filt_rare$map_loaded$Depth <- as.factor(input_filt_rare$map_loaded$Depth)
input_filt_rare$map_loaded$TrtDep <- factor(paste(input_filt_rare$map_loaded$Treatment,
                                                  input_filt_rare$map_loaded$Depth,
                                                  sep = "_"),
                                            levels = c("Control_10", "Control_20", "Control_30",
                                                       "+Basalt_10", "+Basalt_20", "+Basalt_30",
                                                       "+MP1_10", "+MP1_20", "+MP1_30",
                                                       "+MP1+Basalt_10", "+MP1+Basalt_20",
                                                       "+MP1+Basalt_30"))

# Check beta diversity
bc <- calc_dm(input_filt_rare$data_loaded)
set.seed(1150)
ado1 <- adonis2(bc ~ input_filt_rare$map_loaded$MP1 +
                  input_filt_rare$map_loaded$Basalt +
                  input_filt_rare$map_loaded$MP1:input_filt_rare$map_loaded$Basalt +
                  input_filt_rare$map_loaded$Depth,
                by = "terms")
ado1
pcoa <- cmdscale(bc, k = nrow(input_filt_rare$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input_filt_rare$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input_filt_rare$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
df_pcoa <- data.frame(SampleID = rownames(pcoa$points),
                      PCoA1 = pcoa$points[, 1],
                      PCoA2 = pcoa$points[, 2],
                      Treatment = input_filt_rare$map_loaded$Treatment)
centroids <- aggregate(cbind(PCoA1, PCoA2) ~ Treatment, data = df_pcoa, FUN = mean)
eta_sq_m3 <- eta_sq_adonis(ado1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("MP1", "Basalt***", "Depth**", "I*", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g4 <- ggplot(eta_sq_m3, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA, alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m3,
                  aes(y = ypos, label = group), color = "black", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#31688EFF", "grey30", "#FDE725FF", "#35B779FF", "grey70")) +
  theme_void() + 
  theme(legend.position = "none")
micro.hulls <- ddply(input_filt_rare$map_loaded, c("Treatment"), find_hull)
g5 <- ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02)) +
  # stat_ellipse(mapping = aes(x = Axis01, y = Axis02, color = Treatment, fill = Treatment), 
  #              geom = "polygon", alpha = 0.1, linewidth = NA) +
  geom_polygon(data = micro.hulls, aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, aes(fill = Treatment, shape = Depth)) +
  geom_point(data = centroids, aes(PCoA1, PCoA2, fill = Treatment),
             size = 5, shape = 23, stroke = 1) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       fill = "Treatment") +
  scale_shape_manual(values = c(21,24,22)) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 23, size = 4))) +
  theme_bw() +  
  theme(legend.position = "right",
        legend.box.margin = margin(0,0,100,0, unit = "pt"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())
plot.with.inset <-
  ggdraw() +
  draw_plot(g5) +
  draw_plot(g4, x = 0.75, y = 0.1, width = 0.3, height = 0.3)
plot.with.inset
pdf("Beta_Fungi.pdf", width = 7, height = 5)
plot.with.inset
dev.off()

# Depth
micro.hulls <- ddply(input_filt_rare$map_loaded, c("Depth"), find_hull)
g5 <- ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, aes(colour = Depth, fill = Depth),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, aes(shape = Depth, fill = Depth)) +
  labs(x = pcoaA1, 
       y = pcoaA2) +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_colour_manual(values = c("#DDAA33", "#BB5566", "#004488")) +
  scale_fill_manual(values = c("#DDAA33", "#BB5566", "#004488")) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(0, 40, 0, 5))
g5

plot.with.inset <-
  ggdraw() +
  draw_plot(g5) +
  draw_plot(g4, x = 0.75, y = 0.1, width = 0.3, height = 0.3)
plot.with.inset
pdf("Beta_Depth_Fungi.pdf", width = 7, height = 5)
plot.with.inset
dev.off()

# Taxa
cliffplot_taxa_bars(input_filt_rare, 6, variable = "sampleID")

weather_tax <- c("Knufia", "Aspergillus", "Talaromyces", "Trichoderma",
                 "Paxillus")
length(weather_tax)
tax_sum_genus <- summarize_taxonomy(input = input_filt_rare, 
                                    level = 6, 
                                    report_higher_tax = F)
wt <- tax_sum_genus %>%
  filter(rownames(.) %in% weather_tax) # 4/5 of those genera present (missing Paxillus)
bars_genus <- plot_taxa_bars(wt,
                             input_filt_rare$map_loaded,
                             "sampleID",
                             num_taxa = 19,
                             data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare$map_loaded, by = c("group_by" = "sampleID"))
top_genus <- bars_genus %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_genus <- bars_genus %>%
  mutate(taxon = factor(taxon,
                        levels = c(rev(top_genus$taxon))))
nb.cols <- 5
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
pdf("Taxa_Weathering_Fungi.pdf", width = 8, height = 6)
ggplot(bars_genus, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative fungal abundance", fill = "Genus") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(0.001, 0.001)) + 
  facet_nested_wrap(~ Treatment + Depth, scales = "free_x", ncol = 12) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 12),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "right",
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = unit(0.3, "cm"))
dev.off()


#### ___10 ####
# Rerun alpha and beta just on 10 cm samples
top <- filter_data(input_filt_rare,
                   filter_cat = "Depth",
                   keep_vals = 10) # 24 samples remaining

# Check beta diversity
bc <- calc_dm(top$data_loaded)
set.seed(1150)
ado1 <- adonis2(bc ~ top$map_loaded$MP1 * top$map_loaded$Basalt, by = "terms")
ado1 # Basalt***, I*
pcoa <- cmdscale(bc, k = nrow(top$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
top$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
top$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
df_pcoa <- data.frame(SampleID = rownames(pcoa$points),
                      PCoA1 = pcoa$points[, 1],
                      PCoA2 = pcoa$points[, 2],
                      Treatment = top$map_loaded$Treatment)
centroids <- aggregate(cbind(PCoA1, PCoA2) ~ Treatment, data = df_pcoa, FUN = mean)
eta_sq_m3 <- eta_sq_adonis(ado1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("MP1", "Basalt***", "I*", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g4 <- ggplot(eta_sq_m3, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA, alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m3,
                  aes(y = ypos, label = group), color = "black", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#31688EFF", "#FDE725FF", "#35B779FF", "grey70")) +
  theme_void() + 
  theme(legend.position = "none")
micro.hulls <- ddply(top$map_loaded, c("Treatment"), find_hull)
g5 <- ggplot(top$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, shape = 21, aes(fill = Treatment)) +
  geom_point(data = centroids, aes(PCoA1, PCoA2, fill = Treatment),
             size = 5, shape = 23, stroke = 1) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       fill = "Treatment") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 23, size = 4))) +
  theme_bw() +  
  theme(legend.position = "right",
        legend.box.margin = margin(0,0,100,0, unit = "pt"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())
plot.with.inset <-
  ggdraw() +
  draw_plot(g5) +
  draw_plot(g4, x = 0.73, y = 0.1, width = 0.33, height = 0.33)
plot.with.inset
pdf("Beta_10_Fungi.pdf", width = 7, height = 5)
plot.with.inset
dev.off()



#### ___20 ####
mid <- filter_data(input_filt_rare,
                   filter_cat = "Depth",
                   keep_vals = 20) # 23 samples remaining
bc <- calc_dm(mid$data_loaded)
set.seed(1150)
ado1 <- adonis2(bc ~ mid$map_loaded$MP1 * mid$map_loaded$Basalt, by = "terms")
ado1 # No interaction
set.seed(1150)
ado1 <- adonis2(bc ~ mid$map_loaded$MP1 + mid$map_loaded$Basalt, by = "terms")
ado1
pcoa <- cmdscale(bc, k = nrow(mid$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
mid$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
mid$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
df_pcoa <- data.frame(SampleID = rownames(pcoa$points),
                      PCoA1 = pcoa$points[, 1],
                      PCoA2 = pcoa$points[, 2],
                      Treatment = mid$map_loaded$Treatment)
centroids <- aggregate(cbind(PCoA1, PCoA2) ~ Treatment, data = df_pcoa, FUN = mean)
eta_sq_m3 <- eta_sq_adonis(ado1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("MP1", "Basalt", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g4 <- ggplot(eta_sq_m3, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA, alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m3,
                  aes(y = ypos, label = group), color = "black", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#31688EFF", "#35B779FF", "grey70")) +
  theme_void() + 
  theme(legend.position = "none")
micro.hulls <- ddply(mid$map_loaded, c("Treatment"), find_hull)
g5 <- ggplot(mid$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, shape = 24, aes(fill = Treatment)) +
  geom_point(data = centroids, aes(PCoA1, PCoA2, fill = Treatment),
             size = 5, shape = 23, stroke = 1) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       fill = "Treatment") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 23, size = 4))) +
  theme_bw() +  
  theme(legend.position = "right",
        legend.box.margin = margin(0,0,100,0, unit = "pt"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())
plot.with.inset <-
  ggdraw() +
  draw_plot(g5) +
  draw_plot(g4, x = 0.75, y = 0.1, width = 0.3, height = 0.3)
plot.with.inset
pdf("Beta_20_Fungi.pdf", width = 7, height = 5)
plot.with.inset
dev.off()



#### ___30 ####
bot <- filter_data(input_filt_rare,
                   filter_cat = "Depth",
                   keep_vals = 30) # 23 samples remaining
bc <- calc_dm(bot$data_loaded)
set.seed(1150)
ado1 <- adonis2(bc ~ bot$map_loaded$MP1 * bot$map_loaded$Basalt, by = "terms")
ado1 # No interaction
set.seed(1150)
ado1 <- adonis2(bc ~ bot$map_loaded$MP1 + bot$map_loaded$Basalt, by = "terms")
ado1
pcoa <- cmdscale(bc, k = nrow(bot$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
bot$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
bot$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]
df_pcoa <- data.frame(SampleID = rownames(pcoa$points),
                      PCoA1 = pcoa$points[, 1],
                      PCoA2 = pcoa$points[, 2],
                      Treatment = bot$map_loaded$Treatment)
centroids <- aggregate(cbind(PCoA1, PCoA2) ~ Treatment, data = df_pcoa, FUN = mean)
eta_sq_m3 <- eta_sq_adonis(ado1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("MP1", "Basalt*", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g4 <- ggplot(eta_sq_m3, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA, alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m3,
                  aes(y = ypos, label = group), color = "black", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#31688EFF", "#35B779FF", "grey70")) +
  theme_void() + 
  theme(legend.position = "none")
micro.hulls <- ddply(bot$map_loaded, c("Treatment"), find_hull)
g5 <- ggplot(bot$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, shape = 22, aes(fill = Treatment)) +
  geom_point(data = centroids, aes(PCoA1, PCoA2, fill = Treatment),
             size = 5, shape = 23, stroke = 1) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       fill = "Treatment") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 23, size = 4))) +
  theme_bw() +  
  theme(legend.position = "right",
        legend.box.margin = margin(0,0,100,0, unit = "pt"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())
plot.with.inset <-
  ggdraw() +
  draw_plot(g5) +
  draw_plot(g4, x = 0.75, y = 0.1, width = 0.3, height = 0.3)
plot.with.inset
pdf("Beta_30_Fungi.pdf", width = 7, height = 5)
plot.with.inset
dev.off()



#### ___DA ####
# Need to do differential abundance analysis
# Many methods - SIMPER, Multipatt, ANCOM, DESeq2, Aldex2, LFse

# SIMPER (simple, give top drivers of variation)
sim <- simper(t(top$data_loaded), top$map_loaded$Treatment)
sim_sum <- summary(sim)
sim_df1 <- head(sim_sum$`+MP1_Control`, n = 10) %>%
  mutate(Comparison = "Control_+MP1",
         OTU = rownames(.)) %>%
  rename("MeanControl" = avb,
         "MeanTrt" = ava)
sim_df2 <- head(sim_sum$`+Basalt_Control`, n = 10) %>%
  mutate(Comparison = "Control_+Basalt",
         OTU = rownames(.)) %>%
  rename("MeanControl" = avb,
         "MeanTrt" = ava)
sim_df3 <- head(sim_sum$`+MP1+Basalt_Control`, n = 10) %>%
  mutate(Comparison = "Control_+MP1+Basalt",
         OTU = rownames(.)) %>%
  rename("MeanControl" = avb,
         "MeanTrt" = ava)
simper_results <- rbind(sim_df1, sim_df2, sim_df3) %>%
  left_join(., top$taxonomy_loaded, by = c("OTU" = "taxonomy7")) %>%
  dplyr::select(-c(taxonomy1, taxonomy2, taxonomy3, taxonomy4, taxonomy5)) %>%
  mutate(Response = ifelse(MeanTrt > MeanControl, "Positive", "Negative")) %>%
  rename(c("Genus" = "taxonomy6",
           "CumulativeContribution" = "cumsum")) %>%
  mutate(MeanTrt = round((MeanTrt/184666)*100, digits = 2),
         MeanControl = round((MeanControl/184666)*100, digits = 2),
         CumulativeContribution = round(CumulativeContribution*100, digits = 2)) %>%
  mutate(Taxonomy = make.unique(Genus)) %>%
  dplyr::select(Comparison, Taxonomy, Response, MeanTrt,
                MeanControl, CumulativeContribution)
simper_meta <- simper_results %>%
  column_to_rownames(var = "Taxonomy") %>%
  dplyr::select(Response, Comparison)
simper_mat <- simper_results %>%
  column_to_rownames(var = "Taxonomy") %>%
  dplyr::select(MeanControl, MeanTrt, CumulativeContribution) %>%
  as.matrix()
ann_rows <- data.frame(row.names = rownames(simper_meta),
                       "Response" = simper_meta$Response,
                       "Comparison" = simper_meta$Comparison)
ann_colors <- list(Response = c(Positive = "red",
                                Negative = "blue"),
                   Comparison = c(`Control_+MP1` = "#31688EFF", 
                                  `Control_+Basalt` = "#35B779FF", 
                                  `Control_+MP1+Basalt` = "#FDE725FF"))

pheatmap(simper_mat,
         legend = T,
         legend_breaks = c(0.38, 5, 10, 15, max(na.omit(simper_mat))),
         legend_labels = c("0   ", "5   ", "10   ", "15   ", "%\n"),
         main = "",
         border_color = NA,
         scale = "none",
         angle_col = 315,
         fontsize = 6,
         fontsize_row = 6,
         fontsize_col = 8,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         gaps_row = c(10, 20),
         filename = "Simper_Fungi.png",
         width = 5,
         height = 5)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# ANCOMBC (use input_filt, BC for small sample size)
input_filt$map_loaded$Treatment <- factor(input_filt$map_loaded$Treatment,
                                          levels = c("Control",
                                                     "+Basalt",
                                                     "+MP1",
                                                     "+MP1+Basalt"))
otu <- phyloseq::otu_table(input_filt$data_loaded, taxa_are_rows = T)
tax <- phyloseq::tax_table(as.matrix(input_filt$taxonomy_loaded))
map <- phyloseq::sample_data(input_filt$map_loaded)
input.phy <- phyloseq::phyloseq(otu, tax, map)
input.phy
levels(input.phy@sam_data$Treatment)
set.seed(123)
out <- ancombc2(data = input.phy,
               assay_name = "counts",
               tax_level = NULL,
               fix_formula = "Treatment",
               rand_formula = NULL,
               p_adj_method = "holm",
               pseudo = 0,
               pseudo_sens = FALSE,
               prv_cut = 0.10,
               lib_cut = 0,
               s0_perc = 0.05,
               group = "Treatment", 
               struc_zero = TRUE,
               neg_lb = TRUE, 
               alpha = 0.05, 
               n_cl = 1, 
               verbose = TRUE)
res <- out$res
# None significant in any treatment, not even raw p-values.



#### ___Class ####
fun_cla <- read.delim("kaiju_fungi_class.tsv") %>%
  filter(taxon_name != "unclassified") %>%
  filter(taxon_name != "Viruses") %>%
  filter(taxon_name != "cannot be assigned to a (non-viral) class") %>%
  separate(file, remove = T, sep = "_", into = c("sampleID", "Junk")) %>%
  mutate(sampleID = gsub("-", ".", sampleID)) %>%
  dplyr::select(-taxon_id, -Junk, -percent) %>%
  pivot_wider(names_from = sampleID, values_from = reads) %>%
  dplyr::rename(taxonomy = taxon_name) %>%
  dplyr::select(-Ext.Blank.1, -Ext.Blank.2, -Ext.Blank.3) %>%
  mutate(OTU_ID = paste0("OTU_", row_number())) %>%
  mutate(C.R1.20 = 0) %>%
  mutate(B.R1.30 = 0) %>%
  dplyr::select(OTU_ID, all_of(metadata$Sample), taxonomy) %>%
  as.data.frame()
sum(names(fun_cla) %in% metadata$Sample)
out_fp <- "seqtab_wTax_mctoolsr_kaijuCla.txt"
names(fun_cla)[1] = "#OTU_ID"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(fun_cla,
                             out_fp,
                             sep = "\t",
                             row.names = FALSE,
                             append = TRUE))

tax_table_fp <- "seqtab_wTax_mctoolsr_kaijuCla.txt"
map_fp <- "metadata.txt"
input = load_taxa_table(tax_table_fp, map_fp) # 72 samples loaded
input = list()
input$map_loaded <- metadata
rownames(input$map_loaded) <- input$map_loaded$Sample
input$data_loaded <- as.data.frame(fun_cla[,2:73])
rownames(input$data_loaded) <- fun_cla$OTU_ID
input$taxonomy_loaded <- data.frame(taxonomy1 = "Fungi",
                                    taxonomy2 = "Phylum",
                                    taxonomy3 = "Class",
                                    taxonomy4 = "Order",
                                    taxonomy5 = "Family",
                                    taxonomy6 = fun_cla$taxonomy,
                                    taxonomy7 = fun_cla$OTU_ID)
rownames(input$taxonomy_loaded) <- fun_cla$OTU_ID

# Remove singletons and doubletons - none
singdoub <- data.frame("count" = rowSums(input$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.))

# Remove C.R1.20 and B.R1.30, very few reads, (already knew this, few raw reads too)
input_filt <- input

# Rarefaction
rarecurve(t(input_filt$data_loaded), step = 500, label = F)
sort(colSums(input_filt$data_loaded))
mean(colSums(input_filt$data_loaded)) # 290832.3
se(colSums(input_filt$data_loaded)) # 6996.947
set.seed(530)
input_filt_rare <- single_rarefy(input_filt, 203412) # n = 70 still
sort(colSums(input_filt_rare$data_loaded))

# Add rarefied richness and Shannon
input_filt_rare$map_loaded$rich <- specnumber(input_filt_rare$data_loaded, MARGIN = 2)
input_filt_rare$map_loaded$shannon <- diversity(input_filt_rare$data_loaded, 
                                                index = "shannon", MARGIN = 2)
input_filt_rare$map_loaded$rich # All 32 (artifact...)
input_filt_rare$map_loaded$shannon

# Save
#saveRDS(input_filt_rare, "input_filt_rare_kaijuCla.rds")

input_filt_rare <- readRDS("input_filt_rare_kaijuCla.rds")
input_filt_rare$map_loaded$Treatment <- factor(input_filt_rare$map_loaded$Treatment,
                                               levels = c("Control",
                                                          "+Basalt",
                                                          "+MP1",
                                                          "+MP1+Basalt"))
input_filt_rare$map_loaded$Depth <- as.factor(input_filt_rare$map_loaded$Depth)
input_filt_rare$map_loaded$TrtDep <- factor(paste(input_filt_rare$map_loaded$Treatment,
                                                  input_filt_rare$map_loaded$Depth,
                                                  sep = "_"),
                                            levels = c("Control_10", "Control_20", "Control_30",
                                                       "+Basalt_10", "+Basalt_20", "+Basalt_30",
                                                       "+MP1_10", "+MP1_20", "+MP1_30",
                                                       "+MP1+Basalt_10", "+MP1+Basalt_20",
                                                       "+MP1+Basalt_30"))

# Check beta diversity
bc <- calc_dm(input_filt_rare$data_loaded)
set.seed(1150)
ado1 <- adonis2(bc ~ input_filt_rare$map_loaded$MP1 +
                  input_filt_rare$map_loaded$Basalt +
                  input_filt_rare$map_loaded$MP1:input_filt_rare$map_loaded$Basalt +
                  input_filt_rare$map_loaded$Depth,
                by = "terms")
ado1
pcoa <- cmdscale(bc, k = nrow(input_filt_rare$map_loaded) - 1, eig = T)
pcoaA1 <- paste("PC1: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, 1), "%")
pcoaA2 <- paste("PC2: ", round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, 1), "%")
input_filt_rare$map_loaded$Axis01 <- vegan::scores(pcoa)[,1]
input_filt_rare$map_loaded$Axis02 <- vegan::scores(pcoa)[,2]

df_pcoa <- data.frame(SampleID = rownames(pcoa$points),
                      PCoA1 = pcoa$points[, 1],
                      PCoA2 = pcoa$points[, 2],
                      Treatment = input_filt_rare$map_loaded$Treatment)
centroids <- aggregate(cbind(PCoA1, PCoA2) ~ Treatment, data = df_pcoa, FUN = mean)

eta_sq_m3 <- eta_sq_adonis(ado1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("MP1", "Basalt***", "Depth***", "I*", "R")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g4 <- ggplot(eta_sq_m3, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA, alpha = 0.5) +
  coord_polar("y", start = 0) +
  geom_text_repel(data = eta_sq_m3,
                  aes(y = ypos, label = group), color = "black", size = 3,
                  box.padding = 0.1) +
  scale_fill_manual(values = c("#31688EFF", "grey30", "#FDE725FF", "#35B779FF", "grey70")) +
  theme_void() + 
  theme(legend.position = "none")
g4

micro.hulls <- ddply(input_filt_rare$map_loaded, c("Treatment"), find_hull)
g5 <- ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02)) +
  # stat_ellipse(mapping = aes(x = Axis01, y = Axis02, color = Treatment, fill = Treatment), 
  #              geom = "polygon", alpha = 0.1, linewidth = NA) +
  geom_polygon(data = micro.hulls, aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F, linewidth = NA) +
  geom_point(size = 3, alpha = 1, aes(fill = Treatment, shape = Depth)) +
  geom_point(data = centroids, aes(PCoA1, PCoA2, fill = Treatment),
             size = 5, shape = 23, stroke = 1) +
  labs(x = pcoaA1, 
       y = pcoaA2,
       fill = "Treatment") +
  scale_shape_manual(values = c(21,24,22)) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, shape = 23, size = 4))) +
  theme_bw() +  
  theme(legend.position = "right",
        legend.box.margin = margin(0,0,100,0, unit = "pt"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        panel.grid = element_blank())
g5

# Combine
plot.with.inset <-
  ggdraw() +
  draw_plot(g5) +
  draw_plot(g4, x = 0.75, y = 0.1, width = 0.3, height = 0.3)
plot.with.inset
pdf("Beta_Fungi_Class.pdf", width = 7, height = 5)
plot.with.inset
dev.off()

cliffplot_taxa_bars(input_filt_rare, 6, variable = "sampleID")



#### _Sylph ####

#### Lane 1 ####
#### __All ####
d_sort <- input_filt_rare$map_loaded %>%
  arrange(sampleID)
d_sort$sampleID
bacGT <- read.delim("~/Desktop/Fierer/Strains/bac120_metadata_r220.tsv")
arcGT <- read.delim("~/Desktop/Fierer/Strains/ar53_metadata_r220.tsv")
fullGT <- rbind(bacGT, arcGT)
strains_sylph_raw <- read.delim("sylph_profile_andes.tsv") %>%
  mutate(GenomeID = gsub("gtdb_genomes_reps_r220/database/", "", Genome_file)) %>%
  mutate(GenomeID = gsub("_genomic.fna.gz", "", GenomeID)) %>%
  mutate(GenomeID = ifelse(nchar(GenomeID) > 9,
                           substr(GenomeID, start = 17, stop = nchar(GenomeID)),
                           GenomeID))
strains_sylph <- read.delim("sylph_profile_andes.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/AndesAg/02_FiltData/", "", Sample_file)) %>%
  mutate(SampleID = substr(SampleID, start = 1, stop = 7)) %>%
  mutate(SampleID = gsub("-", ".", SampleID)) %>%
  mutate(GenomeID = gsub("gtdb_genomes_reps_r220/database/", "", Genome_file)) %>%
  mutate(GenomeID = gsub("_genomic.fna.gz", "", GenomeID)) %>%
  mutate(GenomeID = ifelse(nchar(GenomeID) > 9,
                           substr(GenomeID, start = 17, stop = nchar(GenomeID)),
                           GenomeID)) %>%
  filter(SampleID %in% d_sort$sampleID)
length(unique(strains_sylph$SampleID)) # 70 samples
length(unique(strains_sylph$GenomeID)) # 200/113,104 + 1

sylph_strains_70 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  pivot_wider(names_from = SampleID, values_from = Sequence_abundance) %>%
  mutate(Type = ifelse(GenomeID == "MP1.fasta", "MP1", "GTDB")) %>%
  column_to_rownames(var = "GenomeID") %>%
  dplyr::select(all_of(d_sort$sampleID)) %>%
  replace(is.na(.), 0)

gtdb_sylph <- fullGT %>%
  filter(grepl(paste(rownames(sylph_strains_70), collapse="|"), accession)) %>% # 181
  mutate(ID = substr(accession, start = 4, stop = nchar(accession))) %>%
  mutate(gtdb_taxonomy = gsub("d__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("p__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("c__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("o__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("f__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("g__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("s__", "", gtdb_taxonomy)) %>%
  separate(gtdb_taxonomy, sep = ";", into = c("Domain", "Phylum", "Class", "Order",
                                              "Family", "Genus", "Species")) %>%
  dplyr::select(accession, ID, checkm2_completeness, checkm2_contamination,
                Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  add_row(accession = "MP1", ID = "MP1.fasta",
          checkm2_completeness = 100, checkm2_contamination = 0.01, Domain = "Bacteria",
          Phylum = "Bacillota", Class = "Bacilli", Order = "Bacillales",
          Family = "Bacillaceae", Genus = "Bacillus", Species = "Bacillus MP1") %>%
  arrange(Species)

sylph_strains_70 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  pivot_wider(names_from = SampleID, values_from = Sequence_abundance) %>%
  mutate(Type = ifelse(GenomeID == "MP1.fasta", "MP1", "GTDB")) %>%
  left_join(., gtdb_sylph, by = c("GenomeID" = "ID")) %>%
  arrange(Species) %>%
  column_to_rownames(var = "Species") %>%
  dplyr::select(all_of(d_sort$sampleID)) %>%
  replace(is.na(.), 0)

ann_cols <- data.frame(row.names = colnames(sylph_strains_70), 
                       Depth = d_sort$Depth,
                       Treatment = d_sort$Treatment)
ann_rows <- data.frame(row.names = rownames(sylph_strains_70)) %>%
  mutate(Type = ifelse(rownames(.) == "Bacillus MP1", "MP1", "GTDB"),
         Domain = gtdb_sylph$Domain,
         #Completeness = gtdb_sylph$checkm2_completeness,
         #Contamination = gtdb_sylph$checkm2_contamination
         )
table(ann_rows$Type)
ann_colors <- list(Type = c("GTDB" = "grey90",
                              "MP1" = "red"),
                   Domain = c("Archaea" = "#F8766D",
                              "Bacteria" = "#619CFF"),
                   Treatment = c("Control" = "#440154FF",
                                 "+Basalt" = "#31688EFF",
                                 "+MP1" = "#35B779FF",
                                 "+MP1+Basalt" = "#FDE725FF"),
                   Depth = c("10" = "grey80",
                             "20" = "grey50",
                             "30" = "grey20"))
max(sylph_strains_70)
breaks <- c(0, seq(0.0466, 27.4631, length.out = 100))
colors <- c("white", colorRampPalette(brewer.pal(9, "Reds"))(99))
pheatmap(sylph_strains_70,
         color = colors,
         breaks = breaks,
         legend = T,
         cluster_rows = T,
         cluster_cols = T,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         fontsize_row = 4,
         fontsize_col = 4,
         na_col = "white",
         border_color = "white",
         filename = "Sylph.png",
         width = 10,
         height = 12)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### __MP1 ####
mp1 <- as.data.frame(t(sylph_strains_70)) %>%
  dplyr::select(`Bacillus MP1`) %>%
  rownames_to_column(var = "sampleID")
sum(d_sort$sampleID != mp1$sampleID)
d_sort$MP1_sylph <- mp1$`Bacillus MP1`

leveneTest(d_sort$MP1_sylph ~ d_sort$Treatment) # Homogeneous
leveneTest(d_sort$MP1_sylph ~ d_sort$Depth) # Homogeneous
m <- aov(MP1_sylph ~ MP1 + Basalt + MP1:Basalt + Depth, data = d_sort)
summary(m)
Anova(m, type = "III", singular.ok = TRUE) # Depth
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
pdf("MP1_sylph.pdf", width = 7, height = 5)
ggplot(d_sort, aes(Treatment, MP1_sylph, fill = Treatment)) +
  geom_boxplot(outliers = F) +
  geom_jitter(size = 3, alpha = 1, width = 0.2, height = 0, aes(shape = Depth)) +
  labs(x = "Treatment", y = "MP1 % abundance (Sylph)") +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24, 22)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)),
         shape = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 14))
dev.off()
colSums(sylph_strains_70)



#### All Lanes ####
# Reran sylph on all lanes concatenated
strains_sylph_raw <- read.delim("sylph_profile_andes_comb.tsv") %>%
  mutate(GenomeID = gsub("gtdb_genomes_reps_r220/database/", "", Genome_file)) %>%
  mutate(GenomeID = gsub("_genomic.fna.gz", "", GenomeID)) %>%
  mutate(GenomeID = ifelse(nchar(GenomeID) > 9,
                           substr(GenomeID, start = 17, stop = nchar(GenomeID)),
                           GenomeID))
strains_sylph <- read.delim("sylph_profile_andes_comb.tsv") %>%
  mutate(SampleID = gsub("/scratch/alpine/clbd1748/AndesAg/Comb/", "", Sample_file)) %>%
  mutate(SampleID = substr(SampleID, start = 1, stop = 7)) %>%
  mutate(SampleID = gsub("-", ".", SampleID)) %>%
  mutate(GenomeID = gsub("gtdb_genomes_reps_r220/database/", "", Genome_file)) %>%
  mutate(GenomeID = gsub("_genomic.fna.gz", "", GenomeID)) %>%
  mutate(GenomeID = ifelse(nchar(GenomeID) > 9,
                           substr(GenomeID, start = 17, stop = nchar(GenomeID)),
                           GenomeID)) %>%
  filter(SampleID %in% d_sort$sampleID)
length(unique(strains_sylph$SampleID)) # 70 samples
length(unique(strains_sylph$GenomeID)) # 398/113,104 + 1
# Note: Using all data detected 398 genomes instead of 200 using Lane 1 only

sylph_strains_70 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  pivot_wider(names_from = SampleID, values_from = Sequence_abundance) %>%
  mutate(Type = ifelse(GenomeID == "MP1.fasta", "MP1", "GTDB")) %>%
  column_to_rownames(var = "GenomeID") %>%
  dplyr::select(all_of(d_sort$sampleID)) %>%
  replace(is.na(.), 0)

ubiq <- sylph_strains_70 %>%
  mutate(Ubiquity = rowSums(across(everything(), ~ . > 0))) %>%
  mutate(MeanAbund = rowMeans(across(everything()))) %>%
  rownames_to_column(var = "GenomeID") %>%
  dplyr::select(GenomeID, Ubiquity, MeanAbund) %>%
  arrange(desc(Ubiquity), desc(MeanAbund))
sum(ubiq$Ubiquity == 70)

df_10 <- sylph_strains_70 %>% select(contains("10"))
df_20 <- sylph_strains_70 %>% select(contains("20"))
df_30 <- sylph_strains_70 %>% select(contains("30"))
ubiq_10 <- df_10 %>%
  mutate(Ubiquity10 = rowSums(across(everything(), ~ . > 0))) %>%
  mutate(MeanAbund10 = rowMeans(across(everything()))) %>%
  rownames_to_column(var = "GenomeID") %>%
  dplyr::select(GenomeID, Ubiquity10, MeanAbund10)
ubiq_20 <- df_20 %>%
  mutate(Ubiquity20 = rowSums(across(everything(), ~ . > 0))) %>%
  mutate(MeanAbund20 = rowMeans(across(everything()))) %>%
  rownames_to_column(var = "GenomeID") %>%
  dplyr::select(GenomeID, Ubiquity20, MeanAbund20)
ubiq_30 <- df_30 %>%
  mutate(Ubiquity30 = rowSums(across(everything(), ~ . > 0))) %>%
  mutate(MeanAbund30 = rowMeans(across(everything()))) %>%
  rownames_to_column(var = "GenomeID") %>%
  dplyr::select(GenomeID, Ubiquity30, MeanAbund30)
ubiq_depth <- ubiq_10 %>%
  left_join(., ubiq_20, by = "GenomeID") %>%
  left_join(., ubiq_30, by = "GenomeID") %>%
  dplyr::select(GenomeID, Ubiquity10, Ubiquity20, Ubiquity30) %>%
  mutate(PresAll = Ubiquity10 > 0 & Ubiquity20 > 0 & Ubiquity30 > 0) %>%
  dplyr::select(GenomeID, Ubiquity10, Ubiquity20, Ubiquity30, PresAll)
sum(ubiq_depth$PresAll == TRUE)

gtdb_sylph <- fullGT %>%
  filter(grepl(paste(rownames(sylph_strains_70), collapse="|"), accession)) %>% # 181
  mutate(ID = substr(accession, start = 4, stop = nchar(accession))) %>%
  mutate(gtdb_taxonomy = gsub("d__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("p__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("c__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("o__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("f__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("g__", "", gtdb_taxonomy)) %>%
  mutate(gtdb_taxonomy = gsub("s__", "", gtdb_taxonomy)) %>%
  separate(gtdb_taxonomy, sep = ";", into = c("Domain", "Phylum", "Class", "Order",
                                              "Family", "Genus", "Species")) %>%
  dplyr::select(accession, ID, checkm2_completeness, checkm2_contamination,
                Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  add_row(accession = "MP1", ID = "MP1.fasta",
          checkm2_completeness = 100, checkm2_contamination = 0.01, Domain = "Bacteria",
          Phylum = "Bacillota", Class = "Bacilli", Order = "Bacillales",
          Family = "Bacillaceae", Genus = "Bacillus", Species = "Bacillus MP1") %>%
  arrange(Species)

sylph_strains_70 <- strains_sylph %>%
  dplyr::select(SampleID, GenomeID, Sequence_abundance) %>%
  pivot_wider(names_from = SampleID, values_from = Sequence_abundance) %>%
  mutate(Type = ifelse(GenomeID == "MP1.fasta", "MP1", "GTDB")) %>%
  left_join(., gtdb_sylph, by = c("GenomeID" = "ID")) %>%
  arrange(Species) %>%
  column_to_rownames(var = "Species") %>%
  dplyr::select(all_of(d_sort$sampleID)) %>%
  replace(is.na(.), 0)

ann_cols <- data.frame(row.names = colnames(sylph_strains_70), 
                       Depth = d_sort$Depth,
                       Treatment = d_sort$Treatment)
ann_rows <- data.frame(row.names = rownames(sylph_strains_70)) %>%
  mutate(Type = ifelse(rownames(.) == "Bacillus MP1", "MP1", "GTDB"),
         Domain = gtdb_sylph$Domain,
         #Completeness = gtdb_sylph$checkm2_completeness,
         #Contamination = gtdb_sylph$checkm2_contamination
  )
table(ann_rows$Type)
ann_colors <- list(Type = c("GTDB" = "grey90",
                            "MP1" = "red"),
                   Domain = c("Archaea" = "#F8766D",
                              "Bacteria" = "#619CFF"),
                   Treatment = c("Control" = "#440154FF",
                                 "+Basalt" = "#31688EFF",
                                 "+MP1" = "#35B779FF",
                                 "+MP1+Basalt" = "#FDE725FF"),
                   Depth = c("10" = "grey80",
                             "20" = "grey50",
                             "30" = "grey20"))
max(sylph_strains_70)
min(sylph_strains_70[sylph_strains_70 > 0], na.rm = TRUE)
breaks <- c(0, seq(0.0074, 21.9267, length.out = 100))
colors <- c("white", colorRampPalette(brewer.pal(9, "Reds"))(99))
pheatmap(sylph_strains_70,
         color = colors,
         breaks = breaks,
         legend = T,
         cluster_rows = T,
         cluster_cols = T,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         fontsize_row = 3,
         fontsize_col = 4,
         na_col = "white",
         border_color = "white",
         filename = "Sylph_Comb.png",
         width = 10,
         height = 13)
dev.off()
dev.set(dev.next())
dev.set(dev.next())


# Plot heatmap just with MP1 and 4 Bacillus
gtdb_sylph_bac <- gtdb_sylph %>%
  filter(Genus == "Bacillus")
sylph_strains_70_Bac <- sylph_strains_70 %>%
  filter(rownames(.) %in% gtdb_sylph_bac$Species)
ann_cols <- data.frame(row.names = colnames(sylph_strains_70_Bac), 
                       Depth = d_sort$Depth,
                       Treatment = d_sort$Treatment)
ann_rows <- data.frame(row.names = rownames(sylph_strains_70_Bac)) %>%
  mutate(Type = ifelse(rownames(.) == "Bacillus MP1", "MP1", "GTDB"),
         #Domain = gtdb_sylph_bac$Domain,
         #Completeness = gtdb_sylph$checkm2_completeness,
         #Contamination = gtdb_sylph$checkm2_contamination
  )
table(ann_rows$Type)
ann_colors <- list(Type = c("GTDB" = "grey90",
                            "MP1" = "red"),
                   # Domain = c("Archaea" = "#F8766D",
                   #            "Bacteria" = "#619CFF"),
                   Treatment = c("UTC" = "#440154FF",
                                 "UTC+B" = "#31688EFF",
                                 "MP1" = "#35B779FF",
                                 "MP1+B" = "#FDE725FF"),
                   Depth = c("10" = "grey80",
                             "20" = "grey50",
                             "30" = "grey20"))
max(sylph_strains_70_Bac)
min(sylph_strains_70_Bac[sylph_strains_70_Bac > 0], na.rm = TRUE)
breaks <- c(0, seq(0.0171, 0.3211, length.out = 100))
colors <- c("white", colorRampPalette(brewer.pal(9, "Reds"))(99))
pheatmap(sylph_strains_70_Bac,
         color = colors,
         breaks = breaks,
         legend = T,
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 90,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         fontsize_row = 10,
         fontsize_col = 4,
         na_col = "white",
         border_color = "white",
         show_colnames = F,
         filename = "FinalFigs/FigureS5.png",
         width = 10,
         height = 8)
dev.off()
dev.set(dev.next())
dev.set(dev.next())


d_sort$Depth <- as.factor(d_sort$Depth)
mp1 <- as.data.frame(t(sylph_strains_70)) %>%
  dplyr::select(`Bacillus MP1`) %>%
  rownames_to_column(var = "sampleID")
sum(d_sort$sampleID != mp1$sampleID)
d_sort$MP1_sylph <- mp1$`Bacillus MP1`
leveneTest(d_sort$MP1_sylph ~ d_sort$Treatment) # Homogeneous
leveneTest(d_sort$MP1_sylph ~ d_sort$Depth) # Homogeneous
m <- aov(MP1_sylph ~ MP1 + Basalt + MP1:Basalt + Depth, data = d_sort)
summary(m)
Anova(m, type = "III", singular.ok = TRUE) # Depth
hist(m$residuals)
shapiro.test(m$residuals)
plot(m$fitted.values, m$residuals)
mp1.sylph <- ggplot(d_sort, aes(Treatment, MP1_sylph, shape = Depth)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment)) +
  geom_point(size = 3, alpha = 1, position = position_jitterdodge(),
             aes(fill = Treatment)) +
  labs(x = NULL, y = "MP1 % abundance (metagenome)", shape = "Depth (cm)") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24, 22)) +
  guides(fill = "none") +
  theme_bw() +
  theme(legend.position = "inside",,
        legend.position.inside = c(0.82, 0.78),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        panel.grid = element_blank())
mp1.sylph
pdf("MP1_sylph_comb.pdf", width = 7, height = 5)
mp1.sylph
dev.off()
png("MP1_sylph_comb.png", width = 7, height = 5, units = "in", res = 300)
mp1.sylph
dev.off()
colSums(sylph_strains_70)
# Wow, now MP1 is detected in all 6 surface samples from +MP1 and +MP1+Basalt
# This shows that sequencing depth does affect Sylph's ability to detect things

# qPCR
qPCR <- read_xlsx("phili, SYBR Green, 08-13-2025, 07Hr 47Min - Text Report.xlsx",
                  sheet = 2) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("UTC", "MP1", "UTC+B", "MP1+B")))
mp1.qPCR <- ggplot(qPCR, aes(Treatment, CFU)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment)) +
  geom_jitter(size = 3, alpha = 1, width = 0.2, height = 0, 
              aes(fill = Treatment), colour = "black", pch = 21) +
  labs(x = NULL, y = "CFU/g soil") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(fill = "none") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid = element_blank())
mp1.qPCR

# leg_s <- get_legend(mp1.sylph)
# mp1.sylph <- mp1.sylph +
#   theme(legend.position = "none")
pdf("FinalFigs/Figure5.pdf", width = 7, height = 3)
plot_grid(mp1.sylph, mp1.qPCR, labels = "AUTO")
dev.off()
png("FinalFigs/Figure5.png", width = 7, height = 3, units = "in", res = 300)
plot_grid(mp1.sylph, mp1.qPCR, labels = "AUTO")
dev.off()



#### 3. Function ####
# Iron cycling genes, IRcyc-A (Epihov and Bryce 2024) or FeGenie
# Carbonic anhydrases (Vienne et al. 2022, Vicca et al. 2022)
# Urease (Vicca)
# See Napieralksi et al. 2019 for a metagenomic study



#### _Import ####
# Import anvi'o functional profile
f1 <- read.delim("~/Desktop/Fierer/AndesAg/functions/A_R1_10_functions.txt")
length(unique(f1$gene_callers_id)) # 321571 gene IDs in the contigs.db
f1_ko <- f1 %>%
  filter(source == "KOfam") # 70405
length(unique(f1_ko$gene_callers_id)) # 69015
70405 - 69015 # 1390 duplicate cases across 1287 gene ids
length(unique(f1_ko$accession)) # 4511 KOs
f1_ko_dup <- f1_ko %>%
  group_by(gene_callers_id) %>%
  summarise(count = n())
f1_count <- f1_ko %>%
  group_by(accession) %>%
  summarise(count = n())

c1 <- read.delim("~/Desktop/Fierer/AndesAg/coverage/A_R1_10-GENE-COVERAGES.txt") # 105835
length(unique(c1$key)) # 459512
c1 <- c1 %>%
  left_join(., f1_ko, by = c("key" = "gene_callers_id")) # 460902
460902-459512 # 1390 - matches duplicate cases

c1_ko <- c1 %>%
  filter(source == "KOfam") # 70405! Matches contigs.db
length(unique(c1_ko$key)) # 69015
length(unique(c1_ko$accession)) # 4511
c1_ko_dup <- c1_ko %>%
  group_by(key) %>%
  summarise(count = n())

# Need to group by gene ID, choose lowest e value
# Then group by KO and sum the coverages
c1_ko_tab <- c1_ko %>%
  group_by(key) %>%
  slice_min(e_value, with_ties = F) %>%
  ungroup() %>%
  group_by(accession) %>%
  summarise(A_R1_10 = sum(A_R1_10))

# Okay this is looking good. Now need to make 24 of these, then merge by KO, then MUSiCC
# Do not ignore!:
# Warning messages:
#   1: In scan(file = file, what = what, sep = sep, quote = quote,  ... :
#                EOF within quoted string
# Very important to include quote = "" or rows are missed!

fs <- list()
cs <- list()
for (i in 1:24) {
  setwd("~/Desktop/Fierer/AndesAg/functions/")
  f_files <- list.files()
  fs[[i]] <- read.delim(f_files[i], quote = "") %>%
    filter(source == "KOfam")
  setwd("~/Desktop/Fierer/AndesAg/coverage/")
  c_files <- list.files()
  sampleID <- substr(c_files[i], start = 1, stop = 7)
  c_tmp <- read.delim(c_files[i], quote = "")
  names(c_tmp)[2] <- "sample"
  cs[[i]] <- c_tmp %>%
    left_join(., fs[[i]], by = c("key" = "gene_callers_id")) %>%
    filter(source == "KOfam") %>%
    group_by(key) %>%
    slice_min(e_value, with_ties = F) %>%
    ungroup() %>%
    group_by(accession) %>%
    summarise(sample = sum(sample))
  names(cs[[i]]) <- c("KO", sampleID)
}
setwd("~/Desktop/Fierer/AndesAg/")

# Now merge
ko_merged <- cs[[1]]
for (i in 2:24) {
  ko_merged <- ko_merged %>%
    full_join(., cs[[i]], by = "KO") # N.B.! Use full join!
}
ko_merged[is.na(ko_merged)] <- 0
#write.table(ko_merged, "coverage_for_musicc.txt", row.names = F, sep = "\t")




# reads_per_gene <- coverage * gene_length / read_length
l1 <- read.delim("genecalls/A_R1_10_genecalls.txt") %>%
  mutate(geneLength = stop - start) %>%
  dplyr::select(gene_callers_id, geneLength)
length(unique(l1$gene_callers_id)) # 455958

c1_ko_tab <- c1_ko %>%
  left_join(., l1, by = c("key" = "gene_callers_id")) %>%
  mutate(geneReads = A_R1_10*geneLength/150) %>%
  group_by(key) %>%
  slice_min(e_value, with_ties = F) %>%
  ungroup() %>%
  group_by(accession) %>%
  summarise(A_R1_10 = sum(geneReads))

# Loop through all samples
ls <- list()
fs <- list()
cs <- list()
for (i in 1:24) {
  setwd("~/Desktop/Fierer/AndesAg/genecalls/")
  l_files <- list.files()
  ls[[i]] <- read.delim(l_files[i], quote = "") %>%
    mutate(geneLength = stop - start) %>%
    dplyr::select(gene_callers_id, geneLength)

  setwd("~/Desktop/Fierer/AndesAg/functions/")
  f_files <- list.files()
  fs[[i]] <- read.delim(f_files[i], quote = "") %>%
    filter(source == "KOfam")
  
  setwd("~/Desktop/Fierer/AndesAg/coverage/")
  c_files <- list.files()
  sampleID <- substr(c_files[i], start = 1, stop = 7)
  c_tmp <- read.delim(c_files[i], quote = "")
  names(c_tmp)[2] <- "sample"
  cs[[i]] <- c_tmp %>%
    left_join(., fs[[i]], by = c("key" = "gene_callers_id")) %>%
    filter(source == "KOfam") %>%
    left_join(., ls[[i]], by = c("key" = "gene_callers_id")) %>%
    mutate(geneReads = sample*geneLength/150) %>%
    group_by(key) %>%
    slice_min(e_value, with_ties = F) %>%
    ungroup() %>%
    group_by(accession) %>%
    summarise(geneReads = sum(geneReads))
  names(cs[[i]]) <- c("KO", sampleID)
}
setwd("~/Desktop/Fierer/AndesAg/")

# Now merge
ko_merged <- cs[[1]]
for (i in 2:24) {
  ko_merged <- ko_merged %>%
    full_join(., cs[[i]], by = "KO") # N.B.! Use full join!
}
ko_merged[is.na(ko_merged)] <- 0
#write.table(ko_merged, "ko_counts.txt", row.names = F, sep = "\t")



#### _MUSiCC ####
# List of genes of interest sent by Andes Ag team
# Removed duplicate lysT, lysS. moved to pyruvate category
goi <- read_excel("GenesList.xlsx")
goi_ko <- goi %>%
  filter(KO != "NA") %>%
  arrange(Pathway, `Gene name`) # 94
length(unique(goi_ko$`Gene name`)) # 94
table(goi_ko$Pathway)

# MUSiCC inter-sample corrected coverages
#ko_musicc <- read.delim("ko_coverage_musicc.tab")
# Run with -n and -c flags
#ko_musicc <- read.delim("ko_coverage_musicc_nc.tab")
# Run on counts instead of coverage. with -n and -c (use this)
ko_musicc <- read.delim("ko_count_musicc_nc.tab")

# MUSiCC goi
ko_musicc_goi <- ko_musicc %>%
  filter(KO %in% goi_ko$KO) %>% 
  column_to_rownames(var = "KO") %>%
  t() %>%
  as.data.frame() # Only 71 of 94 present!
goi_ko_present <- goi_ko %>%
  filter(KO %in% ko_musicc$KO)

# Which ones are absent?
goi_ko_absent <- goi_ko %>%
  filter(KO %notin% ko_musicc$KO)
  
# At this point, you really just need to run ANOVA on those!
# But first need to also make sure assumptions are met
# Need to loop through those 71 KOs and run tests and store outputs
# Then plot as heatmap
meta10 <- read.delim("metadata.txt") %>%
  filter(Depth == 10) %>%
  mutate(SampleID = gsub("\\.", "_", sampleID))
sum(rownames(ko_musicc_goi) != meta10$SampleID)

stat_loop <- as.data.frame(matrix(NA, 71, 7)) %>%
  set_names("KO", "Levene", "Shapiro", "Basalt", "MP1", "Int", "Trt")
for (i in 1:ncol(ko_musicc_goi)) {
  # KO name
  stat_loop$KO[i] <- names(ko_musicc_goi)[i]
  
  # Levene Test
  l1 <- leveneTest(ko_musicc_goi[,i] ~ meta10$Treatment)
  stat_loop$Levene[i] <- l1$`Pr(>F)`[1]
  
  # Models
  m <- aov(ko_musicc_goi[,i] ~ meta10$Basalt * meta10$MP1)
  a <- Anova(m, type = "III", singular.ok = TRUE)
  
  t <- aov(ko_musicc_goi[,i] ~ meta10$Treatment)
  ta <- Anova(t)
  
  stat_loop$Basalt[i] <- a$`Pr(>F)`[2]
  stat_loop$MP1[i] <- a$`Pr(>F)`[3]
  stat_loop$Int[i] <- a$`Pr(>F)`[4]
  stat_loop$Trt[i] <- ta$`Pr(>F)`[1]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  stat_loop$Shapiro[i] <- s1$p.value
}

stat_loop <- stat_loop %>%
  mutate(BasaltPadj = p.adjust(Basalt, method = "fdr"),
         MP1Padj = p.adjust(MP1, method = "fdr"),
         IntPadj = p.adjust(Int, method = "fdr"),
         TrtPadj = p.adjust(Trt, method = "fdr"))

# Plot, with annotation column for sample, annotation rows for sig and pathway
ko_musicc_goi_t <- ko_musicc_goi %>%
  dplyr::select(all_of(goi_ko_present$KO)) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(D_R1_10, D_R2_10, D_R3_10, D_R4_10, D_R5_10, D_R6_10,
                B_R1_10, B_R2_10, B_R3_10, B_R4_10, B_R5_10, B_R6_10,
                A_R1_10, A_R2_10, A_R3_10, A_R4_10, A_R5_10, A_R6_10,
                C_R1_10, C_R2_10, C_R3_10, C_R4_10, C_R5_10, C_R6_10)
meta10 <- meta10 %>%
  arrange(SampleID = colnames(ko_musicc_goi_t)) %>%
  arrange(SampleID = colnames(ko_musicc_goi_t))
sum(colnames(ko_musicc_goi_t) != meta10$SampleID)
sum(rownames(ko_musicc_goi_t) != goi_ko_present$KO)
goi_ko_present <- goi_ko_present %>%
  left_join(., stat_loop, by = "KO") %>%
  mutate(Basalt_Sig = ifelse(Basalt < 0.05, "P < 0.05", " "),
         MP1_Sig = ifelse(MP1 < 0.05, "P < 0.05", " "),
         Int_Sig = ifelse(Int < 0.05, "P < 0.05", " ")) %>%
  mutate(label = paste(KO, `Gene name`, sep = " "))
ann_cols <- data.frame(row.names = colnames(ko_musicc_goi_t), 
                       Treatment = meta10$Treatment)
table(ann_cols$Treatment)
ann_rows <- data.frame(row.names = rownames(ko_musicc_goi_t)) %>%
  mutate(`Pathway` = goi_ko_present$Pathway,
         `Basalt` = goi_ko_present$Basalt_Sig,
         `MP1` = goi_ko_present$MP1_Sig,
         `Int` = goi_ko_present$Int_Sig)
table(ann_rows$Pathway)
table(ann_rows$Basalt)
brewer_pal(palette = "Paired")(12)
viridis_pal()(4)
ann_colors <- list(`Pathway` = c("Acetic acid" = "#A6CEE3",
                                 "Ammonium" = "#1F78B4",
                                 "Biofilm" = "#B2DF8A",
                                 "Fatty Acid Metabolism" = "#33A02C",
                                 "Formic Acid" = "#FB9A99",
                                 "Gluconic acid" = "#E31A1C",
                                 "Glycolic acid" = "#FDBF6F",
                                 "Lactic acid" = "#FF7F00",
                                 "Proprionic acid" = "#CAB2D6",
                                 "Pulcherrimin" = "#6A3D9A",
                                 "Pyruvic acid" = "#FFFF99",
                                 "Siderophore" = "#B15928",
                                 "Succinic acid" = "grey90",
                                 "Urea" = "grey40"),
                   `Basalt` = c("P < 0.05" = "#31688EFF",
                                " " = "white"),
                   `MP1` = c("P < 0.05" = "#35B779FF",
                             " " = "white"),
                   `Int` = c("P < 0.05" = "#FDE725FF",
                             " " = "white"),
                   Treatment = c("Control" = "#440154FF",
                                 "+Basalt" = "#31688EFF",
                                 "+MP1" = "#35B779FF",
                                 "+MP1+Basalt" = "#FDE725FF"))
pheatmap(ko_musicc_goi_t,
         legend = T,
         legend_labels = c("-3  ", "-2  ", "-1  ", "0  ", "1  ", "2  ", "3  "),
         cluster_rows = F,
         cluster_cols = F,
         scale = "row",
         angle_col = 315,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         annotation_names_row = T,
         fontsize = 8,
         fontsize_row = 6,
         na_col = "white",
         border_color = "white",
         gaps_col = c(6, 12, 18),
         labels_row = goi_ko_present$label,
         show_colnames = F,
         filename = "KO_Heatmap_MUSiCC.png",
         width = 6,
         height = 8)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# A bit vague to tell direction - subset to the sig ones and plot
sig <- stat_loop %>%
  filter(Basalt < 0.05 | MP1 < 0.05 | Int < 0.05)
ko_musicc_goi_sig <- ko_musicc_goi %>%
  dplyr::select(sig$KO) %>%
  rownames_to_column(var = "SampleID")
meta10 <- meta10 %>%
  left_join(., ko_musicc_goi_sig, by = "SampleID")
gene_plot <- meta10 %>%
  pivot_longer(cols = c(16:(16 + nrow(sig) - 1))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "+Basalt", "+MP1", "+MP1+Basalt"))) %>%
  left_join(., goi_ko_present, by = c("name" = "KO"))
pdf("SigKO_MUSiCC.pdf", width = 8, height = 6)
ggplot(gene_plot, aes(Treatment, value)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment)) +
  geom_point(size = 3, alpha = 1, pch = 21, aes(fill = Treatment)) +
  labs(x = "Treatment",
       y = NULL) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(~ label, scales = "free_y") +
  guides(colour = "none") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Subset to the 5 with treatment effect
sig2 <- stat_loop %>%
  filter(Trt < 0.05)
ko_musicc_goi_sig2 <- ko_musicc_goi %>%
  dplyr::select(sig2$KO) %>%
  rownames_to_column(var = "SampleID")
meta102 <- meta10 %>%
  left_join(., ko_musicc_goi_sig2, by = "SampleID")
gene_plot2 <- meta102 %>%
  pivot_longer(cols = c(16:19)) %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "+Basalt", "+MP1", "+MP1+Basalt")))
ggplot(gene_plot2, aes(Treatment, value)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment)) +
  geom_jitter(size = 3, alpha = 1, pch = 21, aes(fill = Treatment)) +
  labs(x = "Treatment",
       y = NULL) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(~ name, scales = "free_y") +
  guides(colour = "none") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))




#### _MUSiCC All ####
# Look at all genes not just genes of interest
# Run on counts instead of coverage. with -n and -c (use this)
ko_musicc <- read.delim("ko_count_musicc_nc.tab") %>%
  column_to_rownames(var = "KO") %>%
  t() %>%
  as.data.frame()

# Check how many samples each KO was in?
# Need to remove those in < 3 samples...
prev <- data.frame(KO = colnames(ko_musicc),
                   n = colSums(ko_musicc != 0))
sum(prev$n == 1)
sum(prev$n == 2)
sum(prev$n >= 3)
sum(prev$n == 24)
prev3 <- prev %>%
  filter(n >= 3)

ko_musicc <- ko_musicc %>%
  dplyr::select(prev3$KO)

# Run ANOVA on 5850 KOs present in >= 3 samples
# Check if assumptions are met
# Loop through the KOs and run tests and store outputs
# Then plot as heatmap
meta10 <- read.delim("metadata.txt") %>%
  filter(Depth == 10) %>%
  mutate(SampleID = gsub("\\.", "_", sampleID))
sum(rownames(ko_musicc) != meta10$SampleID)

stat_loop <- as.data.frame(matrix(NA, ncol(ko_musicc), 7)) %>%
  set_names("KO", "Levene", "Shapiro", "Basalt", "MP1", "Int", "Trt")
for (i in 1:ncol(ko_musicc)) {
  # KO name
  stat_loop$KO[i] <- names(ko_musicc)[i]
  
  # Levene Test
  l1 <- leveneTest(ko_musicc[,i] ~ meta10$Treatment)
  stat_loop$Levene[i] <- l1$`Pr(>F)`[1]
  
  # Models
  m <- aov(ko_musicc[,i] ~ meta10$Basalt * meta10$MP1)
  a <- Anova(m, type = "III", singular.ok = TRUE)
  
  t <- aov(ko_musicc[,i] ~ meta10$Treatment)
  ta <- Anova(t)
  
  stat_loop$Basalt[i] <- a$`Pr(>F)`[2]
  stat_loop$MP1[i] <- a$`Pr(>F)`[3]
  stat_loop$Int[i] <- a$`Pr(>F)`[4]
  stat_loop$Trt[i] <- ta$`Pr(>F)`[1]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  stat_loop$Shapiro[i] <- s1$p.value
}

# Count how many passed Levene Test
sum(stat_loop$Levene > 0.05) # 5626/5850

# Count how many passed Shapiro Test
sum(stat_loop$Shapiro > 0.05) # 4662/5850

# Add Pfdr
stat_loop <- stat_loop %>%
  mutate(BasaltPadj = p.adjust(Basalt, method = "fdr"),
         MP1Padj = p.adjust(MP1, method = "fdr"),
         IntPadj = p.adjust(Int, method = "fdr"),
         TrtPadj = p.adjust(Trt, method = "fdr"))

# Sig
stat_loop_sig <- stat_loop %>%
  filter(Basalt < 0.01 | MP1 < 0.01 | Int < 0.01)
stat_loop_sig$Definition <- NA

# Need to confer info to these, use keggFind
for (i in 1:nrow(stat_loop_sig)) {
  stat_loop_sig$Definition[i] <- keggFind(database = "ko", query = stat_loop_sig$KO[i])
}

# Export send to Andes Ag
stat_loop_sig_export <- stat_loop_sig %>%
  dplyr::select(KO, Definition)
write.csv(stat_loop_sig_export, "sig_kos_p01.csv", row.names = F)

# Plot, with annotation column for sample, annotation rows for sig
ko_musicc_t <- ko_musicc %>%
  dplyr::select(all_of(stat_loop_sig$KO)) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(D_R1_10, D_R2_10, D_R3_10, D_R4_10, D_R5_10, D_R6_10,
                B_R1_10, B_R2_10, B_R3_10, B_R4_10, B_R5_10, B_R6_10,
                A_R1_10, A_R2_10, A_R3_10, A_R4_10, A_R5_10, A_R6_10,
                C_R1_10, C_R2_10, C_R3_10, C_R4_10, C_R5_10, C_R6_10)
meta10 <- meta10 %>%
  arrange(SampleID = colnames(ko_musicc_t)) %>%
  arrange(SampleID = colnames(ko_musicc_t))
sum(colnames(ko_musicc_t) != meta10$SampleID)
stat_loop_sig <- stat_loop_sig %>%
  mutate(Basalt_Sig = ifelse(Basalt < 0.01, "P < 0.01", " "),
         MP1_Sig = ifelse(MP1 < 0.01, "P < 0.01", " "),
         Int_Sig = ifelse(Int < 0.01, "P < 0.01", " "))
ann_cols <- data.frame(row.names = colnames(ko_musicc_t), 
                       Treatment = meta10$Treatment)
table(ann_cols$Treatment)
ann_rows <- data.frame(row.names = rownames(ko_musicc_t)) %>%
  mutate(`Basalt` = stat_loop_sig$Basalt_Sig,
         `MP1` = stat_loop_sig$MP1_Sig,
         `Int` = stat_loop_sig$Int_Sig)
table(ann_rows$Basalt)
table(ann_rows$MP1)
table(ann_rows$Int)
ann_colors <- list(`Basalt` = c("P < 0.01" = "#31688EFF",
                                " " = "white"),
                   `MP1` = c("P < 0.01" = "#35B779FF",
                             " " = "white"),
                   `Int` = c("P < 0.01" = "#FDE725FF",
                             " " = "white"),
                   Treatment = c("Control" = "#440154FF",
                                 "+Basalt" = "#31688EFF",
                                 "+MP1" = "#35B779FF",
                                 "+MP1+Basalt" = "#FDE725FF"))
pheatmap(ko_musicc_t,
         legend = T,
         legend_labels = c("-3  ", "-2  ", "-1  ", "0  ", "1  ", "2  ", "3  "),
         cluster_rows = T,
         cluster_cols = F,
         scale = "row",
         angle_col = 315,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         annotation_names_row = T,
         fontsize = 8,
         fontsize_row = 6,
         na_col = "white",
         border_color = "white",
         gaps_col = c(6, 12, 18),
         labels_row = stat_loop_sig$Definition,
         show_colnames = F,
         filename = "KO_Heatmap_MUSiCC_All.png",
         width = 12,
         height = 8)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# A bit vague to tell direction - subset to the sig ones and plot
sig <- stat_loop %>%
  filter(Basalt < 0.05 | MP1 < 0.05 | Int < 0.05)
ko_musicc_goi_sig <- ko_musicc_goi %>%
  dplyr::select(sig$KO) %>%
  rownames_to_column(var = "SampleID")
meta10 <- meta10 %>%
  left_join(., ko_musicc_goi_sig, by = "SampleID")
gene_plot <- meta10 %>%
  pivot_longer(cols = c(16:(16 + nrow(sig) - 1))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "+Basalt", "+MP1", "+MP1+Basalt"))) %>%
  left_join(., goi_ko_present, by = c("name" = "KO"))
pdf("SigKO_MUSiCC.pdf", width = 8, height = 6)
ggplot(gene_plot, aes(Treatment, value)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment)) +
  geom_point(size = 3, alpha = 1, pch = 21, aes(fill = Treatment)) +
  labs(x = "Treatment",
       y = NULL) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(~ label, scales = "free_y") +
  guides(colour = "none") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Subset to the 5 with treatment effect
sig2 <- stat_loop %>%
  filter(Trt < 0.05)
ko_musicc_goi_sig2 <- ko_musicc_goi %>%
  dplyr::select(sig2$KO) %>%
  rownames_to_column(var = "SampleID")
meta102 <- meta10 %>%
  left_join(., ko_musicc_goi_sig2, by = "SampleID")
gene_plot2 <- meta102 %>%
  pivot_longer(cols = c(16:19)) %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "+Basalt", "+MP1", "+MP1+Basalt")))
ggplot(gene_plot2, aes(Treatment, value)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment)) +
  geom_jitter(size = 3, alpha = 1, pch = 21, aes(fill = Treatment)) +
  labs(x = "Treatment",
       y = NULL) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(~ name, scales = "free_y") +
  guides(colour = "none") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))



#### _DESeq2 ####
goi <- read_excel("GenesList.xlsx")
goi_ko <- goi %>%
  filter(KO != "NA") %>%
  arrange(Pathway, `Gene name`) # 94
length(unique(goi_ko$`Gene name`)) # 94
table(goi_ko$Pathway)

# KO read counts
ko_count <- read.delim("ko_counts.txt")

# KO goi
ko_goi <- ko_count %>%
  filter(KO %in% goi_ko$KO) %>% 
  column_to_rownames(var = "KO") %>%
  t() %>%
  as.data.frame() # Only 71 of 94 present!
goi_ko_present <- goi_ko %>%
  filter(KO %in% ko_count$KO)

# Which ones are absent?
goi_ko_absent <- goi_ko %>%
  filter(KO %notin% ko_count$KO)

# DESeq normalization, with no design
ko_count <- ko_count %>%
  column_to_rownames(var = "KO")
meta10 <- read.delim("metadata.txt") %>%
  filter(Depth == 10) %>%
  mutate(SampleID = gsub("\\.", "_", sampleID))
sum(rownames(ko_goi) != meta10$SampleID)

dds_input_1 <- DESeqDataSetFromMatrix(countData = round(ko_count),
                                      colData = meta10,
                                      design = ~ Treatment)
dds_input_SF <- estimateSizeFactors(dds_input_1)
dds_input_D <- estimateDispersions(dds_input_SF)
ko_DESeq <- as.data.frame((counts(dds_input_D, normalized = T))) %>%
  t() %>%
  as.data.frame()
sum(rownames(ko_DESeq) != meta10$SampleID)

# Wald Differential abundance test
dds_input_da <- DESeqDataSetFromMatrix(countData = round(ko_count),
                                       colData = meta10,
                                       design = ~ Treatment)
wald <- DESeq(object = dds_input_da, test = "Wald", fitType = "parametric")
res_b <- results(wald,
                 contrast = c("Treatment", "+Basalt", "Control"),
                 independentFiltering = FALSE,
                 pAdjustMethod = "fdr")
sum(res_b$pvalue < 0.05, na.rm = T)
summary(res_b)
plotMA(res_b)
res_m <- results(wald,
                 contrast = c("Treatment", "+MP1", "Control"),
                 independentFiltering = FALSE,
                 pAdjustMethod = "fdr")
sum(res_m$pvalue < 0.05, na.rm = T)
summary(res_m)
plotMA(res_m)
res_i <- results(wald,
                 contrast = c("Treatment", "+MP1+Basalt", "Control"),
                 independentFiltering = FALSE,
                 pAdjustMethod = "fdr")
sum(res_i$pvalue < 0.05, na.rm = T)
summary(res_i)
plotMA(res_i)

# LRT
lrt <- DESeq(object = dds_input_da, test = "LRT", reduced = ~1)
res_bL <- results(lrt,
                  contrast = c("Treatment", "+Basalt", "Control"),
                  independentFiltering = FALSE,
                  pAdjustMethod = "fdr")
sum(res_bL$pvalue < 0.05, na.rm = T)
summary(res_bL)
plotMA(res_bL)
res_mL <- results(lrt,
                  contrast = c("Treatment", "+MP1", "Control"),
                  independentFiltering = FALSE,
                  pAdjustMethod = "fdr")
sum(res_mL$pvalue < 0.05, na.rm = T)
summary(res_mL)
plotMA(res_mL)
res_iL <- results(lrt,
                  contrast = c("Treatment", "+MP1+Basalt", "Control"),
                  independentFiltering = FALSE,
                  pAdjustMethod = "fdr")
sum(res_iL$pvalue < 0.05, na.rm = T)
summary(res_iL)
plotMA(res_iL)

# Manipulate results - collate to goi, recalc padj, note sig Basalt or MP1 or interaction
wald_results <- data.frame(KO = res_b@rownames,
                           BasaltFC = res_b@listData$log2FoldChange,
                           BasaltP = res_b@listData$pvalue,
                           MP1FC= res_m@listData$log2FoldChange,
                           MP1P = res_m@listData$pvalue,
                           IFC = res_i@listData$log2FoldChange,
                           IP = res_i@listData$pvalue) %>%
  filter(KO %in% goi_ko_present$KO) %>%
  mutate(BasaltPadj = p.adjust(BasaltP, method = "fdr"),
         MP1Padj = p.adjust(MP1P, method = "fdr"),
         IntPadj = p.adjust(IP, method = "fdr"))

lrt_results <- data.frame(KO = res_bL@rownames,
                          BasaltFC = res_bL@listData$log2FoldChange,
                          BasaltP = res_bL@listData$pvalue,
                          MP1FC= res_mL@listData$log2FoldChange,
                          MP1P = res_mL@listData$pvalue,
                          IFC = res_iL@listData$log2FoldChange,
                          IP = res_iL@listData$pvalue) %>%
  filter(KO %in% goi_ko_present$KO) %>%
  mutate(BasaltPadj = p.adjust(BasaltP, method = "fdr"),
         MP1Padj = p.adjust(MP1P, method = "fdr"),
         IntPadj = p.adjust(IP, method = "fdr"))

# Use Wald results

# Plot, with annotation column for sample, annotation rows for sig and pathway
ko_goi_t <- ko_goi %>%
  dplyr::select(all_of(goi_ko_present$KO)) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(D_R1_10, D_R2_10, D_R3_10, D_R4_10, D_R5_10, D_R6_10,
                B_R1_10, B_R2_10, B_R3_10, B_R4_10, B_R5_10, B_R6_10,
                A_R1_10, A_R2_10, A_R3_10, A_R4_10, A_R5_10, A_R6_10,
                C_R1_10, C_R2_10, C_R3_10, C_R4_10, C_R5_10, C_R6_10)
meta10 <- meta10 %>%
  arrange(SampleID = colnames(ko_goi_t)) %>%
  arrange(SampleID = colnames(ko_goi_t))
sum(colnames(ko_goi_t) != meta10$SampleID)
sum(rownames(ko_goi_t) != goi_ko_present$KO)
goi_ko_present <- goi_ko_present %>%
  left_join(., wald_results, by = "KO") %>%
  mutate(Basalt_Sig = ifelse(BasaltP < 0.05, "P < 0.05", " "),
         MP1_Sig = ifelse(MP1P < 0.05, "P < 0.05", " "),
         Int_Sig = ifelse(IP < 0.05, "P < 0.05", " ")) %>%
  mutate(label = paste(KO, `Gene name`, sep = " "))
ann_cols <- data.frame(row.names = colnames(ko_musicc_goi_t), 
                       Treatment = meta10$Treatment)
table(ann_cols$Treatment)
ann_rows <- data.frame(row.names = rownames(ko_musicc_goi_t)) %>%
  mutate(`Pathway` = goi_ko_present$Pathway,
         `Basalt` = goi_ko_present$Basalt_Sig,
         `MP1` = goi_ko_present$MP1_Sig,
         `Int` = goi_ko_present$Int_Sig)
table(ann_rows$Pathway)
table(ann_rows$Basalt)
brewer_pal(palette = "Paired")(12)
viridis_pal()(4)
ann_colors <- list(`Pathway` = c("Acetic acid" = "#A6CEE3",
                                 "Ammonium" = "#1F78B4",
                                 "Biofilm" = "#B2DF8A",
                                 "Fatty Acid Metabolism" = "#33A02C",
                                 "Formic Acid" = "#FB9A99",
                                 "Gluconic acid" = "#E31A1C",
                                 "Glycolic acid" = "#FDBF6F",
                                 "Lactic acid" = "#FF7F00",
                                 "Proprionic acid" = "#CAB2D6",
                                 "Pulcherrimin" = "#6A3D9A",
                                 "Pyruvic acid" = "#FFFF99",
                                 "Siderophore" = "#B15928",
                                 "Succinic acid" = "grey90",
                                 "Urea" = "grey40"),
                   `Basalt` = c("P < 0.05" = "#31688EFF",
                                " " = "white"),
                   `MP1` = c("P < 0.05" = "#35B779FF",
                                " " = "white"),
                   `Int` = c("P < 0.05" = "#FDE725FF",
                                " " = "white"),
                   Treatment = c("Control" = "#440154FF",
                                 "+Basalt" = "#31688EFF",
                                 "+MP1" = "#35B779FF",
                                 "+MP1+Basalt" = "#FDE725FF"))
pheatmap(ko_goi_t,
         legend = T,
         #legend_labels = c("-3  ", "-2  ", "-1  ", "0  ", "1  ", "2  ", "3  "),
         cluster_rows = F,
         cluster_cols = F,
         scale = "row",
         angle_col = 315,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         annotation_names_row = T,
         fontsize = 8,
         fontsize_row = 6,
         na_col = "white",
         border_color = "white",
         gaps_col = c(6, 12, 18),
         labels_row = goi_ko_present$label,
         show_colnames = F,
         filename = "KO_Heatmap_DESeq.png",
         width = 6,
         height = 8)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### _CPM ####
# List of genes of interest sent by Andes Ag team
# Removed duplicate lysT, lysS. moved to pyruvate category
goi <- read_excel("GenesList.xlsx")
goi_ko <- goi %>%
  filter(KO != "NA") %>%
  arrange(Pathway, `Gene name`) # 94
length(unique(goi_ko$`Gene name`)) # 94
table(goi_ko$Pathway)

# KO read counts
ko_count <- read.delim("ko_counts.txt") %>%
  column_to_rownames(var = "KO")

# Normalize the KO Counts by calculating counts per million total mapped reads!
mapped_reads <- read.delim("combined_defaults.txt") %>%
  left_join(., rc, by = c("layers" = "SampleID")) %>%
  mutate(PerMapped = total_reads_mapped/ReadsAll*100)
rownames(mapped_reads) <- mapped_reads$layers
ko_cpm <- sweep(ko_count, 2, mapped_reads$total_reads_mapped, FUN = "/") * 1e6
ko_cpm <- ko_cpm %>%
  rownames_to_column(var = "KO")

# CPM goi
ko_cpm_goi <- ko_cpm %>%
  filter(KO %in% goi_ko$KO) %>% 
  column_to_rownames(var = "KO") %>%
  t() %>%
  as.data.frame() # Only 71 of 94 present!
goi_ko_present <- goi_ko %>%
  filter(KO %in% ko_cpm$KO)

# Which ones are absent?
goi_ko_absent <- goi_ko %>%
  filter(KO %notin% ko_cpm$KO)

# At this point, you really just need to run ANOVA on those!
# But first need to also make sure assumptions are met
# Need to loop through those 71 KOs and run tests and store outputs
# Then plot as heatmap
meta10 <- read.delim("metadata.txt") %>%
  filter(Depth == 10) %>%
  mutate(SampleID = gsub("\\.", "_", sampleID))
sum(rownames(ko_cpm_goi) != meta10$SampleID)

stat_loop <- as.data.frame(matrix(NA, 71, 7)) %>%
  set_names("KO", "Levene", "Shapiro", "Basalt", "MP1", "Int", "Trt")
for (i in 1:ncol(ko_cpm_goi)) {
  # KO name
  stat_loop$KO[i] <- names(ko_cpm_goi)[i]
  
  # Levene Test
  l1 <- leveneTest(ko_cpm_goi[,i] ~ meta10$Treatment)
  stat_loop$Levene[i] <- l1$`Pr(>F)`[1]
  
  # Models
  m <- aov(ko_cpm_goi[,i] ~ meta10$Basalt * meta10$MP1)
  a <- Anova(m, type = "III", singular.ok = TRUE)
  
  t <- aov(ko_cpm_goi[,i] ~ meta10$Treatment)
  ta <- Anova(t)
  
  stat_loop$Basalt[i] <- a$`Pr(>F)`[2]
  stat_loop$MP1[i] <- a$`Pr(>F)`[3]
  stat_loop$Int[i] <- a$`Pr(>F)`[4]
  stat_loop$Trt[i] <- ta$`Pr(>F)`[1]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  stat_loop$Shapiro[i] <- s1$p.value
}

stat_loop <- stat_loop %>%
  mutate(BasaltPadj = p.adjust(Basalt, method = "fdr"),
         MP1Padj = p.adjust(MP1, method = "fdr"),
         IntPadj = p.adjust(Int, method = "fdr"),
         TrtPadj = p.adjust(Trt, method = "fdr"))

# Plot, with annotation column for sample, annotation rows for sig and pathway
ko_cpm_goi_t <- ko_cpm_goi %>%
  dplyr::select(all_of(goi_ko_present$KO)) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(D_R1_10, D_R2_10, D_R3_10, D_R4_10, D_R5_10, D_R6_10,
                B_R1_10, B_R2_10, B_R3_10, B_R4_10, B_R5_10, B_R6_10,
                A_R1_10, A_R2_10, A_R3_10, A_R4_10, A_R5_10, A_R6_10,
                C_R1_10, C_R2_10, C_R3_10, C_R4_10, C_R5_10, C_R6_10)
meta10 <- meta10 %>%
  arrange(SampleID = colnames(ko_cpm_goi_t)) %>%
  arrange(SampleID = colnames(ko_cpm_goi_t))
sum(colnames(ko_cpm_goi_t) != meta10$SampleID)
sum(rownames(ko_cpm_goi_t) != goi_ko_present$KO)
goi_ko_present <- goi_ko_present %>%
  left_join(., stat_loop, by = "KO") %>%
  mutate(Basalt_Sig = ifelse(Basalt < 0.05, "P < 0.05", " "),
         MP1_Sig = ifelse(MP1 < 0.05, "P < 0.05", " "),
         Int_Sig = ifelse(Int < 0.05, "P < 0.05", " ")) %>%
  mutate(label = paste(KO, `Gene name`, sep = " "))
ann_cols <- data.frame(row.names = colnames(ko_cpm_goi_t), 
                       Treatment = meta10$Treatment)
table(ann_cols$Treatment)
ann_rows <- data.frame(row.names = rownames(ko_cpm_goi_t)) %>%
  mutate(`Pathway` = goi_ko_present$Pathway,
         `Basalt` = goi_ko_present$Basalt_Sig,
         `MP1` = goi_ko_present$MP1_Sig,
         `Int` = goi_ko_present$Int_Sig)
table(ann_rows$Pathway)
table(ann_rows$Basalt)
brewer_pal(palette = "Paired")(12)
viridis_pal()(4)
ann_colors <- list(`Pathway` = c("Acetic acid" = "#A6CEE3",
                                 "Ammonium" = "#1F78B4",
                                 "Biofilm" = "#B2DF8A",
                                 "Fatty Acid Metabolism" = "#33A02C",
                                 "Formic Acid" = "#FB9A99",
                                 "Gluconic acid" = "#E31A1C",
                                 "Glycolic acid" = "#FDBF6F",
                                 "Lactic acid" = "#FF7F00",
                                 "Proprionic acid" = "#CAB2D6",
                                 "Pulcherrimin" = "#6A3D9A",
                                 "Pyruvic acid" = "#FFFF99",
                                 "Siderophore" = "#B15928",
                                 "Succinic acid" = "grey90",
                                 "Urea" = "grey40"),
                   `Basalt` = c("P < 0.05" = "#31688EFF",
                                " " = "white"),
                   `MP1` = c("P < 0.05" = "#35B779FF",
                             " " = "white"),
                   `Int` = c("P < 0.05" = "#FDE725FF",
                             " " = "white"),
                   Treatment = c("Control" = "#440154FF",
                                 "+Basalt" = "#31688EFF",
                                 "+MP1" = "#35B779FF",
                                 "+MP1+Basalt" = "#FDE725FF"))
pheatmap(ko_cpm_goi_t,
         legend = T,
         #legend_labels = c("-3  ", "-2  ", "-1  ", "0  ", "1  ", "2  ", "3  "),
         cluster_rows = F,
         cluster_cols = F,
         scale = "row",
         angle_col = 315,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         annotation_names_row = T,
         fontsize = 8,
         fontsize_row = 6,
         na_col = "white",
         border_color = "white",
         gaps_col = c(6, 12, 18),
         labels_row = goi_ko_present$label,
         show_colnames = F,
         filename = "KO_Heatmap_CPM.png",
         width = 6,
         height = 8)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# A bit vague to tell direction - subset to the sig ones and plot
sig <- stat_loop %>%
  filter(Basalt < 0.05 | MP1 < 0.05 | Int < 0.05)
ko_cpm_goi_sig <- ko_cpm_goi %>%
  dplyr::select(sig$KO) %>%
  rownames_to_column(var = "SampleID")
meta10 <- meta10 %>%
  left_join(., ko_cpm_goi_sig, by = "SampleID")
gene_plot <- meta10 %>%
  pivot_longer(cols = c(16:(16 + nrow(sig) - 1))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "+Basalt", "+MP1", "+MP1+Basalt"))) %>%
  left_join(., goi_ko_present, by = c("name" = "KO"))
pdf("SigKO_CPM.pdf", width = 8, height = 6)
ggplot(gene_plot, aes(Treatment, value)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment)) +
  geom_point(size = 3, alpha = 1, pch = 21, aes(fill = Treatment)) +
  labs(x = "Treatment",
       y = NULL) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(~ label, scales = "free_y") +
  guides(colour = "none") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Subset to the 4 with treatment effect
sig2 <- stat_loop %>%
  filter(Trt < 0.05)
ko_cpm_goi_sig2 <- ko_cpm_goi %>%
  dplyr::select(sig2$KO) %>%
  rownames_to_column(var = "SampleID")
meta102 <- meta10 %>%
  left_join(., ko_cpm_goi_sig2, by = "SampleID")
gene_plot2 <- meta102 %>%
  pivot_longer(cols = c(16:19)) %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "+Basalt", "+MP1", "+MP1+Basalt")))
ggplot(gene_plot2, aes(Treatment, value)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment)) +
  geom_jitter(size = 3, alpha = 1, pch = 21, aes(fill = Treatment)) +
  labs(x = "Treatment",
       y = NULL) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(~ name, scales = "free_y") +
  guides(colour = "none") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))



#### _FeGenie ####
fe <- read.csv("FeGenie-heatmap-data.csv") %>%
  column_to_rownames(var = "X") %>%
  set_names(substr(names(.), start = 1, stop = 7)) %>%
  dplyr::select(-X.1)
fe <- fe[rowSums(fe) != 0, ]
pheatmap(fe, 
         cluster_cols = F,
         cluster_rows = F,
         scale = "row")

# CPM - per total mapped reads!
mapped_reads <- read.delim("combined_defaults.txt") %>%
  left_join(., rc, by = c("layers" = "SampleID")) %>%
  mutate(PerMapped = total_reads_mapped/ReadsAll*100)
rownames(mapped_reads) <- mapped_reads$layers
range(mapped_reads$PerMapped) # 19.62% to 34.68%
plot(mapped_reads$ReadsAll, mapped_reads$total_reads_mapped)
all(colnames(fe) == rownames(mapped_reads))

fe_cpm <- sweep(fe, 2, mapped_reads$total_reads_mapped, FUN = "/") * 1e6
pheatmap(fe_cpm, 
         cluster_cols = F,
         cluster_rows = F,
         scale = "row")
fe_cpm_t <- as.data.frame(t(fe_cpm))

# Stats
meta10 <- read.delim("metadata.txt") %>%
  filter(Depth == 10) %>%
  mutate(SampleID = gsub("\\.", "_", sampleID))
sum(rownames(fe_cpm_t) != meta10$SampleID)
stat_loop <- as.data.frame(matrix(NA, 11, 7)) %>%
  set_names("KO", "Levene", "Shapiro", "Basalt", "MP1", "Int", "Trt")
for (i in 1:ncol(fe_cpm_t)) {
  # KO name
  stat_loop$KO[i] <- names(fe_cpm_t)[i]
  
  # Levene Test
  l1 <- leveneTest(fe_cpm_t[,i] ~ meta10$Treatment)
  stat_loop$Levene[i] <- l1$`Pr(>F)`[1]
  
  # Models
  m <- aov(fe_cpm_t[,i] ~ meta10$Basalt * meta10$MP1)
  a <- Anova(m, type = "III", singular.ok = TRUE)
  
  t <- aov(fe_cpm_t[,i] ~ meta10$Treatment)
  ta <- Anova(t)
  
  stat_loop$Basalt[i] <- a$`Pr(>F)`[2]
  stat_loop$MP1[i] <- a$`Pr(>F)`[3]
  stat_loop$Int[i] <- a$`Pr(>F)`[4]
  stat_loop$Trt[i] <- ta$`Pr(>F)`[1]
  
  # Shapiro Test
  s1 <- shapiro.test(m$residuals)
  stat_loop$Shapiro[i] <- s1$p.value
}

stat_loop <- stat_loop %>%
  mutate(BasaltPadj = p.adjust(Basalt, method = "fdr"),
         MP1Padj = p.adjust(MP1, method = "fdr"),
         IntPadj = p.adjust(Int, method = "fdr"),
         TrtPadj = p.adjust(Trt, method = "fdr"))

# Looks like levene passed but shapiro did not.
hist(fe_cpm)
fe_cpm_switch <- as.data.frame(t(fe_cpm))
hist(fe_cpm_switch$`iron_aquisition-iron_transport`)
hist(fe_cpm_switch$`iron_aquisition-heme_transport`)
hist(fe_cpm_switch$`iron_aquisition-siderophore_synthesis`)
hist(fe_cpm_switch$`iron_aquisition-siderophore_transport`)
hist(fe_cpm_switch$`iron_aquisition-siderophore_transport_potential`)
hist(fe_cpm_switch$iron_gene_regulation)
hist(fe_cpm_switch$iron_oxidation)
hist(fe_cpm_switch$possible_iron_oxidation_and_possible_iron_reduction)
hist(fe_cpm_switch$probable_iron_reduction)
hist(fe_cpm_switch$iron_reduction)
hist(fe_cpm_switch$iron_storage)

# Plot, with annotation column for sample, annotation rows for sig and pathway
fe_cpm_t <- fe_cpm_t %>%
  dplyr::select(D_R1_10, D_R2_10, D_R3_10, D_R4_10, D_R5_10, D_R6_10,
                B_R1_10, B_R2_10, B_R3_10, B_R4_10, B_R5_10, B_R6_10,
                A_R1_10, A_R2_10, A_R3_10, A_R4_10, A_R5_10, A_R6_10,
                C_R1_10, C_R2_10, C_R3_10, C_R4_10, C_R5_10, C_R6_10)
meta10 <- meta10 %>%
  arrange(SampleID = colnames(fe_cpm_t)) %>%
  arrange(SampleID = colnames(fe_cpm_t))
sum(colnames(fe_cpm_t) != meta10$SampleID)
goi_ko_present <- fe_cpm %>%
  rownames_to_column(var = "KO") %>%
  left_join(., stat_loop, by = "KO") %>%
  mutate(Basalt_Sig = ifelse(Basalt < 0.05, "P < 0.05", " "),
         MP1_Sig = ifelse(MP1 < 0.05, "P < 0.05", " "),
         Int_Sig = ifelse(Int < 0.05, "P < 0.05", " "))
ann_cols <- data.frame(row.names = colnames(fe_cpm_t), 
                       Treatment = meta10$Treatment)
table(ann_cols$Treatment)
ann_rows <- data.frame(row.names = rownames(fe_cpm_t)) %>%
  mutate(`Basalt` = goi_ko_present$Basalt_Sig,
         `MP1` = goi_ko_present$MP1_Sig,
         `Int` = goi_ko_present$Int_Sig)
table(ann_rows$Basalt)
ann_colors <- list(`Basalt` = c("P < 0.05" = "#31688EFF",
                                " " = "white"),
                   `MP1` = c("P < 0.05" = "#35B779FF",
                             " " = "white"),
                   `Int` = c("P < 0.05" = "#FDE725FF",
                             " " = "white"),
                   Treatment = c("Control" = "#440154FF",
                                 "+Basalt" = "#31688EFF",
                                 "+MP1" = "#35B779FF",
                                 "+MP1+Basalt" = "#FDE725FF"))
pheatmap(fe_cpm_t,
         legend = T,
         legend_labels = c("-4  ", "-2  ", "0  ", "2  ", "4  "),
         cluster_rows = F,
         cluster_cols = F,
         scale = "row",
         angle_col = 315,
         display_numbers = F,
         number_color = "black",
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         annotation_names_row = T,
         fontsize = 8,
         fontsize_row = 6,
         na_col = "white",
         border_color = "white",
         gaps_col = c(6, 12, 18),
         labels_row = goi_ko_present$label,
         show_colnames = F,
         filename = "KO_Heatmap_Fe_CPM.png",
         width = 5,
         height = 3)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# A bit vague to tell direction - subset to the sig ones and plot
sig <- stat_loop %>%
  filter(Basalt < 0.05 | MP1 < 0.05 | Int < 0.05)
ko_cpm_goi_sig <- fe_cpm_switch %>%
  dplyr::select(sig$KO) %>%
  rownames_to_column(var = "SampleID")
meta10 <- meta10 %>%
  left_join(., ko_cpm_goi_sig, by = "SampleID")
gene_plot <- meta10 %>%
  pivot_longer(cols = c(16:(16 + nrow(sig) - 1))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "+Basalt", "+MP1", "+MP1+Basalt"))) %>%
  left_join(., goi_ko_present, by = c("name" = "KO"))
pdf("SigKO_Fe_CPM.pdf", width = 8, height = 6)
ggplot(gene_plot, aes(Treatment, value)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment)) +
  geom_point(size = 3, alpha = 1, pch = 21, aes(fill = Treatment)) +
  labs(x = "Treatment",
       y = NULL) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(~ name, scales = "free_y") +
  guides(colour = "none") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



#### _KO Comp ####
# Look at richness and composition
meta10 <- read.delim("metadata.txt") %>%
  filter(Depth == 10) %>%
  mutate(SampleID = gsub("\\.", "_", sampleID)) %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "+Basalt", "+MP1", "+MP1+Basalt")))
ko_musicc <- read.delim("ko_count_musicc_nc.tab") %>%
  column_to_rownames(var = "KO")
sum(colnames(ko_musicc) != meta10$SampleID)
input <- list()
input$map_loaded <- meta10
rownames(input$map_loaded) <- meta10$SampleID
input$data_loaded <- ko_musicc
input$taxonomy_loaded <- data.frame(KO = rownames(ko_musicc),
                                    taxonomy1 = NA,
                                    taxonomy2 = NA,
                                    taxonomy3 = NA,
                                    taxonomy4 = NA,
                                    taxonomy5 = NA,
                                    taxonomy6 = NA,
                                    taxonomy7 = NA) %>%
  column_to_rownames(var = "KO")
input$map_loaded$rich <- specnumber(input$data_loaded, MARGIN = 2)
range(input$map_loaded$rich)
bc <- calc_dm(input$data_loaded)
plot_nmds(bc, metadata_map = input$map_loaded, color_cat = "Treatment") +
  scale_color_viridis_d() +
  labs(color = "Treatment")
adonis2(bc ~ input$map_loaded$Treatment)
ggplot(input$map_loaded, aes(Treatment, rich)) +
  geom_boxplot(outliers = F, show.legend = F, aes(colour = Treatment)) +
  geom_jitter(size = 3, width = 0.2, alpha = 1, pch = 21, aes(fill = Treatment)) +
  labs(x = "Treatment",
       y = NULL) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = "none") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))



#### _Compare Methods ####
# Compare MUSiCC normalized coverage vs DESeq normalized counts vs CPM KO-mapped reads
# Could also explore other CPM options - # reads, # assembled bases, KO-mapped reads
# These are remarkably similar...
cor(ko_goi_t, ko_musicc_goi_t)
plot(ko_goi_t$A_R1_10, ko_musicc_goi_t$A_R1_10)
plot(ko_goi_t$A_R2_10, ko_musicc_goi_t$A_R2_10)
plot(ko_goi_t$A_R3_10, ko_musicc_goi_t$A_R3_10)
plot(ko_goi_t$A_R4_10, ko_musicc_goi_t$A_R4_10)
plot(ko_goi_t$A_R5_10, ko_musicc_goi_t$A_R5_10)
plot(ko_goi_t$A_R6_10, ko_musicc_goi_t$A_R6_10)
plot(ko_goi_t$B_R1_10, ko_musicc_goi_t$B_R1_10)
plot(ko_goi_t$B_R2_10, ko_musicc_goi_t$B_R2_10)
plot(ko_goi_t$B_R3_10, ko_musicc_goi_t$B_R3_10)
plot(ko_goi_t$B_R4_10, ko_musicc_goi_t$B_R4_10)
plot(ko_goi_t$B_R5_10, ko_musicc_goi_t$B_R5_10)
plot(ko_goi_t$B_R6_10, ko_musicc_goi_t$B_R6_10)

plot(ko_count$A_R1_10, ko_DESeq[2,])

