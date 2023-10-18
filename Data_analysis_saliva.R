# 1. Clearing workspace, loading libraries, setting seed ----

# Clear environment and set working directory
rm(list = ls())
setwd("/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project")

# Load libraries
library(Phenoflow)
library(flowViz)
library(flowFDA)
library(flowAI)
library(vegan)
library(ggplot2)
library(ggtext)
library(RColorBrewer)
#library(ggrepel)
#library(ape)
library(gridExtra)
#library(grid)
library(scales)
library(cowplot)
library(reshape2)
library(dplyr)
library(tidyverse)
library(ggcyto)
library(doParallel)
library(randomForest)
library(caret)
library(pairwiseAdonis)
#library(ROCR)
library(mltools)
library(data.table)

# Make R use multiple cores
# nCores <- detectCores()
# useCores <- nCores - 10 # Keep at least 10 cores free for other users on server
# cl <- makeCluster(useCores, type = "PSOCK") # type argument: PSOCK ~ all systems, FORK ~ Unix/Mac
# registerDoParallel(cl)
#stopCluster(cl)

# Set seed
seed <- 777
set.seed(seed = seed)

# Load data
Datapath <- "/media/cmetfcm/Fabian/Oral_Microbiology/MOA_Antimicrobials/20210428_MOA_Antibiotics_Saliva_Hanna"
fcsfiles <- list.files(path = Datapath, recursive = TRUE, pattern = ".fcs", full.names = TRUE)
flowData <- flowCore::read.flowSet(files = fcsfiles, transformation = FALSE, emptyValue = F)

# Load metadata
metadata <- read.table("/media/cmetfcm/Fabian/Oral_Microbiology/MOA_Antimicrobials/20210428_MOA_Antibiotics_Saliva_Hanna/20210428_Metadata_MOA_Saliva.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE) 


# 2. Transformation of data ----

flowData_transformed <- transform(flowData,
                                  `FSC-H` = asinh(`FSC-H`),
                                  `SSC-H` = asinh(`SSC-H`),
                                  `BL1-H` = asinh(`BL1-H`),
                                  `BL3-H` = asinh(`BL3-H`),
                                  `FSC-A` = asinh(`FSC-A`),
                                  `SSC-A` = asinh(`SSC-A`),
                                  `BL1-A` = asinh(`BL1-A`),
                                  `BL3-A` = asinh(`BL3-A`),
                                  `FSC-W` = asinh(`FSC-W`),
                                  `SSC-W` = asinh(`SSC-W`),
                                  `BL1-W` = asinh(`BL1-W`),
                                  `BL3-W` = asinh(`BL3-W`))

param = c("FSC-H", "SSC-H", "BL1-H", "BL3-H", "FSC-A", "SSC-A", "BL1-A", "BL3-A", "FSC-W", "SSC-W", "BL1-W", "BL3-W")


# 3. Exploration of data and quality control of data ----

# Create plot for green vs red fluorescence
xyplot(`BL3-H`~`BL1-H`, data = flowData_transformed,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(5, 15))),
       axis = axis.default, nbin = 125, main = "Green vs red fluorescence (BL1-H vs BL3-H)", xlab = "BL1-H", ylab = "BL3-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

# Create plot for forward vs side scatter
xyplot(`FSC-H`~`SSC-H`, data = flowData_transformed,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(0, 15))),
       axis = axis.default, nbin = 125, main = "Side vs forward scatter (SSC-H vs FSC-H)", xlab = "SSC-H", ylab = "FSC-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

# Create plot for green fluorescence vs side scatter
xyplot(`SSC-H`~`BL1-H`, data = flowData_transformed,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(5, 15))),
       axis = axis.default, nbin = 125, main = "Green fluorescence vs side scatter (BL1-H vs SSC-H)", xlab = "BL1-H", ylab = "SSC-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)


# 4. Gating of cells ----

## 4.1. Gating all cells based on BL1-BL3 ----
# Construction of gate
sqrcut_total <- matrix(c(6.7, 7.9, 7.9, 7, 7.3, 14.5, 14.5, 6.7,
                         5.8, 7.9, 8.9, 8.9, 14.7, 14.7, 7.5, 4), ncol = 2, nrow = 8)
colnames(sqrcut_total) <- c("BL1-H", "BL3-H")
polyGateTotal <- polygonGate(.gate = sqrcut_total, filterId = "Cells")

# Gating quality check
# Check whether the gate is in the right position
xyplot(`BL3-H`~`BL1-H`, data = flowData_transformed[c(44:55)], filter = polyGateTotal,
       scales = list(y = list(limits = c(2, 15)),
                     x = list(limits = c(5, 15))),
       axis = axis.default, nbin = 125, main = "Gating quality check", xlab = "BL1-H", ylab = "BL3-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

# Check whether the gate is in the right position (all data)
p_QC_gating <- xyplot(`BL3-H`~`BL1-H`, data = flowData_transformed, filter = polyGateTotal,
                      scales = list(y = list(limits = c(2, 15)),
                                    x = list(limits = c(5, 15))),
                      axis = axis.default, nbin = 125, main = "Gating quality check", xlab = "BL1-H", ylab = "BL3-H",
                      par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_QC_gating)

p_gating_strat_total <- xyplot(`BL3-H`~`BL1-H`, data = flowData_transformed[c(80, 84, 87)], filter = polyGateTotal,
                               scales = list(y = list(limits = c(4, 15), cex = 1),
                                             x = list(limits = c(6, 15), cex = 1)),
                               axis = axis.default, nbin = 125, main = NULL, xlab = list(label = "BL1-H", cex = 1.5), ylab = list(label = "BL3-H", cex = 1.5),
                               strip = strip.custom(factor.levels = c("Heat Killed Control", "Positive Control", "Sterile PBS")),
                               par.strip.text = list(col = "white", font = 1, cex = 1.5), smooth = FALSE)
print(p_gating_strat_total)

# Subset the data to the cells only
flowData_transformed_gated <- Subset(flowData_transformed, polyGateTotal)

## 4.2. Additional quality control ----
# Check singlets by showing height vs width of the same parameter
p_singlets <- xyplot(`BL1-W`~`BL1-H`, data = flowData_transformed_gated,
                     scales = list(y = list(limits = c(0, 9)),
                                   x = list(limits = c(6, 15))),
                     axis = axis.default, nbin = 125, main = "Singlet analysis (BL1-H vs BL1-W)", xlab = "BL1-H", ylab = "BL1-W",
                     par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_singlets)

p_singlets_BL1_AH <- xyplot(`BL1-A`~`BL1-H`, data = flowData_transformed_gated,
                            scales = list(y = list(limits = c(5, 15)),
                                          x = list(limits = c(5, 15))),
                            axis = axis.default, nbin = 125, main = "Singlet analysis (BL1-H vs BL1-A)", xlab = "BL1-H", ylab = "BL1-A",
                            par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_singlets_BL1_AH)

p_singlets_SSC <- xyplot(`SSC-W`~`SSC-H`, data = flowData_transformed_gated,
                         scales = list(y = list(limits = c(1, 9)),
                                       x = list(limits = c(4, 15))),
                         axis = axis.default, nbin = 125, main = "Singlet analysis (SSC-H vs SSC-W)", xlab = "SSC-H", ylab = "SSC-W",
                         par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_singlets_SSC)

p_singlets_FSC <- xyplot(`FSC-W`~`FSC-H`, data = flowData_transformed_gated,
                         scales = list(y = list(limits = c(1, 9)),
                                       x = list(limits = c(4, 15))),
                         axis = axis.default, nbin = 125, main = "Singlet analysis (FSC-H vs FSC-W)", xlab = "FSC-H", ylab = "FSC-W",
                         par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_singlets_FSC)

# Density plots
p_density_SSC_W <- autoplot(flowData_transformed_gated, "SSC-W")+
  labs(title = "Density plot SSC-W",
       x = "SSC-W",
       y = "Density")
print(p_density_SSC_W)

p_density_BL1_H <- autoplot(flowData_transformed_gated, "BL1-H")+
  labs(title = "Density plot BL1-H",
       x = "BL1-H",
       y = "Density")
print(p_density_BL1_H)

# Side scatter vs green fluorescence
xyplot(`BL1-H`~`SSC-H`, data = flowData_transformed_gated,
       scales = list(y = list(limits = c(5, 15)),
                     x = list(limits = c(4, 16))),
       axis = axis.default, nbin = 125, main = "Side scatter vs green fluorescence (SSC-H vs BL1-H) gated saliva", xlab = "SSC-H", ylab = "BL1-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)


# 5. Phenotypic diversity analysis raw data ----
# All cells will be considered for the diversity analysis

# Analysis for D2 only
flowData_transformed_D2_gated <- flowData_transformed_gated[c(44:86)]
metadata_D2 <- metadata[c(44:86), ]
row.names(metadata_D2) <- NULL

## 5.1. Normalization of data ----
# Normalize over the highest value of all parameters (as this will result in a value between -1 and 1). Use same transformation for height and area (same range). Leave out width for now, as using different values for normalization will increase or decrease the weight of certain measured parameters. If one uses a Gaussian model, normalization is not needed, but then optimization of the band width is necessary.
summary <- fsApply(x = flowData_transformed_D2_gated, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
max = max(summary[, "BL1-H"])
mytrans<- function(x) x/max

# Actual normalization of data
flowData_transformed_D2_norm <- transform(flowData_transformed_D2_gated,
                                          `FSC-H` = mytrans(`FSC-H`),
                                          `SSC-H` = mytrans(`SSC-H`),
                                          `BL1-H` = mytrans(`BL1-H`),
                                          `BL3-H` = mytrans(`BL3-H`),
                                          `FSC-A` = mytrans(`FSC-A`),
                                          `SSC-A` = mytrans(`SSC-A`),
                                          `BL1-A` = mytrans(`BL1-A`),
                                          `BL3-A` = mytrans(`BL3-A`),
                                          `FSC-W` = mytrans(`FSC-W`),
                                          `SSC-W` = mytrans(`SSC-W`),
                                          `BL1-W` = mytrans(`BL1-W`),
                                          `BL3-W` = mytrans(`BL3-W`))

## 5.2. Diversity analysis ----
param_div <- c("FSC-H", "SSC-H", "BL1-H", "BL3-H", "FSC-A", "SSC-A", "BL1-A", "BL3-A")

#### 5.2.1. Fingerprinting ----
# Randomly resample to a fixed number of cells
flowData_transformed_D2_resample <- FCS_resample(flowData_transformed_D2_norm, sample = 3400, replace = TRUE)

# Make metadata file for resampled data
names_flowData_transformed_D2_resample <- sampleNames(flowData_transformed_D2_resample)
metadata_D2_resample <- metadata_D2[match(names_flowData_transformed_D2_resample, metadata_D2$Filename), ]
row.names(metadata_D2_resample) <- NULL

# Select data without MembraneDisruption class
metadata_D2_resample_NoMDNoH <- metadata_D2_resample[c(1:9, 11:16, 32:40), ]
rownames(metadata_D2_resample_NoMDNoH) <- NULL

# Calculating fingerprint with bandwidth = 0.01
fbasis_D2 <- flowBasis(flowData_transformed_D2_resample, param = param_div, nbin = 128, 
                       bw = 0.01, normalize = function(x) x)

fbasis_D2_NoMDNoH <- flowBasis(flowData_transformed_D2_resample[c(1:9, 11:16, 32:40)], param = param_div, nbin = 128,
                               bw = 0.01, normalize = function(x) x)

#### 5.2.2. Ordination ----
# Beta-diversity assessment of fingerprint (Bray-Curtis)
beta.div_D2 <- beta_div_fcm(fbasis_D2, ord.type = "PCoA")
beta.div_D2_NMDS <- beta_div_fcm(fbasis_D2, ord.type = "NMDS")

beta.div_D2_NoMDNoH <- beta_div_fcm(fbasis_D2_NoMDNoH, ord.type = "PCoA")
beta.div_D2_NoMDNoH_NMDS <- beta_div_fcm(fbasis_D2_NoMDNoH, ord.type = "NMDS")

# Plot ordination
pbetadiv_D2_PCoA <- plot_beta_fcm(beta.div_D2, color = as.factor(metadata_D2_resample$MOA), shape = as.factor(metadata_D2_resample$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL)+
  scale_color_manual(values = c("darkorchid3", "#1919ff", "#ffae19", "red", "#4daf4a", "#a65628"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "Heat", "Membrane Disruption", "Control", "Protein Synthesis: 50S Inhibition"))+
  scale_shape_manual(values=c(1:14))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(pbetadiv_D2_PCoA)

pbetadiv_D2_NMDS <- plot_beta_fcm(beta.div_D2_NMDS, color = as.factor(metadata_D2_resample$MOA), shape = as.factor(metadata_D2_resample$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL, x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("darkorchid3", "#1919ff", "#ffae19", "red", "#4daf4a", "#a65628"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "Heat", "Membrane Disruption", "Control", "Protein Synthesis: 50S Inhibition"))+
  scale_shape_manual(values=c(1:14))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(pbetadiv_D2_NMDS)

pbetadiv_D2_NoMDNoH_PCoA <- plot_beta_fcm(beta.div_D2_NoMDNoH, color = as.factor(metadata_D2_resample_NoMDNoH$MOA), shape = as.factor(metadata_D2_resample_NoMDNoH$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL)+
  scale_color_manual(values = c("red", "#1919ff", "#4daf4a", "darkgoldenrod3"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "Control", "Protein Synthesis: 50S Inhibition"))+
  scale_shape_manual(values=c(1, 2, 5:8, 10, 11))+ # same symbols as other plots: 2, 3, 4, 7, 8, 9, 13, 35
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 26),
        plot.title = element_text(size = 30, face = "bold"),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 26),
        legend.text.align = 0,
        legend.position = "bottom")+
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 2, title.position = "top"), size = guide_legend(ncol = 1, title.position = "top"))
print(pbetadiv_D2_NoMDNoH_PCoA)

pbetadiv_D2_NoMDNoH_NMDS <- plot_beta_fcm(beta.div_D2_NoMDNoH_NMDS, color = as.factor(metadata_D2_resample_NoMDNoH$MOA), shape = as.factor(metadata_D2_resample_NoMDNoH$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL, x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("red", "#1919ff", "#4daf4a", "darkgoldenrod3"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "Control", "Protein Synthesis: 50S Inhibition"))+
  scale_shape_manual(values=c(1, 2, 5:8, 10, 11))+ # same symbols as other plots: 2, 3, 4, 7, 8, 9, 13, 35
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 26),
        plot.title = element_text(size = 30, face = "bold"),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 26),
        legend.text.align = 0,
        legend.position = "bottom")+
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 2, title.position = "top"), size = guide_legend(ncol = 1, title.position = "top"))
print(pbetadiv_D2_NoMDNoH_NMDS)


# 6. PhenoGMM ----

# Parameters to base the model on
paramGMM <- c("FSC-H", "SSC-H", "BL1-H", "BL3-H", "FSC-A", "SSC-A", "BL1-A", "BL3-A")

## 6.1. Optimization model ----
# Select data exlcuding MD and Heat
flowData_transformed_D2_GMM <- flowData_transformed_D2_gated[c(1:9, 11:16, 32:40)]

# Optimization of amount of phenotypes (Bayesian Information Criterion, BIC)
# Retain useful information only and pool all samples by MOA
fcs_PGMM_D2_NoMDNoH <- Phenoflow::FCS_pool(flowData_transformed_D2_GMM,
                                           stub = c("20210428_MOA_Antibiotics_Saliva_Hanna_Donor2_CellWallSynthesis",
                                                    "20210428_MOA_Antibiotics_Saliva_Hanna_Donor2_ProteinSynthesis50SInhibition",
                                                    "20210428_MOA_Antibiotics_Saliva_Hanna_Donor2_Control",
                                                    "20210428_MOA_Antibiotics_Saliva_Hanna_Donor2_DNAReplication"))
fcs_PGMM_D2_NoMDNoH <- FCS_resample(fcs_PGMM_D2_NoMDNoH, replace = TRUE, sample = 15000)
fcs_PGMM_D2_NoMDNoH <- fcs_PGMM_D2_NoMDNoH[, paramGMM]
fcs_PGMM_D2_NoMDNoH <- Phenoflow::FCS_pool(fcs_PGMM_D2_NoMDNoH, stub = "*")
fcs_PGMM_D2_NoMDNoH <- fcs_PGMM_D2_NoMDNoH[, paramGMM]
PhenoGMM_D2_NoMDNoH <- Phenoflow::PhenoGMM(fcs_x = fcs_PGMM_D2_NoMDNoH, param = paramGMM, downsample = FALSE, nG = 80, auto_nG = TRUE, nG_interval = 10, fcs_scale = FALSE, diagnostic_plot = TRUE)
saveRDS(object = PhenoGMM_D2_NoMDNoH, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/D2_24h_PhenoGMM_8param_10to80per10_NoMDNoH.rds")

# Visualization of BIC values
#PhenoGMM_D2_NoMDNoH <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/D2_24h_PhenoGMM_8param_10to80per10_NoMDNoH.rds") # Call for results that are to be plotted
NoClusters_PGMM <- dimnames(PhenoGMM_D2_NoMDNoH[[2]]$BIC)[[1]]
BICValues_PGMM <- data.frame(as.matrix(PhenoGMM_D2_NoMDNoH[[2]]$BIC)[1:length(NoClusters_PGMM), ])
BICValues_PGMM$NoClusters <- rownames(BICValues_PGMM)
BICValues_PGMM <- reshape2::melt(BICValues_PGMM, id.vars = "NoClusters")
colnames(BICValues_PGMM) <- c("NoClusters", "ModelType", "BIC")
BICValues_PGMM$NoClusters <- as.numeric(BICValues_PGMM$NoClusters)
BICValues_PGMM <- BICValues_PGMM[!is.na(BICValues_PGMM$BIC), ] # Remove NA values
BICValues_PGMM$ModelType <- droplevels(BICValues_PGMM$ModelType, except = unique(BICValues_PGMM$ModelType)) # Remove levels that are not being used

p_BIC_PGMM_D2_NoMDNoH <- BICValues_PGMM %>% 
  ggplot(data = ., aes(x = NoClusters, y = BIC))+
  geom_line(alpha = 1, aes (color = ModelType), show.legend = FALSE)+
  geom_point(shape = 21, size = 3, alpha = 1, aes(fill = ModelType))+
  labs(title = NULL, x = "Number of clusters", y = "BIC", fill = "Model type")+
  theme_bw()+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        legend.text.align = 0)
print(p_BIC_PGMM_D2_NoMDNoH)

## 6.2. Allocation data to model ----
NoClusters_PGMM_D2_NoMDNoH <- 50
PhenoGMM_D2_fixedclust_NoMDNoH <- Phenoflow::PhenoGMM(fcs_x = fcs_PGMM_D2_NoMDNoH, param = paramGMM, downsample = FALSE, nG = NoClusters_PGMM_D2_NoMDNoH, fcs_scale = FALSE, diagnostic_plot = TRUE)
saveRDS(object = PhenoGMM_D2_fixedclust_NoMDNoH, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/D2_24h_PhenoGMM_8param_50clust_NoMDNoH.rds")
#PhenoGMM_D2_fixedclust_NoMDNoH <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/D2_24h_PhenoGMM_8param_50clust_NoMDNoH.rds") # Call for results that are to be plotted

# Applying mask to data
testPred_D2_NoMDNoH <- PhenoMaskGMM(fcs_x = flowData_transformed_D2_GMM, gmm = PhenoGMM_D2_fixedclust_NoMDNoH, fcs_scale = FALSE)

results_D2_NoMDNoH <- testPred_D2_NoMDNoH[[1]]
rownames(results_D2_NoMDNoH) <- sampleNames(flowData_transformed_D2_GMM)
results_D2_NoMDNoH <- select(results_D2_NoMDNoH, -c(Sample_names))
results_D2_NoMDNoH[is.na(results_D2_NoMDNoH)] <- 0 # Replace NA values by 0
results_D2_NoMDNoH[1:ncol(results_D2_NoMDNoH)] <- lapply(results_D2_NoMDNoH[1:ncol(results_D2_NoMDNoH)], as.numeric)

results_D2_NoMDNoH_rel <- sweep(results_D2_NoMDNoH, MARGIN = 1, rowSums(results_D2_NoMDNoH), `/`)


## 6.3. Ordination of GMM output ----
distance_GMM_D2_NoMDNoH <- vegan::vegdist(x = results_D2_NoMDNoH_rel, method = "bray", binary = FALSE)
mds_GMM_D2_NoMDNoH <- stats::cmdscale(distance_GMM_D2_NoMDNoH, k = 2, eig = TRUE, add = TRUE)
NMDS_GMM_D2_NoMDNoH <- vegan::metaMDS(distance_GMM_D2_NoMDNoH, autotransform = FALSE, k = 2, trymax = 100)

plot_PCoA_GMM_D2_NoMDNoH <- plot_beta_fcm(mds_GMM_D2_NoMDNoH, color = as.factor(metadata_D2_resample_NoMDNoH$MOA), shape = as.factor(metadata_D2_resample_NoMDNoH$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 14, alpha = 1, stroke = 1)+
  labs(title = NULL)+
  scale_color_manual(values = c("red", "#1919ff", "#4daf4a", "darkgoldenrod3"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "Control", "Protein Synthesis: 50S Inhibition"))+
  scale_shape_manual(values=c(1, 2, 5:8, 10, 11))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 26),
        plot.title = element_text(size = 30, face = "bold"),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 26),
        legend.text.align = 0,
        legend.position = "bottom")+
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 2, title.position = "top"), size = guide_legend(ncol = 1, title.position = "top"))
print(plot_PCoA_GMM_D2_NoMDNoH)

plot_NMDS_GMM_D2_NoMDNoH <- plot_beta_fcm(NMDS_GMM_D2_NoMDNoH, color = as.factor(metadata_D2_resample_NoMDNoH$MOA), shape = as.factor(metadata_D2_resample_NoMDNoH$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 14, alpha = 1, stroke = 1)+
  labs(title = NULL, x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("red", "#1919ff", "#4daf4a", "darkgoldenrod3"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "Control", "Protein Synthesis: 50S Inhibition"))+
  scale_shape_manual(values=c(1, 2, 5:8, 10, 11))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 26),
        plot.title = element_text(size = 30, face = "bold"),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 26),
        legend.text.align = 0,
        legend.position = "bottom")+
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 2, title.position = "top"), size = guide_legend(ncol = 1, title.position = "top"))
print(plot_NMDS_GMM_D2_NoMDNoH)

# Adding ellipses to the ordinations
NMDS_GMM_D2_NoMDNoH_visual <- as.data.frame(vegan::scores(NMDS_GMM_D2_NoMDNoH))
NMDS_GMM_D2_NoMDNoH_visual$Compound <- metadata_D2_resample_NoMDNoH$Compound
NMDS_GMM_D2_NoMDNoH_visual$MOA <- metadata_D2_resample_NoMDNoH$MOA
saveRDS(object = NMDS_GMM_D2_NoMDNoH_visual, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/Objects_plots/NMDS_GMM_D2_visual_24h.rds")

plot_NMDS_GMM_D2_NoMDNoH_ellipse <- ggplot(NMDS_GMM_D2_NoMDNoH_visual, aes(x = NMDS1, y = NMDS2))+
  geom_point(size = 14, alpha = 1, stroke = 1, aes(shape = Compound, color = MOA))+
  labs(title = NULL, color = "Mechanism of Action", shape = "Compound", x = "NMDS1", y = "NMDS2")+
  theme_bw() +
  scale_color_manual(values = c("red", "#1919ff", "#4daf4a", "darkgoldenrod3"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "Control", "Protein Synthesis: 50S Inhibition"))+
  scale_shape_manual(values=c(1, 2, 5:8, 10, 11))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 26),
        plot.title = element_text(size = 30, face = "bold"),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 26),
        legend.text.align = 0,
        legend.position = "bottom")+
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 2, title.position = "top"), size = guide_legend(ncol = 1, title.position = "top"))+
  stat_ellipse(aes(color = MOA), level = 0.95)
print(plot_NMDS_GMM_D2_NoMDNoH_ellipse)

var_pcoa_D2_NoMDNoH <- vegan::eigenvals(mds_GMM_D2_NoMDNoH)/sum(vegan::eigenvals(mds_GMM_D2_NoMDNoH))
PcoA_D2_NoMDNoH <- as.data.frame(mds_GMM_D2_NoMDNoH$points)
names(PcoA_D2_NoMDNoH)[1:2] <- c("PCoA1", "PCoA2")
PcoA_D2_NoMDNoH$Compound <- metadata_D2_resample_NoMDNoH$Compound
PcoA_D2_NoMDNoH$MOA <- metadata_D2_resample_NoMDNoH$MOA
saveRDS(object = PcoA_D2_NoMDNoH, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/Objects_plots/PCoA_GMM_D2_visual_24h.rds")
saveRDS(object = var_pcoa_D2_NoMDNoH, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/Objects_plots/PCoA_GMM_D2_var_24h.rds")

plot_PCoA_GMM_D2_NoMDNoH_ellipse <- ggplot(PcoA_D2_NoMDNoH, aes(x = PCoA1, y = PCoA2))+
  geom_point(size = 14, alpha = 1, stroke = 1, aes(shape = Compound, color = MOA))+
  labs(title = NULL, color = "Mechanism of Action", shape = "Compound", x = paste0("PCoA1 (", round(100 * var_pcoa_D2_NoMDNoH[1], 1), "%)"), y = paste0("PCoA2 (", round(100 * var_pcoa_D2_NoMDNoH[2], 1), "%)")) +
  theme_bw() +
  scale_color_manual(values = c("red", "#1919ff", "#4daf4a", "darkgoldenrod3"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "Control", "Protein Synthesis: 50S Inhibition"))+
  scale_shape_manual(values=c(1, 2, 5:8, 10, 11))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 26),
        plot.title = element_text(size = 30, face = "bold"),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 26),
        legend.text.align = 0,
        legend.position = "bottom")+
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 2, title.position = "top"), size = guide_legend(ncol = 1, title.position = "top"))+
  stat_ellipse(aes(color = MOA), level = 0.95)
print(plot_PCoA_GMM_D2_NoMDNoH_ellipse)


## 6.4. Cluster analysis ----
labels_compounds_D2 <- metadata_D2_resample_NoMDNoH$Compound

clust_D2 <- stats::hclust(d = distance_GMM_D2_NoMDNoH)
clust_D2$labels <- labels_compounds_D2

dendro_data_D2 <- clust_D2 %>% 
  as.dendrogram() %>% 
  dendro_data()

dendro_data_D2$labels$MOA <- dendro_data_D2$labels$label
dendro_data_D2$labels$MOA <- gsub('Amoxicillin', 'Cell Wall Synthesis', dendro_data_D2$labels$MOA)
dendro_data_D2$labels$MOA <- gsub('Amoxiclav', 'Cell Wall Synthesis', dendro_data_D2$labels$MOA)
dendro_data_D2$labels$MOA <- gsub('Control', 'Control', dendro_data_D2$labels$MOA)
dendro_data_D2$labels$MOA <- gsub('Ciprofloxacin', 'DNA Replication', dendro_data_D2$labels$MOA)
dendro_data_D2$labels$MOA <- gsub('Metronidazole', 'DNA Replication', dendro_data_D2$labels$MOA)
dendro_data_D2$labels$MOA <- gsub('Clindamycin', 'Protein Synthesis: 50S Inhibition', dendro_data_D2$labels$MOA)
dendro_data_D2$labels$MOA <- gsub('Erythromycin', 'Protein Synthesis: 50S Inhibition', dendro_data_D2$labels$MOA)
dendro_data_D2$labels$MOA <- gsub('Azithromycin', 'Protein Synthesis: 50S Inhibition', dendro_data_D2$labels$MOA)

saveRDS(object = dendro_data_D2, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/Objects_plots/dendro_D2_24h.rds")

plot_clust_D2 <- ggplot(segment(dendro_data_D2))+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dendro_data_D2$labels, aes(x = x, y = y, label = label, color = MOA), angle = 90, hjust = 1.2, fontface = 'bold', size = 8)+
  labs(x = "Compound", y = "Height", color = "Mechanism of Action")+
  theme_bw()+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 24),
        legend.position = "right",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_color_manual(values = c("red", "#1919ff", "#4daf4a", "darkgoldenrod3"))+
  guides(color = guide_legend(ncol = 1, title.position = "top", override.aes = list(size = 12)))+
  ylim(c(-0.2, 1))
print(plot_clust_D2)