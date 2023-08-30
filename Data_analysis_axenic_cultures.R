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
Datapath_Av <- "/media/cmetfcm/Fabian/Oral_Microbiology/MOA_Antimicrobials/20210902_MOA_Antibiotics/Av"
fcsfiles_Av <- list.files(path = Datapath_Av, recursive = TRUE, pattern = ".fcs", full.names = TRUE)
flowData_Av <- flowCore::read.flowSet(files = fcsfiles_Av, transformation = FALSE, emptyValue = F)

Datapath_Fn <- "/media/cmetfcm/Fabian/Oral_Microbiology/MOA_Antimicrobials/20210902_MOA_Antibiotics/Fn"
fcsfiles_Fn <- list.files(path = Datapath_Fn, recursive = TRUE, pattern = ".fcs", full.names = TRUE)
flowData_Fn <- flowCore::read.flowSet(files = fcsfiles_Fn, transformation = FALSE, emptyValue = F)

# Load metadata
metadata_Av <- read.table("/media/cmetfcm/Fabian/Oral_Microbiology/MOA_Antimicrobials/20210902_MOA_Antibiotics/Av/Metadata_20210902_MOA_Antibiotics_Av.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE) 
metadata_Fn <- read.table("/media/cmetfcm/Fabian/Oral_Microbiology/MOA_Antimicrobials/20210902_MOA_Antibiotics/Fn/Metadata_20210902_MOA_Antibiotics_Fn.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE) 


# 2. Transformation of data ----

flowData_transformed_Av <- transform(flowData_Av,
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

flowData_transformed_Fn <- transform(flowData_Fn,
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

## 3.1. Av ----
# Create plot for green vs red fluorescence
xyplot(`BL3-H`~`BL1-H`, data = flowData_transformed_Av,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(5, 15))),
       axis = axis.default, nbin = 125, main = "Green vs red fluorescence (BL1-H vs BL3-H) Av", xlab = "BL1-H", ylab = "BL3-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

# Create plot for forward vs side scatter
xyplot(`FSC-H`~`SSC-H`, data = flowData_transformed_Av,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(0, 15))),
       axis = axis.default, nbin = 125, main = "Side vs forward scatter (SSC-H vs FSC-H) Av", xlab = "SSC-H", ylab = "FSC-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

# Create plot for green fluorescence vs side scatter
xyplot(`SSC-H`~`BL1-H`, data = flowData_transformed_Av,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(5, 15))),
       axis = axis.default, nbin = 125, main = "Green fluorescence vs side scatter (BL1-H vs SSC-H) Av", xlab = "BL1-H", ylab = "SSC-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

## 3.2. Fn ----
# Create plot for green vs red fluorescence
xyplot(`BL3-H`~`BL1-H`, data = flowData_transformed_Fn,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(5, 15))),
       axis = axis.default, nbin = 125, main = "Green vs red fluorescence (BL1-H vs BL3-H) Fn", xlab = "BL1-H", ylab = "BL3-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

# Create plot for forward vs side scatter
xyplot(`FSC-H`~`SSC-H`, data = flowData_transformed_Fn,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(0, 15))),
       axis = axis.default, nbin = 125, main = "Side vs forward scatter (SSC-H vs FSC-H) Fn", xlab = "SSC-H", ylab = "FSC-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

# Create plot for green fluorescence vs side scatter
xyplot(`SSC-H`~`BL1-H`, data = flowData_transformed_Fn,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(5, 15))),
       axis = axis.default, nbin = 125, main = "Green fluorescence vs side scatter (BL1-H vs SSC-H) Fn", xlab = "BL1-H", ylab = "SSC-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)


# 4. Gating of cells ----

## 4.1. Gating all cells based on BL1-BL3 ----
# Construction of gate
sqrcut_total <- matrix(c(6.7, 7.9, 7.9, 7, 7.3, 14.5, 14.5, 6.7,
                         5.8, 7.9, 8.9, 8.9, 14.7, 14.7, 7.5, 4), ncol = 2, nrow = 8)
colnames(sqrcut_total) <- c("BL1-H", "BL3-H")
polyGateTotal <- polygonGate(.gate = sqrcut_total, filterId = "Cells")

# Gating quality check
# Av
# Check whether the gate is in the right position
xyplot(`BL3-H`~`BL1-H`, data = flowData_transformed_Av[c(1, 13:16, 19, 22, 28, 36, 40, 48, 49, 57, 68, 75)], filter = polyGateTotal,
       scales = list(y = list(limits = c(2, 15)),
                     x = list(limits = c(5, 15))),
       axis = axis.default, nbin = 125, main = "Gating quality check Av", xlab = "BL1-H", ylab = "BL3-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

# Check whether the gate is in the right position (all data)
p_QC_gating_Av <- xyplot(`BL3-H`~`BL1-H`, data = flowData_transformed_Av, filter = polyGateTotal,
                         scales = list(y = list(limits = c(2, 15)),
                                       x = list(limits = c(5, 15))),
                         axis = axis.default, nbin = 125, main = "Gating quality check Av", xlab = "BL1-H", ylab = "BL3-H",
                         par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_QC_gating_Av)

p_gating_strat_Av_total <- xyplot(`BL3-H`~`BL1-H`, data = flowData_transformed_Av[c(13, 14, 19, 40)], filter = polyGateTotal,
                                  scales = list(y = list(limits = c(4, 15), cex = 1),
                                                x = list(limits = c(6, 15), cex = 1)),
                                  axis = axis.default, nbin = 125, main = NULL, xlab = list(label = "BL1-H", cex = 1.5), ylab = list(label = "BL3-H", cex = 1.5),
                                  strip = strip.custom(factor.levels = c("Sterile PBS", "Sterile BHI", "Positive Control", "Heat Killed Control")),
                                  par.strip.text = list(col = "white", font = 1, cex = 1.5), smooth = FALSE)
print(p_gating_strat_Av_total)

# Subset the data to the cells only
flowData_transformed_Av_gated <- Subset(flowData_transformed_Av, polyGateTotal)

# Fn
# Check whether the gate is in the right position
xyplot(`BL3-H`~`BL1-H`, data = flowData_transformed_Fn[c(1, 13:16, 19, 22, 28, 36, 40, 48, 49, 57, 68, 75)], filter = polyGateTotal,
       scales = list(y = list(limits = c(2, 15)),
                     x = list(limits = c(5, 15))),
       axis = axis.default, nbin = 125, main = "Gating quality check Fn", xlab = "BL1-H", ylab = "BL3-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

# Check whether the gate is in the right position (all data)
p_QC_gating_Fn <- xyplot(`BL3-H`~`BL1-H`, data = flowData_transformed_Fn, filter = polyGateTotal,
                         scales = list(y = list(limits = c(2, 15)),
                                       x = list(limits = c(5, 15))),
                         axis = axis.default, nbin = 125, main = "Gating quality check Fn", xlab = "BL1-H", ylab = "BL3-H",
                         par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_QC_gating_Fn)

p_gating_strat_Fn_total <- xyplot(`BL3-H`~`BL1-H`, data = flowData_transformed_Fn[c(13, 14, 19, 40)], filter = polyGateTotal,
                                  scales = list(y = list(limits = c(4, 15), cex = 1),
                                                x = list(limits = c(6, 15), cex = 1)),
                                  axis = axis.default, nbin = 125, main = NULL, xlab = list(label = "BL1-H", cex = 1.5), ylab = list(label = "BL3-H", cex = 1.5),
                                  strip = strip.custom(factor.levels = c("Sterile PBS", "Sterile BHI", "Positive Control", "Heat Killed Control")),
                                  par.strip.text = list(col = "white", font = 1, cex = 1.5), smooth = FALSE)
print(p_gating_strat_Fn_total)

# Subset the data to the cells only
flowData_transformed_Fn_gated <- Subset(flowData_transformed_Fn, polyGateTotal)

## 4.2. Additional quality control ----
# Check singlets by showing height vs width of the same parameter
p_singlets_Av <- xyplot(`BL1-W`~`BL1-H`, data = flowData_transformed_Av_gated,
                        scales = list(y = list(limits = c(0, 9)),
                                      x = list(limits = c(6, 15))),
                        axis = axis.default, nbin = 125, main = "Singlet analysis Av (BL1-H vs BL1-W)", xlab = "BL1-H", ylab = "BL1-W",
                        par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_singlets_Av)

p_singlets_Av_BL1_AH <- xyplot(`BL1-A`~`BL1-H`, data = flowData_transformed_Av_gated,
                               scales = list(y = list(limits = c(5, 15)),
                                             x = list(limits = c(5, 15))),
                               axis = axis.default, nbin = 125, main = "Singlet analysis Av (BL1-H vs BL1-A)", xlab = "BL1-H", ylab = "BL1-A",
                               par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_singlets_Av_BL1_AH)

p_singlets_Av_SSC <- xyplot(`SSC-W`~`SSC-H`, data = flowData_transformed_Av_gated,
                            scales = list(y = list(limits = c(1, 9)),
                                          x = list(limits = c(4, 15))),
                            axis = axis.default, nbin = 125, main = "Singlet analysis Av (SSC-H vs SSC-W)", xlab = "SSC-H", ylab = "SSC-W",
                            par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_singlets_Av_SSC)

p_singlets_Av_FSC <- xyplot(`FSC-W`~`FSC-H`, data = flowData_transformed_Av_gated,
                            scales = list(y = list(limits = c(1, 9)),
                                          x = list(limits = c(4, 15))),
                            axis = axis.default, nbin = 125, main = "Singlet analysis Av (FSC-H vs FSC-W)", xlab = "FSC-H", ylab = "FSC-W",
                            par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_singlets_Av_FSC)

p_singlets_Fn <- xyplot(`BL1-W`~`BL1-H`, data = flowData_transformed_Fn_gated,
                        scales = list(y = list(limits = c(0, 9)),
                                      x = list(limits = c(6, 15))),
                        axis = axis.default, nbin = 125, main = "Singlet analysis Fn (BL1-H vs BL1-W)", xlab = "BL1-H", ylab = "BL1-W",
                        par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_singlets_Fn)

p_singlets_Fn_BL1_AH <- xyplot(`BL1-A`~`BL1-H`, data = flowData_transformed_Fn_gated,
                               scales = list(y = list(limits = c(5, 15)),
                                             x = list(limits = c(5, 15))),
                               axis = axis.default, nbin = 125, main = "Singlet analysis Fn (BL1-H vs BL1-A)", xlab = "BL1-H", ylab = "BL1-A",
                               par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_singlets_Fn_BL1_AH)

p_singlets_Fn_SSC <- xyplot(`SSC-W`~`SSC-H`, data = flowData_transformed_Fn_gated,
                            scales = list(y = list(limits = c(1, 9)),
                                          x = list(limits = c(4, 15))),
                            axis = axis.default, nbin = 125, main = "Singlet analysis Fn (SSC-H vs SSC-W)", xlab = "SSC-H", ylab = "SSC-W",
                            par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_singlets_Fn_SSC)

p_singlets_Fn_FSC <- xyplot(`FSC-W`~`FSC-H`, data = flowData_transformed_Fn_gated,
                            scales = list(y = list(limits = c(1, 9)),
                                          x = list(limits = c(4, 15))),
                            axis = axis.default, nbin = 125, main = "Singlet analysis Fn (FSC-H vs FSC-W)", xlab = "FSC-H", ylab = "FSC-W",
                            par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_singlets_Fn_FSC)

# Density plots
p_density_SSC_W_Av <- autoplot(flowData_transformed_Av_gated, "SSC-W")+
  labs(title = "Density plot SSC-W Av",
       x = "SSC-W",
       y = "Density")
print(p_density_SSC_W_Av)

p_density_BL1_H_Av <- autoplot(flowData_transformed_Av_gated, "BL1-H")+
  labs(title = "Density plot BL1-H Av",
       x = "BL1-H",
       y = "Density")
print(p_density_BL1_H_Av)

p_density_SSC_W_Fn <- autoplot(flowData_transformed_Fn_gated, "SSC-W")+
  labs(title = "Density plot SSC-W Fn",
       x = "SSC-W",
       y = "Density")
print(p_density_SSC_W_Fn)

p_density_BL1_H_Fn <- autoplot(flowData_transformed_Fn_gated, "BL1-H")+
  labs(title = "Density plot BL1-H Fn",
       x = "BL1-H",
       y = "Density")
print(p_density_BL1_H_Fn)

# Side scatter vs green fluorescence
xyplot(`BL1-H`~`SSC-H`, data = flowData_transformed_Av_gated,
       scales = list(y = list(limits = c(5, 15)),
                     x = list(limits = c(4, 16))),
       axis = axis.default, nbin = 125, main = "Side scatter vs green fluorescence (SSC-H vs BL1-H) gated Av", xlab = "SSC-H", ylab = "BL1-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL1-H`~`SSC-H`, data = flowData_transformed_Fn_gated,
       scales = list(y = list(limits = c(5, 15)),
                     x = list(limits = c(4, 16))),
       axis = axis.default, nbin = 125, main = "Side scatter vs green fluorescence (SSC-H vs BL1-H) gated Fn", xlab = "SSC-H", ylab = "BL1-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)


# 5. Phenotypic diversity analysis raw data ----

# Load custom plot_beta_fcm() function
source("/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/plot_beta_fcm_custom.R")

# All cells will be considered for the diversity analysis

## 5.1. Normalization of data ----
# Normalize over the highest value of all parameters (as this will result in a value between -1 and 1). Use same transformation for height and area (same range). Leave out width for now, as using different values for normalization will increase or decrease the weight of certain measured parameters. If one uses a Gaussian model, normalization is not needed, but then optimization of the band width is necessary.
# Av
summary_Av <- fsApply(x = flowData_transformed_Av_gated, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
max_Av = max(summary_Av[, "BL1-H"])
mytrans_Av <- function(x) x/max_Av

# Actual normalization of data
flowData_transformed_Av_norm <- transform(flowData_transformed_Av_gated,
                                          `FSC-H` = mytrans_Av(`FSC-H`),
                                          `SSC-H` = mytrans_Av(`SSC-H`),
                                          `BL1-H` = mytrans_Av(`BL1-H`),
                                          `BL3-H` = mytrans_Av(`BL3-H`),
                                          `FSC-A` = mytrans_Av(`FSC-A`),
                                          `SSC-A` = mytrans_Av(`SSC-A`),
                                          `BL1-A` = mytrans_Av(`BL1-A`),
                                          `BL3-A` = mytrans_Av(`BL3-A`),
                                          `FSC-W` = mytrans_Av(`FSC-W`),
                                          `SSC-W` = mytrans_Av(`SSC-W`),
                                          `BL1-W` = mytrans_Av(`BL1-W`),
                                          `BL3-W` = mytrans_Av(`BL3-W`))

# Select data without negative controls
flowData_transformed_Av_norm_noneg <- flowData_transformed_Av_norm[c(1:12, 16:75)]
metadata_Av_noneg <- metadata_Av[c(1:12, 16:75), ]
row.names(metadata_Av_noneg) <- NULL

# Fn
summary_Fn <- fsApply(x = flowData_transformed_Fn_gated, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
max_Fn = max(summary_Fn[, "BL1-H"])
mytrans_Fn <- function(x) x/max_Fn

# Actual normalization of data
flowData_transformed_Fn_norm <- transform(flowData_transformed_Fn_gated,
                                          `FSC-H` = mytrans_Fn(`FSC-H`),
                                          `SSC-H` = mytrans_Fn(`SSC-H`),
                                          `BL1-H` = mytrans_Fn(`BL1-H`),
                                          `BL3-H` = mytrans_Fn(`BL3-H`),
                                          `FSC-A` = mytrans_Fn(`FSC-A`),
                                          `SSC-A` = mytrans_Fn(`SSC-A`),
                                          `BL1-A` = mytrans_Fn(`BL1-A`),
                                          `BL3-A` = mytrans_Fn(`BL3-A`),
                                          `FSC-W` = mytrans_Fn(`FSC-W`),
                                          `SSC-W` = mytrans_Fn(`SSC-W`),
                                          `BL1-W` = mytrans_Fn(`BL1-W`),
                                          `BL3-W` = mytrans_Fn(`BL3-W`))

# Select data without negative controls
flowData_transformed_Fn_norm_noneg <- flowData_transformed_Fn_norm[c(1:12, 16:75)]
metadata_Fn_noneg <- metadata_Fn[c(1:12, 16:75), ]
row.names(metadata_Fn_noneg) <- NULL

## 5.2. Diversity analysis ----
param_div <- c("FSC-H", "SSC-H", "BL1-H", "BL3-H", "FSC-A", "SSC-A", "BL1-A", "BL3-A")

### 5.2.1. Av ----
#### 5.2.1.1. Fingerprinting ----
# Randomly resample to a fixed number of cells
flowData_transformed_Av_resample <- FCS_resample(flowData_transformed_Av_norm_noneg, sample = 10000, replace = TRUE)

# Make metadata file for resampled data
names_flowData_transformed_Av_resample <- sampleNames(flowData_transformed_Av_resample)
metadata_Av_resample <- metadata_Av_noneg[match(names_flowData_transformed_Av_resample, metadata_Av_noneg$Filename), ]
row.names(metadata_Av_resample) <- NULL

metadata_Av_resample_novel_compounds <- metadata_Av_resample[c(1, 5:72), ]
rownames(metadata_Av_resample_novel_compounds) <- NULL

# Calculating fingerprint with bandwidth = 0.01
fbasis_Av <- flowBasis(flowData_transformed_Av_resample, param = param_div, nbin = 128, 
                       bw = 0.01, normalize = function(x) x)

fbasis_Av_novel_compounds <- flowBasis(flowData_transformed_Av_resample[c(1, 5:72)], param = param_div, nbin = 128,
                                       bw = 0.01, normalize = function(x) x)

#### 5.2.1.2. Ordination ----
# Beta-diversity assessment of fingerprint (Bray-Curtis)
beta.div_Av <- beta_div_fcm(fbasis_Av, ord.type = "PCoA")
beta.div_Av_NMDS <- beta_div_fcm(fbasis_Av, ord.type = "NMDS")

# Beta-diversity assessment of fingerprints excluding cephalothin (Bray-Curtis)
beta.div_Av_novel_compounds <- beta_div_fcm(fbasis_Av_novel_compounds, ord.type = "PCoA")
beta.div_Av_novel_compounds_NMDS <- beta_div_fcm(fbasis_Av_novel_compounds, ord.type = "NMDS")

# Plot ordination
pbetadiv_Av_PCoA <- plot_beta_fcm(beta.div_Av, color = as.factor(metadata_Av_resample$MOA), shape = as.factor(metadata_Av_resample$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL)+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(pbetadiv_Av_PCoA)

pbetadiv_Av_NMDS <- plot_beta_fcm(beta.div_Av_NMDS, color = as.factor(metadata_Av_resample$MOA), shape = as.factor(metadata_Av_resample$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL, x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(pbetadiv_Av_NMDS)

# Plot ordination excluding cephalothin
pbetadiv_Av_novel_compounds_PCoA <- plot_beta_fcm(beta.div_Av_novel_compounds, color = as.factor(metadata_Av_resample_novel_compounds$MOA), shape = as.factor(metadata_Av_resample_novel_compounds$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL)+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:4, 6:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(pbetadiv_Av_novel_compounds_PCoA)

pbetadiv_Av_novel_compounds_NMDS <- plot_beta_fcm(beta.div_Av_novel_compounds_NMDS, color = as.factor(metadata_Av_resample_novel_compounds$MOA), shape = as.factor(metadata_Av_resample_novel_compounds$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL, x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:4, 6:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(pbetadiv_Av_novel_compounds_NMDS)

#### 5.2.1.3. Statistics ----
### Bray-Curtis dissimilarity matrix
# All samples
fbasis_Av_2 <- fbasis_Av@basis/apply(fbasis_Av@basis, 1, max)
fbasis_Av_2 <- base::round(fbasis_Av_2, 4)
dist_Av <- vegan::vegdist(fbasis_Av_2, method = "bray", binary = FALSE)

metadata_Av_dist <- metadata_Av_noneg
rownames(metadata_Av_dist) <- metadata_Av_dist$Filename

# Excluding cephalothin
fbasis_Av_3 <- fbasis_Av_2[c(1, 5:72),]
dist_Av_novel_compounds <- vegan::vegdist(fbasis_Av_3, method = "bray", binary = FALSE)

metadata_Av_dist_novel_compounds <- metadata_Av_noneg[c(1, 5:72),]
rownames(metadata_Av_dist_novel_compounds) <- metadata_Av_dist_novel_compounds$Filename

### ANOSIM
permutations = 10000

# All treatments
ano_Av <- vegan::anosim(dist_Av, metadata_Av_dist$MOA, distance = "bray", permutations = permutations)
ano_Av_export <- data.frame("ANOSIM Statistic R" = ano_Av$statistic, "P" = ano_Av$signif, "permutations" = ano_Av$permutations)
write.csv2(file = "ANOSIM_Av_all_treatments.csv", ano_Av_export) # Export results as csv file
plot_dissimilarity_anosim_Av <- plot(ano_Av, xaxt = "n", main = "Dissimilarity ranks Av (ANOSIM)", xlab = "Class", ylab = "Dissimilarity rank")
axis(1, at = 1:11, labels = c("Between", "CWS", "DNARepl", "DNATrans", "FAM", "Heat", "MD", "Control", "PS-30SInh", "PS-50SInh", "PS-tRNAInt"))

# Excluding cephalothin
ano_Av_novel_compounds <- vegan::anosim(dist_Av_novel_compounds, metadata_Av_dist_novel_compounds$MOA, distance = "bray", permutations = permutations)
ano_Av_novel_compounds_export <- data.frame("ANOSIM Statistic R" = ano_Av_novel_compounds$statistic, "P" = ano_Av_novel_compounds$signif, "permutations" = ano_Av_novel_compounds$permutations)
write.csv2(file = "ANOSIM_Av_novel_compounds.csv", ano_Av_novel_compounds_export) # Export results as csv file
plot_dissimilarity_anosim_Av_novel_compounds <- plot(ano_Av_novel_compounds, xaxt = "n", main = "Dissimilarity ranks Av (ANOSIM)", xlab = "Class", ylab = "Dissimilarity rank")
axis(1, at = 1:11, labels = c("Between", "CWS", "DNARepl", "DNATrans", "FAM", "Heat", "MD", "Control", "PS-30SInh", "PS-50SInh", "PS-tRNAInt"))

## Pairwise comparison different MOA
# All treatments
pairwise_anosim_Av <- numeric()

# CellWallSynthesis vs Control
dist_Av_CWS_PC <- dist_Av %>% 
  as.matrix() %>% 
  .[c(1:18), c(1:18)] %>% 
  as.dist()

metadata_Av_dist_CWS_PC <- metadata_Av_dist[c(1:18), ]
ano_Av_CWS_PC <- vegan::anosim(dist_Av_CWS_PC, metadata_Av_dist_CWS_PC$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["CellWallSynthesis_Control"] <- ano_Av_CWS_PC$signif

# CellWallSynthesis vs DNAReplication
dist_Av_CWS_DR <- dist_Av %>% 
  as.matrix() %>% 
  .[c(1:12, 19:27), c(1:12, 19:27)] %>% 
  as.dist()

metadata_Av_dist_CWS_DR <- metadata_Av_dist[c(1:12, 19:27), ]
ano_Av_CWS_DR <- vegan::anosim(dist_Av_CWS_DR, metadata_Av_dist_CWS_DR$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["CellWallSynthesis_DNAReplication"] <- ano_Av_CWS_DR$signif

# CellWallSynthesis vs DNATranscription
dist_Av_CWS_DT <- dist_Av %>% 
  as.matrix() %>% 
  .[c(1:12, 28:33), c(1:12, 28:33)] %>% 
  as.dist()

metadata_Av_dist_CWS_DT <- metadata_Av_dist[c(1:12, 28:33), ]
ano_Av_CWS_DT <- vegan::anosim(dist_Av_CWS_DT, metadata_Av_dist_CWS_DT$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["CellWallSynthesis_DNATranscription"] <- ano_Av_CWS_DT$signif

# CellWallSynthesis vs FolicAcidMetabolism
dist_Av_CWS_FAM <- dist_Av %>% 
  as.matrix() %>% 
  .[c(1:12, 34:36), c(1:12, 34:36)] %>% 
  as.dist()

metadata_Av_dist_CWS_FAM <- metadata_Av_dist[c(1:12, 34:36), ]
ano_Av_CWS_FAM <- vegan::anosim(dist_Av_CWS_FAM, metadata_Av_dist_CWS_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["CellWallSynthesis_FolicAcidMetabolism"] <- ano_Av_CWS_FAM$signif

# CellWallSynthesis vs Heat
dist_Av_CWS_H <- dist_Av %>% 
  as.matrix() %>% 
  .[c(1:12, 37:39), c(1:12, 37:39)] %>% 
  as.dist()

metadata_Av_dist_CWS_H <- metadata_Av_dist[c(1:12, 37:39), ]
ano_Av_CWS_H <- vegan::anosim(dist_Av_CWS_H, metadata_Av_dist_CWS_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["CellWallSynthesis_Heat"] <- ano_Av_CWS_H$signif

# CellWallSynthesis vs MembraneDisruption
dist_Av_CWS_MD <- dist_Av %>% 
  as.matrix() %>% 
  .[c(1:12, 40:51), c(1:12, 40:51)] %>% 
  as.dist()

metadata_Av_dist_CWS_MD <- metadata_Av_dist[c(1:12, 40:51), ]
ano_Av_CWS_MD <- vegan::anosim(dist_Av_CWS_MD, metadata_Av_dist_CWS_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["CellWallSynthesis_MembraneDisruption"] <- ano_Av_CWS_MD$signif

# CellWallSynthesis vs ProteinSynthesis30SInhibition
dist_Av_CWS_PS30SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(1:12, 52:60), c(1:12, 52:60)] %>% 
  as.dist()

metadata_Av_dist_CWS_PS30SI <- metadata_Av_dist[c(1:12, 52:60), ]
ano_Av_CWS_PS30SI <- vegan::anosim(dist_Av_CWS_PS30SI, metadata_Av_dist_CWS_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["CellWallSynthesis_ProteinSynthesis30SInhibition"] <- ano_Av_CWS_PS30SI$signif

# CellWallSynthesis vs ProteinSynthesis50SInhibition
dist_Av_CWS_PS50SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(1:12, 61:69), c(1:12, 61:69)] %>% 
  as.dist()

metadata_Av_dist_CWS_PS50SI <- metadata_Av_dist[c(1:12, 61:69), ]
ano_Av_CWS_PS50SI <- vegan::anosim(dist_Av_CWS_PS50SI, metadata_Av_dist_CWS_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["CellWallSynthesis_ProteinSynthesis50SInhibition"] <- ano_Av_CWS_PS50SI$signif

# CellWallSynthesis vs ProteinSynthesistRNAInterference
dist_Av_CWS_PStRI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(1:12, 70:72), c(1:12, 70:72)] %>% 
  as.dist()

metadata_Av_dist_CWS_PStRI <- metadata_Av_dist[c(1:12, 70:72), ]
ano_Av_CWS_PStRI <- vegan::anosim(dist_Av_CWS_PStRI, metadata_Av_dist_CWS_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["CellWallSynthesis_ProteinSynthesistRNAInterference"] <- ano_Av_CWS_PStRI$signif

# Control vs DNAReplication
dist_Av_PC_DR <- dist_Av %>% 
  as.matrix() %>% 
  .[c(13:18, 19:27), c(13:18, 19:27)] %>% 
  as.dist()

metadata_Av_dist_PC_DR <- metadata_Av_dist[c(13:18, 19:27), ]
ano_Av_PC_DR <- vegan::anosim(dist_Av_PC_DR, metadata_Av_dist_PC_DR$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["Control_DNAReplication"] <- ano_Av_PC_DR$signif

# Control vs DNATranscription
dist_Av_PC_DT <- dist_Av %>% 
  as.matrix() %>% 
  .[c(13:18, 28:33), c(13:18, 28:33)] %>% 
  as.dist()

metadata_Av_dist_PC_DT <- metadata_Av_dist[c(13:18, 28:33), ]
ano_Av_PC_DT <- vegan::anosim(dist_Av_PC_DT, metadata_Av_dist_PC_DT$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["Control_DNATranscription"] <- ano_Av_PC_DT$signif

# Control vs FolicAcidMetabolism
dist_Av_PC_FAM <- dist_Av %>% 
  as.matrix() %>% 
  .[c(13:18, 34:36), c(13:18, 34:36)] %>% 
  as.dist()

metadata_Av_dist_PC_FAM <- metadata_Av_dist[c(13:18, 34:36), ]
ano_Av_PC_FAM <- vegan::anosim(dist_Av_PC_FAM, metadata_Av_dist_PC_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["Control_FolicAcidMetabolism"] <- ano_Av_PC_FAM$signif

# Control vs Heat
dist_Av_PC_H <- dist_Av %>% 
  as.matrix() %>% 
  .[c(13:18, 37:39), c(13:18, 37:39)] %>% 
  as.dist()

metadata_Av_dist_PC_H <- metadata_Av_dist[c(13:18, 37:39), ]
ano_Av_PC_H <- vegan::anosim(dist_Av_PC_H, metadata_Av_dist_PC_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["Control_Heat"] <- ano_Av_PC_H$signif

# Control vs MembraneDisruption
dist_Av_PC_MD <- dist_Av %>% 
  as.matrix() %>% 
  .[c(13:18, 40:51), c(13:18, 40:51)] %>% 
  as.dist()

metadata_Av_dist_PC_MD <- metadata_Av_dist[c(13:18, 40:51), ]
ano_Av_PC_MD <- vegan::anosim(dist_Av_PC_MD, metadata_Av_dist_PC_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["Control_MembraneDisruption"] <- ano_Av_PC_MD$signif

# Control vs ProteinSynthesis30SInhibition
dist_Av_PC_PS30SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(13:18, 52:60), c(13:18, 52:60)] %>% 
  as.dist()

metadata_Av_dist_PC_PS30SI <- metadata_Av_dist[c(13:18, 52:60), ]
ano_Av_PC_PS30SI <- vegan::anosim(dist_Av_PC_PS30SI, metadata_Av_dist_PC_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["Control_ProteinSynthesis30SInhibition"] <- ano_Av_PC_PS30SI$signif

# Control vs ProteinSynthesis50SInhibition
dist_Av_PC_PS50SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(13:18, 61:69), c(13:18, 61:69)] %>% 
  as.dist()

metadata_Av_dist_PC_PS50SI <- metadata_Av_dist[c(13:18, 61:69), ]
ano_Av_PC_PS50SI <- vegan::anosim(dist_Av_PC_PS50SI, metadata_Av_dist_PC_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["Control_ProteinSynthesis50SInhibition"] <- ano_Av_PC_PS50SI$signif

# Control vs ProteinSynthesistRNAInterference
dist_Av_PC_PStRI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(13:18, 70:72), c(13:18, 70:72)] %>% 
  as.dist()

metadata_Av_dist_PC_PStRI <- metadata_Av_dist[c(13:18, 70:72), ]
ano_Av_PC_PStRI <- vegan::anosim(dist_Av_PC_PStRI, metadata_Av_dist_PC_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["Control_ProteinSynthesistRNAInterference"] <- ano_Av_PC_PStRI$signif

# DNAReplication vs DNATranscription
dist_Av_DR_DT <- dist_Av %>% 
  as.matrix() %>% 
  .[c(19:27, 28:33), c(19:27, 28:33)] %>% 
  as.dist()

metadata_Av_dist_DR_DT <- metadata_Av_dist[c(19:27, 28:33), ]
ano_Av_DR_DT <- vegan::anosim(dist_Av_DR_DT, metadata_Av_dist_DR_DT$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["DNAReplication_DNATranscription"] <- ano_Av_DR_DT$signif

# DNAReplication vs FolicAcidMetabolism
dist_Av_DR_FAM <- dist_Av %>% 
  as.matrix() %>% 
  .[c(19:27, 34:36), c(19:27, 34:36)] %>% 
  as.dist()

metadata_Av_dist_DR_FAM <- metadata_Av_dist[c(19:27, 34:36), ]
ano_Av_DR_FAM <- vegan::anosim(dist_Av_DR_FAM, metadata_Av_dist_DR_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["DNAReplication_FolicAcidMetabolism"] <- ano_Av_DR_FAM$signif

# DNAReplication vs Heat
dist_Av_DR_H <- dist_Av %>% 
  as.matrix() %>% 
  .[c(19:27, 37:39), c(19:27, 37:39)] %>% 
  as.dist()

metadata_Av_dist_DR_H <- metadata_Av_dist[c(19:27, 37:39), ]
ano_Av_DR_H <- vegan::anosim(dist_Av_DR_H, metadata_Av_dist_DR_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["DNAReplication_Heat"] <- ano_Av_DR_H$signif

# DNAReplication vs MembraneDisruption
dist_Av_DR_MD <- dist_Av %>% 
  as.matrix() %>% 
  .[c(19:27, 40:51), c(19:27, 40:51)] %>% 
  as.dist()

metadata_Av_dist_DR_MD <- metadata_Av_dist[c(19:27, 40:51), ]
ano_Av_DR_MD <- vegan::anosim(dist_Av_DR_MD, metadata_Av_dist_DR_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["DNAReplication_MembraneDisruption"] <- ano_Av_DR_MD$signif

# DNAReplication vs ProteinSynthesis30SInhibition
dist_Av_DR_PS30SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(19:27, 52:60), c(19:27, 52:60)] %>% 
  as.dist()

metadata_Av_dist_DR_PS30SI <- metadata_Av_dist[c(19:27, 52:60), ]
ano_Av_DR_PS30SI <- vegan::anosim(dist_Av_DR_PS30SI, metadata_Av_dist_DR_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["DNAReplication_ProteinSynthesis30SInhibition"] <- ano_Av_DR_PS30SI$signif

# DNAReplication vs ProteinSynthesis50SInhibition
dist_Av_DR_PS50SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(19:27, 61:69), c(19:27, 61:69)] %>% 
  as.dist()

metadata_Av_dist_DR_PS50SI <- metadata_Av_dist[c(19:27, 61:69), ]
ano_Av_DR_PS50SI <- vegan::anosim(dist_Av_DR_PS50SI, metadata_Av_dist_DR_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["DNAReplication_ProteinSynthesis50SInhibition"] <- ano_Av_DR_PS50SI$signif

# DNAReplication vs ProteinSynthesistRNAInterference
dist_Av_DR_PStRI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(19:27, 70:72), c(19:27, 70:72)] %>% 
  as.dist()

metadata_Av_dist_DR_PStRI <- metadata_Av_dist[c(19:27, 70:72), ]
ano_Av_DR_PStRI <- vegan::anosim(dist_Av_DR_PStRI, metadata_Av_dist_DR_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["DNAReplication_ProteinSynthesistRNAInterference"] <- ano_Av_DR_PStRI$signif

# DNATranscription vs FolicAcidMetabolism
dist_Av_DT_FAM <- dist_Av %>% 
  as.matrix() %>% 
  .[c(28:33, 34:36), c(28:33, 34:36)] %>% 
  as.dist()

metadata_Av_dist_DT_FAM <- metadata_Av_dist[c(28:33, 34:36), ]
ano_Av_DT_FAM <- vegan::anosim(dist_Av_DT_FAM, metadata_Av_dist_DT_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["DNATranscription_FolicAcidMetabolism"] <- ano_Av_DT_FAM$signif

# DNATranscription vs Heat
dist_Av_DT_H <- dist_Av %>% 
  as.matrix() %>% 
  .[c(28:33, 37:39), c(28:33, 37:39)] %>% 
  as.dist()

metadata_Av_dist_DT_H <- metadata_Av_dist[c(28:33, 37:39), ]
ano_Av_DT_H <- vegan::anosim(dist_Av_DT_H, metadata_Av_dist_DT_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["DNATranscription_Heat"] <- ano_Av_DT_H$signif

# DNATranscription vs MembraneDisruption
dist_Av_DT_MD <- dist_Av %>% 
  as.matrix() %>% 
  .[c(28:33, 40:51), c(28:33, 40:51)] %>% 
  as.dist()

metadata_Av_dist_DT_MD <- metadata_Av_dist[c(28:33, 40:51), ]
ano_Av_DT_MD <- vegan::anosim(dist_Av_DT_MD, metadata_Av_dist_DT_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["DNATranscription_MembraneDisruption"] <- ano_Av_DT_MD$signif

# DNATranscription vs ProteinSynthesis30SInhibition
dist_Av_DT_PS30SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(28:33, 52:60), c(28:33, 52:60)] %>% 
  as.dist()

metadata_Av_dist_DT_PS30SI <- metadata_Av_dist[c(28:33, 52:60), ]
ano_Av_DT_PS30SI <- vegan::anosim(dist_Av_DT_PS30SI, metadata_Av_dist_DT_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["DNATranscription_ProteinSynthesis30SInhibition"] <- ano_Av_DT_PS30SI$signif

# DNATranscription vs ProteinSynthesis50SInhibition
dist_Av_DT_PS50SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(28:33, 61:69), c(28:33, 61:69)] %>% 
  as.dist()

metadata_Av_dist_DT_PS50SI <- metadata_Av_dist[c(28:33, 61:69), ]
ano_Av_DT_PS50SI <- vegan::anosim(dist_Av_DT_PS50SI, metadata_Av_dist_DT_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["DNATranscription_ProteinSynthesis50SInhibition"] <- ano_Av_DT_PS50SI$signif

# DNATranscription vs ProteinSynthesistRNAInterference
dist_Av_DT_PStRI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(28:33, 70:72), c(28:33, 70:72)] %>% 
  as.dist()

metadata_Av_dist_DT_PStRI <- metadata_Av_dist[c(28:33, 70:72), ]
ano_Av_DT_PStRI <- vegan::anosim(dist_Av_DT_PStRI, metadata_Av_dist_DT_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["DNATranscription_ProteinSynthesistRNAInterference"] <- ano_Av_DT_PStRI$signif

# FolicAcidMetabolism vs Heat
dist_Av_FAM_H <- dist_Av %>% 
  as.matrix() %>% 
  .[c(34:36, 37:39), c(34:36, 37:39)] %>% 
  as.dist()

metadata_Av_dist_FAM_H <- metadata_Av_dist[c(34:36, 37:39), ]
ano_Av_FAM_H <- vegan::anosim(dist_Av_FAM_H, metadata_Av_dist_FAM_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["FolicAcidMetabolism_Heat"] <- ano_Av_FAM_H$signif

# FolicAcidMetabolism vs MembraneDisruption
dist_Av_FAM_MD <- dist_Av %>% 
  as.matrix() %>% 
  .[c(34:36, 40:51), c(34:36, 40:51)] %>% 
  as.dist()

metadata_Av_dist_FAM_MD <- metadata_Av_dist[c(34:36, 40:51), ]
ano_Av_FAM_MD <- vegan::anosim(dist_Av_FAM_MD, metadata_Av_dist_FAM_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["FolicAcidMetabolism_MembraneDisruption"] <- ano_Av_FAM_MD$signif

# FolicAcidMetabolism vs ProteinSynthesis30SInhibition
dist_Av_FAM_PS30SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(34:36, 52:60), c(34:36, 52:60)] %>% 
  as.dist()

metadata_Av_dist_FAM_PS30SI <- metadata_Av_dist[c(34:36, 52:60), ]
ano_Av_FAM_PS30SI <- vegan::anosim(dist_Av_FAM_PS30SI, metadata_Av_dist_FAM_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["FolicAcidMetabolism_ProteinSynthesis30SInhibition"] <- ano_Av_FAM_PS30SI$signif

# FolicAcidMetabolism vs ProteinSynthesis50SInhibition
dist_Av_FAM_PS50SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(34:36, 61:69), c(34:36, 61:69)] %>% 
  as.dist()

metadata_Av_dist_FAM_PS50SI <- metadata_Av_dist[c(34:36, 61:69), ]
ano_Av_FAM_PS50SI <- vegan::anosim(dist_Av_FAM_PS50SI, metadata_Av_dist_FAM_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["FolicAcidMetabolism_ProteinSynthesis50SInhibition"] <- ano_Av_FAM_PS50SI$signif

# FolicAcidMetabolism vs ProteinSynthesistRNAInterference
dist_Av_FAM_PStRI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(34:36, 70:72), c(34:36, 70:72)] %>% 
  as.dist()

metadata_Av_dist_FAM_PStRI <- metadata_Av_dist[c(34:36, 70:72), ]
ano_Av_FAM_PStRI <- vegan::anosim(dist_Av_FAM_PStRI, metadata_Av_dist_FAM_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["FolicAcidMetabolism_ProteinSynthesistRNAInterference"] <- ano_Av_FAM_PStRI$signif

# Heat vs MembraneDisruption
dist_Av_H_MD <- dist_Av %>% 
  as.matrix() %>% 
  .[c(37:39, 40:51), c(37:39, 40:51)] %>% 
  as.dist()

metadata_Av_dist_H_MD <- metadata_Av_dist[c(37:39, 40:51), ]
ano_Av_H_MD <- vegan::anosim(dist_Av_H_MD, metadata_Av_dist_H_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["Heat_MembraneDisruption"] <- ano_Av_H_MD$signif

# Heat vs ProteinSynthesis30SInhibition
dist_Av_H_PS30SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(37:39, 52:60), c(37:39, 52:60)] %>% 
  as.dist()

metadata_Av_dist_H_PS30SI <- metadata_Av_dist[c(37:39, 52:60), ]
ano_Av_H_PS30SI <- vegan::anosim(dist_Av_H_PS30SI, metadata_Av_dist_H_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["Heat_ProteinSynthesis30SInhibition"] <- ano_Av_H_PS30SI$signif

# Heat vs ProteinSynthesis50SInhibition
dist_Av_H_PS50SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(37:39, 61:69), c(37:39, 61:69)] %>% 
  as.dist()

metadata_Av_dist_H_PS50SI <- metadata_Av_dist[c(37:39, 61:69), ]
ano_Av_H_PS50SI <- vegan::anosim(dist_Av_H_PS50SI, metadata_Av_dist_H_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["Heat_ProteinSynthesis50SInhibition"] <- ano_Av_H_PS50SI$signif

# Heat vs ProteinSynthesistRNAInterference
dist_Av_H_PStRI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(37:39, 70:72), c(37:39, 70:72)] %>% 
  as.dist()

metadata_Av_dist_H_PStRI <- metadata_Av_dist[c(37:39, 70:72), ]
ano_Av_H_PStRI <- vegan::anosim(dist_Av_H_PStRI, metadata_Av_dist_H_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["Heat_ProteinSynthesistRNAInterference"] <- ano_Av_H_PStRI$signif

# MembraneDisruption vs ProteinSynthesis30SInhibition
dist_Av_MD_PS30SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(40:51, 52:60), c(40:51, 52:60)] %>% 
  as.dist()

metadata_Av_dist_MD_PS30SI <- metadata_Av_dist[c(40:51, 52:60), ]
ano_Av_MD_PS30SI <- vegan::anosim(dist_Av_MD_PS30SI, metadata_Av_dist_MD_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["MembraneDisruption_ProteinSynthesis30SInhibition"] <- ano_Av_MD_PS30SI$signif

# MembraneDisruption vs ProteinSynthesis50SInhibition
dist_Av_MD_PS50SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(40:51, 61:69), c(40:51, 61:69)] %>% 
  as.dist()

metadata_Av_dist_MD_PS50SI <- metadata_Av_dist[c(40:51, 61:69), ]
ano_Av_MD_PS50SI <- vegan::anosim(dist_Av_MD_PS50SI, metadata_Av_dist_MD_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["MembraneDisruption_ProteinSynthesis50SInhibition"] <- ano_Av_MD_PS50SI$signif

# MembraneDisruption vs ProteinSynthesistRNAInterference
dist_Av_MD_PStRI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(40:51, 70:72), c(40:51, 70:72)] %>% 
  as.dist()

metadata_Av_dist_MD_PStRI <- metadata_Av_dist[c(40:51, 70:72), ]
ano_Av_MD_PStRI <- vegan::anosim(dist_Av_MD_PStRI, metadata_Av_dist_MD_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["MembraneDisruption_ProteinSynthesistRNAInterference"] <- ano_Av_MD_PStRI$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesis50SInhibition
dist_Av_PS30SI_PS50SI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(52:60, 61:69), c(52:60, 61:69)] %>% 
  as.dist()

metadata_Av_dist_PS30SI_PS50SI <- metadata_Av_dist[c(52:60, 61:69), ]
ano_Av_PS30SI_PS50SI <- vegan::anosim(dist_Av_PS30SI_PS50SI, metadata_Av_dist_PS30SI_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["ProteinSynthesis30SInhibition_ProteinSynthesis50SInhibition"] <- ano_Av_PS30SI_PS50SI$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesistRNAInterference
dist_Av_PS30SI_PStRI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(52:60, 70:72), c(52:60, 70:72)] %>% 
  as.dist()

metadata_Av_dist_PS30SI_PStRI <- metadata_Av_dist[c(52:60, 70:72), ]
ano_Av_PS30SI_PStRI <- vegan::anosim(dist_Av_PS30SI_PStRI, metadata_Av_dist_PS30SI_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["ProteinSynthesis30SInhibition_ProteinSynthesistRNAInterference"] <- ano_Av_PS30SI_PStRI$signif

# ProteinSynthesis50SInhibition vs ProteinSynthesistRNAInterference
dist_Av_PS50SI_PStRI <- dist_Av %>% 
  as.matrix() %>% 
  .[c(61:69, 70:72), c(61:69, 70:72)] %>% 
  as.dist()

metadata_Av_dist_PS50SI_PStRI <- metadata_Av_dist[c(61:69, 70:72), ]
ano_Av_PS50SI_PStRI <- vegan::anosim(dist_Av_PS50SI_PStRI, metadata_Av_dist_PS50SI_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av["ProteinSynthesis50SInhibition_ProteinSynthesistRNAInterference"] <- ano_Av_PS50SI_PStRI$signif

# Correction for pairwise comparisons
pairwise_anosim_Av_bonferroni <- p.adjust(pairwise_anosim_Av, method = "bonferroni")
pairwise_anosim_Av_BH <- p.adjust(pairwise_anosim_Av, method = "BH")

pairwise_anosim_Av_df <- as.data.frame(pairwise_anosim_Av)
pairwise_anosim_Av_bonferroni_df <- as.data.frame(pairwise_anosim_Av_bonferroni)
pairwise_anosim_Av_BH_df <- as.data.frame(pairwise_anosim_Av_BH)

pairwise_anosim_Av_export <- merge(pairwise_anosim_Av_df, pairwise_anosim_Av_bonferroni_df, by = "row.names") %>% 
  column_to_rownames(var = "Row.names") %>% 
  merge(pairwise_anosim_Av_BH_df, by = "row.names") %>% 
  rename("classes" = "Row.names", "p" = "pairwise_anosim_Av", "p.adjust.bonferroni" = "pairwise_anosim_Av_bonferroni", "p.adjust.BH" = "pairwise_anosim_Av_BH")
write.csv2(file = "ANOSIM_Av_all_treatments_pairwise.csv", pairwise_anosim_Av_export) # Export results as csv file

# Novel compounds (exlcuding cephalothin)
pairwise_anosim_Av_novel_compounds <- numeric()

# Indices MOA classes
ind_CWS <- c(1:9)
ind_PC <- c(10:15)
ind_DR <- c(16:24)
ind_DT <- c(25:30)
ind_FAM <- c(31:33)
ind_H <- c(34:36)
ind_MD <- c(37:48)
ind_PS30SI <- c(49:57)
ind_PS50SI <- c(58:66)
ind_PStRI <- c(67:69)

# CellWallSynthesis vs Control
dist_Av_novel_compounds_CWS_PC <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PC), c(ind_CWS, ind_PC)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_CWS_PC <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_PC), ]
ano_Av_novel_compounds_CWS_PC <- vegan::anosim(dist_Av_novel_compounds_CWS_PC, metadata_Av_dist_novel_compounds_CWS_PC$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["CellWallSynthesis_Control"] <- ano_Av_novel_compounds_CWS_PC$signif

# CellWallSynthesis vs DNAReplication
dist_Av_novel_compounds_CWS_DR <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_DR), c(ind_CWS, ind_DR)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_CWS_DR <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_DR), ]
ano_Av_novel_compounds_CWS_DR <- vegan::anosim(dist_Av_novel_compounds_CWS_DR, metadata_Av_dist_novel_compounds_CWS_DR$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["CellWallSynthesis_DNAReplication"] <- ano_Av_novel_compounds_CWS_DR$signif

# CellWallSynthesis vs DNATranscription
dist_Av_novel_compounds_CWS_DT <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_DT), c(ind_CWS, ind_DT)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_CWS_DT <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_DT), ]
ano_Av_novel_compounds_CWS_DT <- vegan::anosim(dist_Av_novel_compounds_CWS_DT, metadata_Av_dist_novel_compounds_CWS_DT$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["CellWallSynthesis_DNATranscription"] <- ano_Av_novel_compounds_CWS_DT$signif

# CellWallSynthesis vs FolicAcidMetabolism
dist_Av_novel_compounds_CWS_FAM <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_FAM), c(ind_CWS, ind_FAM)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_CWS_FAM <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_FAM), ]
ano_Av_novel_compounds_CWS_FAM <- vegan::anosim(dist_Av_novel_compounds_CWS_FAM, metadata_Av_dist_novel_compounds_CWS_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["CellWallSynthesis_FolicAcidMetabolism"] <- ano_Av_novel_compounds_CWS_FAM$signif

# CellWallSynthesis vs Heat
dist_Av_novel_compounds_CWS_H <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_H), c(ind_CWS, ind_H)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_CWS_H <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_H), ]
ano_Av_novel_compounds_CWS_H <- vegan::anosim(dist_Av_novel_compounds_CWS_H, metadata_Av_dist_novel_compounds_CWS_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["CellWallSynthesis_Heat"] <- ano_Av_novel_compounds_CWS_H$signif

# CellWallSynthesis vs MembraneDisruption
dist_Av_novel_compounds_CWS_MD <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_MD), c(ind_CWS, ind_MD)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_CWS_MD <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_MD), ]
ano_Av_novel_compounds_CWS_MD <- vegan::anosim(dist_Av_novel_compounds_CWS_MD, metadata_Av_dist_novel_compounds_CWS_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["CellWallSynthesis_MembraneDisruption"] <- ano_Av_novel_compounds_CWS_MD$signif

# CellWallSynthesis vs ProteinSynthesis30SInhibition
dist_Av_novel_compounds_CWS_PS30SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PS30SI), c(ind_CWS, ind_PS30SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_CWS_PS30SI <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_PS30SI), ]
ano_Av_novel_compounds_CWS_PS30SI <- vegan::anosim(dist_Av_novel_compounds_CWS_PS30SI, metadata_Av_dist_novel_compounds_CWS_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["CellWallSynthesis_ProteinSynthesis30SInhibition"] <- ano_Av_novel_compounds_CWS_PS30SI$signif

# CellWallSynthesis vs ProteinSynthesis50SInhibition
dist_Av_novel_compounds_CWS_PS50SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PS50SI), c(ind_CWS, ind_PS50SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_CWS_PS50SI <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_PS50SI), ]
ano_Av_novel_compounds_CWS_PS50SI <- vegan::anosim(dist_Av_novel_compounds_CWS_PS50SI, metadata_Av_dist_novel_compounds_CWS_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["CellWallSynthesis_ProteinSynthesis50SInhibition"] <- ano_Av_novel_compounds_CWS_PS50SI$signif

# CellWallSynthesis vs ProteinSynthesistRNAInterference
dist_Av_novel_compounds_CWS_PStRI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PStRI), c(ind_CWS, ind_PStRI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_CWS_PStRI <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_PStRI), ]
ano_Av_novel_compounds_CWS_PStRI <- vegan::anosim(dist_Av_novel_compounds_CWS_PStRI, metadata_Av_dist_novel_compounds_CWS_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["CellWallSynthesis_ProteinSynthesistRNAInterference"] <- ano_Av_novel_compounds_CWS_PStRI$signif

# Control vs DNAReplication
dist_Av_novel_compounds_PC_DR <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_DR), c(ind_PC, ind_DR)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_PC_DR <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_DR), ]
ano_Av_novel_compounds_PC_DR <- vegan::anosim(dist_Av_novel_compounds_PC_DR, metadata_Av_dist_novel_compounds_PC_DR$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["Control_DNAReplication"] <- ano_Av_novel_compounds_PC_DR$signif

# Control vs DNATranscription
dist_Av_novel_compounds_PC_DT <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_DT), c(ind_PC, ind_DT)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_PC_DT <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_DT), ]
ano_Av_novel_compounds_PC_DT <- vegan::anosim(dist_Av_novel_compounds_PC_DT, metadata_Av_dist_novel_compounds_PC_DT$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["Control_DNATranscription"] <- ano_Av_novel_compounds_PC_DT$signif

# Control vs FolicAcidMetabolism
dist_Av_novel_compounds_PC_FAM <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_FAM), c(ind_PC, ind_FAM)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_PC_FAM <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_FAM), ]
ano_Av_novel_compounds_PC_FAM <- vegan::anosim(dist_Av_novel_compounds_PC_FAM, metadata_Av_dist_novel_compounds_PC_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["Control_FolicAcidMetabolism"] <- ano_Av_novel_compounds_PC_FAM$signif

# Control vs Heat
dist_Av_novel_compounds_PC_H <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_H), c(ind_PC, ind_H)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_PC_H <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_H), ]
ano_Av_novel_compounds_PC_H <- vegan::anosim(dist_Av_novel_compounds_PC_H, metadata_Av_dist_novel_compounds_PC_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["Control_Heat"] <- ano_Av_novel_compounds_PC_H$signif

# Control vs MembraneDisruption
dist_Av_novel_compounds_PC_MD <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_MD), c(ind_PC, ind_MD)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_PC_MD <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_MD), ]
ano_Av_novel_compounds_PC_MD <- vegan::anosim(dist_Av_novel_compounds_PC_MD, metadata_Av_dist_novel_compounds_PC_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["Control_MembraneDisruption"] <- ano_Av_novel_compounds_PC_MD$signif

# Control vs ProteinSynthesis30SInhibition
dist_Av_novel_compounds_PC_PS30SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_PS30SI), c(ind_PC, ind_PS30SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_PC_PS30SI <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_PS30SI), ]
ano_Av_novel_compounds_PC_PS30SI <- vegan::anosim(dist_Av_novel_compounds_PC_PS30SI, metadata_Av_dist_novel_compounds_PC_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["Control_ProteinSynthesis30SInhibition"] <- ano_Av_novel_compounds_PC_PS30SI$signif

# Control vs ProteinSynthesis50SInhibition
dist_Av_novel_compounds_PC_PS50SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_PS50SI), c(ind_PC, ind_PS50SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_PC_PS50SI <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_PS50SI), ]
ano_Av_novel_compounds_PC_PS50SI <- vegan::anosim(dist_Av_novel_compounds_PC_PS50SI, metadata_Av_dist_novel_compounds_PC_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["Control_ProteinSynthesis50SInhibition"] <- ano_Av_novel_compounds_PC_PS50SI$signif

# Control vs ProteinSynthesistRNAInterference
dist_Av_novel_compounds_PC_PStRI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_PStRI), c(ind_PC, ind_PStRI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_PC_PStRI <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_PStRI), ]
ano_Av_novel_compounds_PC_PStRI <- vegan::anosim(dist_Av_novel_compounds_PC_PStRI, metadata_Av_dist_novel_compounds_PC_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["Control_ProteinSynthesistRNAInterference"] <- ano_Av_novel_compounds_PC_PStRI$signif

# DNAReplication vs DNATranscription
dist_Av_novel_compounds_DR_DT <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_DT), c(ind_DR, ind_DT)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_DR_DT <- metadata_Av_dist_novel_compounds[c(ind_DR, ind_DT), ]
ano_Av_novel_compounds_DR_DT <- vegan::anosim(dist_Av_novel_compounds_DR_DT, metadata_Av_dist_novel_compounds_DR_DT$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["DNAReplication_DNATranscription"] <- ano_Av_novel_compounds_DR_DT$signif

# DNAReplication vs FolicAcidMetabolism
dist_Av_novel_compounds_DR_FAM <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_FAM), c(ind_DR, ind_FAM)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_DR_FAM <- metadata_Av_dist_novel_compounds[c(ind_DR, ind_FAM), ]
ano_Av_novel_compounds_DR_FAM <- vegan::anosim(dist_Av_novel_compounds_DR_FAM, metadata_Av_dist_novel_compounds_DR_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["DNAReplication_FolicAcidMetabolism"] <- ano_Av_novel_compounds_DR_FAM$signif

# DNAReplication vs Heat
dist_Av_novel_compounds_DR_H <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_H), c(ind_DR, ind_H)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_DR_H <- metadata_Av_dist_novel_compounds[c(ind_DR, ind_H), ]
ano_Av_novel_compounds_DR_H <- vegan::anosim(dist_Av_novel_compounds_DR_H, metadata_Av_dist_novel_compounds_DR_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["DNAReplication_Heat"] <- ano_Av_novel_compounds_DR_H$signif

# DNAReplication vs MembraneDisruption
dist_Av_novel_compounds_DR_MD <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_MD), c(ind_DR, ind_MD)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_DR_MD <- metadata_Av_dist_novel_compounds[c(ind_DR, ind_MD), ]
ano_Av_novel_compounds_DR_MD <- vegan::anosim(dist_Av_novel_compounds_DR_MD, metadata_Av_dist_novel_compounds_DR_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["DNAReplication_MembraneDisruption"] <- ano_Av_novel_compounds_DR_MD$signif

# DNAReplication vs ProteinSynthesis30SInhibition
dist_Av_novel_compounds_DR_PS30SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_PS30SI), c(ind_DR, ind_PS30SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_DR_PS30SI <- metadata_Av_dist_novel_compounds[c(ind_DR, ind_PS30SI), ]
ano_Av_novel_compounds_DR_PS30SI <- vegan::anosim(dist_Av_novel_compounds_DR_PS30SI, metadata_Av_dist_novel_compounds_DR_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["DNAReplication_ProteinSynthesis30SInhibition"] <- ano_Av_novel_compounds_DR_PS30SI$signif

# DNAReplication vs ProteinSynthesis50SInhibition
dist_Av_novel_compounds_DR_PS50SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_PS50SI), c(ind_DR, ind_PS50SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_DR_PS50SI <- metadata_Av_dist_novel_compounds[c(ind_DR, ind_PS50SI), ]
ano_Av_novel_compounds_DR_PS50SI <- vegan::anosim(dist_Av_novel_compounds_DR_PS50SI, metadata_Av_dist_novel_compounds_DR_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["DNAReplication_ProteinSynthesis50SInhibition"] <- ano_Av_novel_compounds_DR_PS50SI$signif

# DNAReplication vs ProteinSynthesistRNAInterference
dist_Av_novel_compounds_DR_PStRI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_PStRI), c(ind_DR, ind_PStRI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_DR_PStRI <- metadata_Av_dist_novel_compounds[c(ind_DR, ind_PStRI), ]
ano_Av_novel_compounds_DR_PStRI <- vegan::anosim(dist_Av_novel_compounds_DR_PStRI, metadata_Av_dist_novel_compounds_DR_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["DNAReplication_ProteinSynthesistRNAInterference"] <- ano_Av_novel_compounds_DR_PStRI$signif

# DNATranscription vs FolicAcidMetabolism
dist_Av_novel_compounds_DT_FAM <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_FAM), c(ind_DT, ind_FAM)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_DT_FAM <- metadata_Av_dist_novel_compounds[c(ind_DT, ind_FAM), ]
ano_Av_novel_compounds_DT_FAM <- vegan::anosim(dist_Av_novel_compounds_DT_FAM, metadata_Av_dist_novel_compounds_DT_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["DNATranscription_FolicAcidMetabolism"] <- ano_Av_novel_compounds_DT_FAM$signif

# DNATranscription vs Heat
dist_Av_novel_compounds_DT_H <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_H), c(ind_DT, ind_H)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_DT_H <- metadata_Av_dist_novel_compounds[c(ind_DT, ind_H), ]
ano_Av_novel_compounds_DT_H <- vegan::anosim(dist_Av_novel_compounds_DT_H, metadata_Av_dist_novel_compounds_DT_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["DNATranscription_Heat"] <- ano_Av_novel_compounds_DT_H$signif

# DNATranscription vs MembraneDisruption
dist_Av_novel_compounds_DT_MD <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_MD), c(ind_DT, ind_MD)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_DT_MD <- metadata_Av_dist_novel_compounds[c(ind_DT, ind_MD), ]
ano_Av_novel_compounds_DT_MD <- vegan::anosim(dist_Av_novel_compounds_DT_MD, metadata_Av_dist_novel_compounds_DT_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["DNATranscription_MembraneDisruption"] <- ano_Av_novel_compounds_DT_MD$signif

# DNATranscription vs ProteinSynthesis30SInhibition
dist_Av_novel_compounds_DT_PS30SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_PS30SI), c(ind_DT, ind_PS30SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_DT_PS30SI <- metadata_Av_dist_novel_compounds[c(ind_DT, ind_PS30SI), ]
ano_Av_novel_compounds_DT_PS30SI <- vegan::anosim(dist_Av_novel_compounds_DT_PS30SI, metadata_Av_dist_novel_compounds_DT_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["DNATranscription_ProteinSynthesis30SInhibition"] <- ano_Av_novel_compounds_DT_PS30SI$signif

# DNATranscription vs ProteinSynthesis50SInhibition
dist_Av_novel_compounds_DT_PS50SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_PS50SI), c(ind_DT, ind_PS50SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_DT_PS50SI <- metadata_Av_dist_novel_compounds[c(ind_DT, ind_PS50SI), ]
ano_Av_novel_compounds_DT_PS50SI <- vegan::anosim(dist_Av_novel_compounds_DT_PS50SI, metadata_Av_dist_novel_compounds_DT_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["DNATranscription_ProteinSynthesis50SInhibition"] <- ano_Av_novel_compounds_DT_PS50SI$signif

# DNATranscription vs ProteinSynthesistRNAInterference
dist_Av_novel_compounds_DT_PStRI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_PStRI), c(ind_DT, ind_PStRI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_DT_PStRI <- metadata_Av_dist_novel_compounds[c(ind_DT, ind_PStRI), ]
ano_Av_novel_compounds_DT_PStRI <- vegan::anosim(dist_Av_novel_compounds_DT_PStRI, metadata_Av_dist_novel_compounds_DT_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["DNATranscription_ProteinSynthesistRNAInterference"] <- ano_Av_novel_compounds_DT_PStRI$signif

# FolicAcidMetabolism vs Heat
dist_Av_novel_compounds_FAM_H <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_H), c(ind_FAM, ind_H)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_FAM_H <- metadata_Av_dist_novel_compounds[c(ind_FAM, ind_H), ]
ano_Av_novel_compounds_FAM_H <- vegan::anosim(dist_Av_novel_compounds_FAM_H, metadata_Av_dist_novel_compounds_FAM_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["FolicAcidMetabolism_Heat"] <- ano_Av_novel_compounds_FAM_H$signif

# FolicAcidMetabolism vs MembraneDisruption
dist_Av_novel_compounds_FAM_MD <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_MD), c(ind_FAM, ind_MD)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_FAM_MD <- metadata_Av_dist_novel_compounds[c(ind_FAM, ind_MD), ]
ano_Av_novel_compounds_FAM_MD <- vegan::anosim(dist_Av_novel_compounds_FAM_MD, metadata_Av_dist_novel_compounds_FAM_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["FolicAcidMetabolism_MembraneDisruption"] <- ano_Av_novel_compounds_FAM_MD$signif

# FolicAcidMetabolism vs ProteinSynthesis30SInhibition
dist_Av_novel_compounds_FAM_PS30SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_PS30SI), c(ind_FAM, ind_PS30SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_FAM_PS30SI <- metadata_Av_dist_novel_compounds[c(ind_FAM, ind_PS30SI), ]
ano_Av_novel_compounds_FAM_PS30SI <- vegan::anosim(dist_Av_novel_compounds_FAM_PS30SI, metadata_Av_dist_novel_compounds_FAM_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["FolicAcidMetabolism_ProteinSynthesis30SInhibition"] <- ano_Av_novel_compounds_FAM_PS30SI$signif

# FolicAcidMetabolism vs ProteinSynthesis50SInhibition
dist_Av_novel_compounds_FAM_PS50SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_PS50SI), c(ind_FAM, ind_PS50SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_FAM_PS50SI <- metadata_Av_dist_novel_compounds[c(ind_FAM, ind_PS50SI), ]
ano_Av_novel_compounds_FAM_PS50SI <- vegan::anosim(dist_Av_novel_compounds_FAM_PS50SI, metadata_Av_dist_novel_compounds_FAM_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["FolicAcidMetabolism_ProteinSynthesis50SInhibition"] <- ano_Av_novel_compounds_FAM_PS50SI$signif

# FolicAcidMetabolism vs ProteinSynthesistRNAInterference
dist_Av_novel_compounds_FAM_PStRI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_PStRI), c(ind_FAM, ind_PStRI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_FAM_PStRI <- metadata_Av_dist_novel_compounds[c(ind_FAM, ind_PStRI), ]
ano_Av_novel_compounds_FAM_PStRI <- vegan::anosim(dist_Av_novel_compounds_FAM_PStRI, metadata_Av_dist_novel_compounds_FAM_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["FolicAcidMetabolism_ProteinSynthesistRNAInterference"] <- ano_Av_novel_compounds_FAM_PStRI$signif

# Heat vs MembraneDisruption
dist_Av_novel_compounds_H_MD <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_MD), c(ind_H, ind_MD)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_H_MD <- metadata_Av_dist_novel_compounds[c(ind_H, ind_MD), ]
ano_Av_novel_compounds_H_MD <- vegan::anosim(dist_Av_novel_compounds_H_MD, metadata_Av_dist_novel_compounds_H_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["Heat_MembraneDisruption"] <- ano_Av_novel_compounds_H_MD$signif

# Heat vs ProteinSynthesis30SInhibition
dist_Av_novel_compounds_H_PS30SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_PS30SI), c(ind_H, ind_PS30SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_H_PS30SI <- metadata_Av_dist_novel_compounds[c(ind_H, ind_PS30SI), ]
ano_Av_novel_compounds_H_PS30SI <- vegan::anosim(dist_Av_novel_compounds_H_PS30SI, metadata_Av_dist_novel_compounds_H_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["Heat_ProteinSynthesis30SInhibition"] <- ano_Av_novel_compounds_H_PS30SI$signif

# Heat vs ProteinSynthesis50SInhibition
dist_Av_novel_compounds_H_PS50SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_PS50SI), c(ind_H, ind_PS50SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_H_PS50SI <- metadata_Av_dist_novel_compounds[c(ind_H, ind_PS50SI), ]
ano_Av_novel_compounds_H_PS50SI <- vegan::anosim(dist_Av_novel_compounds_H_PS50SI, metadata_Av_dist_novel_compounds_H_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["Heat_ProteinSynthesis50SInhibition"] <- ano_Av_novel_compounds_H_PS50SI$signif

# Heat vs ProteinSynthesistRNAInterference
dist_Av_novel_compounds_H_PStRI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_PStRI), c(ind_H, ind_PStRI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_H_PStRI <- metadata_Av_dist_novel_compounds[c(ind_H, ind_PStRI), ]
ano_Av_novel_compounds_H_PStRI <- vegan::anosim(dist_Av_novel_compounds_H_PStRI, metadata_Av_dist_novel_compounds_H_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["Heat_ProteinSynthesistRNAInterference"] <- ano_Av_novel_compounds_H_PStRI$signif

# MembraneDisruption vs ProteinSynthesis30SInhibition
dist_Av_novel_compounds_MD_PS30SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_MD, ind_PS30SI), c(ind_MD, ind_PS30SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_MD_PS30SI <- metadata_Av_dist_novel_compounds[c(ind_MD, ind_PS30SI), ]
ano_Av_novel_compounds_MD_PS30SI <- vegan::anosim(dist_Av_novel_compounds_MD_PS30SI, metadata_Av_dist_novel_compounds_MD_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["MembraneDisruption_ProteinSynthesis30SInhibition"] <- ano_Av_novel_compounds_MD_PS30SI$signif

# MembraneDisruption vs ProteinSynthesis50SInhibition
dist_Av_novel_compounds_MD_PS50SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_MD, ind_PS50SI), c(ind_MD, ind_PS50SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_MD_PS50SI <- metadata_Av_dist_novel_compounds[c(ind_MD, ind_PS50SI), ]
ano_Av_novel_compounds_MD_PS50SI <- vegan::anosim(dist_Av_novel_compounds_MD_PS50SI, metadata_Av_dist_novel_compounds_MD_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["MembraneDisruption_ProteinSynthesis50SInhibition"] <- ano_Av_novel_compounds_MD_PS50SI$signif

# MembraneDisruption vs ProteinSynthesistRNAInterference
dist_Av_novel_compounds_MD_PStRI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_MD, ind_PStRI), c(ind_MD, ind_PStRI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_MD_PStRI <- metadata_Av_dist_novel_compounds[c(ind_MD, ind_PStRI), ]
ano_Av_novel_compounds_MD_PStRI <- vegan::anosim(dist_Av_novel_compounds_MD_PStRI, metadata_Av_dist_novel_compounds_MD_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["MembraneDisruption_ProteinSynthesistRNAInterference"] <- ano_Av_novel_compounds_MD_PStRI$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesis50SInhibition
dist_Av_novel_compounds_PS30SI_PS50SI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PS30SI, ind_PS50SI), c(ind_PS30SI, ind_PS50SI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_PS30SI_PS50SI <- metadata_Av_dist_novel_compounds[c(ind_PS30SI, ind_PS50SI), ]
ano_Av_novel_compounds_PS30SI_PS50SI <- vegan::anosim(dist_Av_novel_compounds_PS30SI_PS50SI, metadata_Av_dist_novel_compounds_PS30SI_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["ProteinSynthesis30SInhibition_ProteinSynthesis50SInhibition"] <- ano_Av_novel_compounds_PS30SI_PS50SI$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesistRNAInterference
dist_Av_novel_compounds_PS30SI_PStRI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PS30SI, ind_PStRI), c(ind_PS30SI, ind_PStRI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_PS30SI_PStRI <- metadata_Av_dist_novel_compounds[c(ind_PS30SI, ind_PStRI), ]
ano_Av_novel_compounds_PS30SI_PStRI <- vegan::anosim(dist_Av_novel_compounds_PS30SI_PStRI, metadata_Av_dist_novel_compounds_PS30SI_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["ProteinSynthesis30SInhibition_ProteinSynthesistRNAInterference"] <- ano_Av_novel_compounds_PS30SI_PStRI$signif

# ProteinSynthesis50SInhibition vs ProteinSynthesistRNAInterference
dist_Av_novel_compounds_PS50SI_PStRI <- dist_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PS50SI, ind_PStRI), c(ind_PS50SI, ind_PStRI)] %>% 
  as.dist()

metadata_Av_dist_novel_compounds_PS50SI_PStRI <- metadata_Av_dist_novel_compounds[c(ind_PS50SI, ind_PStRI), ]
ano_Av_novel_compounds_PS50SI_PStRI <- vegan::anosim(dist_Av_novel_compounds_PS50SI_PStRI, metadata_Av_dist_novel_compounds_PS50SI_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Av_novel_compounds["ProteinSynthesis50SInhibition_ProteinSynthesistRNAInterference"] <- ano_Av_novel_compounds_PS50SI_PStRI$signif

# Correction for pairwise comparisons
pairwise_anosim_Av_novel_compounds_bonferroni <- p.adjust(pairwise_anosim_Av_novel_compounds, method = "bonferroni")
pairwise_anosim_Av_novel_compounds_BH <- p.adjust(pairwise_anosim_Av_novel_compounds, method = "BH")

pairwise_anosim_Av_novel_compounds_df <- as.data.frame(pairwise_anosim_Av_novel_compounds)
pairwise_anosim_Av_novel_compounds_bonferroni_df <- as.data.frame(pairwise_anosim_Av_novel_compounds_bonferroni)
pairwise_anosim_Av_novel_compounds_BH_df <- as.data.frame(pairwise_anosim_Av_novel_compounds_BH)

pairwise_anosim_Av_novel_compounds_export <- merge(pairwise_anosim_Av_novel_compounds_df, pairwise_anosim_Av_novel_compounds_bonferroni_df, by = "row.names") %>% 
  column_to_rownames(var = "Row.names") %>% 
  merge(pairwise_anosim_Av_novel_compounds_BH_df, by = "row.names") %>% 
  rename("classes" = "Row.names", "p" = "pairwise_anosim_Av_novel_compounds", "p.adjust.bonferroni" = "pairwise_anosim_Av_novel_compounds_bonferroni", "p.adjust.BH" = "pairwise_anosim_Av_novel_compounds_BH")
write.csv2(file = "ANOSIM_Av_novel_compounds_pairwise.csv", pairwise_anosim_Av_novel_compounds_export) # Export results as csv file

### 5.2.2. Fn ----
#### 5.2.2.1. Fingerprinting ----
# Randomly resample to a fixed number of cells
flowData_transformed_Fn_resample <- FCS_resample(flowData_transformed_Fn_norm_noneg, sample = 9775, replace = TRUE)

# Make metadata file for resampled data
names_flowData_transformed_Fn_resample <- sampleNames(flowData_transformed_Fn_resample)
metadata_Fn_resample <- metadata_Fn_noneg[match(names_flowData_transformed_Fn_resample, metadata_Fn_noneg$Filename), ]
row.names(metadata_Fn_resample) <- NULL

metadata_Fn_resample_novel_compounds <- metadata_Fn_resample[c(1, 5:72), ]
rownames(metadata_Fn_resample_novel_compounds) <- NULL

# Calculating fingerprint with bandwidth = 0.01
fbasis_Fn <- flowBasis(flowData_transformed_Fn_resample, param = param_div, nbin = 128, 
                       bw = 0.01, normalize = function(x) x)

fbasis_Fn_novel_compounds <- flowBasis(flowData_transformed_Fn_resample[c(1, 5:72)], param = param_div, nbin = 128,
                                       bw = 0.01, normalize = function(x) x)

#### 5.2.2.2. Ordination ----
# Beta-diversity assessment of fingerprint (Bray-Curtis)
beta.div_Fn <- beta_div_fcm(fbasis_Fn, ord.type = "PCoA")
beta.div_Fn_NMDS <- beta_div_fcm(fbasis_Fn, ord.type = "NMDS")

# Beta-diversity assessment of fingerprints excluding cephalothin (Bray-Curtis)
beta.div_Fn_novel_compounds <- beta_div_fcm(fbasis_Fn_novel_compounds, ord.type = "PCoA")
beta.div_Fn_novel_compounds_NMDS <- beta_div_fcm(fbasis_Fn_novel_compounds, ord.type = "NMDS")

# Plot ordination
pbetadiv_Fn_PCoA <- plot_beta_fcm(beta.div_Fn, color = as.factor(metadata_Fn_resample$MOA), shape = as.factor(metadata_Fn_resample$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL)+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(pbetadiv_Fn_PCoA)

pbetadiv_Fn_NMDS <- plot_beta_fcm(beta.div_Fn_NMDS, color = as.factor(metadata_Fn_resample$MOA), shape = as.factor(metadata_Fn_resample$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL, x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(pbetadiv_Fn_NMDS)

# Plot ordination excluding cephalothin
pbetadiv_Fn_novel_compounds_PCoA <- plot_beta_fcm(beta.div_Fn_novel_compounds, color = as.factor(metadata_Fn_resample_novel_compounds$MOA), shape = as.factor(metadata_Fn_resample_novel_compounds$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL)+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:4, 6:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(pbetadiv_Fn_novel_compounds_PCoA)

pbetadiv_Fn_novel_compounds_NMDS <- plot_beta_fcm(beta.div_Fn_novel_compounds_NMDS, color = as.factor(metadata_Fn_resample_novel_compounds$MOA), shape = as.factor(metadata_Fn_resample_novel_compounds$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL, x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:4, 6:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(pbetadiv_Fn_novel_compounds_NMDS)

#### 5.2.2.3. Statistics ----
### Bray-Curtis dissimilarity matrix
# All samples
fbasis_Fn_2 <- fbasis_Fn@basis/apply(fbasis_Fn@basis, 1, max)
fbasis_Fn_2 <- base::round(fbasis_Fn_2, 4)
dist_Fn <- vegan::vegdist(fbasis_Fn_2, method = "bray", binary = FALSE)

metadata_Fn_dist <- metadata_Fn_noneg
rownames(metadata_Fn_dist) <- metadata_Fn_dist$Filename

# Excluding cephalothin
fbasis_Fn_3 <- fbasis_Fn_2[c(1, 5:72),]
dist_Fn_novel_compounds <- vegan::vegdist(fbasis_Fn_3, method = "bray", binary = FALSE)

metadata_Fn_dist_novel_compounds <- metadata_Fn_noneg[c(1, 5:72),]
rownames(metadata_Fn_dist_novel_compounds) <- metadata_Fn_dist_novel_compounds$Filename

### ANOSIM
permutations = 10000

# All treatments
ano_Fn <- vegan::anosim(dist_Fn, metadata_Fn_dist$MOA, distance = "bray", permutations = permutations)
ano_Fn_export <- data.frame("ANOSIM Statistic R" = ano_Fn$statistic, "P" = ano_Fn$signif, "permutations" = ano_Fn$permutations)
write.csv2(file = "ANOSIM_Fn_all_treatments.csv", ano_Fn_export) # Export results as csv file
plot_dissimilarity_anosim_Fn <- plot(ano_Fn, xaxt = "n", main = "Dissimilarity ranks Fn (ANOSIM)", xlab = "Class", ylab = "Dissimilarity rank")
axis(1, at = 1:11, labels = c("Between", "CWS", "DNARepl", "DNATrans", "FAM", "Heat", "MD", "Control", "PS-30SInh", "PS-50SInh", "PS-tRNAInt"))

# Excluding cephalothin
ano_Fn_novel_compounds <- vegan::anosim(dist_Fn_novel_compounds, metadata_Fn_dist_novel_compounds$MOA, distance = "bray", permutations = permutations)
ano_Fn_novel_compounds_export <- data.frame("ANOSIM Statistic R" = ano_Fn_novel_compounds$statistic, "P" = ano_Fn_novel_compounds$signif, "permutations" = ano_Fn_novel_compounds$permutations)
write.csv2(file = "ANOSIM_Fn_novel_compounds.csv", ano_Fn_novel_compounds_export) # Export results as csv file
plot_dissimilarity_anosim_Fn_novel_compounds <- plot(ano_Fn_novel_compounds, xaxt = "n", main = "Dissimilarity ranks Fn (ANOSIM)", xlab = "Class", ylab = "Dissimilarity rank")
axis(1, at = 1:11, labels = c("Between", "CWS", "DNARepl", "DNATrans", "FAM", "Heat", "MD", "Control", "PS-30SInh", "PS-50SInh", "PS-tRNAInt"))

## Pairwise comparison different MOA
# All treatments
pairwise_anosim_Fn <- numeric()

# CellWallSynthesis vs Control
dist_Fn_CWS_PC <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(1:18), c(1:18)] %>% 
  as.dist()

metadata_Fn_dist_CWS_PC <- metadata_Fn_dist[c(1:18), ]
ano_Fn_CWS_PC <- vegan::anosim(dist_Fn_CWS_PC, metadata_Fn_dist_CWS_PC$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["CellWallSynthesis_Control"] <- ano_Fn_CWS_PC$signif

# CellWallSynthesis vs DNAReplication
dist_Fn_CWS_DR <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(1:12, 19:27), c(1:12, 19:27)] %>% 
  as.dist()

metadata_Fn_dist_CWS_DR <- metadata_Fn_dist[c(1:12, 19:27), ]
ano_Fn_CWS_DR <- vegan::anosim(dist_Fn_CWS_DR, metadata_Fn_dist_CWS_DR$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["CellWallSynthesis_DNAReplication"] <- ano_Fn_CWS_DR$signif

# CellWallSynthesis vs DNATranscription
dist_Fn_CWS_DT <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(1:12, 28:33), c(1:12, 28:33)] %>% 
  as.dist()

metadata_Fn_dist_CWS_DT <- metadata_Fn_dist[c(1:12, 28:33), ]
ano_Fn_CWS_DT <- vegan::anosim(dist_Fn_CWS_DT, metadata_Fn_dist_CWS_DT$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["CellWallSynthesis_DNATranscription"] <- ano_Fn_CWS_DT$signif

# CellWallSynthesis vs FolicAcidMetabolism
dist_Fn_CWS_FAM <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(1:12, 34:36), c(1:12, 34:36)] %>% 
  as.dist()

metadata_Fn_dist_CWS_FAM <- metadata_Fn_dist[c(1:12, 34:36), ]
ano_Fn_CWS_FAM <- vegan::anosim(dist_Fn_CWS_FAM, metadata_Fn_dist_CWS_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["CellWallSynthesis_FolicAcidMetabolism"] <- ano_Fn_CWS_FAM$signif

# CellWallSynthesis vs Heat
dist_Fn_CWS_H <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(1:12, 37:39), c(1:12, 37:39)] %>% 
  as.dist()

metadata_Fn_dist_CWS_H <- metadata_Fn_dist[c(1:12, 37:39), ]
ano_Fn_CWS_H <- vegan::anosim(dist_Fn_CWS_H, metadata_Fn_dist_CWS_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["CellWallSynthesis_Heat"] <- ano_Fn_CWS_H$signif

# CellWallSynthesis vs MembraneDisruption
dist_Fn_CWS_MD <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(1:12, 40:51), c(1:12, 40:51)] %>% 
  as.dist()

metadata_Fn_dist_CWS_MD <- metadata_Fn_dist[c(1:12, 40:51), ]
ano_Fn_CWS_MD <- vegan::anosim(dist_Fn_CWS_MD, metadata_Fn_dist_CWS_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["CellWallSynthesis_MembraneDisruption"] <- ano_Fn_CWS_MD$signif

# CellWallSynthesis vs ProteinSynthesis30SInhibition
dist_Fn_CWS_PS30SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(1:12, 52:60), c(1:12, 52:60)] %>% 
  as.dist()

metadata_Fn_dist_CWS_PS30SI <- metadata_Fn_dist[c(1:12, 52:60), ]
ano_Fn_CWS_PS30SI <- vegan::anosim(dist_Fn_CWS_PS30SI, metadata_Fn_dist_CWS_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["CellWallSynthesis_ProteinSynthesis30SInhibition"] <- ano_Fn_CWS_PS30SI$signif

# CellWallSynthesis vs ProteinSynthesis50SInhibition
dist_Fn_CWS_PS50SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(1:12, 61:69), c(1:12, 61:69)] %>% 
  as.dist()

metadata_Fn_dist_CWS_PS50SI <- metadata_Fn_dist[c(1:12, 61:69), ]
ano_Fn_CWS_PS50SI <- vegan::anosim(dist_Fn_CWS_PS50SI, metadata_Fn_dist_CWS_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["CellWallSynthesis_ProteinSynthesis50SInhibition"] <- ano_Fn_CWS_PS50SI$signif

# CellWallSynthesis vs ProteinSynthesistRNAInterference
dist_Fn_CWS_PStRI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(1:12, 70:72), c(1:12, 70:72)] %>% 
  as.dist()

metadata_Fn_dist_CWS_PStRI <- metadata_Fn_dist[c(1:12, 70:72), ]
ano_Fn_CWS_PStRI <- vegan::anosim(dist_Fn_CWS_PStRI, metadata_Fn_dist_CWS_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["CellWallSynthesis_ProteinSynthesistRNAInterference"] <- ano_Fn_CWS_PStRI$signif

# Control vs DNAReplication
dist_Fn_PC_DR <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(13:18, 19:27), c(13:18, 19:27)] %>% 
  as.dist()

metadata_Fn_dist_PC_DR <- metadata_Fn_dist[c(13:18, 19:27), ]
ano_Fn_PC_DR <- vegan::anosim(dist_Fn_PC_DR, metadata_Fn_dist_PC_DR$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["Control_DNAReplication"] <- ano_Fn_PC_DR$signif

# Control vs DNATranscription
dist_Fn_PC_DT <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(13:18, 28:33), c(13:18, 28:33)] %>% 
  as.dist()

metadata_Fn_dist_PC_DT <- metadata_Fn_dist[c(13:18, 28:33), ]
ano_Fn_PC_DT <- vegan::anosim(dist_Fn_PC_DT, metadata_Fn_dist_PC_DT$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["Control_DNATranscription"] <- ano_Fn_PC_DT$signif

# Control vs FolicAcidMetabolism
dist_Fn_PC_FAM <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(13:18, 34:36), c(13:18, 34:36)] %>% 
  as.dist()

metadata_Fn_dist_PC_FAM <- metadata_Fn_dist[c(13:18, 34:36), ]
ano_Fn_PC_FAM <- vegan::anosim(dist_Fn_PC_FAM, metadata_Fn_dist_PC_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["Control_FolicAcidMetabolism"] <- ano_Fn_PC_FAM$signif

# Control vs Heat
dist_Fn_PC_H <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(13:18, 37:39), c(13:18, 37:39)] %>% 
  as.dist()

metadata_Fn_dist_PC_H <- metadata_Fn_dist[c(13:18, 37:39), ]
ano_Fn_PC_H <- vegan::anosim(dist_Fn_PC_H, metadata_Fn_dist_PC_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["Control_Heat"] <- ano_Fn_PC_H$signif

# Control vs MembraneDisruption
dist_Fn_PC_MD <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(13:18, 40:51), c(13:18, 40:51)] %>% 
  as.dist()

metadata_Fn_dist_PC_MD <- metadata_Fn_dist[c(13:18, 40:51), ]
ano_Fn_PC_MD <- vegan::anosim(dist_Fn_PC_MD, metadata_Fn_dist_PC_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["Control_MembraneDisruption"] <- ano_Fn_PC_MD$signif

# Control vs ProteinSynthesis30SInhibition
dist_Fn_PC_PS30SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(13:18, 52:60), c(13:18, 52:60)] %>% 
  as.dist()

metadata_Fn_dist_PC_PS30SI <- metadata_Fn_dist[c(13:18, 52:60), ]
ano_Fn_PC_PS30SI <- vegan::anosim(dist_Fn_PC_PS30SI, metadata_Fn_dist_PC_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["Control_ProteinSynthesis30SInhibition"] <- ano_Fn_PC_PS30SI$signif

# Control vs ProteinSynthesis50SInhibition
dist_Fn_PC_PS50SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(13:18, 61:69), c(13:18, 61:69)] %>% 
  as.dist()

metadata_Fn_dist_PC_PS50SI <- metadata_Fn_dist[c(13:18, 61:69), ]
ano_Fn_PC_PS50SI <- vegan::anosim(dist_Fn_PC_PS50SI, metadata_Fn_dist_PC_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["Control_ProteinSynthesis50SInhibition"] <- ano_Fn_PC_PS50SI$signif

# Control vs ProteinSynthesistRNAInterference
dist_Fn_PC_PStRI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(13:18, 70:72), c(13:18, 70:72)] %>% 
  as.dist()

metadata_Fn_dist_PC_PStRI <- metadata_Fn_dist[c(13:18, 70:72), ]
ano_Fn_PC_PStRI <- vegan::anosim(dist_Fn_PC_PStRI, metadata_Fn_dist_PC_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["Control_ProteinSynthesistRNAInterference"] <- ano_Fn_PC_PStRI$signif

# DNAReplication vs DNATranscription
dist_Fn_DR_DT <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(19:27, 28:33), c(19:27, 28:33)] %>% 
  as.dist()

metadata_Fn_dist_DR_DT <- metadata_Fn_dist[c(19:27, 28:33), ]
ano_Fn_DR_DT <- vegan::anosim(dist_Fn_DR_DT, metadata_Fn_dist_DR_DT$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["DNAReplication_DNATranscription"] <- ano_Fn_DR_DT$signif

# DNAReplication vs FolicAcidMetabolism
dist_Fn_DR_FAM <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(19:27, 34:36), c(19:27, 34:36)] %>% 
  as.dist()

metadata_Fn_dist_DR_FAM <- metadata_Fn_dist[c(19:27, 34:36), ]
ano_Fn_DR_FAM <- vegan::anosim(dist_Fn_DR_FAM, metadata_Fn_dist_DR_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["DNAReplication_FolicAcidMetabolism"] <- ano_Fn_DR_FAM$signif

# DNAReplication vs Heat
dist_Fn_DR_H <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(19:27, 37:39), c(19:27, 37:39)] %>% 
  as.dist()

metadata_Fn_dist_DR_H <- metadata_Fn_dist[c(19:27, 37:39), ]
ano_Fn_DR_H <- vegan::anosim(dist_Fn_DR_H, metadata_Fn_dist_DR_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["DNAReplication_Heat"] <- ano_Fn_DR_H$signif

# DNAReplication vs MembraneDisruption
dist_Fn_DR_MD <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(19:27, 40:51), c(19:27, 40:51)] %>% 
  as.dist()

metadata_Fn_dist_DR_MD <- metadata_Fn_dist[c(19:27, 40:51), ]
ano_Fn_DR_MD <- vegan::anosim(dist_Fn_DR_MD, metadata_Fn_dist_DR_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["DNAReplication_MembraneDisruption"] <- ano_Fn_DR_MD$signif

# DNAReplication vs ProteinSynthesis30SInhibition
dist_Fn_DR_PS30SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(19:27, 52:60), c(19:27, 52:60)] %>% 
  as.dist()

metadata_Fn_dist_DR_PS30SI <- metadata_Fn_dist[c(19:27, 52:60), ]
ano_Fn_DR_PS30SI <- vegan::anosim(dist_Fn_DR_PS30SI, metadata_Fn_dist_DR_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["DNAReplication_ProteinSynthesis30SInhibition"] <- ano_Fn_DR_PS30SI$signif

# DNAReplication vs ProteinSynthesis50SInhibition
dist_Fn_DR_PS50SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(19:27, 61:69), c(19:27, 61:69)] %>% 
  as.dist()

metadata_Fn_dist_DR_PS50SI <- metadata_Fn_dist[c(19:27, 61:69), ]
ano_Fn_DR_PS50SI <- vegan::anosim(dist_Fn_DR_PS50SI, metadata_Fn_dist_DR_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["DNAReplication_ProteinSynthesis50SInhibition"] <- ano_Fn_DR_PS50SI$signif

# DNAReplication vs ProteinSynthesistRNAInterference
dist_Fn_DR_PStRI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(19:27, 70:72), c(19:27, 70:72)] %>% 
  as.dist()

metadata_Fn_dist_DR_PStRI <- metadata_Fn_dist[c(19:27, 70:72), ]
ano_Fn_DR_PStRI <- vegan::anosim(dist_Fn_DR_PStRI, metadata_Fn_dist_DR_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["DNAReplication_ProteinSynthesistRNAInterference"] <- ano_Fn_DR_PStRI$signif

# DNATranscription vs FolicAcidMetabolism
dist_Fn_DT_FAM <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(28:33, 34:36), c(28:33, 34:36)] %>% 
  as.dist()

metadata_Fn_dist_DT_FAM <- metadata_Fn_dist[c(28:33, 34:36), ]
ano_Fn_DT_FAM <- vegan::anosim(dist_Fn_DT_FAM, metadata_Fn_dist_DT_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["DNATranscription_FolicAcidMetabolism"] <- ano_Fn_DT_FAM$signif

# DNATranscription vs Heat
dist_Fn_DT_H <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(28:33, 37:39), c(28:33, 37:39)] %>% 
  as.dist()

metadata_Fn_dist_DT_H <- metadata_Fn_dist[c(28:33, 37:39), ]
ano_Fn_DT_H <- vegan::anosim(dist_Fn_DT_H, metadata_Fn_dist_DT_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["DNATranscription_Heat"] <- ano_Fn_DT_H$signif

# DNATranscription vs MembraneDisruption
dist_Fn_DT_MD <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(28:33, 40:51), c(28:33, 40:51)] %>% 
  as.dist()

metadata_Fn_dist_DT_MD <- metadata_Fn_dist[c(28:33, 40:51), ]
ano_Fn_DT_MD <- vegan::anosim(dist_Fn_DT_MD, metadata_Fn_dist_DT_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["DNATranscription_MembraneDisruption"] <- ano_Fn_DT_MD$signif

# DNATranscription vs ProteinSynthesis30SInhibition
dist_Fn_DT_PS30SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(28:33, 52:60), c(28:33, 52:60)] %>% 
  as.dist()

metadata_Fn_dist_DT_PS30SI <- metadata_Fn_dist[c(28:33, 52:60), ]
ano_Fn_DT_PS30SI <- vegan::anosim(dist_Fn_DT_PS30SI, metadata_Fn_dist_DT_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["DNATranscription_ProteinSynthesis30SInhibition"] <- ano_Fn_DT_PS30SI$signif

# DNATranscription vs ProteinSynthesis50SInhibition
dist_Fn_DT_PS50SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(28:33, 61:69), c(28:33, 61:69)] %>% 
  as.dist()

metadata_Fn_dist_DT_PS50SI <- metadata_Fn_dist[c(28:33, 61:69), ]
ano_Fn_DT_PS50SI <- vegan::anosim(dist_Fn_DT_PS50SI, metadata_Fn_dist_DT_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["DNATranscription_ProteinSynthesis50SInhibition"] <- ano_Fn_DT_PS50SI$signif

# DNATranscription vs ProteinSynthesistRNAInterference
dist_Fn_DT_PStRI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(28:33, 70:72), c(28:33, 70:72)] %>% 
  as.dist()

metadata_Fn_dist_DT_PStRI <- metadata_Fn_dist[c(28:33, 70:72), ]
ano_Fn_DT_PStRI <- vegan::anosim(dist_Fn_DT_PStRI, metadata_Fn_dist_DT_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["DNATranscription_ProteinSynthesistRNAInterference"] <- ano_Fn_DT_PStRI$signif

# FolicAcidMetabolism vs Heat
dist_Fn_FAM_H <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(34:36, 37:39), c(34:36, 37:39)] %>% 
  as.dist()

metadata_Fn_dist_FAM_H <- metadata_Fn_dist[c(34:36, 37:39), ]
ano_Fn_FAM_H <- vegan::anosim(dist_Fn_FAM_H, metadata_Fn_dist_FAM_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["FolicAcidMetabolism_Heat"] <- ano_Fn_FAM_H$signif

# FolicAcidMetabolism vs MembraneDisruption
dist_Fn_FAM_MD <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(34:36, 40:51), c(34:36, 40:51)] %>% 
  as.dist()

metadata_Fn_dist_FAM_MD <- metadata_Fn_dist[c(34:36, 40:51), ]
ano_Fn_FAM_MD <- vegan::anosim(dist_Fn_FAM_MD, metadata_Fn_dist_FAM_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["FolicAcidMetabolism_MembraneDisruption"] <- ano_Fn_FAM_MD$signif

# FolicAcidMetabolism vs ProteinSynthesis30SInhibition
dist_Fn_FAM_PS30SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(34:36, 52:60), c(34:36, 52:60)] %>% 
  as.dist()

metadata_Fn_dist_FAM_PS30SI <- metadata_Fn_dist[c(34:36, 52:60), ]
ano_Fn_FAM_PS30SI <- vegan::anosim(dist_Fn_FAM_PS30SI, metadata_Fn_dist_FAM_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["FolicAcidMetabolism_ProteinSynthesis30SInhibition"] <- ano_Fn_FAM_PS30SI$signif

# FolicAcidMetabolism vs ProteinSynthesis50SInhibition
dist_Fn_FAM_PS50SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(34:36, 61:69), c(34:36, 61:69)] %>% 
  as.dist()

metadata_Fn_dist_FAM_PS50SI <- metadata_Fn_dist[c(34:36, 61:69), ]
ano_Fn_FAM_PS50SI <- vegan::anosim(dist_Fn_FAM_PS50SI, metadata_Fn_dist_FAM_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["FolicAcidMetabolism_ProteinSynthesis50SInhibition"] <- ano_Fn_FAM_PS50SI$signif

# FolicAcidMetabolism vs ProteinSynthesistRNAInterference
dist_Fn_FAM_PStRI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(34:36, 70:72), c(34:36, 70:72)] %>% 
  as.dist()

metadata_Fn_dist_FAM_PStRI <- metadata_Fn_dist[c(34:36, 70:72), ]
ano_Fn_FAM_PStRI <- vegan::anosim(dist_Fn_FAM_PStRI, metadata_Fn_dist_FAM_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["FolicAcidMetabolism_ProteinSynthesistRNAInterference"] <- ano_Fn_FAM_PStRI$signif

# Heat vs MembraneDisruption
dist_Fn_H_MD <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(37:39, 40:51), c(37:39, 40:51)] %>% 
  as.dist()

metadata_Fn_dist_H_MD <- metadata_Fn_dist[c(37:39, 40:51), ]
ano_Fn_H_MD <- vegan::anosim(dist_Fn_H_MD, metadata_Fn_dist_H_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["Heat_MembraneDisruption"] <- ano_Fn_H_MD$signif

# Heat vs ProteinSynthesis30SInhibition
dist_Fn_H_PS30SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(37:39, 52:60), c(37:39, 52:60)] %>% 
  as.dist()

metadata_Fn_dist_H_PS30SI <- metadata_Fn_dist[c(37:39, 52:60), ]
ano_Fn_H_PS30SI <- vegan::anosim(dist_Fn_H_PS30SI, metadata_Fn_dist_H_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["Heat_ProteinSynthesis30SInhibition"] <- ano_Fn_H_PS30SI$signif

# Heat vs ProteinSynthesis50SInhibition
dist_Fn_H_PS50SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(37:39, 61:69), c(37:39, 61:69)] %>% 
  as.dist()

metadata_Fn_dist_H_PS50SI <- metadata_Fn_dist[c(37:39, 61:69), ]
ano_Fn_H_PS50SI <- vegan::anosim(dist_Fn_H_PS50SI, metadata_Fn_dist_H_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["Heat_ProteinSynthesis50SInhibition"] <- ano_Fn_H_PS50SI$signif

# Heat vs ProteinSynthesistRNAInterference
dist_Fn_H_PStRI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(37:39, 70:72), c(37:39, 70:72)] %>% 
  as.dist()

metadata_Fn_dist_H_PStRI <- metadata_Fn_dist[c(37:39, 70:72), ]
ano_Fn_H_PStRI <- vegan::anosim(dist_Fn_H_PStRI, metadata_Fn_dist_H_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["Heat_ProteinSynthesistRNAInterference"] <- ano_Fn_H_PStRI$signif

# MembraneDisruption vs ProteinSynthesis30SInhibition
dist_Fn_MD_PS30SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(40:51, 52:60), c(40:51, 52:60)] %>% 
  as.dist()

metadata_Fn_dist_MD_PS30SI <- metadata_Fn_dist[c(40:51, 52:60), ]
ano_Fn_MD_PS30SI <- vegan::anosim(dist_Fn_MD_PS30SI, metadata_Fn_dist_MD_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["MembraneDisruption_ProteinSynthesis30SInhibition"] <- ano_Fn_MD_PS30SI$signif

# MembraneDisruption vs ProteinSynthesis50SInhibition
dist_Fn_MD_PS50SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(40:51, 61:69), c(40:51, 61:69)] %>% 
  as.dist()

metadata_Fn_dist_MD_PS50SI <- metadata_Fn_dist[c(40:51, 61:69), ]
ano_Fn_MD_PS50SI <- vegan::anosim(dist_Fn_MD_PS50SI, metadata_Fn_dist_MD_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["MembraneDisruption_ProteinSynthesis50SInhibition"] <- ano_Fn_MD_PS50SI$signif

# MembraneDisruption vs ProteinSynthesistRNAInterference
dist_Fn_MD_PStRI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(40:51, 70:72), c(40:51, 70:72)] %>% 
  as.dist()

metadata_Fn_dist_MD_PStRI <- metadata_Fn_dist[c(40:51, 70:72), ]
ano_Fn_MD_PStRI <- vegan::anosim(dist_Fn_MD_PStRI, metadata_Fn_dist_MD_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["MembraneDisruption_ProteinSynthesistRNAInterference"] <- ano_Fn_MD_PStRI$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesis50SInhibition
dist_Fn_PS30SI_PS50SI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(52:60, 61:69), c(52:60, 61:69)] %>% 
  as.dist()

metadata_Fn_dist_PS30SI_PS50SI <- metadata_Fn_dist[c(52:60, 61:69), ]
ano_Fn_PS30SI_PS50SI <- vegan::anosim(dist_Fn_PS30SI_PS50SI, metadata_Fn_dist_PS30SI_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["ProteinSynthesis30SInhibition_ProteinSynthesis50SInhibition"] <- ano_Fn_PS30SI_PS50SI$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesistRNAInterference
dist_Fn_PS30SI_PStRI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(52:60, 70:72), c(52:60, 70:72)] %>% 
  as.dist()

metadata_Fn_dist_PS30SI_PStRI <- metadata_Fn_dist[c(52:60, 70:72), ]
ano_Fn_PS30SI_PStRI <- vegan::anosim(dist_Fn_PS30SI_PStRI, metadata_Fn_dist_PS30SI_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["ProteinSynthesis30SInhibition_ProteinSynthesistRNAInterference"] <- ano_Fn_PS30SI_PStRI$signif

# ProteinSynthesis50SInhibition vs ProteinSynthesistRNAInterference
dist_Fn_PS50SI_PStRI <- dist_Fn %>% 
  as.matrix() %>% 
  .[c(61:69, 70:72), c(61:69, 70:72)] %>% 
  as.dist()

metadata_Fn_dist_PS50SI_PStRI <- metadata_Fn_dist[c(61:69, 70:72), ]
ano_Fn_PS50SI_PStRI <- vegan::anosim(dist_Fn_PS50SI_PStRI, metadata_Fn_dist_PS50SI_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn["ProteinSynthesis50SInhibition_ProteinSynthesistRNAInterference"] <- ano_Fn_PS50SI_PStRI$signif

# Correction for pairwise comparisons
pairwise_anosim_Fn_bonferroni <- p.adjust(pairwise_anosim_Fn, method = "bonferroni")
pairwise_anosim_Fn_BH <- p.adjust(pairwise_anosim_Fn, method = "BH")

pairwise_anosim_Fn_df <- as.data.frame(pairwise_anosim_Fn)
pairwise_anosim_Fn_bonferroni_df <- as.data.frame(pairwise_anosim_Fn_bonferroni)
pairwise_anosim_Fn_BH_df <- as.data.frame(pairwise_anosim_Fn_BH)

pairwise_anosim_Fn_export <- merge(pairwise_anosim_Fn_df, pairwise_anosim_Fn_bonferroni_df, by = "row.names") %>% 
  column_to_rownames(var = "Row.names") %>% 
  merge(pairwise_anosim_Fn_BH_df, by = "row.names") %>% 
  rename("classes" = "Row.names", "p" = "pairwise_anosim_Fn", "p.adjust.bonferroni" = "pairwise_anosim_Fn_bonferroni", "p.adjust.BH" = "pairwise_anosim_Fn_BH")
write.csv2(file = "ANOSIM_Fn_all_treatments_pairwise.csv", pairwise_anosim_Fn_export) # Export results as csv file

# Novel compounds (exlcuding cephalothin)
pairwise_anosim_Fn_novel_compounds <- numeric()

# Indices MOA classes
ind_CWS <- c(1:9)
ind_PC <- c(10:15)
ind_DR <- c(16:24)
ind_DT <- c(25:30)
ind_FAM <- c(31:33)
ind_H <- c(34:36)
ind_MD <- c(37:48)
ind_PS30SI <- c(49:57)
ind_PS50SI <- c(58:66)
ind_PStRI <- c(67:69)

# CellWallSynthesis vs Control
dist_Fn_novel_compounds_CWS_PC <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PC), c(ind_CWS, ind_PC)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_CWS_PC <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_PC), ]
ano_Fn_novel_compounds_CWS_PC <- vegan::anosim(dist_Fn_novel_compounds_CWS_PC, metadata_Fn_dist_novel_compounds_CWS_PC$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["CellWallSynthesis_Control"] <- ano_Fn_novel_compounds_CWS_PC$signif

# CellWallSynthesis vs DNAReplication
dist_Fn_novel_compounds_CWS_DR <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_DR), c(ind_CWS, ind_DR)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_CWS_DR <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_DR), ]
ano_Fn_novel_compounds_CWS_DR <- vegan::anosim(dist_Fn_novel_compounds_CWS_DR, metadata_Fn_dist_novel_compounds_CWS_DR$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["CellWallSynthesis_DNAReplication"] <- ano_Fn_novel_compounds_CWS_DR$signif

# CellWallSynthesis vs DNATranscription
dist_Fn_novel_compounds_CWS_DT <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_DT), c(ind_CWS, ind_DT)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_CWS_DT <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_DT), ]
ano_Fn_novel_compounds_CWS_DT <- vegan::anosim(dist_Fn_novel_compounds_CWS_DT, metadata_Fn_dist_novel_compounds_CWS_DT$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["CellWallSynthesis_DNATranscription"] <- ano_Fn_novel_compounds_CWS_DT$signif

# CellWallSynthesis vs FolicAcidMetabolism
dist_Fn_novel_compounds_CWS_FAM <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_FAM), c(ind_CWS, ind_FAM)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_CWS_FAM <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_FAM), ]
ano_Fn_novel_compounds_CWS_FAM <- vegan::anosim(dist_Fn_novel_compounds_CWS_FAM, metadata_Fn_dist_novel_compounds_CWS_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["CellWallSynthesis_FolicAcidMetabolism"] <- ano_Fn_novel_compounds_CWS_FAM$signif

# CellWallSynthesis vs Heat
dist_Fn_novel_compounds_CWS_H <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_H), c(ind_CWS, ind_H)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_CWS_H <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_H), ]
ano_Fn_novel_compounds_CWS_H <- vegan::anosim(dist_Fn_novel_compounds_CWS_H, metadata_Fn_dist_novel_compounds_CWS_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["CellWallSynthesis_Heat"] <- ano_Fn_novel_compounds_CWS_H$signif

# CellWallSynthesis vs MembraneDisruption
dist_Fn_novel_compounds_CWS_MD <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_MD), c(ind_CWS, ind_MD)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_CWS_MD <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_MD), ]
ano_Fn_novel_compounds_CWS_MD <- vegan::anosim(dist_Fn_novel_compounds_CWS_MD, metadata_Fn_dist_novel_compounds_CWS_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["CellWallSynthesis_MembraneDisruption"] <- ano_Fn_novel_compounds_CWS_MD$signif

# CellWallSynthesis vs ProteinSynthesis30SInhibition
dist_Fn_novel_compounds_CWS_PS30SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PS30SI), c(ind_CWS, ind_PS30SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_CWS_PS30SI <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_PS30SI), ]
ano_Fn_novel_compounds_CWS_PS30SI <- vegan::anosim(dist_Fn_novel_compounds_CWS_PS30SI, metadata_Fn_dist_novel_compounds_CWS_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["CellWallSynthesis_ProteinSynthesis30SInhibition"] <- ano_Fn_novel_compounds_CWS_PS30SI$signif

# CellWallSynthesis vs ProteinSynthesis50SInhibition
dist_Fn_novel_compounds_CWS_PS50SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PS50SI), c(ind_CWS, ind_PS50SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_CWS_PS50SI <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_PS50SI), ]
ano_Fn_novel_compounds_CWS_PS50SI <- vegan::anosim(dist_Fn_novel_compounds_CWS_PS50SI, metadata_Fn_dist_novel_compounds_CWS_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["CellWallSynthesis_ProteinSynthesis50SInhibition"] <- ano_Fn_novel_compounds_CWS_PS50SI$signif

# CellWallSynthesis vs ProteinSynthesistRNAInterference
dist_Fn_novel_compounds_CWS_PStRI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PStRI), c(ind_CWS, ind_PStRI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_CWS_PStRI <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_PStRI), ]
ano_Fn_novel_compounds_CWS_PStRI <- vegan::anosim(dist_Fn_novel_compounds_CWS_PStRI, metadata_Fn_dist_novel_compounds_CWS_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["CellWallSynthesis_ProteinSynthesistRNAInterference"] <- ano_Fn_novel_compounds_CWS_PStRI$signif

# Control vs DNAReplication
dist_Fn_novel_compounds_PC_DR <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_DR), c(ind_PC, ind_DR)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_PC_DR <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_DR), ]
ano_Fn_novel_compounds_PC_DR <- vegan::anosim(dist_Fn_novel_compounds_PC_DR, metadata_Fn_dist_novel_compounds_PC_DR$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["Control_DNAReplication"] <- ano_Fn_novel_compounds_PC_DR$signif

# Control vs DNATranscription
dist_Fn_novel_compounds_PC_DT <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_DT), c(ind_PC, ind_DT)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_PC_DT <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_DT), ]
ano_Fn_novel_compounds_PC_DT <- vegan::anosim(dist_Fn_novel_compounds_PC_DT, metadata_Fn_dist_novel_compounds_PC_DT$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["Control_DNATranscription"] <- ano_Fn_novel_compounds_PC_DT$signif

# Control vs FolicAcidMetabolism
dist_Fn_novel_compounds_PC_FAM <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_FAM), c(ind_PC, ind_FAM)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_PC_FAM <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_FAM), ]
ano_Fn_novel_compounds_PC_FAM <- vegan::anosim(dist_Fn_novel_compounds_PC_FAM, metadata_Fn_dist_novel_compounds_PC_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["Control_FolicAcidMetabolism"] <- ano_Fn_novel_compounds_PC_FAM$signif

# Control vs Heat
dist_Fn_novel_compounds_PC_H <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_H), c(ind_PC, ind_H)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_PC_H <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_H), ]
ano_Fn_novel_compounds_PC_H <- vegan::anosim(dist_Fn_novel_compounds_PC_H, metadata_Fn_dist_novel_compounds_PC_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["Control_Heat"] <- ano_Fn_novel_compounds_PC_H$signif

# Control vs MembraneDisruption
dist_Fn_novel_compounds_PC_MD <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_MD), c(ind_PC, ind_MD)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_PC_MD <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_MD), ]
ano_Fn_novel_compounds_PC_MD <- vegan::anosim(dist_Fn_novel_compounds_PC_MD, metadata_Fn_dist_novel_compounds_PC_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["Control_MembraneDisruption"] <- ano_Fn_novel_compounds_PC_MD$signif

# Control vs ProteinSynthesis30SInhibition
dist_Fn_novel_compounds_PC_PS30SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_PS30SI), c(ind_PC, ind_PS30SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_PC_PS30SI <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_PS30SI), ]
ano_Fn_novel_compounds_PC_PS30SI <- vegan::anosim(dist_Fn_novel_compounds_PC_PS30SI, metadata_Fn_dist_novel_compounds_PC_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["Control_ProteinSynthesis30SInhibition"] <- ano_Fn_novel_compounds_PC_PS30SI$signif

# Control vs ProteinSynthesis50SInhibition
dist_Fn_novel_compounds_PC_PS50SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_PS50SI), c(ind_PC, ind_PS50SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_PC_PS50SI <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_PS50SI), ]
ano_Fn_novel_compounds_PC_PS50SI <- vegan::anosim(dist_Fn_novel_compounds_PC_PS50SI, metadata_Fn_dist_novel_compounds_PC_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["Control_ProteinSynthesis50SInhibition"] <- ano_Fn_novel_compounds_PC_PS50SI$signif

# Control vs ProteinSynthesistRNAInterference
dist_Fn_novel_compounds_PC_PStRI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_PStRI), c(ind_PC, ind_PStRI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_PC_PStRI <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_PStRI), ]
ano_Fn_novel_compounds_PC_PStRI <- vegan::anosim(dist_Fn_novel_compounds_PC_PStRI, metadata_Fn_dist_novel_compounds_PC_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["Control_ProteinSynthesistRNAInterference"] <- ano_Fn_novel_compounds_PC_PStRI$signif

# DNAReplication vs DNATranscription
dist_Fn_novel_compounds_DR_DT <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_DT), c(ind_DR, ind_DT)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_DR_DT <- metadata_Fn_dist_novel_compounds[c(ind_DR, ind_DT), ]
ano_Fn_novel_compounds_DR_DT <- vegan::anosim(dist_Fn_novel_compounds_DR_DT, metadata_Fn_dist_novel_compounds_DR_DT$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["DNAReplication_DNATranscription"] <- ano_Fn_novel_compounds_DR_DT$signif

# DNAReplication vs FolicAcidMetabolism
dist_Fn_novel_compounds_DR_FAM <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_FAM), c(ind_DR, ind_FAM)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_DR_FAM <- metadata_Fn_dist_novel_compounds[c(ind_DR, ind_FAM), ]
ano_Fn_novel_compounds_DR_FAM <- vegan::anosim(dist_Fn_novel_compounds_DR_FAM, metadata_Fn_dist_novel_compounds_DR_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["DNAReplication_FolicAcidMetabolism"] <- ano_Fn_novel_compounds_DR_FAM$signif

# DNAReplication vs Heat
dist_Fn_novel_compounds_DR_H <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_H), c(ind_DR, ind_H)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_DR_H <- metadata_Fn_dist_novel_compounds[c(ind_DR, ind_H), ]
ano_Fn_novel_compounds_DR_H <- vegan::anosim(dist_Fn_novel_compounds_DR_H, metadata_Fn_dist_novel_compounds_DR_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["DNAReplication_Heat"] <- ano_Fn_novel_compounds_DR_H$signif

# DNAReplication vs MembraneDisruption
dist_Fn_novel_compounds_DR_MD <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_MD), c(ind_DR, ind_MD)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_DR_MD <- metadata_Fn_dist_novel_compounds[c(ind_DR, ind_MD), ]
ano_Fn_novel_compounds_DR_MD <- vegan::anosim(dist_Fn_novel_compounds_DR_MD, metadata_Fn_dist_novel_compounds_DR_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["DNAReplication_MembraneDisruption"] <- ano_Fn_novel_compounds_DR_MD$signif

# DNAReplication vs ProteinSynthesis30SInhibition
dist_Fn_novel_compounds_DR_PS30SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_PS30SI), c(ind_DR, ind_PS30SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_DR_PS30SI <- metadata_Fn_dist_novel_compounds[c(ind_DR, ind_PS30SI), ]
ano_Fn_novel_compounds_DR_PS30SI <- vegan::anosim(dist_Fn_novel_compounds_DR_PS30SI, metadata_Fn_dist_novel_compounds_DR_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["DNAReplication_ProteinSynthesis30SInhibition"] <- ano_Fn_novel_compounds_DR_PS30SI$signif

# DNAReplication vs ProteinSynthesis50SInhibition
dist_Fn_novel_compounds_DR_PS50SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_PS50SI), c(ind_DR, ind_PS50SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_DR_PS50SI <- metadata_Fn_dist_novel_compounds[c(ind_DR, ind_PS50SI), ]
ano_Fn_novel_compounds_DR_PS50SI <- vegan::anosim(dist_Fn_novel_compounds_DR_PS50SI, metadata_Fn_dist_novel_compounds_DR_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["DNAReplication_ProteinSynthesis50SInhibition"] <- ano_Fn_novel_compounds_DR_PS50SI$signif

# DNAReplication vs ProteinSynthesistRNAInterference
dist_Fn_novel_compounds_DR_PStRI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_PStRI), c(ind_DR, ind_PStRI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_DR_PStRI <- metadata_Fn_dist_novel_compounds[c(ind_DR, ind_PStRI), ]
ano_Fn_novel_compounds_DR_PStRI <- vegan::anosim(dist_Fn_novel_compounds_DR_PStRI, metadata_Fn_dist_novel_compounds_DR_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["DNAReplication_ProteinSynthesistRNAInterference"] <- ano_Fn_novel_compounds_DR_PStRI$signif

# DNATranscription vs FolicAcidMetabolism
dist_Fn_novel_compounds_DT_FAM <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_FAM), c(ind_DT, ind_FAM)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_DT_FAM <- metadata_Fn_dist_novel_compounds[c(ind_DT, ind_FAM), ]
ano_Fn_novel_compounds_DT_FAM <- vegan::anosim(dist_Fn_novel_compounds_DT_FAM, metadata_Fn_dist_novel_compounds_DT_FAM$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["DNATranscription_FolicAcidMetabolism"] <- ano_Fn_novel_compounds_DT_FAM$signif

# DNATranscription vs Heat
dist_Fn_novel_compounds_DT_H <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_H), c(ind_DT, ind_H)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_DT_H <- metadata_Fn_dist_novel_compounds[c(ind_DT, ind_H), ]
ano_Fn_novel_compounds_DT_H <- vegan::anosim(dist_Fn_novel_compounds_DT_H, metadata_Fn_dist_novel_compounds_DT_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["DNATranscription_Heat"] <- ano_Fn_novel_compounds_DT_H$signif

# DNATranscription vs MembraneDisruption
dist_Fn_novel_compounds_DT_MD <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_MD), c(ind_DT, ind_MD)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_DT_MD <- metadata_Fn_dist_novel_compounds[c(ind_DT, ind_MD), ]
ano_Fn_novel_compounds_DT_MD <- vegan::anosim(dist_Fn_novel_compounds_DT_MD, metadata_Fn_dist_novel_compounds_DT_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["DNATranscription_MembraneDisruption"] <- ano_Fn_novel_compounds_DT_MD$signif

# DNATranscription vs ProteinSynthesis30SInhibition
dist_Fn_novel_compounds_DT_PS30SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_PS30SI), c(ind_DT, ind_PS30SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_DT_PS30SI <- metadata_Fn_dist_novel_compounds[c(ind_DT, ind_PS30SI), ]
ano_Fn_novel_compounds_DT_PS30SI <- vegan::anosim(dist_Fn_novel_compounds_DT_PS30SI, metadata_Fn_dist_novel_compounds_DT_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["DNATranscription_ProteinSynthesis30SInhibition"] <- ano_Fn_novel_compounds_DT_PS30SI$signif

# DNATranscription vs ProteinSynthesis50SInhibition
dist_Fn_novel_compounds_DT_PS50SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_PS50SI), c(ind_DT, ind_PS50SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_DT_PS50SI <- metadata_Fn_dist_novel_compounds[c(ind_DT, ind_PS50SI), ]
ano_Fn_novel_compounds_DT_PS50SI <- vegan::anosim(dist_Fn_novel_compounds_DT_PS50SI, metadata_Fn_dist_novel_compounds_DT_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["DNATranscription_ProteinSynthesis50SInhibition"] <- ano_Fn_novel_compounds_DT_PS50SI$signif

# DNATranscription vs ProteinSynthesistRNAInterference
dist_Fn_novel_compounds_DT_PStRI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_PStRI), c(ind_DT, ind_PStRI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_DT_PStRI <- metadata_Fn_dist_novel_compounds[c(ind_DT, ind_PStRI), ]
ano_Fn_novel_compounds_DT_PStRI <- vegan::anosim(dist_Fn_novel_compounds_DT_PStRI, metadata_Fn_dist_novel_compounds_DT_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["DNATranscription_ProteinSynthesistRNAInterference"] <- ano_Fn_novel_compounds_DT_PStRI$signif

# FolicAcidMetabolism vs Heat
dist_Fn_novel_compounds_FAM_H <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_H), c(ind_FAM, ind_H)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_FAM_H <- metadata_Fn_dist_novel_compounds[c(ind_FAM, ind_H), ]
ano_Fn_novel_compounds_FAM_H <- vegan::anosim(dist_Fn_novel_compounds_FAM_H, metadata_Fn_dist_novel_compounds_FAM_H$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["FolicAcidMetabolism_Heat"] <- ano_Fn_novel_compounds_FAM_H$signif

# FolicAcidMetabolism vs MembraneDisruption
dist_Fn_novel_compounds_FAM_MD <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_MD), c(ind_FAM, ind_MD)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_FAM_MD <- metadata_Fn_dist_novel_compounds[c(ind_FAM, ind_MD), ]
ano_Fn_novel_compounds_FAM_MD <- vegan::anosim(dist_Fn_novel_compounds_FAM_MD, metadata_Fn_dist_novel_compounds_FAM_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["FolicAcidMetabolism_MembraneDisruption"] <- ano_Fn_novel_compounds_FAM_MD$signif

# FolicAcidMetabolism vs ProteinSynthesis30SInhibition
dist_Fn_novel_compounds_FAM_PS30SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_PS30SI), c(ind_FAM, ind_PS30SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_FAM_PS30SI <- metadata_Fn_dist_novel_compounds[c(ind_FAM, ind_PS30SI), ]
ano_Fn_novel_compounds_FAM_PS30SI <- vegan::anosim(dist_Fn_novel_compounds_FAM_PS30SI, metadata_Fn_dist_novel_compounds_FAM_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["FolicAcidMetabolism_ProteinSynthesis30SInhibition"] <- ano_Fn_novel_compounds_FAM_PS30SI$signif

# FolicAcidMetabolism vs ProteinSynthesis50SInhibition
dist_Fn_novel_compounds_FAM_PS50SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_PS50SI), c(ind_FAM, ind_PS50SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_FAM_PS50SI <- metadata_Fn_dist_novel_compounds[c(ind_FAM, ind_PS50SI), ]
ano_Fn_novel_compounds_FAM_PS50SI <- vegan::anosim(dist_Fn_novel_compounds_FAM_PS50SI, metadata_Fn_dist_novel_compounds_FAM_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["FolicAcidMetabolism_ProteinSynthesis50SInhibition"] <- ano_Fn_novel_compounds_FAM_PS50SI$signif

# FolicAcidMetabolism vs ProteinSynthesistRNAInterference
dist_Fn_novel_compounds_FAM_PStRI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_PStRI), c(ind_FAM, ind_PStRI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_FAM_PStRI <- metadata_Fn_dist_novel_compounds[c(ind_FAM, ind_PStRI), ]
ano_Fn_novel_compounds_FAM_PStRI <- vegan::anosim(dist_Fn_novel_compounds_FAM_PStRI, metadata_Fn_dist_novel_compounds_FAM_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["FolicAcidMetabolism_ProteinSynthesistRNAInterference"] <- ano_Fn_novel_compounds_FAM_PStRI$signif

# Heat vs MembraneDisruption
dist_Fn_novel_compounds_H_MD <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_MD), c(ind_H, ind_MD)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_H_MD <- metadata_Fn_dist_novel_compounds[c(ind_H, ind_MD), ]
ano_Fn_novel_compounds_H_MD <- vegan::anosim(dist_Fn_novel_compounds_H_MD, metadata_Fn_dist_novel_compounds_H_MD$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["Heat_MembraneDisruption"] <- ano_Fn_novel_compounds_H_MD$signif

# Heat vs ProteinSynthesis30SInhibition
dist_Fn_novel_compounds_H_PS30SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_PS30SI), c(ind_H, ind_PS30SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_H_PS30SI <- metadata_Fn_dist_novel_compounds[c(ind_H, ind_PS30SI), ]
ano_Fn_novel_compounds_H_PS30SI <- vegan::anosim(dist_Fn_novel_compounds_H_PS30SI, metadata_Fn_dist_novel_compounds_H_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["Heat_ProteinSynthesis30SInhibition"] <- ano_Fn_novel_compounds_H_PS30SI$signif

# Heat vs ProteinSynthesis50SInhibition
dist_Fn_novel_compounds_H_PS50SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_PS50SI), c(ind_H, ind_PS50SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_H_PS50SI <- metadata_Fn_dist_novel_compounds[c(ind_H, ind_PS50SI), ]
ano_Fn_novel_compounds_H_PS50SI <- vegan::anosim(dist_Fn_novel_compounds_H_PS50SI, metadata_Fn_dist_novel_compounds_H_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["Heat_ProteinSynthesis50SInhibition"] <- ano_Fn_novel_compounds_H_PS50SI$signif

# Heat vs ProteinSynthesistRNAInterference
dist_Fn_novel_compounds_H_PStRI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_PStRI), c(ind_H, ind_PStRI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_H_PStRI <- metadata_Fn_dist_novel_compounds[c(ind_H, ind_PStRI), ]
ano_Fn_novel_compounds_H_PStRI <- vegan::anosim(dist_Fn_novel_compounds_H_PStRI, metadata_Fn_dist_novel_compounds_H_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["Heat_ProteinSynthesistRNAInterference"] <- ano_Fn_novel_compounds_H_PStRI$signif

# MembraneDisruption vs ProteinSynthesis30SInhibition
dist_Fn_novel_compounds_MD_PS30SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_MD, ind_PS30SI), c(ind_MD, ind_PS30SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_MD_PS30SI <- metadata_Fn_dist_novel_compounds[c(ind_MD, ind_PS30SI), ]
ano_Fn_novel_compounds_MD_PS30SI <- vegan::anosim(dist_Fn_novel_compounds_MD_PS30SI, metadata_Fn_dist_novel_compounds_MD_PS30SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["MembraneDisruption_ProteinSynthesis30SInhibition"] <- ano_Fn_novel_compounds_MD_PS30SI$signif

# MembraneDisruption vs ProteinSynthesis50SInhibition
dist_Fn_novel_compounds_MD_PS50SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_MD, ind_PS50SI), c(ind_MD, ind_PS50SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_MD_PS50SI <- metadata_Fn_dist_novel_compounds[c(ind_MD, ind_PS50SI), ]
ano_Fn_novel_compounds_MD_PS50SI <- vegan::anosim(dist_Fn_novel_compounds_MD_PS50SI, metadata_Fn_dist_novel_compounds_MD_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["MembraneDisruption_ProteinSynthesis50SInhibition"] <- ano_Fn_novel_compounds_MD_PS50SI$signif

# MembraneDisruption vs ProteinSynthesistRNAInterference
dist_Fn_novel_compounds_MD_PStRI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_MD, ind_PStRI), c(ind_MD, ind_PStRI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_MD_PStRI <- metadata_Fn_dist_novel_compounds[c(ind_MD, ind_PStRI), ]
ano_Fn_novel_compounds_MD_PStRI <- vegan::anosim(dist_Fn_novel_compounds_MD_PStRI, metadata_Fn_dist_novel_compounds_MD_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["MembraneDisruption_ProteinSynthesistRNAInterference"] <- ano_Fn_novel_compounds_MD_PStRI$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesis50SInhibition
dist_Fn_novel_compounds_PS30SI_PS50SI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PS30SI, ind_PS50SI), c(ind_PS30SI, ind_PS50SI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_PS30SI_PS50SI <- metadata_Fn_dist_novel_compounds[c(ind_PS30SI, ind_PS50SI), ]
ano_Fn_novel_compounds_PS30SI_PS50SI <- vegan::anosim(dist_Fn_novel_compounds_PS30SI_PS50SI, metadata_Fn_dist_novel_compounds_PS30SI_PS50SI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["ProteinSynthesis30SInhibition_ProteinSynthesis50SInhibition"] <- ano_Fn_novel_compounds_PS30SI_PS50SI$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesistRNAInterference
dist_Fn_novel_compounds_PS30SI_PStRI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PS30SI, ind_PStRI), c(ind_PS30SI, ind_PStRI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_PS30SI_PStRI <- metadata_Fn_dist_novel_compounds[c(ind_PS30SI, ind_PStRI), ]
ano_Fn_novel_compounds_PS30SI_PStRI <- vegan::anosim(dist_Fn_novel_compounds_PS30SI_PStRI, metadata_Fn_dist_novel_compounds_PS30SI_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["ProteinSynthesis30SInhibition_ProteinSynthesistRNAInterference"] <- ano_Fn_novel_compounds_PS30SI_PStRI$signif

# ProteinSynthesis50SInhibition vs ProteinSynthesistRNAInterference
dist_Fn_novel_compounds_PS50SI_PStRI <- dist_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(ind_PS50SI, ind_PStRI), c(ind_PS50SI, ind_PStRI)] %>% 
  as.dist()

metadata_Fn_dist_novel_compounds_PS50SI_PStRI <- metadata_Fn_dist_novel_compounds[c(ind_PS50SI, ind_PStRI), ]
ano_Fn_novel_compounds_PS50SI_PStRI <- vegan::anosim(dist_Fn_novel_compounds_PS50SI_PStRI, metadata_Fn_dist_novel_compounds_PS50SI_PStRI$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_Fn_novel_compounds["ProteinSynthesis50SInhibition_ProteinSynthesistRNAInterference"] <- ano_Fn_novel_compounds_PS50SI_PStRI$signif

# Correction for pairwise comparisons
pairwise_anosim_Fn_novel_compounds_bonferroni <- p.adjust(pairwise_anosim_Fn_novel_compounds, method = "bonferroni")
pairwise_anosim_Fn_novel_compounds_BH <- p.adjust(pairwise_anosim_Fn_novel_compounds, method = "BH")

pairwise_anosim_Fn_novel_compounds_df <- as.data.frame(pairwise_anosim_Fn_novel_compounds)
pairwise_anosim_Fn_novel_compounds_bonferroni_df <- as.data.frame(pairwise_anosim_Fn_novel_compounds_bonferroni)
pairwise_anosim_Fn_novel_compounds_BH_df <- as.data.frame(pairwise_anosim_Fn_novel_compounds_BH)

pairwise_anosim_Fn_novel_compounds_export <- merge(pairwise_anosim_Fn_novel_compounds_df, pairwise_anosim_Fn_novel_compounds_bonferroni_df, by = "row.names") %>% 
  column_to_rownames(var = "Row.names") %>% 
  merge(pairwise_anosim_Fn_novel_compounds_BH_df, by = "row.names") %>% 
  rename("classes" = "Row.names", "p" = "pairwise_anosim_Fn_novel_compounds", "p.adjust.bonferroni" = "pairwise_anosim_Fn_novel_compounds_bonferroni", "p.adjust.BH" = "pairwise_anosim_Fn_novel_compounds_BH")
write.csv2(file = "ANOSIM_Fn_novel_compounds_pairwise.csv", pairwise_anosim_Fn_novel_compounds_export) # Export results as csv file


### 5.2.3. Both strains combined ----
# Diversity analysis for both strains combined (indicate whether strain or treatment is dominant)
flowData_transformed_gated <- flowCore::rbind2(flowData_transformed_Av_gated, flowData_transformed_Fn_gated)
metadata <- rbind(metadata_Av, metadata_Fn)

# Normalization of data
summary <- fsApply(x = flowData_transformed_gated, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
max = max(summary[, "BL1-H"])
mytrans <- function(x) x/max

# Actual normalization of data
flowData_transformed_norm <- transform(flowData_transformed_gated,
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

# Select data without negative controls
flowData_transformed_norm_noneg <- flowData_transformed_norm[c(1:12, 16:87, 91:150)]
metadata_noneg <- metadata[c(1:12, 16:87, 91:150), ]
row.names(metadata_noneg) <- NULL

#### 5.2.3.1. Fingerprinting ----
# Randomly resample to a fixed number of cells
flowData_transformed_resample <- FCS_resample(flowData_transformed_norm_noneg, sample = 9775, replace = TRUE)

# Make metadata file for resampled data
names_flowData_transformed_resample <- sampleNames(flowData_transformed_resample)
metadata_resample <- metadata_noneg[match(names_flowData_transformed_resample, metadata_noneg$Filename), ]
row.names(metadata_resample) <- NULL

# For visualization purposes
size_metadata <- c(rep(3, 72), rep(8, 72))
metadata_resample <- cbind(metadata_resample, size_metadata)
colnames(metadata_resample)[colnames(metadata_resample) == "size_metadata"] <- "Size"

# Calculating fingerprint with bandwidth = 0.01
fbasis <- flowBasis(flowData_transformed_resample, param_div, nbin = 128, 
                    bw = 0.01, normalize = function(x) x)

#### 5.2.3.2. Ordination ----
# Beta-diversity assessment of fingerprint
beta.div <- beta_div_fcm(fbasis, ord.type = "PCoA")
beta.div_NMDS <- beta_div_fcm(fbasis, ord.type = "NMDS")

# Plot ordination
pbetadiv <- plot_beta_fcm_custom(beta.div, color = as.factor(metadata_resample$MOA), shape = as.factor(metadata_resample$Compound), size = as.factor(metadata_resample$Size))+ 
  theme_bw() +
  geom_point(alpha = 0.7)+
  labs(title = NULL, color = "Mode of Action", shape = "Compound", size = "Strain")+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_size_manual(values = c(3, 8),
                    labels = c(expression(italic("Actinomyces viscosus")), expression(italic("Fusobacterium nucleatum"))))+
  scale_shape_manual(values = c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        legend.text.align = 0)
print(pbetadiv)

pbetadiv_color <- plot_beta_fcm(beta.div, color = as.factor(metadata_resample$Strain), shape = as.factor(metadata_resample$Compound))+ 
  theme_bw() +
  geom_point(alpha = 0.7)+
  labs(title = NULL, color = "Strain", shape = "Compound")+
  scale_color_manual(values = c("blue", "red"),
                     labels = c(expression(italic("Actinomyces viscosus")), expression(italic("Fusobacterium nucleatum"))))+
  scale_shape_manual(values = c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        legend.text.align = 0)
print(pbetadiv_color)

pbetadiv_NMDS <- plot_beta_fcm_custom(beta.div_NMDS, color = as.factor(metadata_resample$MOA), shape = as.factor(metadata_resample$Compound), size = as.factor(metadata_resample$Size))+ 
  theme_bw() +
  geom_point(alpha = 0.7)+
  labs(title = NULL, color = "Mode of Action", shape = "Compound", size = "Strain", x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_size_manual(values = c(3, 8),
                    labels = c(expression(italic("Actinomyces viscosus")), expression(italic("Fusobacterium nucleatum"))))+
  scale_shape_manual(values = c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        legend.text.align = 0)
print(pbetadiv_NMDS)

pbetadiv_color_NMDS <- plot_beta_fcm(beta.div_NMDS, color = as.factor(metadata_resample$Strain), shape = as.factor(metadata_resample$Compound))+ 
  theme_bw() +
  geom_point(alpha = 0.7)+
  labs(title = NULL, color = "Strain", shape = "Compound", x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("blue", "red"),
                     labels = c(expression(italic("Actinomyces viscosus")), expression(italic("Fusobacterium nucleatum"))))+
  scale_shape_manual(values = c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        legend.text.align = 0)
print(pbetadiv_color_NMDS)


# 6. PhenoGMM ----

# Parameters to base the model on
paramGMM <- c("FSC-H", "SSC-H", "BL1-H", "BL3-H", "FSC-A", "SSC-A", "BL1-A", "BL3-A")

## 6.1. Av ----
### 6.1.1. Optimization model ----
# Select data exlcuding cephalothin
flowData_transformed_Av_GMM_novel_compounds <- flowData_transformed_Av_norm_noneg[c(1, 5:72)]

# Optimization of amount of phenotypes (Bayesian Information Criterion, BIC)
# Retain useful information only and pool all samples by MOA
fcs_PGMM_Av_novel_compounds <- Phenoflow::FCS_pool(flowData_transformed_Av_GMM_novel_compounds,
                                                   stub = c("20210902_Fabian_MOA_Av_Experiment_CellWallSynthesis",
                                                            "20210902_Fabian_MOA_Av_Experiment_Control",
                                                            "20210902_Fabian_MOA_Av_Experiment_DNAReplication",
                                                            "20210902_Fabian_MOA_Av_Experiment_DNATranscription",
                                                            "20210902_Fabian_MOA_Av_Experiment_FolicAcidMetabolism",
                                                            "20210902_Fabian_MOA_Av_Experiment_Heat",
                                                            "20210902_Fabian_MOA_Av_Experiment_MembraneDisruption",
                                                            "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis30SInhibition",
                                                            "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis50SInhibition",
                                                            "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesistRNAInterference"))
fcs_PGMM_Av_novel_compounds <- FCS_resample(fcs_PGMM_Av_novel_compounds, replace = TRUE, sample = 100000)
fcs_PGMM_Av_novel_compounds <- fcs_PGMM_Av_novel_compounds[, paramGMM]
fcs_PGMM_Av_novel_compounds <- Phenoflow::FCS_pool(fcs_PGMM_Av_novel_compounds, stub = "*")
fcs_PGMM_Av_novel_compounds <- fcs_PGMM_Av_novel_compounds[, paramGMM]
PhenoGMM_Av_noneg_novel_compounds <- Phenoflow::PhenoGMM(fcs_x = fcs_PGMM_Av_novel_compounds, param = paramGMM, downsample = FALSE, nG = 80, auto_nG = TRUE, nG_interval = 10, fcs_scale = FALSE, diagnostic_plot = TRUE)
saveRDS(object = PhenoGMM_Av_noneg_novel_compounds, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/Av_24h_PhenoGMM_8param_10to80per10_novel_compounds.rds")

# Visualization of BIC values
#PhenoGMM_Av_noneg_novel_compounds <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/Av_24h_PhenoGMM_8param_10to80per10_novel_compounds.rds") # Call for results that are to be plotted
NoClusters_PGMM <- dimnames(PhenoGMM_Av_noneg_novel_compounds[[2]]$BIC)[[1]]
BICValues_PGMM <- data.frame(as.matrix(PhenoGMM_Av_noneg_novel_compounds[[2]]$BIC)[1:length(NoClusters_PGMM), ])
BICValues_PGMM$NoClusters <- rownames(BICValues_PGMM)
BICValues_PGMM <- reshape2::melt(BICValues_PGMM, id.vars = "NoClusters")
colnames(BICValues_PGMM) <- c("NoClusters", "ModelType", "BIC")
BICValues_PGMM$NoClusters <- as.numeric(BICValues_PGMM$NoClusters)
BICValues_PGMM <- BICValues_PGMM[!is.na(BICValues_PGMM$BIC), ] # Remove NA values
BICValues_PGMM$ModelType <- droplevels(BICValues_PGMM$ModelType, except = unique(BICValues_PGMM$ModelType)) # Remove levels that are not being used

p_BIC_PGMM_Av_novel_compounds <- BICValues_PGMM %>% 
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
print(p_BIC_PGMM_Av_novel_compounds)

### 6.1.2. Allocation data to model ----
## Make model for optimum clusters determined in previous step (excluding Cephalothin)
NoClusters_PGMM_Av_novel_compounds <- 50
PhenoGMM_Av_noneg_fixedclust_novel_compounds <- Phenoflow::PhenoGMM(fcs_x = fcs_PGMM_Av_novel_compounds, param = paramGMM, downsample = FALSE, nG = NoClusters_PGMM_Av_novel_compounds, fcs_scale = FALSE, diagnostic_plot = TRUE)
saveRDS(object = PhenoGMM_Av_noneg_fixedclust_novel_compounds, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/Av_24h_PhenoGMM_8param_50clust_novel_compounds.rds")
#PhenoGMM_Av_noneg_fixedclust_novel_compounds <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/Av_24h_PhenoGMM_8param_50clust_novel_compounds.rds") # Call for results that are to be plotted

# Applying mask to data
testPred_Av_noneg_novel_compounds <- PhenoMaskGMM(fcs_x = flowData_transformed_Av_norm_noneg, gmm = PhenoGMM_Av_noneg_fixedclust_novel_compounds, fcs_scale = FALSE)

results_Av_novel_compounds <- testPred_Av_noneg_novel_compounds[[1]]
rownames(results_Av_novel_compounds) <- sampleNames(flowData_transformed_Av_norm_noneg)
results_Av_novel_compounds <- select(results_Av_novel_compounds, -c(Sample_names))
results_Av_novel_compounds[is.na(results_Av_novel_compounds)] <- 0 # Replace NA values by 0
results_Av_novel_compounds[1:ncol(results_Av_novel_compounds)] <- lapply(results_Av_novel_compounds[1:ncol(results_Av_novel_compounds)], as.numeric)
# Normalize abundances to sum = 1
results_Av_rel_novel_compounds <- sweep(results_Av_novel_compounds, MARGIN = 1, rowSums(results_Av_novel_compounds), `/`)

### 6.1.3. Visualization PhenoGMM models ----
source("/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/gmm_identifier.R")

## Relative abundance
# Visualization PhenoGMM model (excluding Cephalothin in construction of GMM model)
results_Av_rel_novel_compounds <- results_Av_rel_novel_compounds %>% 
  add_column(Sample_name = metadata_Av_noneg$Sample_complete)
results_Av_rel_novel_compounds_melted <- melt(results_Av_rel_novel_compounds, id.vars = "Sample_name")

p_PhenoGMM_Av_novel_compounds <- results_Av_rel_novel_compounds_melted %>% 
  ggplot(data = ., aes(x = Sample_name, y = value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(y = "Relative abundance", x = "Sample", title = NULL, fill = "Cluster GMM model") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        legend.text.align = 0)
print(p_PhenoGMM_Av_novel_compounds)

## Flow cytometric space
PGMM_Av_visual_pooled <- gmm_identifier(flow_set = fcs_PGMM_Av_novel_compounds, PhenoGMM_model = PhenoGMM_Av_noneg_fixedclust_novel_compounds)
PGMM_Av_visual <- gmm_identifier(flow_set = flowData_transformed_Av_norm_noneg, PhenoGMM_model = PhenoGMM_Av_noneg_fixedclust_novel_compounds)

p_PGMM_Av_cytometric_pooled <- ggplot(data.frame(PGMM_Av_visual_pooled[[1]]))+
  geom_point(aes(x = BL1.H, y = BL3.H, color = as.factor(gmm_id)), size = 0.1)+
  xlim(0.4, 1)+
  ylim(0.2, 1)+
  labs(x = 'BL1-H', y = 'BL3-H', title = NULL, color = 'Cluster GMM model')+
  theme_bw()
print(p_PGMM_Av_cytometric_pooled)

# Flow cytometric space for one sample
PGMM_Av_visual_select <- PGMM_Av_visual[37]

p_PGMM_Av_cytometric <- ggplot(data.frame(PGMM_Av_visual_select[[1]]))+
  geom_point(aes(x = BL1.H, y = BL3.H, color = as.factor(gmm_id)), size = 0.1)+
  xlim(0.4, 1)+
  ylim(0.2, 1)+
  labs(x = 'BL1-H', y = 'BL3-H', title = metadata_Av_noneg$Sample_complete[37], color = 'Cluster GMM model')+
  theme_bw()+
  theme(plot.title = element_text(size = 26))
print(p_PGMM_Av_cytometric)

# Flow cytometric space for All samples
plist_PGMM_Av <- list()

for (i in 1:length(PGMM_Av_visual)) {
  PGMM_Av_visual_loop <- PGMM_Av_visual[i]
  plist_PGMM_Av[[i]] <- ggplot(data.frame(PGMM_Av_visual_loop[[1]]))+
    geom_point(aes(x = BL1.H, y = BL3.H, color = as.factor(gmm_id)), size = 0.1)+
    xlim(0.4, 1)+
    ylim(0.2, 1)+
    labs(x = 'BL1-H', y = 'BL3-H', title = metadata_Av_noneg$Sample_complete[i], color = 'Cluster GMM model')+
    theme_bw()+
    theme(legend.position = "none",
          plot.title = element_text(size = 12))
}

grobs_PGMM_Av <- ggplotGrob(p_PGMM_Av_cytometric_pooled)$grobs
legend_PGMM_grobs_Av <- grobs_PGMM_Av[[which(sapply(grobs_PGMM_Av, function(x) x$name) == "guide-box")]]

nCol_PGMM_Av <- floor(sqrt(length(plist_PGMM_Av)))
p_PGMM_grid_Av <- cowplot::plot_grid(plotlist = plist_PGMM_Av, ncol = nCol_PGMM_Av)
p_PGMM_sample_Av <- cowplot::plot_grid(p_PGMM_grid_Av, legend_PGMM_grobs_Av, ncol = 2, rel_widths = c(9, 0.4))
print(p_PGMM_sample_Av)

### 6.1.4. Ordination of GMM output ----
# Model excluding cephalothin
rownames_GMM_Av_novel_compounds <- row.names(results_Av_novel_compounds)
new_rownames_GMM_Av_novel_compounds <- gsub("^.*?nt_", "", rownames_GMM_Av_novel_compounds)
new_rownames_GMM_Av_novel_compounds <- gsub("_100_SGPI", "", new_rownames_GMM_Av_novel_compounds)
new_rownames_GMM_Av_novel_compounds <- gsub(".fcs", "", new_rownames_GMM_Av_novel_compounds)
results_Av_ordination_novel_compounds <- results_Av_novel_compounds
rownames(results_Av_ordination_novel_compounds) <- new_rownames_GMM_Av_novel_compounds
results_Av_ordination_novel_compounds <- as.matrix(results_Av_ordination_novel_compounds)

distance_GMM_Av_novel_compounds <- vegan::vegdist(x = results_Av_ordination_novel_compounds, method = "bray", binary = FALSE)
mds_GMM_Av_novel_compounds <- stats::cmdscale(distance_GMM_Av_novel_compounds, k = 2, eig = TRUE, add = TRUE)
NMDS_GMM_Av_novel_compounds <- vegan::metaMDS(distance_GMM_Av_novel_compounds, autotransform = FALSE, k = 2, trymax = 100)

plot_PCoA_GMM_Av_novel_compounds <- plot_beta_fcm(mds_GMM_Av_novel_compounds, color = as.factor(metadata_Av_dist$MOA), shape = as.factor(metadata_Av_dist$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL)+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(plot_PCoA_GMM_Av_novel_compounds)

plot_NMDS_GMM_Av_novel_compounds <- plot_beta_fcm(NMDS_GMM_Av_novel_compounds, color = as.factor(metadata_Av_dist$MOA), shape = as.factor(metadata_Av_dist$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL, x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(plot_NMDS_GMM_Av_novel_compounds)

# Model excluding cephalothin and ordination excluding cephalothin
results_Av_ordination_novel_compounds_noceph <- results_Av_ordination_novel_compounds[c(1, 5:72), ]
distance_GMM_Av_novel_compounds_noceph <- vegan::vegdist(x = results_Av_ordination_novel_compounds_noceph, method = "bray", binary = FALSE)
mds_GMM_Av_novel_compounds_noceph <- stats::cmdscale(distance_GMM_Av_novel_compounds_noceph, k = 2, eig = TRUE, add = TRUE)
NMDS_GMM_Av_novel_compounds_noceph <- vegan::metaMDS(distance_GMM_Av_novel_compounds_noceph, autotransform = FALSE, k = 2, trymax = 100)

metadata_Av_dist_noceph <- metadata_Av_dist[c(1, 5:72), ]

plot_PCoA_GMM_Av_novel_compounds_noceph <- plot_beta_fcm(mds_GMM_Av_novel_compounds_noceph, color = as.factor(metadata_Av_dist_noceph$MOA), shape = as.factor(metadata_Av_dist_noceph$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL)+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:4, 6:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(plot_PCoA_GMM_Av_novel_compounds_noceph)

plot_NMDS_GMM_Av_novel_compounds_noceph <- plot_beta_fcm(NMDS_GMM_Av_novel_compounds_noceph, color = as.factor(metadata_Av_dist_noceph$MOA), shape = as.factor(metadata_Av_dist_noceph$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL, x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:4, 6:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(plot_NMDS_GMM_Av_novel_compounds_noceph)


### 6.1.5. Statistics ----
### ANOSIM
ano_GMM_Av_novel_compounds <- vegan::anosim(distance_GMM_Av_novel_compounds, metadata_Av_dist$MOA, distance = "bray", permutations = permutations)
ano_GMM_Av_novel_compounds_export <- data.frame("ANOSIM Statistic R" = ano_GMM_Av_novel_compounds$statistic, "P" = ano_GMM_Av_novel_compounds$signif, "permutations" = ano_GMM_Av_novel_compounds$permutations)
write.csv2(file = "ANOSIM_GMM_Av_novel_compounds.csv", ano_GMM_Av_novel_compounds_export) # Export results as csv file
plot_pairwise_anosim_GMM_Av_novel_compounds <- plot(ano_GMM_Av_novel_compounds, xaxt = "n", main = "Dissimilarity ranks GMM Av (ANOSIM)", xlab = "Class", ylab = "Dissimilarity rank")
axis(1, at = 1:11, labels = c("Between", "CWS", "DNARepl", "DNATrans", "FAM", "Heat", "MD", "Control", "PS-30SInh", "PS-50SInh", "PS-tRNAInt"))

# Exclude cephalothin from calculation
ano_GMM_Av_novel_compounds_no_ceph <- vegan::anosim(distance_GMM_Av_novel_compounds_noceph, metadata_Av_dist_noceph$MOA, distance = "bray", permutations = permutations)
ano_GMM_Av_novel_compounds_no_ceph_export <- data.frame("ANOSIM Statistic R" = ano_GMM_Av_novel_compounds_no_ceph$statistic, "P" = ano_GMM_Av_novel_compounds_no_ceph$signif, "permutations" = ano_GMM_Av_novel_compounds_no_ceph$permutations)
write.csv2(file = "ANOSIM_GMM_Av_novel_compounds_no_ceph.csv", ano_GMM_Av_novel_compounds_no_ceph_export) # Export results as csv file
plot_pairwise_anosim_GMM_Av_novel_compounds_no_ceph <- plot(ano_GMM_Av_novel_compounds_no_ceph, xaxt = "n", main = "Dissimilarity ranks GMM Av (ANOSIM)", xlab = "Class", ylab = "Dissimilarity rank")
axis(1, at = 1:11, labels = c("Between", "CWS", "DNARepl", "DNATrans", "FAM", "Heat", "MD", "Control", "PS-30SInh", "PS-50SInh", "PS-tRNAInt"))

# Pairwise comparison different MOA
pairwise_anosim_GMM_Av <- numeric()

# CellWallSynthesis vs Control
distance_GMM_Av_CWS_PC_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 13:18), c(1:12, 13:18)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_PC_inclCeph <- metadata_Av_dist[c(1:12, 13:18), ]
ano_GMM_Av_CWS_PC_inclCeph <- vegan::anosim(distance_GMM_Av_CWS_PC_inclCeph, metadata_GMM_Av_dist_CWS_PC_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_Control_inclCephalothin"] <- ano_GMM_Av_CWS_PC_inclCeph$signif

# CellWallSynthesis vs DNAReplication
distance_GMM_Av_CWS_DR_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 19:27), c(1:12, 19:27)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_DR_inclCeph <- metadata_Av_dist[c(1:12, 19:27), ]
ano_GMM_Av_CWS_DR_inclCeph <- vegan::anosim(distance_GMM_Av_CWS_DR_inclCeph, metadata_GMM_Av_dist_CWS_DR_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_DNAReplication_inclCephalothin"] <- ano_GMM_Av_CWS_DR_inclCeph$signif

# CellWallSynthesis vs DNATranscription
distance_GMM_Av_CWS_DT_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 28:33), c(1:12, 28:33)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_DT_inclCeph <- metadata_Av_dist[c(1:12, 28:33), ]
ano_GMM_Av_CWS_DT_inclCeph <- vegan::anosim(distance_GMM_Av_CWS_DT_inclCeph, metadata_GMM_Av_dist_CWS_DT_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_DNATranscription_inclCephalothin"] <- ano_GMM_Av_CWS_DT_inclCeph$signif

# CellWallSynthesis vs FolicAcidMetabolism
distance_GMM_Av_CWS_FAM_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 34:36), c(1:12, 34:36)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_FAM_inclCeph <- metadata_Av_dist[c(1:12, 34:36), ]
ano_GMM_Av_CWS_FAM_inclCeph <- vegan::anosim(distance_GMM_Av_CWS_FAM_inclCeph, metadata_GMM_Av_dist_CWS_FAM_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_FolicAcidMetabolism_inclCephalothin"] <- ano_GMM_Av_CWS_FAM_inclCeph$signif

# CellWallSynthesis vs Heat
distance_GMM_Av_CWS_H_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 37:39), c(1:12, 37:39)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_H_inclCeph <- metadata_Av_dist[c(1:12, 37:39), ]
ano_GMM_Av_CWS_H_inclCeph <- vegan::anosim(distance_GMM_Av_CWS_H_inclCeph, metadata_GMM_Av_dist_CWS_H_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_Heat_inclCephalothin"] <- ano_GMM_Av_CWS_H_inclCeph$signif

# CellWallSynthesis vs MembraneDisruption
distance_GMM_Av_CWS_MD_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 40:51), c(1:12, 40:51)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_MD_inclCeph <- metadata_Av_dist[c(1:12, 40:51), ]
ano_GMM_Av_CWS_MD_inclCeph <- vegan::anosim(distance_GMM_Av_CWS_MD_inclCeph, metadata_GMM_Av_dist_CWS_MD_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_MembraneDisruption_inclCephalothin"] <- ano_GMM_Av_CWS_MD_inclCeph$signif

# CellWallSynthesis vs ProteinSynthesis30SInhibition
distance_GMM_Av_CWS_PS30SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 52:60), c(1:12, 52:60)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_PS30SI_inclCeph <- metadata_Av_dist[c(1:12, 52:60), ]
ano_GMM_Av_CWS_PS30SI_inclCeph <- vegan::anosim(distance_GMM_Av_CWS_PS30SI_inclCeph, metadata_GMM_Av_dist_CWS_PS30SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_ProteinSynthesis30SInhibition_inclCephalothin"] <- ano_GMM_Av_CWS_PS30SI_inclCeph$signif

# CellWallSynthesis vs ProteinSynthesis50SInhibition
distance_GMM_Av_CWS_PS50SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 61:69), c(1:12, 61:69)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_PS50SI_inclCeph <- metadata_Av_dist[c(1:12, 61:69), ]
ano_GMM_Av_CWS_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Av_CWS_PS50SI_inclCeph, metadata_GMM_Av_dist_CWS_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Av_CWS_PS50SI_inclCeph$signif

# CellWallSynthesis vs ProteinSynthesistRNAInterference
distance_GMM_Av_CWS_PStRI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 70:72), c(1:12, 70:72)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_PStRI_inclCeph <- metadata_Av_dist[c(1:12, 70:72), ]
ano_GMM_Av_CWS_PStRI_inclCeph <- vegan::anosim(distance_GMM_Av_CWS_PStRI_inclCeph, metadata_GMM_Av_dist_CWS_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Av_CWS_PStRI_inclCeph$signif

# Control vs DNAReplication
distance_GMM_Av_PC_DR_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 19:27), c(13:18, 19:27)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_DR_inclCeph <- metadata_Av_dist[c(13:18, 19:27), ]
ano_GMM_Av_PC_DR_inclCeph <- vegan::anosim(distance_GMM_Av_PC_DR_inclCeph, metadata_GMM_Av_dist_PC_DR_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_DNAReplication_inclCephalothin"] <- ano_GMM_Av_PC_DR_inclCeph$signif

# Control vs DNATranscription
distance_GMM_Av_PC_DT_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 28:33), c(13:18, 28:33)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_DT_inclCeph <- metadata_Av_dist[c(13:18, 28:33), ]
ano_GMM_Av_PC_DT_inclCeph <- vegan::anosim(distance_GMM_Av_PC_DT_inclCeph, metadata_GMM_Av_dist_PC_DT_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_DNATranscription_inclCephalothin"] <- ano_GMM_Av_PC_DT_inclCeph$signif

# Control vs FolicAcidMetabolism
distance_GMM_Av_PC_FAM_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 34:36), c(13:18, 34:36)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_FAM_inclCeph <- metadata_Av_dist[c(13:18, 34:36), ]
ano_GMM_Av_PC_FAM_inclCeph <- vegan::anosim(distance_GMM_Av_PC_FAM_inclCeph, metadata_GMM_Av_dist_PC_FAM_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_FolicAcidMetabolism_inclCephalothin"] <- ano_GMM_Av_PC_FAM_inclCeph$signif

# Control vs Heat
distance_GMM_Av_PC_H_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 37:39), c(13:18, 37:39)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_H_inclCeph <- metadata_Av_dist[c(13:18, 37:39), ]
ano_GMM_Av_PC_H_inclCeph <- vegan::anosim(distance_GMM_Av_PC_H_inclCeph, metadata_GMM_Av_dist_PC_H_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_Heat_inclCephalothin"] <- ano_GMM_Av_PC_H_inclCeph$signif

# Control vs MembraneDisruption
distance_GMM_Av_PC_MD_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 40:51), c(13:18, 40:51)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_MD_inclCeph <- metadata_Av_dist[c(13:18, 40:51), ]
ano_GMM_Av_PC_MD_inclCeph <- vegan::anosim(distance_GMM_Av_PC_MD_inclCeph, metadata_GMM_Av_dist_PC_MD_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_MembraneDisruption_inclCephalothin"] <- ano_GMM_Av_PC_MD_inclCeph$signif

# Control vs ProteinSynthesis30SInhibition
distance_GMM_Av_PC_PS30SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 52:60), c(13:18, 52:60)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_PS30SI_inclCeph <- metadata_Av_dist[c(13:18, 52:60), ]
ano_GMM_Av_PC_PS30SI_inclCeph <- vegan::anosim(distance_GMM_Av_PC_PS30SI_inclCeph, metadata_GMM_Av_dist_PC_PS30SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_ProteinSynthesis30SInhibition_inclCephalothin"] <- ano_GMM_Av_PC_PS30SI_inclCeph$signif

# Control vs ProteinSynthesis50SInhibition
distance_GMM_Av_PC_PS50SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 61:69), c(13:18, 61:69)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_PS50SI_inclCeph <- metadata_Av_dist[c(13:18, 61:69), ]
ano_GMM_Av_PC_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Av_PC_PS50SI_inclCeph, metadata_GMM_Av_dist_PC_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Av_PC_PS50SI_inclCeph$signif

# Control vs ProteinSynthesistRNAInterference
distance_GMM_Av_PC_PStRI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 70:72), c(13:18, 70:72)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_PStRI_inclCeph <- metadata_Av_dist[c(13:18, 70:72), ]
ano_GMM_Av_PC_PStRI_inclCeph <- vegan::anosim(distance_GMM_Av_PC_PStRI_inclCeph, metadata_GMM_Av_dist_PC_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Av_PC_PStRI_inclCeph$signif

# DNAReplication vs DNATranscription
distance_GMM_Av_DR_DT_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(19:27, 28:33), c(19:27, 28:33)] %>% 
  as.dist()

metadata_GMM_Av_dist_DR_DT_inclCeph <- metadata_Av_dist[c(19:27, 28:33), ]
ano_GMM_Av_DR_DT_inclCeph <- vegan::anosim(distance_GMM_Av_DR_DT_inclCeph, metadata_GMM_Av_dist_DR_DT_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNAReplication_DNATranscription_inclCephalothin"] <- ano_GMM_Av_DR_DT_inclCeph$signif

# DNAReplication vs FolicAcidMetabolism
distance_GMM_Av_DR_FAM_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(19:27, 34:36), c(19:27, 34:36)] %>% 
  as.dist()

metadata_GMM_Av_dist_DR_FAM_inclCeph <- metadata_Av_dist[c(19:27, 34:36), ]
ano_GMM_Av_DR_FAM_inclCeph <- vegan::anosim(distance_GMM_Av_DR_FAM_inclCeph, metadata_GMM_Av_dist_DR_FAM_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNAReplication_FolicAcidMetabolism_inclCephalothin"] <- ano_GMM_Av_DR_FAM_inclCeph$signif

# DNAReplication vs Heat
distance_GMM_Av_DR_H_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(19:27, 37:39), c(19:27, 37:39)] %>% 
  as.dist()

metadata_GMM_Av_dist_DR_H_inclCeph <- metadata_Av_dist[c(19:27, 37:39), ]
ano_GMM_Av_DR_H_inclCeph <- vegan::anosim(distance_GMM_Av_DR_H_inclCeph, metadata_GMM_Av_dist_DR_H_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNAReplication_Heat_inclCephalothin"] <- ano_GMM_Av_DR_H_inclCeph$signif

# DNAReplication vs MembraneDisruption
distance_GMM_Av_DR_MD_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(19:27, 40:51), c(19:27, 40:51)] %>% 
  as.dist()

metadata_GMM_Av_dist_DR_MD_inclCeph <- metadata_Av_dist[c(19:27, 40:51), ]
ano_GMM_Av_DR_MD_inclCeph <- vegan::anosim(distance_GMM_Av_DR_MD_inclCeph, metadata_GMM_Av_dist_DR_MD_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNAReplication_MembraneDisruption_inclCephalothin"] <- ano_GMM_Av_DR_MD_inclCeph$signif

# DNAReplication vs ProteinSynthesis30SInhibition
distance_GMM_Av_DR_PS30SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(19:27, 52:60), c(19:27, 52:60)] %>% 
  as.dist()

metadata_GMM_Av_dist_DR_PS30SI_inclCeph <- metadata_Av_dist[c(19:27, 52:60), ]
ano_GMM_Av_DR_PS30SI_inclCeph <- vegan::anosim(distance_GMM_Av_DR_PS30SI_inclCeph, metadata_GMM_Av_dist_DR_PS30SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNAReplication_ProteinSynthesis30SInhibition_inclCephalothin"] <- ano_GMM_Av_DR_PS30SI_inclCeph$signif

# DNAReplication vs ProteinSynthesis50SInhibition
distance_GMM_Av_DR_PS50SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(19:27, 61:69), c(19:27, 61:69)] %>% 
  as.dist()

metadata_GMM_Av_dist_DR_PS50SI_inclCeph <- metadata_Av_dist[c(19:27, 61:69), ]
ano_GMM_Av_DR_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Av_DR_PS50SI_inclCeph, metadata_GMM_Av_dist_DR_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNAReplication_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Av_DR_PS50SI_inclCeph$signif

# DNAReplication vs ProteinSynthesistRNAInterference
distance_GMM_Av_DR_PStRI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(19:27, 70:72), c(19:27, 70:72)] %>% 
  as.dist()

metadata_GMM_Av_dist_DR_PStRI_inclCeph <- metadata_Av_dist[c(19:27, 70:72), ]
ano_GMM_Av_DR_PStRI_inclCeph <- vegan::anosim(distance_GMM_Av_DR_PStRI_inclCeph, metadata_GMM_Av_dist_DR_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNAReplication_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Av_DR_PStRI_inclCeph$signif

# DNATranscription vs FolicAcidMetabolism
distance_GMM_Av_DT_FAM_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(28:33, 34:36), c(28:33, 34:36)] %>% 
  as.dist()

metadata_GMM_Av_dist_DT_FAM_inclCeph <- metadata_Av_dist[c(28:33, 34:36), ]
ano_GMM_Av_DT_FAM_inclCeph <- vegan::anosim(distance_GMM_Av_DT_FAM_inclCeph, metadata_GMM_Av_dist_DT_FAM_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNATranscription_FolicAcidMetabolism_inclCephalothin"] <- ano_GMM_Av_DT_FAM_inclCeph$signif

# DNATranscription vs Heat
distance_GMM_Av_DT_H_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(28:33, 37:39), c(28:33, 37:39)] %>% 
  as.dist()

metadata_GMM_Av_dist_DT_H_inclCeph <- metadata_Av_dist[c(28:33, 37:39), ]
ano_GMM_Av_DT_H_inclCeph <- vegan::anosim(distance_GMM_Av_DT_H_inclCeph, metadata_GMM_Av_dist_DT_H_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNATranscription_Heat_inclCephalothin"] <- ano_GMM_Av_DT_H_inclCeph$signif

# DNATranscription vs MembraneDisruption
distance_GMM_Av_DT_MD_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(28:33, 40:51), c(28:33, 40:51)] %>% 
  as.dist()

metadata_GMM_Av_dist_DT_MD_inclCeph <- metadata_Av_dist[c(28:33, 40:51), ]
ano_GMM_Av_DT_MD_inclCeph <- vegan::anosim(distance_GMM_Av_DT_MD_inclCeph, metadata_GMM_Av_dist_DT_MD_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNATranscription_MembraneDisruption_inclCephalothin"] <- ano_GMM_Av_DT_MD_inclCeph$signif

# DNATranscription vs ProteinSynthesis30SInhibition
distance_GMM_Av_DT_PS30SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(28:33, 52:60), c(28:33, 52:60)] %>% 
  as.dist()

metadata_GMM_Av_dist_DT_PS30SI_inclCeph <- metadata_Av_dist[c(28:33, 52:60), ]
ano_GMM_Av_DT_PS30SI_inclCeph <- vegan::anosim(distance_GMM_Av_DT_PS30SI_inclCeph, metadata_GMM_Av_dist_DT_PS30SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNATranscription_ProteinSynthesis30SInhibition_inclCephalothin"] <- ano_GMM_Av_DT_PS30SI_inclCeph$signif

# DNATranscription vs ProteinSynthesis50SInhibition
distance_GMM_Av_DT_PS50SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(28:33, 61:69), c(28:33, 61:69)] %>% 
  as.dist()

metadata_GMM_Av_dist_DT_PS50SI_inclCeph <- metadata_Av_dist[c(28:33, 61:69), ]
ano_GMM_Av_DT_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Av_DT_PS50SI_inclCeph, metadata_GMM_Av_dist_DT_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNATranscription_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Av_DT_PS50SI_inclCeph$signif

# DNATranscription vs ProteinSynthesistRNAInterference
distance_GMM_Av_DT_PStRI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(28:33, 70:72), c(28:33, 70:72)] %>% 
  as.dist()

metadata_GMM_Av_dist_DT_PStRI_inclCeph <- metadata_Av_dist[c(28:33, 70:72), ]
ano_GMM_Av_DT_PStRI_inclCeph <- vegan::anosim(distance_GMM_Av_DT_PStRI_inclCeph, metadata_GMM_Av_dist_DT_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNATranscription_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Av_DT_PStRI_inclCeph$signif

# FolicAcidMetabolism vs Heat
distance_GMM_Av_FAM_H_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(34:36, 37:39), c(34:36, 37:39)] %>% 
  as.dist()

metadata_GMM_Av_dist_FAM_H_inclCeph <- metadata_Av_dist[c(34:36, 37:39), ]
ano_GMM_Av_FAM_H_inclCeph <- vegan::anosim(distance_GMM_Av_FAM_H_inclCeph, metadata_GMM_Av_dist_FAM_H_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["FolicAcidMetabolism_Heat_inclCephalothin"] <- ano_GMM_Av_FAM_H_inclCeph$signif

# FolicAcidMetabolism vs MembraneDisruption
distance_GMM_Av_FAM_MD_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(34:36, 40:51), c(34:36, 40:51)] %>% 
  as.dist()

metadata_GMM_Av_dist_FAM_MD_inclCeph <- metadata_Av_dist[c(34:36, 40:51), ]
ano_GMM_Av_FAM_MD_inclCeph <- vegan::anosim(distance_GMM_Av_FAM_MD_inclCeph, metadata_GMM_Av_dist_FAM_MD_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["FolicAcidMetabolism_MembraneDisruption_inclCephalothin"] <- ano_GMM_Av_FAM_MD_inclCeph$signif

# FolicAcidMetabolism vs ProteinSynthesis30SInhibition
distance_GMM_Av_FAM_PS30SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(34:36, 52:60), c(34:36, 52:60)] %>% 
  as.dist()

metadata_GMM_Av_dist_FAM_PS30SI_inclCeph <- metadata_Av_dist[c(34:36, 52:60), ]
ano_GMM_Av_FAM_PS30SI_inclCeph <- vegan::anosim(distance_GMM_Av_FAM_PS30SI_inclCeph, metadata_GMM_Av_dist_FAM_PS30SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["FolicAcidMetabolism_ProteinSynthesis30SInhibition_inclCephalothin"] <- ano_GMM_Av_FAM_PS30SI_inclCeph$signif

# FolicAcidMetabolism vs ProteinSynthesis50SInhibition
distance_GMM_Av_FAM_PS50SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(34:36, 61:69), c(34:36, 61:69)] %>% 
  as.dist()

metadata_GMM_Av_dist_FAM_PS50SI_inclCeph <- metadata_Av_dist[c(34:36, 61:69), ]
ano_GMM_Av_FAM_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Av_FAM_PS50SI_inclCeph, metadata_GMM_Av_dist_FAM_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["FolicAcidMetabolism_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Av_FAM_PS50SI_inclCeph$signif

# FolicAcidMetabolism vs ProteinSynthesistRNAInterference
distance_GMM_Av_FAM_PStRI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(34:36, 70:72), c(34:36, 70:72)] %>% 
  as.dist()

metadata_GMM_Av_dist_FAM_PStRI_inclCeph <- metadata_Av_dist[c(34:36, 70:72), ]
ano_GMM_Av_FAM_PStRI_inclCeph <- vegan::anosim(distance_GMM_Av_FAM_PStRI_inclCeph, metadata_GMM_Av_dist_FAM_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["FolicAcidMetabolism_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Av_FAM_PStRI_inclCeph$signif

# Heat vs MembraneDisruption
distance_GMM_Av_H_MD_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(37:39, 40:51), c(37:39, 40:51)] %>% 
  as.dist()

metadata_GMM_Av_dist_H_MD_inclCeph <- metadata_Av_dist[c(37:39, 40:51), ]
ano_GMM_Av_H_MD_inclCeph <- vegan::anosim(distance_GMM_Av_H_MD_inclCeph, metadata_GMM_Av_dist_H_MD_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Heat_MembraneDisruption_inclCephalothin"] <- ano_GMM_Av_H_MD_inclCeph$signif

# Heat vs ProteinSynthesis30SInhibition
distance_GMM_Av_H_PS30SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(37:39, 52:60), c(37:39, 52:60)] %>% 
  as.dist()

metadata_GMM_Av_dist_H_PS30SI_inclCeph <- metadata_Av_dist[c(37:39, 52:60), ]
ano_GMM_Av_H_PS30SI_inclCeph <- vegan::anosim(distance_GMM_Av_H_PS30SI_inclCeph, metadata_GMM_Av_dist_H_PS30SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Heat_ProteinSynthesis30SInhibition_inclCephalothin"] <- ano_GMM_Av_H_PS30SI_inclCeph$signif

# Heat vs ProteinSynthesis50SInhibition
distance_GMM_Av_H_PS50SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(37:39, 61:69), c(37:39, 61:69)] %>% 
  as.dist()

metadata_GMM_Av_dist_H_PS50SI_inclCeph <- metadata_Av_dist[c(37:39, 61:69), ]
ano_GMM_Av_H_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Av_H_PS50SI_inclCeph, metadata_GMM_Av_dist_H_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Heat_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Av_H_PS50SI_inclCeph$signif

# Heat vs ProteinSynthesistRNAInterference
distance_GMM_Av_H_PStRI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(37:39, 70:72), c(37:39, 70:72)] %>% 
  as.dist()

metadata_GMM_Av_dist_H_PStRI_inclCeph <- metadata_Av_dist[c(37:39, 70:72), ]
ano_GMM_Av_H_PStRI_inclCeph <- vegan::anosim(distance_GMM_Av_H_PStRI_inclCeph, metadata_GMM_Av_dist_H_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Heat_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Av_H_PStRI_inclCeph$signif

# MembraneDisruption vs ProteinSynthesis30SInhibition
distance_GMM_Av_MD_PS30SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(40:51, 52:60), c(40:51, 52:60)] %>% 
  as.dist()

metadata_GMM_Av_dist_MD_PS30SI_inclCeph <- metadata_Av_dist[c(40:51, 52:60), ]
ano_GMM_Av_MD_PS30SI_inclCeph <- vegan::anosim(distance_GMM_Av_MD_PS30SI_inclCeph, metadata_GMM_Av_dist_MD_PS30SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["MembraneDisruption_ProteinSynthesis30SInhibition_inclCephalothin"] <- ano_GMM_Av_MD_PS30SI_inclCeph$signif

# MembraneDisruption vs ProteinSynthesis50SInhibition
distance_GMM_Av_MD_PS50SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(40:51, 61:69), c(40:51, 61:69)] %>% 
  as.dist()

metadata_GMM_Av_dist_MD_PS50SI_inclCeph <- metadata_Av_dist[c(40:51, 61:69), ]
ano_GMM_Av_MD_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Av_MD_PS50SI_inclCeph, metadata_GMM_Av_dist_MD_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["MembraneDisruption_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Av_MD_PS50SI_inclCeph$signif

# MembraneDisruption vs ProteinSynthesistRNAInterference
distance_GMM_Av_MD_PStRI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(40:51, 70:72), c(40:51, 70:72)] %>% 
  as.dist()

metadata_GMM_Av_dist_MD_PStRI_inclCeph <- metadata_Av_dist[c(40:51, 70:72), ]
ano_GMM_Av_MD_PStRI_inclCeph <- vegan::anosim(distance_GMM_Av_MD_PStRI_inclCeph, metadata_GMM_Av_dist_MD_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["MembraneDisruption_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Av_MD_PStRI_inclCeph$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesis50SInhibition
distance_GMM_Av_PS30SI_PS50SI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(52:60, 61:69), c(52:60, 61:69)] %>% 
  as.dist()

metadata_GMM_Av_dist_PS30SI_PS50SI_inclCeph <- metadata_Av_dist[c(52:60, 61:69), ]
ano_GMM_Av_PS30SI_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Av_PS30SI_PS50SI_inclCeph, metadata_GMM_Av_dist_PS30SI_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["ProteinSynthesis30SInhibition_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Av_PS30SI_PS50SI_inclCeph$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesistRNAInterference
distance_GMM_Av_PS30SI_PStRI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(52:60, 70:72), c(52:60, 70:72)] %>% 
  as.dist()

metadata_GMM_Av_dist_PS30SI_PStRI_inclCeph <- metadata_Av_dist[c(52:60, 70:72), ]
ano_GMM_Av_PS30SI_PStRI_inclCeph <- vegan::anosim(distance_GMM_Av_PS30SI_PStRI_inclCeph, metadata_GMM_Av_dist_PS30SI_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["ProteinSynthesis30SInhibition_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Av_PS30SI_PStRI_inclCeph$signif

# ProteinSynthesis50SInhibition vs ProteinSynthesistRNAInterference
distance_GMM_Av_PS50SI_PStRI_inclCeph <- distance_GMM_Av_novel_compounds %>% 
  as.matrix() %>% 
  .[c(61:69, 70:72), c(61:69, 70:72)] %>% 
  as.dist()

metadata_GMM_Av_dist_PS50SI_PStRI_inclCeph <- metadata_Av_dist[c(61:69, 70:72), ]
ano_GMM_Av_PS50SI_PStRI_inclCeph <- vegan::anosim(distance_GMM_Av_PS50SI_PStRI_inclCeph, metadata_GMM_Av_dist_PS50SI_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["ProteinSynthesis50SInhibition_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Av_PS50SI_PStRI_inclCeph$signif

## Pairwise comparison excluding cephalothin

# CellWallSynthesis vs Control
distance_GMM_Av_CWS_PC_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PC), c(ind_CWS, ind_PC)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_PC_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_PC), ]
ano_GMM_Av_CWS_PC_exclCeph <- vegan::anosim(distance_GMM_Av_CWS_PC_exclCeph, metadata_GMM_Av_dist_CWS_PC_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_Control_exclCephalothin"] <- ano_GMM_Av_CWS_PC_exclCeph$signif

# CellWallSynthesis vs DNAReplication
distance_GMM_Av_CWS_DR_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_DR), c(ind_CWS, ind_DR)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_DR_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_DR), ]
ano_GMM_Av_CWS_DR_exclCeph <- vegan::anosim(distance_GMM_Av_CWS_DR_exclCeph, metadata_GMM_Av_dist_CWS_DR_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_DNAReplication_exclCephalothin"] <- ano_GMM_Av_CWS_DR_exclCeph$signif

# CellWallSynthesis vs DNATranscription
distance_GMM_Av_CWS_DT_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_DT), c(ind_CWS, ind_DT)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_DT_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_DT), ]
ano_GMM_Av_CWS_DT_exclCeph <- vegan::anosim(distance_GMM_Av_CWS_DT_exclCeph, metadata_GMM_Av_dist_CWS_DT_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_DNATranscription_exclCephalothin"] <- ano_GMM_Av_CWS_DT_exclCeph$signif

# CellWallSynthesis vs FolicAcidMetabolism
distance_GMM_Av_CWS_FAM_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_FAM), c(ind_CWS, ind_FAM)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_FAM_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_FAM), ]
ano_GMM_Av_CWS_FAM_exclCeph <- vegan::anosim(distance_GMM_Av_CWS_FAM_exclCeph, metadata_GMM_Av_dist_CWS_FAM_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_FolicAcidMetabolism_exclCephalothin"] <- ano_GMM_Av_CWS_FAM_exclCeph$signif

# CellWallSynthesis vs Heat
distance_GMM_Av_CWS_H_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_H), c(ind_CWS, ind_H)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_H_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_H), ]
ano_GMM_Av_CWS_H_exclCeph <- vegan::anosim(distance_GMM_Av_CWS_H_exclCeph, metadata_GMM_Av_dist_CWS_H_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_Heat_exclCephalothin"] <- ano_GMM_Av_CWS_H_exclCeph$signif

# CellWallSynthesis vs MembraneDisruption
distance_GMM_Av_CWS_MD_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_MD), c(ind_CWS, ind_MD)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_MD_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_MD), ]
ano_GMM_Av_CWS_MD_exclCeph <- vegan::anosim(distance_GMM_Av_CWS_MD_exclCeph, metadata_GMM_Av_dist_CWS_MD_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_MembraneDisruption_exclCephalothin"] <- ano_GMM_Av_CWS_MD_exclCeph$signif

# CellWallSynthesis vs ProteinSynthesis30SInhibition
distance_GMM_Av_CWS_PS30SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PS30SI), c(ind_CWS, ind_PS30SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_PS30SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_PS30SI), ]
ano_GMM_Av_CWS_PS30SI_exclCeph <- vegan::anosim(distance_GMM_Av_CWS_PS30SI_exclCeph, metadata_GMM_Av_dist_CWS_PS30SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_ProteinSynthesis30SInhibition_exclCephalothin"] <- ano_GMM_Av_CWS_PS30SI_exclCeph$signif

# CellWallSynthesis vs ProteinSynthesis50SInhibition
distance_GMM_Av_CWS_PS50SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PS50SI), c(ind_CWS, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_PS50SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_PS50SI), ]
ano_GMM_Av_CWS_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Av_CWS_PS50SI_exclCeph, metadata_GMM_Av_dist_CWS_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Av_CWS_PS50SI_exclCeph$signif

# CellWallSynthesis vs ProteinSynthesistRNAInterference
distance_GMM_Av_CWS_PStRI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PStRI), c(ind_CWS, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Av_dist_CWS_PStRI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_CWS, ind_PStRI), ]
ano_GMM_Av_CWS_PStRI_exclCeph <- vegan::anosim(distance_GMM_Av_CWS_PStRI_exclCeph, metadata_GMM_Av_dist_CWS_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["CellWallSynthesis_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Av_CWS_PStRI_exclCeph$signif

# Control vs DNAReplication
distance_GMM_Av_PC_DR_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_DR), c(ind_PC, ind_DR)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_DR_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_DR), ]
ano_GMM_Av_PC_DR_exclCeph <- vegan::anosim(distance_GMM_Av_PC_DR_exclCeph, metadata_GMM_Av_dist_PC_DR_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_DNAReplication_exclCephalothin"] <- ano_GMM_Av_PC_DR_exclCeph$signif

# Control vs DNATranscription
distance_GMM_Av_PC_DT_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_DT), c(ind_PC, ind_DT)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_DT_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_DT), ]
ano_GMM_Av_PC_DT_exclCeph <- vegan::anosim(distance_GMM_Av_PC_DT_exclCeph, metadata_GMM_Av_dist_PC_DT_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_DNATranscription_exclCephalothin"] <- ano_GMM_Av_PC_DT_exclCeph$signif

# Control vs FolicAcidMetabolism
distance_GMM_Av_PC_FAM_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_FAM), c(ind_PC, ind_FAM)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_FAM_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_FAM), ]
ano_GMM_Av_PC_FAM_exclCeph <- vegan::anosim(distance_GMM_Av_PC_FAM_exclCeph, metadata_GMM_Av_dist_PC_FAM_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_FolicAcidMetabolism_exclCephalothin"] <- ano_GMM_Av_PC_FAM_exclCeph$signif

# Control vs Heat
distance_GMM_Av_PC_H_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_H), c(ind_PC, ind_H)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_H_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_H), ]
ano_GMM_Av_PC_H_exclCeph <- vegan::anosim(distance_GMM_Av_PC_H_exclCeph, metadata_GMM_Av_dist_PC_H_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_Heat_exclCephalothin"] <- ano_GMM_Av_PC_H_exclCeph$signif

# Control vs MembraneDisruption
distance_GMM_Av_PC_MD_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_MD), c(ind_PC, ind_MD)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_MD_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_MD), ]
ano_GMM_Av_PC_MD_exclCeph <- vegan::anosim(distance_GMM_Av_PC_MD_exclCeph, metadata_GMM_Av_dist_PC_MD_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_MembraneDisruption_exclCephalothin"] <- ano_GMM_Av_PC_MD_exclCeph$signif

# Control vs ProteinSynthesis30SInhibition
distance_GMM_Av_PC_PS30SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_PS30SI), c(ind_PC, ind_PS30SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_PS30SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_PS30SI), ]
ano_GMM_Av_PC_PS30SI_exclCeph <- vegan::anosim(distance_GMM_Av_PC_PS30SI_exclCeph, metadata_GMM_Av_dist_PC_PS30SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_ProteinSynthesis30SInhibition_exclCephalothin"] <- ano_GMM_Av_PC_PS30SI_exclCeph$signif

# Control vs ProteinSynthesis50SInhibition
distance_GMM_Av_PC_PS50SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_PS50SI), c(ind_PC, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_PS50SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_PS50SI), ]
ano_GMM_Av_PC_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Av_PC_PS50SI_exclCeph, metadata_GMM_Av_dist_PC_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Av_PC_PS50SI_exclCeph$signif

# Control vs ProteinSynthesistRNAInterference
distance_GMM_Av_PC_PStRI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_PStRI), c(ind_PC, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Av_dist_PC_PStRI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_PC, ind_PStRI), ]
ano_GMM_Av_PC_PStRI_exclCeph <- vegan::anosim(distance_GMM_Av_PC_PStRI_exclCeph, metadata_GMM_Av_dist_PC_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Control_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Av_PC_PStRI_exclCeph$signif

# DNAReplication vs DNATranscription
distance_GMM_Av_DR_DT_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_DT), c(ind_DR, ind_DT)] %>% 
  as.dist()

metadata_GMM_Av_dist_DR_DT_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_DR, ind_DT), ]
ano_GMM_Av_DR_DT_exclCeph <- vegan::anosim(distance_GMM_Av_DR_DT_exclCeph, metadata_GMM_Av_dist_DR_DT_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNAReplication_DNATranscription_exclCephalothin"] <- ano_GMM_Av_DR_DT_exclCeph$signif

# DNAReplication vs FolicAcidMetabolism
distance_GMM_Av_DR_FAM_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_FAM), c(ind_DR, ind_FAM)] %>% 
  as.dist()

metadata_GMM_Av_dist_DR_FAM_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_DR, ind_FAM), ]
ano_GMM_Av_DR_FAM_exclCeph <- vegan::anosim(distance_GMM_Av_DR_FAM_exclCeph, metadata_GMM_Av_dist_DR_FAM_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNAReplication_FolicAcidMetabolism_exclCephalothin"] <- ano_GMM_Av_DR_FAM_exclCeph$signif

# DNAReplication vs Heat
distance_GMM_Av_DR_H_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_H), c(ind_DR, ind_H)] %>% 
  as.dist()

metadata_GMM_Av_dist_DR_H_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_DR, ind_H), ]
ano_GMM_Av_DR_H_exclCeph <- vegan::anosim(distance_GMM_Av_DR_H_exclCeph, metadata_GMM_Av_dist_DR_H_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNAReplication_Heat_exclCephalothin"] <- ano_GMM_Av_DR_H_exclCeph$signif

# DNAReplication vs MembraneDisruption
distance_GMM_Av_DR_MD_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_MD), c(ind_DR, ind_MD)] %>% 
  as.dist()

metadata_GMM_Av_dist_DR_MD_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_DR, ind_MD), ]
ano_GMM_Av_DR_MD_exclCeph <- vegan::anosim(distance_GMM_Av_DR_MD_exclCeph, metadata_GMM_Av_dist_DR_MD_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNAReplication_MembraneDisruption_exclCephalothin"] <- ano_GMM_Av_DR_MD_exclCeph$signif

# DNAReplication vs ProteinSynthesis30SInhibition
distance_GMM_Av_DR_PS30SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_PS30SI), c(ind_DR, ind_PS30SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_DR_PS30SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_DR, ind_PS30SI), ]
ano_GMM_Av_DR_PS30SI_exclCeph <- vegan::anosim(distance_GMM_Av_DR_PS30SI_exclCeph, metadata_GMM_Av_dist_DR_PS30SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNAReplication_ProteinSynthesis30SInhibition_exclCephalothin"] <- ano_GMM_Av_DR_PS30SI_exclCeph$signif

# DNAReplication vs ProteinSynthesis50SInhibition
distance_GMM_Av_DR_PS50SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_PS50SI), c(ind_DR, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_DR_PS50SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_DR, ind_PS50SI), ]
ano_GMM_Av_DR_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Av_DR_PS50SI_exclCeph, metadata_GMM_Av_dist_DR_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNAReplication_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Av_DR_PS50SI_exclCeph$signif

# DNAReplication vs ProteinSynthesistRNAInterference
distance_GMM_Av_DR_PStRI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_PStRI), c(ind_DR, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Av_dist_DR_PStRI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_DR, ind_PStRI), ]
ano_GMM_Av_DR_PStRI_exclCeph <- vegan::anosim(distance_GMM_Av_DR_PStRI_exclCeph, metadata_GMM_Av_dist_DR_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNAReplication_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Av_DR_PStRI_exclCeph$signif

# DNATranscription vs FolicAcidMetabolism
distance_GMM_Av_DT_FAM_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_FAM), c(ind_DT, ind_FAM)] %>% 
  as.dist()

metadata_GMM_Av_dist_DT_FAM_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_DT, ind_FAM), ]
ano_GMM_Av_DT_FAM_exclCeph <- vegan::anosim(distance_GMM_Av_DT_FAM_exclCeph, metadata_GMM_Av_dist_DT_FAM_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNATranscription_FolicAcidMetabolism_exclCephalothin"] <- ano_GMM_Av_DT_FAM_exclCeph$signif

# DNATranscription vs Heat
distance_GMM_Av_DT_H_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_H), c(ind_DT, ind_H)] %>% 
  as.dist()

metadata_GMM_Av_dist_DT_H_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_DT, ind_H), ]
ano_GMM_Av_DT_H_exclCeph <- vegan::anosim(distance_GMM_Av_DT_H_exclCeph, metadata_GMM_Av_dist_DT_H_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNATranscription_Heat_exclCephalothin"] <- ano_GMM_Av_DT_H_exclCeph$signif

# DNATranscription vs MembraneDisruption
distance_GMM_Av_DT_MD_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_MD), c(ind_DT, ind_MD)] %>% 
  as.dist()

metadata_GMM_Av_dist_DT_MD_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_DT, ind_MD), ]
ano_GMM_Av_DT_MD_exclCeph <- vegan::anosim(distance_GMM_Av_DT_MD_exclCeph, metadata_GMM_Av_dist_DT_MD_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNATranscription_MembraneDisruption_exclCephalothin"] <- ano_GMM_Av_DT_MD_exclCeph$signif

# DNATranscription vs ProteinSynthesis30SInhibition
distance_GMM_Av_DT_PS30SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_PS30SI), c(ind_DT, ind_PS30SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_DT_PS30SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_DT, ind_PS30SI), ]
ano_GMM_Av_DT_PS30SI_exclCeph <- vegan::anosim(distance_GMM_Av_DT_PS30SI_exclCeph, metadata_GMM_Av_dist_DT_PS30SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNATranscription_ProteinSynthesis30SInhibition_exclCephalothin"] <- ano_GMM_Av_DT_PS30SI_exclCeph$signif

# DNATranscription vs ProteinSynthesis50SInhibition
distance_GMM_Av_DT_PS50SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_PS50SI), c(ind_DT, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_DT_PS50SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_DT, ind_PS50SI), ]
ano_GMM_Av_DT_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Av_DT_PS50SI_exclCeph, metadata_GMM_Av_dist_DT_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNATranscription_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Av_DT_PS50SI_exclCeph$signif

# DNATranscription vs ProteinSynthesistRNAInterference
distance_GMM_Av_DT_PStRI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_PStRI), c(ind_DT, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Av_dist_DT_PStRI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_DT, ind_PStRI), ]
ano_GMM_Av_DT_PStRI_exclCeph <- vegan::anosim(distance_GMM_Av_DT_PStRI_exclCeph, metadata_GMM_Av_dist_DT_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["DNATranscription_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Av_DT_PStRI_exclCeph$signif

# FolicAcidMetabolism vs Heat
distance_GMM_Av_FAM_H_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_H), c(ind_FAM, ind_H)] %>% 
  as.dist()

metadata_GMM_Av_dist_FAM_H_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_FAM, ind_H), ]
ano_GMM_Av_FAM_H_exclCeph <- vegan::anosim(distance_GMM_Av_FAM_H_exclCeph, metadata_GMM_Av_dist_FAM_H_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["FolicAcidMetabolism_Heat_exclCephalothin"] <- ano_GMM_Av_FAM_H_exclCeph$signif

# FolicAcidMetabolism vs MembraneDisruption
distance_GMM_Av_FAM_MD_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_MD), c(ind_FAM, ind_MD)] %>% 
  as.dist()

metadata_GMM_Av_dist_FAM_MD_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_FAM, ind_MD), ]
ano_GMM_Av_FAM_MD_exclCeph <- vegan::anosim(distance_GMM_Av_FAM_MD_exclCeph, metadata_GMM_Av_dist_FAM_MD_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["FolicAcidMetabolism_MembraneDisruption_exclCephalothin"] <- ano_GMM_Av_FAM_MD_exclCeph$signif

# FolicAcidMetabolism vs ProteinSynthesis30SInhibition
distance_GMM_Av_FAM_PS30SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_PS30SI), c(ind_FAM, ind_PS30SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_FAM_PS30SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_FAM, ind_PS30SI), ]
ano_GMM_Av_FAM_PS30SI_exclCeph <- vegan::anosim(distance_GMM_Av_FAM_PS30SI_exclCeph, metadata_GMM_Av_dist_FAM_PS30SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["FolicAcidMetabolism_ProteinSynthesis30SInhibition_exclCephalothin"] <- ano_GMM_Av_FAM_PS30SI_exclCeph$signif

# FolicAcidMetabolism vs ProteinSynthesis50SInhibition
distance_GMM_Av_FAM_PS50SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_PS50SI), c(ind_FAM, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_FAM_PS50SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_FAM, ind_PS50SI), ]
ano_GMM_Av_FAM_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Av_FAM_PS50SI_exclCeph, metadata_GMM_Av_dist_FAM_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["FolicAcidMetabolism_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Av_FAM_PS50SI_exclCeph$signif

# FolicAcidMetabolism vs ProteinSynthesistRNAInterference
distance_GMM_Av_FAM_PStRI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_PStRI), c(ind_FAM, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Av_dist_FAM_PStRI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_FAM, ind_PStRI), ]
ano_GMM_Av_FAM_PStRI_exclCeph <- vegan::anosim(distance_GMM_Av_FAM_PStRI_exclCeph, metadata_GMM_Av_dist_FAM_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["FolicAcidMetabolism_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Av_FAM_PStRI_exclCeph$signif

# Heat vs MembraneDisruption
distance_GMM_Av_H_MD_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_MD), c(ind_H, ind_MD)] %>% 
  as.dist()

metadata_GMM_Av_dist_H_MD_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_H, ind_MD), ]
ano_GMM_Av_H_MD_exclCeph <- vegan::anosim(distance_GMM_Av_H_MD_exclCeph, metadata_GMM_Av_dist_H_MD_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Heat_MembraneDisruption_exclCephalothin"] <- ano_GMM_Av_H_MD_exclCeph$signif

# Heat vs ProteinSynthesis30SInhibition
distance_GMM_Av_H_PS30SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_PS30SI), c(ind_H, ind_PS30SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_H_PS30SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_H, ind_PS30SI), ]
ano_GMM_Av_H_PS30SI_exclCeph <- vegan::anosim(distance_GMM_Av_H_PS30SI_exclCeph, metadata_GMM_Av_dist_H_PS30SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Heat_ProteinSynthesis30SInhibition_exclCephalothin"] <- ano_GMM_Av_H_PS30SI_exclCeph$signif

# Heat vs ProteinSynthesis50SInhibition
distance_GMM_Av_H_PS50SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_PS50SI), c(ind_H, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_H_PS50SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_H, ind_PS50SI), ]
ano_GMM_Av_H_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Av_H_PS50SI_exclCeph, metadata_GMM_Av_dist_H_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Heat_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Av_H_PS50SI_exclCeph$signif

# Heat vs ProteinSynthesistRNAInterference
distance_GMM_Av_H_PStRI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_PStRI), c(ind_H, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Av_dist_H_PStRI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_H, ind_PStRI), ]
ano_GMM_Av_H_PStRI_exclCeph <- vegan::anosim(distance_GMM_Av_H_PStRI_exclCeph, metadata_GMM_Av_dist_H_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["Heat_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Av_H_PStRI_exclCeph$signif

# MembraneDisruption vs ProteinSynthesis30SInhibition
distance_GMM_Av_MD_PS30SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_MD, ind_PS30SI), c(ind_MD, ind_PS30SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_MD_PS30SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_MD, ind_PS30SI), ]
ano_GMM_Av_MD_PS30SI_exclCeph <- vegan::anosim(distance_GMM_Av_MD_PS30SI_exclCeph, metadata_GMM_Av_dist_MD_PS30SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["MembraneDisruption_ProteinSynthesis30SInhibition_exclCephalothin"] <- ano_GMM_Av_MD_PS30SI_exclCeph$signif

# MembraneDisruption vs ProteinSynthesis50SInhibition
distance_GMM_Av_MD_PS50SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_MD, ind_PS50SI), c(ind_MD, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_MD_PS50SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_MD, ind_PS50SI), ]
ano_GMM_Av_MD_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Av_MD_PS50SI_exclCeph, metadata_GMM_Av_dist_MD_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["MembraneDisruption_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Av_MD_PS50SI_exclCeph$signif

# MembraneDisruption vs ProteinSynthesistRNAInterference
distance_GMM_Av_MD_PStRI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_MD, ind_PStRI), c(ind_MD, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Av_dist_MD_PStRI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_MD, ind_PStRI), ]
ano_GMM_Av_MD_PStRI_exclCeph <- vegan::anosim(distance_GMM_Av_MD_PStRI_exclCeph, metadata_GMM_Av_dist_MD_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["MembraneDisruption_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Av_MD_PStRI_exclCeph$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesis50SInhibition
distance_GMM_Av_PS30SI_PS50SI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PS30SI, ind_PS50SI), c(ind_PS30SI, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Av_dist_PS30SI_PS50SI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_PS30SI, ind_PS50SI), ]
ano_GMM_Av_PS30SI_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Av_PS30SI_PS50SI_exclCeph, metadata_GMM_Av_dist_PS30SI_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["ProteinSynthesis30SInhibition_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Av_PS30SI_PS50SI_exclCeph$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesistRNAInterference
distance_GMM_Av_PS30SI_PStRI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PS30SI, ind_PStRI), c(ind_PS30SI, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Av_dist_PS30SI_PStRI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_PS30SI, ind_PStRI), ]
ano_GMM_Av_PS30SI_PStRI_exclCeph <- vegan::anosim(distance_GMM_Av_PS30SI_PStRI_exclCeph, metadata_GMM_Av_dist_PS30SI_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["ProteinSynthesis30SInhibition_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Av_PS30SI_PStRI_exclCeph$signif

# ProteinSynthesis50SInhibition vs ProteinSynthesistRNAInterference
distance_GMM_Av_PS50SI_PStRI_exclCeph <- distance_GMM_Av_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PS50SI, ind_PStRI), c(ind_PS50SI, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Av_dist_PS50SI_PStRI_exclCeph <- metadata_Av_dist_novel_compounds[c(ind_PS50SI, ind_PStRI), ]
ano_GMM_Av_PS50SI_PStRI_exclCeph <- vegan::anosim(distance_GMM_Av_PS50SI_PStRI_exclCeph, metadata_GMM_Av_dist_PS50SI_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Av["ProteinSynthesis50SInhibition_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Av_PS50SI_PStRI_exclCeph$signif

# Correction for pairwise comparisons
pairwise_anosim_GMM_Av_bonferroni <- p.adjust(pairwise_anosim_GMM_Av, method = "bonferroni")
pairwise_anosim_GMM_Av_BH <- p.adjust(pairwise_anosim_GMM_Av, method = "BH")

pairwise_anosim_GMM_Av_df <- as.data.frame(pairwise_anosim_GMM_Av)
pairwise_anosim_GMM_Av_bonferroni_df <- as.data.frame(pairwise_anosim_GMM_Av_bonferroni)
pairwise_anosim_GMM_Av_BH_df <- as.data.frame(pairwise_anosim_GMM_Av_BH)

pairwise_anosim_GMM_Av_export <- merge(pairwise_anosim_GMM_Av_df, pairwise_anosim_GMM_Av_bonferroni_df, by = "row.names") %>% 
  column_to_rownames(var = "Row.names") %>% 
  merge(pairwise_anosim_GMM_Av_BH_df, by = "row.names") %>% 
  rename("classes" = "Row.names", "p" = "pairwise_anosim_GMM_Av", "p.adjust.bonferroni" = "pairwise_anosim_GMM_Av_bonferroni", "p.adjust.BH" = "pairwise_anosim_GMM_Av_BH")
write.csv2(file = "ANOSIM_GMM_Av_pairwise_novel_compounds.csv", pairwise_anosim_GMM_Av_export) # Export results as csv file


## 6.2. Fn ----
### 6.2.1. Optimization model ----
# Select data exlcuding cephalothin
flowData_transformed_Fn_GMM_novel_compounds <- flowData_transformed_Fn_norm_noneg[c(1, 5:72)]

# Optimization of amount of phenotypes (Bayesian Information Criterion, BIC)
# Retain usefull information only and pool all samples by MOA
fcs_PGMM_Fn_novel_compounds <- Phenoflow::FCS_pool(flowData_transformed_Fn_GMM_novel_compounds,
                                                   stub = c("20210902_Fabian_MOA_Fn_Experiment_CellWallSynthesis",
                                                            "20210902_Fabian_MOA_Fn_Experiment_Control",
                                                            "20210902_Fabian_MOA_Fn_Experiment_DNAReplication",
                                                            "20210902_Fabian_MOA_Fn_Experiment_DNATranscription",
                                                            "20210902_Fabian_MOA_Fn_Experiment_FolicAcidMetabolism",
                                                            "20210902_Fabian_MOA_Fn_Experiment_Heat",
                                                            "20210902_Fabian_MOA_Fn_Experiment_MembraneDisruption",
                                                            "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis30SInhibition",
                                                            "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis50SInhibition",
                                                            "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesistRNAInterference"))
fcs_PGMM_Fn_novel_compounds <- FCS_resample(fcs_PGMM_Fn_novel_compounds, replace = TRUE, sample = 100000)
fcs_PGMM_Fn_novel_compounds <- fcs_PGMM_Fn_novel_compounds[, paramGMM]
fcs_PGMM_Fn_novel_compounds <- Phenoflow::FCS_pool(fcs_PGMM_Fn_novel_compounds, stub = "*")
fcs_PGMM_Fn_novel_compounds <- fcs_PGMM_Fn_novel_compounds[, paramGMM]
PhenoGMM_Fn_noneg_novel_compounds <- Phenoflow::PhenoGMM(fcs_x = fcs_PGMM_Fn_novel_compounds, param = paramGMM, downsample = FALSE, nG = 80, auto_nG = TRUE, nG_interval = 10, fcs_scale = FALSE, diagnostic_plot = TRUE)
saveRDS(object = PhenoGMM_Fn_noneg_novel_compounds, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/Fn_24h_PhenoGMM_8param_10to80per10_novel_compounds.rds")

# Visualization of BIC values
#PhenoGMM_Fn_noneg_novel_compounds <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/Fn_PhenoGMM_8param_10to80per10_novel_compounds.rds") # Call for results that are to be plotted
NoClusters_PGMM <- dimnames(PhenoGMM_Fn_noneg_novel_compounds[[2]]$BIC)[[1]]
BICValues_PGMM <- data.frame(as.matrix(PhenoGMM_Fn_noneg_novel_compounds[[2]]$BIC)[1:length(NoClusters_PGMM), ])
BICValues_PGMM$NoClusters <- rownames(BICValues_PGMM)
BICValues_PGMM <- reshape2::melt(BICValues_PGMM, id.vars = "NoClusters")
colnames(BICValues_PGMM) <- c("NoClusters", "ModelType", "BIC")
BICValues_PGMM$NoClusters <- as.numeric(BICValues_PGMM$NoClusters)
BICValues_PGMM <- BICValues_PGMM[!is.na(BICValues_PGMM$BIC), ] # Remove NA values
BICValues_PGMM$ModelType <- droplevels(BICValues_PGMM$ModelType, except = unique(BICValues_PGMM$ModelType)) # Remove levels that are not being used

p_BIC_PGMM_Fn_novel_compounds <- BICValues_PGMM %>% 
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
print(p_BIC_PGMM_Fn_novel_compounds)

### 6.2.2. Allocation data to model ----
## Make model for optimum clusters determined in previous step (excluding Cephalothin)
NoClusters_PGMM_Fn_novel_compounds <- 40
PhenoGMM_Fn_noneg_fixedclust_novel_compounds <- Phenoflow::PhenoGMM(fcs_x = fcs_PGMM_Fn_novel_compounds, param = paramGMM, downsample = FALSE, nG = NoClusters_PGMM_Fn_novel_compounds, fcs_scale = FALSE, diagnostic_plot = TRUE)
saveRDS(object = PhenoGMM_Fn_noneg_fixedclust_novel_compounds, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/Fn_24h_PhenoGMM_8param_40clust_novel_compounds.rds")
#PhenoGMM_Fn_noneg_fixedclust_novel_compounds <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/Fn_24h_PhenoGMM_8param_40clust_novel_compounds.rds") # Call for results that are to be plotted

# Applying mask to data
testPred_Fn_noneg_novel_compounds <- PhenoMaskGMM(fcs_x = flowData_transformed_Fn_norm_noneg, gmm = PhenoGMM_Fn_noneg_fixedclust_novel_compounds, fcs_scale = FALSE)

results_Fn_novel_compounds <- testPred_Fn_noneg_novel_compounds[[1]]
rownames(results_Fn_novel_compounds) <- sampleNames(flowData_transformed_Fn_norm_noneg)
results_Fn_novel_compounds <- select(results_Fn_novel_compounds, -c(Sample_names))
results_Fn_novel_compounds[is.na(results_Fn_novel_compounds)] <- 0 # Replace NA values by 0
results_Fn_novel_compounds[1:ncol(results_Fn_novel_compounds)] <- lapply(results_Fn_novel_compounds[1:ncol(results_Fn_novel_compounds)], as.numeric)
# Normalize abundances to sum = 1
results_Fn_rel_novel_compounds <- sweep(results_Fn_novel_compounds, MARGIN = 1, rowSums(results_Fn_novel_compounds), `/`)

### 6.2.3. Visualization PhenoGMM models ----
source("/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/gmm_identifier.R")

## Relative abundance
# Visualization PhenoGMM model (excluding Cephalothin in construction of GMM model)
results_Fn_rel_novel_compounds <- results_Fn_rel_novel_compounds %>% 
  add_column(Sample_name = metadata_Fn_noneg$Sample_complete)
results_Fn_rel_novel_compounds_melted <- melt(results_Fn_rel_novel_compounds, id.vars = "Sample_name")

p_PhenoGMM_Fn_novel_compounds <- results_Fn_rel_novel_compounds_melted %>% 
  ggplot(data = ., aes(x = Sample_name, y = value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(y = "Relative abundance", x = "Sample", title = NULL, fill = "Cluster GMM model") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.text.align = 0)
print(p_PhenoGMM_Fn_novel_compounds)

## Flow cytometric space
PGMM_Fn_visual_pooled <- gmm_identifier(flow_set = fcs_PGMM_Fn_novel_compounds, PhenoGMM_model = PhenoGMM_Fn_noneg_fixedclust_novel_compounds)
PGMM_Fn_visual <- gmm_identifier(flow_set = flowData_transformed_Fn_norm_noneg, PhenoGMM_model = PhenoGMM_Fn_noneg_fixedclust_novel_compounds)

p_PGMM_Fn_cytometric_pooled <- ggplot(data.frame(PGMM_Fn_visual_pooled[[1]]))+
  geom_point(aes(x = BL1.H, y = BL3.H, color = as.factor(gmm_id)), size = 0.1)+
  xlim(0.4, 1)+
  ylim(0.2, 1)+
  labs(x = 'BL1-H', y = 'BL3-H', title = NULL, color = 'Cluster GMM model')+
  theme_bw()
print(p_PGMM_Fn_cytometric_pooled)

# Flow cytometric space for one sample
PGMM_Fn_visual_select <- PGMM_Fn_visual[37]

p_PGMM_Fn_cytometric <- ggplot(data.frame(PGMM_Fn_visual_select[[1]]))+
  geom_point(aes(x = BL1.H, y = BL3.H, color = as.factor(gmm_id)), size = 0.1)+
  xlim(0.4, 1)+
  ylim(0.2, 1)+
  labs(x = 'BL1-H', y = 'BL3-H', title = metadata_Fn_noneg$Sample_complete[37], color = 'Cluster GMM model')+
  theme_bw()+
  theme(plot.title = element_text(size = 26))
print(p_PGMM_Fn_cytometric)

# Flow cytometric space for All samples
plist_PGMM_Fn <- list()

for (i in 1:length(PGMM_Fn_visual)) {
  PGMM_Fn_visual_loop <- PGMM_Fn_visual[i]
  plist_PGMM_Fn[[i]] <- ggplot(data.frame(PGMM_Fn_visual_loop[[1]]))+
    geom_point(aes(x = BL1.H, y = BL3.H, color = as.factor(gmm_id)), size = 0.1)+
    xlim(0.4, 1)+
    ylim(0.2, 1)+
    labs(x = 'BL1-H', y = 'BL3-H', title = metadata_Fn_noneg$Sample_complete[i], color = 'Cluster GMM model')+
    theme_bw()+
    theme(legend.position = "none",
          plot.title = element_text(size = 12))
}

grobs_PGMM_Fn <- ggplotGrob(p_PGMM_Fn_cytometric_pooled)$grobs
legend_PGMM_grobs_Fn <- grobs_PGMM_Fn[[which(sapply(grobs_PGMM_Fn, function(x) x$name) == "guide-box")]]

nCol_PGMM_Fn <- floor(sqrt(length(plist_PGMM_Fn)))
p_PGMM_grid_Fn <- cowplot::plot_grid(plotlist = plist_PGMM_Fn, ncol = nCol_PGMM_Fn)
p_PGMM_sample_Fn <- cowplot::plot_grid(p_PGMM_grid_Fn, legend_PGMM_grobs_Fn, ncol = 2, rel_widths = c(9, 0.4))
print(p_PGMM_sample_Fn)

### 6.2.4. Ordination of GMM output ----
# Model excluding cephalothin
rownames_GMM_Fn_novel_compounds <- row.names(results_Fn_novel_compounds)
new_rownames_GMM_Fn_novel_compounds <- gsub("^.*?nt_", "", rownames_GMM_Fn_novel_compounds)
new_rownames_GMM_Fn_novel_compounds <- gsub("_100_SGPI", "", new_rownames_GMM_Fn_novel_compounds)
new_rownames_GMM_Fn_novel_compounds <- gsub(".fcs", "", new_rownames_GMM_Fn_novel_compounds)
results_Fn_ordination_novel_compounds <- results_Fn_novel_compounds
rownames(results_Fn_ordination_novel_compounds) <- new_rownames_GMM_Fn_novel_compounds
results_Fn_ordination_novel_compounds <- as.matrix(results_Fn_ordination_novel_compounds)

distance_GMM_Fn_novel_compounds <- vegan::vegdist(x = results_Fn_ordination_novel_compounds, method = "bray", binary = FALSE)
mds_GMM_Fn_novel_compounds <- stats::cmdscale(distance_GMM_Fn_novel_compounds, k = 2, eig = TRUE, add = TRUE)
NMDS_GMM_Fn_novel_compounds <- vegan::metaMDS(distance_GMM_Fn_novel_compounds, autotransform = FALSE, k = 2, trymax = 100)

plot_PCoA_GMM_Fn_novel_compounds <- plot_beta_fcm(mds_GMM_Fn_novel_compounds, color = as.factor(metadata_Fn_dist$MOA), shape = as.factor(metadata_Fn_dist$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL)+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(plot_PCoA_GMM_Fn_novel_compounds)

plot_NMDS_GMM_Fn_novel_compounds <- plot_beta_fcm(NMDS_GMM_Fn_novel_compounds, color = as.factor(metadata_Fn_dist$MOA), shape = as.factor(metadata_Fn_dist$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL, x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(plot_NMDS_GMM_Fn_novel_compounds)

# Model excluding cephalothin and ordination excluding cephalothin
results_Fn_ordination_novel_compounds_noceph <- results_Fn_ordination_novel_compounds[c(1, 5:72), ]
distance_GMM_Fn_novel_compounds_noceph <- vegan::vegdist(x = results_Fn_ordination_novel_compounds_noceph, method = "bray", binary = FALSE)
mds_GMM_Fn_novel_compounds_noceph <- stats::cmdscale(distance_GMM_Fn_novel_compounds_noceph, k = 2, eig = TRUE, add = TRUE)
NMDS_GMM_Fn_novel_compounds_noceph <- vegan::metaMDS(distance_GMM_Fn_novel_compounds_noceph, autotransform = FALSE, k = 2, trymax = 100)

metadata_Fn_dist_noceph <- metadata_Fn_dist[c(1, 5:72), ]

plot_PCoA_GMM_Fn_novel_compounds_noceph <- plot_beta_fcm(mds_GMM_Fn_novel_compounds_noceph, color = as.factor(metadata_Fn_dist_noceph$MOA), shape = as.factor(metadata_Fn_dist_noceph$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL)+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:4, 6:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(plot_PCoA_GMM_Fn_novel_compounds_noceph)

plot_NMDS_GMM_Fn_novel_compounds_noceph <- plot_beta_fcm(NMDS_GMM_Fn_novel_compounds_noceph, color = as.factor(metadata_Fn_dist_noceph$MOA), shape = as.factor(metadata_Fn_dist_noceph$Compound), labels = list("Mode of Action", "Compound")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.7)+
  labs(title = NULL, x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_shape_manual(values=c(1:4, 6:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(plot_NMDS_GMM_Fn_novel_compounds_noceph)

### 6.2.5. Statistics ----
### ANOSIM
ano_GMM_Fn_novel_compounds <- vegan::anosim(distance_GMM_Fn_novel_compounds, metadata_Fn_dist$MOA, distance = "bray", permutations = permutations)
ano_GMM_Fn_novel_compounds_export <- data.frame("ANOSIM Statistic R" = ano_GMM_Fn_novel_compounds$statistic, "P" = ano_GMM_Fn_novel_compounds$signif, "permutations" = ano_GMM_Fn_novel_compounds$permutations)
write.csv2(file = "ANOSIM_GMM_Fn_novel_compounds.csv", ano_GMM_Fn_novel_compounds_export) # Export results as csv file
plot_pairwise_anosim_GMM_Fn_novel_compounds <- plot(ano_GMM_Fn_novel_compounds, xaxt = "n", main = "Dissimilarity ranks GMM Fn (ANOSIM)", xlab = "Class", ylab = "Dissimilarity rank")
axis(1, at = 1:11, labels = c("Between", "CWS", "DNARepl", "DNATrans", "FAM", "Heat", "MD", "Control", "PS-30SInh", "PS-50SInh", "PS-tRNAInt"))

# Exclude cephalothin from calculation
ano_GMM_Fn_novel_compounds_no_ceph <- vegan::anosim(distance_GMM_Fn_novel_compounds_noceph, metadata_Fn_dist_noceph$MOA, distance = "bray", permutations = permutations)
ano_GMM_Fn_novel_compounds_no_ceph_export <- data.frame("ANOSIM Statistic R" = ano_GMM_Fn_novel_compounds_no_ceph$statistic, "P" = ano_GMM_Fn_novel_compounds_no_ceph$signif, "permutations" = ano_GMM_Fn_novel_compounds_no_ceph$permutations)
write.csv2(file = "ANOSIM_GMM_Fn_novel_compounds_no_ceph.csv", ano_GMM_Fn_novel_compounds_no_ceph_export) # Export results as csv file
plot_pairwise_anosim_GMM_Fn_novel_compounds_no_ceph <- plot(ano_GMM_Fn_novel_compounds_no_ceph, xaxt = "n", main = "Dissimilarity ranks GMM Fn (ANOSIM)", xlab = "Class", ylab = "Dissimilarity rank")
axis(1, at = 1:11, labels = c("Between", "CWS", "DNARepl", "DNATrans", "FAM", "Heat", "MD", "Control", "PS-30SInh", "PS-50SInh", "PS-tRNAInt"))

# Pairwise comparison different MOA
pairwise_anosim_GMM_Fn <- numeric()

# CellWallSynthesis vs Control
distance_GMM_Fn_CWS_PC_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 13:18), c(1:12, 13:18)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_PC_inclCeph <- metadata_Fn_dist[c(1:12, 13:18), ]
ano_GMM_Fn_CWS_PC_inclCeph <- vegan::anosim(distance_GMM_Fn_CWS_PC_inclCeph, metadata_GMM_Fn_dist_CWS_PC_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_Control_inclCephalothin"] <- ano_GMM_Fn_CWS_PC_inclCeph$signif

# CellWallSynthesis vs DNAReplication
distance_GMM_Fn_CWS_DR_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 19:27), c(1:12, 19:27)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_DR_inclCeph <- metadata_Fn_dist[c(1:12, 19:27), ]
ano_GMM_Fn_CWS_DR_inclCeph <- vegan::anosim(distance_GMM_Fn_CWS_DR_inclCeph, metadata_GMM_Fn_dist_CWS_DR_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_DNAReplication_inclCephalothin"] <- ano_GMM_Fn_CWS_DR_inclCeph$signif

# CellWallSynthesis vs DNATranscription
distance_GMM_Fn_CWS_DT_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 28:33), c(1:12, 28:33)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_DT_inclCeph <- metadata_Fn_dist[c(1:12, 28:33), ]
ano_GMM_Fn_CWS_DT_inclCeph <- vegan::anosim(distance_GMM_Fn_CWS_DT_inclCeph, metadata_GMM_Fn_dist_CWS_DT_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_DNATranscription_inclCephalothin"] <- ano_GMM_Fn_CWS_DT_inclCeph$signif

# CellWallSynthesis vs FolicAcidMetabolism
distance_GMM_Fn_CWS_FAM_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 34:36), c(1:12, 34:36)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_FAM_inclCeph <- metadata_Fn_dist[c(1:12, 34:36), ]
ano_GMM_Fn_CWS_FAM_inclCeph <- vegan::anosim(distance_GMM_Fn_CWS_FAM_inclCeph, metadata_GMM_Fn_dist_CWS_FAM_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_FolicAcidMetabolism_inclCephalothin"] <- ano_GMM_Fn_CWS_FAM_inclCeph$signif

# CellWallSynthesis vs Heat
distance_GMM_Fn_CWS_H_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 37:39), c(1:12, 37:39)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_H_inclCeph <- metadata_Fn_dist[c(1:12, 37:39), ]
ano_GMM_Fn_CWS_H_inclCeph <- vegan::anosim(distance_GMM_Fn_CWS_H_inclCeph, metadata_GMM_Fn_dist_CWS_H_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_Heat_inclCephalothin"] <- ano_GMM_Fn_CWS_H_inclCeph$signif

# CellWallSynthesis vs MembraneDisruption
distance_GMM_Fn_CWS_MD_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 40:51), c(1:12, 40:51)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_MD_inclCeph <- metadata_Fn_dist[c(1:12, 40:51), ]
ano_GMM_Fn_CWS_MD_inclCeph <- vegan::anosim(distance_GMM_Fn_CWS_MD_inclCeph, metadata_GMM_Fn_dist_CWS_MD_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_MembraneDisruption_inclCephalothin"] <- ano_GMM_Fn_CWS_MD_inclCeph$signif

# CellWallSynthesis vs ProteinSynthesis30SInhibition
distance_GMM_Fn_CWS_PS30SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 52:60), c(1:12, 52:60)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_PS30SI_inclCeph <- metadata_Fn_dist[c(1:12, 52:60), ]
ano_GMM_Fn_CWS_PS30SI_inclCeph <- vegan::anosim(distance_GMM_Fn_CWS_PS30SI_inclCeph, metadata_GMM_Fn_dist_CWS_PS30SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_ProteinSynthesis30SInhibition_inclCephalothin"] <- ano_GMM_Fn_CWS_PS30SI_inclCeph$signif

# CellWallSynthesis vs ProteinSynthesis50SInhibition
distance_GMM_Fn_CWS_PS50SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 61:69), c(1:12, 61:69)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_PS50SI_inclCeph <- metadata_Fn_dist[c(1:12, 61:69), ]
ano_GMM_Fn_CWS_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Fn_CWS_PS50SI_inclCeph, metadata_GMM_Fn_dist_CWS_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Fn_CWS_PS50SI_inclCeph$signif

# CellWallSynthesis vs ProteinSynthesistRNAInterference
distance_GMM_Fn_CWS_PStRI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(1:12, 70:72), c(1:12, 70:72)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_PStRI_inclCeph <- metadata_Fn_dist[c(1:12, 70:72), ]
ano_GMM_Fn_CWS_PStRI_inclCeph <- vegan::anosim(distance_GMM_Fn_CWS_PStRI_inclCeph, metadata_GMM_Fn_dist_CWS_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Fn_CWS_PStRI_inclCeph$signif

# Control vs DNAReplication
distance_GMM_Fn_PC_DR_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 19:27), c(13:18, 19:27)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_DR_inclCeph <- metadata_Fn_dist[c(13:18, 19:27), ]
ano_GMM_Fn_PC_DR_inclCeph <- vegan::anosim(distance_GMM_Fn_PC_DR_inclCeph, metadata_GMM_Fn_dist_PC_DR_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_DNAReplication_inclCephalothin"] <- ano_GMM_Fn_PC_DR_inclCeph$signif

# Control vs DNATranscription
distance_GMM_Fn_PC_DT_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 28:33), c(13:18, 28:33)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_DT_inclCeph <- metadata_Fn_dist[c(13:18, 28:33), ]
ano_GMM_Fn_PC_DT_inclCeph <- vegan::anosim(distance_GMM_Fn_PC_DT_inclCeph, metadata_GMM_Fn_dist_PC_DT_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_DNATranscription_inclCephalothin"] <- ano_GMM_Fn_PC_DT_inclCeph$signif

# Control vs FolicAcidMetabolism
distance_GMM_Fn_PC_FAM_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 34:36), c(13:18, 34:36)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_FAM_inclCeph <- metadata_Fn_dist[c(13:18, 34:36), ]
ano_GMM_Fn_PC_FAM_inclCeph <- vegan::anosim(distance_GMM_Fn_PC_FAM_inclCeph, metadata_GMM_Fn_dist_PC_FAM_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_FolicAcidMetabolism_inclCephalothin"] <- ano_GMM_Fn_PC_FAM_inclCeph$signif

# Control vs Heat
distance_GMM_Fn_PC_H_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 37:39), c(13:18, 37:39)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_H_inclCeph <- metadata_Fn_dist[c(13:18, 37:39), ]
ano_GMM_Fn_PC_H_inclCeph <- vegan::anosim(distance_GMM_Fn_PC_H_inclCeph, metadata_GMM_Fn_dist_PC_H_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_Heat_inclCephalothin"] <- ano_GMM_Fn_PC_H_inclCeph$signif

# Control vs MembraneDisruption
distance_GMM_Fn_PC_MD_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 40:51), c(13:18, 40:51)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_MD_inclCeph <- metadata_Fn_dist[c(13:18, 40:51), ]
ano_GMM_Fn_PC_MD_inclCeph <- vegan::anosim(distance_GMM_Fn_PC_MD_inclCeph, metadata_GMM_Fn_dist_PC_MD_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_MembraneDisruption_inclCephalothin"] <- ano_GMM_Fn_PC_MD_inclCeph$signif

# Control vs ProteinSynthesis30SInhibition
distance_GMM_Fn_PC_PS30SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 52:60), c(13:18, 52:60)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_PS30SI_inclCeph <- metadata_Fn_dist[c(13:18, 52:60), ]
ano_GMM_Fn_PC_PS30SI_inclCeph <- vegan::anosim(distance_GMM_Fn_PC_PS30SI_inclCeph, metadata_GMM_Fn_dist_PC_PS30SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_ProteinSynthesis30SInhibition_inclCephalothin"] <- ano_GMM_Fn_PC_PS30SI_inclCeph$signif

# Control vs ProteinSynthesis50SInhibition
distance_GMM_Fn_PC_PS50SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 61:69), c(13:18, 61:69)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_PS50SI_inclCeph <- metadata_Fn_dist[c(13:18, 61:69), ]
ano_GMM_Fn_PC_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Fn_PC_PS50SI_inclCeph, metadata_GMM_Fn_dist_PC_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Fn_PC_PS50SI_inclCeph$signif

# Control vs ProteinSynthesistRNAInterference
distance_GMM_Fn_PC_PStRI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(13:18, 70:72), c(13:18, 70:72)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_PStRI_inclCeph <- metadata_Fn_dist[c(13:18, 70:72), ]
ano_GMM_Fn_PC_PStRI_inclCeph <- vegan::anosim(distance_GMM_Fn_PC_PStRI_inclCeph, metadata_GMM_Fn_dist_PC_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Fn_PC_PStRI_inclCeph$signif

# DNAReplication vs DNATranscription
distance_GMM_Fn_DR_DT_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(19:27, 28:33), c(19:27, 28:33)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DR_DT_inclCeph <- metadata_Fn_dist[c(19:27, 28:33), ]
ano_GMM_Fn_DR_DT_inclCeph <- vegan::anosim(distance_GMM_Fn_DR_DT_inclCeph, metadata_GMM_Fn_dist_DR_DT_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNAReplication_DNATranscription_inclCephalothin"] <- ano_GMM_Fn_DR_DT_inclCeph$signif

# DNAReplication vs FolicAcidMetabolism
distance_GMM_Fn_DR_FAM_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(19:27, 34:36), c(19:27, 34:36)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DR_FAM_inclCeph <- metadata_Fn_dist[c(19:27, 34:36), ]
ano_GMM_Fn_DR_FAM_inclCeph <- vegan::anosim(distance_GMM_Fn_DR_FAM_inclCeph, metadata_GMM_Fn_dist_DR_FAM_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNAReplication_FolicAcidMetabolism_inclCephalothin"] <- ano_GMM_Fn_DR_FAM_inclCeph$signif

# DNAReplication vs Heat
distance_GMM_Fn_DR_H_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(19:27, 37:39), c(19:27, 37:39)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DR_H_inclCeph <- metadata_Fn_dist[c(19:27, 37:39), ]
ano_GMM_Fn_DR_H_inclCeph <- vegan::anosim(distance_GMM_Fn_DR_H_inclCeph, metadata_GMM_Fn_dist_DR_H_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNAReplication_Heat_inclCephalothin"] <- ano_GMM_Fn_DR_H_inclCeph$signif

# DNAReplication vs MembraneDisruption
distance_GMM_Fn_DR_MD_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(19:27, 40:51), c(19:27, 40:51)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DR_MD_inclCeph <- metadata_Fn_dist[c(19:27, 40:51), ]
ano_GMM_Fn_DR_MD_inclCeph <- vegan::anosim(distance_GMM_Fn_DR_MD_inclCeph, metadata_GMM_Fn_dist_DR_MD_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNAReplication_MembraneDisruption_inclCephalothin"] <- ano_GMM_Fn_DR_MD_inclCeph$signif

# DNAReplication vs ProteinSynthesis30SInhibition
distance_GMM_Fn_DR_PS30SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(19:27, 52:60), c(19:27, 52:60)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DR_PS30SI_inclCeph <- metadata_Fn_dist[c(19:27, 52:60), ]
ano_GMM_Fn_DR_PS30SI_inclCeph <- vegan::anosim(distance_GMM_Fn_DR_PS30SI_inclCeph, metadata_GMM_Fn_dist_DR_PS30SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNAReplication_ProteinSynthesis30SInhibition_inclCephalothin"] <- ano_GMM_Fn_DR_PS30SI_inclCeph$signif

# DNAReplication vs ProteinSynthesis50SInhibition
distance_GMM_Fn_DR_PS50SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(19:27, 61:69), c(19:27, 61:69)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DR_PS50SI_inclCeph <- metadata_Fn_dist[c(19:27, 61:69), ]
ano_GMM_Fn_DR_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Fn_DR_PS50SI_inclCeph, metadata_GMM_Fn_dist_DR_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNAReplication_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Fn_DR_PS50SI_inclCeph$signif

# DNAReplication vs ProteinSynthesistRNAInterference
distance_GMM_Fn_DR_PStRI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(19:27, 70:72), c(19:27, 70:72)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DR_PStRI_inclCeph <- metadata_Fn_dist[c(19:27, 70:72), ]
ano_GMM_Fn_DR_PStRI_inclCeph <- vegan::anosim(distance_GMM_Fn_DR_PStRI_inclCeph, metadata_GMM_Fn_dist_DR_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNAReplication_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Fn_DR_PStRI_inclCeph$signif

# DNATranscription vs FolicAcidMetabolism
distance_GMM_Fn_DT_FAM_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(28:33, 34:36), c(28:33, 34:36)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DT_FAM_inclCeph <- metadata_Fn_dist[c(28:33, 34:36), ]
ano_GMM_Fn_DT_FAM_inclCeph <- vegan::anosim(distance_GMM_Fn_DT_FAM_inclCeph, metadata_GMM_Fn_dist_DT_FAM_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNATranscription_FolicAcidMetabolism_inclCephalothin"] <- ano_GMM_Fn_DT_FAM_inclCeph$signif

# DNATranscription vs Heat
distance_GMM_Fn_DT_H_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(28:33, 37:39), c(28:33, 37:39)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DT_H_inclCeph <- metadata_Fn_dist[c(28:33, 37:39), ]
ano_GMM_Fn_DT_H_inclCeph <- vegan::anosim(distance_GMM_Fn_DT_H_inclCeph, metadata_GMM_Fn_dist_DT_H_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNATranscription_Heat_inclCephalothin"] <- ano_GMM_Fn_DT_H_inclCeph$signif

# DNATranscription vs MembraneDisruption
distance_GMM_Fn_DT_MD_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(28:33, 40:51), c(28:33, 40:51)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DT_MD_inclCeph <- metadata_Fn_dist[c(28:33, 40:51), ]
ano_GMM_Fn_DT_MD_inclCeph <- vegan::anosim(distance_GMM_Fn_DT_MD_inclCeph, metadata_GMM_Fn_dist_DT_MD_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNATranscription_MembraneDisruption_inclCephalothin"] <- ano_GMM_Fn_DT_MD_inclCeph$signif

# DNATranscription vs ProteinSynthesis30SInhibition
distance_GMM_Fn_DT_PS30SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(28:33, 52:60), c(28:33, 52:60)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DT_PS30SI_inclCeph <- metadata_Fn_dist[c(28:33, 52:60), ]
ano_GMM_Fn_DT_PS30SI_inclCeph <- vegan::anosim(distance_GMM_Fn_DT_PS30SI_inclCeph, metadata_GMM_Fn_dist_DT_PS30SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNATranscription_ProteinSynthesis30SInhibition_inclCephalothin"] <- ano_GMM_Fn_DT_PS30SI_inclCeph$signif

# DNATranscription vs ProteinSynthesis50SInhibition
distance_GMM_Fn_DT_PS50SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(28:33, 61:69), c(28:33, 61:69)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DT_PS50SI_inclCeph <- metadata_Fn_dist[c(28:33, 61:69), ]
ano_GMM_Fn_DT_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Fn_DT_PS50SI_inclCeph, metadata_GMM_Fn_dist_DT_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNATranscription_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Fn_DT_PS50SI_inclCeph$signif

# DNATranscription vs ProteinSynthesistRNAInterference
distance_GMM_Fn_DT_PStRI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(28:33, 70:72), c(28:33, 70:72)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DT_PStRI_inclCeph <- metadata_Fn_dist[c(28:33, 70:72), ]
ano_GMM_Fn_DT_PStRI_inclCeph <- vegan::anosim(distance_GMM_Fn_DT_PStRI_inclCeph, metadata_GMM_Fn_dist_DT_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNATranscription_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Fn_DT_PStRI_inclCeph$signif

# FolicAcidMetabolism vs Heat
distance_GMM_Fn_FAM_H_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(34:36, 37:39), c(34:36, 37:39)] %>% 
  as.dist()

metadata_GMM_Fn_dist_FAM_H_inclCeph <- metadata_Fn_dist[c(34:36, 37:39), ]
ano_GMM_Fn_FAM_H_inclCeph <- vegan::anosim(distance_GMM_Fn_FAM_H_inclCeph, metadata_GMM_Fn_dist_FAM_H_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["FolicAcidMetabolism_Heat_inclCephalothin"] <- ano_GMM_Fn_FAM_H_inclCeph$signif

# FolicAcidMetabolism vs MembraneDisruption
distance_GMM_Fn_FAM_MD_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(34:36, 40:51), c(34:36, 40:51)] %>% 
  as.dist()

metadata_GMM_Fn_dist_FAM_MD_inclCeph <- metadata_Fn_dist[c(34:36, 40:51), ]
ano_GMM_Fn_FAM_MD_inclCeph <- vegan::anosim(distance_GMM_Fn_FAM_MD_inclCeph, metadata_GMM_Fn_dist_FAM_MD_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["FolicAcidMetabolism_MembraneDisruption_inclCephalothin"] <- ano_GMM_Fn_FAM_MD_inclCeph$signif

# FolicAcidMetabolism vs ProteinSynthesis30SInhibition
distance_GMM_Fn_FAM_PS30SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(34:36, 52:60), c(34:36, 52:60)] %>% 
  as.dist()

metadata_GMM_Fn_dist_FAM_PS30SI_inclCeph <- metadata_Fn_dist[c(34:36, 52:60), ]
ano_GMM_Fn_FAM_PS30SI_inclCeph <- vegan::anosim(distance_GMM_Fn_FAM_PS30SI_inclCeph, metadata_GMM_Fn_dist_FAM_PS30SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["FolicAcidMetabolism_ProteinSynthesis30SInhibition_inclCephalothin"] <- ano_GMM_Fn_FAM_PS30SI_inclCeph$signif

# FolicAcidMetabolism vs ProteinSynthesis50SInhibition
distance_GMM_Fn_FAM_PS50SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(34:36, 61:69), c(34:36, 61:69)] %>% 
  as.dist()

metadata_GMM_Fn_dist_FAM_PS50SI_inclCeph <- metadata_Fn_dist[c(34:36, 61:69), ]
ano_GMM_Fn_FAM_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Fn_FAM_PS50SI_inclCeph, metadata_GMM_Fn_dist_FAM_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["FolicAcidMetabolism_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Fn_FAM_PS50SI_inclCeph$signif

# FolicAcidMetabolism vs ProteinSynthesistRNAInterference
distance_GMM_Fn_FAM_PStRI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(34:36, 70:72), c(34:36, 70:72)] %>% 
  as.dist()

metadata_GMM_Fn_dist_FAM_PStRI_inclCeph <- metadata_Fn_dist[c(34:36, 70:72), ]
ano_GMM_Fn_FAM_PStRI_inclCeph <- vegan::anosim(distance_GMM_Fn_FAM_PStRI_inclCeph, metadata_GMM_Fn_dist_FAM_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["FolicAcidMetabolism_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Fn_FAM_PStRI_inclCeph$signif

# Heat vs MembraneDisruption
distance_GMM_Fn_H_MD_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(37:39, 40:51), c(37:39, 40:51)] %>% 
  as.dist()

metadata_GMM_Fn_dist_H_MD_inclCeph <- metadata_Fn_dist[c(37:39, 40:51), ]
ano_GMM_Fn_H_MD_inclCeph <- vegan::anosim(distance_GMM_Fn_H_MD_inclCeph, metadata_GMM_Fn_dist_H_MD_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Heat_MembraneDisruption_inclCephalothin"] <- ano_GMM_Fn_H_MD_inclCeph$signif

# Heat vs ProteinSynthesis30SInhibition
distance_GMM_Fn_H_PS30SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(37:39, 52:60), c(37:39, 52:60)] %>% 
  as.dist()

metadata_GMM_Fn_dist_H_PS30SI_inclCeph <- metadata_Fn_dist[c(37:39, 52:60), ]
ano_GMM_Fn_H_PS30SI_inclCeph <- vegan::anosim(distance_GMM_Fn_H_PS30SI_inclCeph, metadata_GMM_Fn_dist_H_PS30SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Heat_ProteinSynthesis30SInhibition_inclCephalothin"] <- ano_GMM_Fn_H_PS30SI_inclCeph$signif

# Heat vs ProteinSynthesis50SInhibition
distance_GMM_Fn_H_PS50SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(37:39, 61:69), c(37:39, 61:69)] %>% 
  as.dist()

metadata_GMM_Fn_dist_H_PS50SI_inclCeph <- metadata_Fn_dist[c(37:39, 61:69), ]
ano_GMM_Fn_H_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Fn_H_PS50SI_inclCeph, metadata_GMM_Fn_dist_H_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Heat_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Fn_H_PS50SI_inclCeph$signif

# Heat vs ProteinSynthesistRNAInterference
distance_GMM_Fn_H_PStRI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(37:39, 70:72), c(37:39, 70:72)] %>% 
  as.dist()

metadata_GMM_Fn_dist_H_PStRI_inclCeph <- metadata_Fn_dist[c(37:39, 70:72), ]
ano_GMM_Fn_H_PStRI_inclCeph <- vegan::anosim(distance_GMM_Fn_H_PStRI_inclCeph, metadata_GMM_Fn_dist_H_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Heat_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Fn_H_PStRI_inclCeph$signif

# MembraneDisruption vs ProteinSynthesis30SInhibition
distance_GMM_Fn_MD_PS30SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(40:51, 52:60), c(40:51, 52:60)] %>% 
  as.dist()

metadata_GMM_Fn_dist_MD_PS30SI_inclCeph <- metadata_Fn_dist[c(40:51, 52:60), ]
ano_GMM_Fn_MD_PS30SI_inclCeph <- vegan::anosim(distance_GMM_Fn_MD_PS30SI_inclCeph, metadata_GMM_Fn_dist_MD_PS30SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["MembraneDisruption_ProteinSynthesis30SInhibition_inclCephalothin"] <- ano_GMM_Fn_MD_PS30SI_inclCeph$signif

# MembraneDisruption vs ProteinSynthesis50SInhibition
distance_GMM_Fn_MD_PS50SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(40:51, 61:69), c(40:51, 61:69)] %>% 
  as.dist()

metadata_GMM_Fn_dist_MD_PS50SI_inclCeph <- metadata_Fn_dist[c(40:51, 61:69), ]
ano_GMM_Fn_MD_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Fn_MD_PS50SI_inclCeph, metadata_GMM_Fn_dist_MD_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["MembraneDisruption_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Fn_MD_PS50SI_inclCeph$signif

# MembraneDisruption vs ProteinSynthesistRNAInterference
distance_GMM_Fn_MD_PStRI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(40:51, 70:72), c(40:51, 70:72)] %>% 
  as.dist()

metadata_GMM_Fn_dist_MD_PStRI_inclCeph <- metadata_Fn_dist[c(40:51, 70:72), ]
ano_GMM_Fn_MD_PStRI_inclCeph <- vegan::anosim(distance_GMM_Fn_MD_PStRI_inclCeph, metadata_GMM_Fn_dist_MD_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["MembraneDisruption_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Fn_MD_PStRI_inclCeph$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesis50SInhibition
distance_GMM_Fn_PS30SI_PS50SI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(52:60, 61:69), c(52:60, 61:69)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PS30SI_PS50SI_inclCeph <- metadata_Fn_dist[c(52:60, 61:69), ]
ano_GMM_Fn_PS30SI_PS50SI_inclCeph <- vegan::anosim(distance_GMM_Fn_PS30SI_PS50SI_inclCeph, metadata_GMM_Fn_dist_PS30SI_PS50SI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["ProteinSynthesis30SInhibition_ProteinSynthesis50SInhibition_inclCephalothin"] <- ano_GMM_Fn_PS30SI_PS50SI_inclCeph$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesistRNAInterference
distance_GMM_Fn_PS30SI_PStRI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(52:60, 70:72), c(52:60, 70:72)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PS30SI_PStRI_inclCeph <- metadata_Fn_dist[c(52:60, 70:72), ]
ano_GMM_Fn_PS30SI_PStRI_inclCeph <- vegan::anosim(distance_GMM_Fn_PS30SI_PStRI_inclCeph, metadata_GMM_Fn_dist_PS30SI_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["ProteinSynthesis30SInhibition_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Fn_PS30SI_PStRI_inclCeph$signif

# ProteinSynthesis50SInhibition vs ProteinSynthesistRNAInterference
distance_GMM_Fn_PS50SI_PStRI_inclCeph <- distance_GMM_Fn_novel_compounds %>% 
  as.matrix() %>% 
  .[c(61:69, 70:72), c(61:69, 70:72)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PS50SI_PStRI_inclCeph <- metadata_Fn_dist[c(61:69, 70:72), ]
ano_GMM_Fn_PS50SI_PStRI_inclCeph <- vegan::anosim(distance_GMM_Fn_PS50SI_PStRI_inclCeph, metadata_GMM_Fn_dist_PS50SI_PStRI_inclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["ProteinSynthesis50SInhibition_ProteinSynthesistRNAInterference_inclCephalothin"] <- ano_GMM_Fn_PS50SI_PStRI_inclCeph$signif

## Pairwise comparison excluding cephalothin

# CellWallSynthesis vs Control
distance_GMM_Fn_CWS_PC_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PC), c(ind_CWS, ind_PC)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_PC_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_PC), ]
ano_GMM_Fn_CWS_PC_exclCeph <- vegan::anosim(distance_GMM_Fn_CWS_PC_exclCeph, metadata_GMM_Fn_dist_CWS_PC_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_Control_exclCephalothin"] <- ano_GMM_Fn_CWS_PC_exclCeph$signif

# CellWallSynthesis vs DNAReplication
distance_GMM_Fn_CWS_DR_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_DR), c(ind_CWS, ind_DR)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_DR_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_DR), ]
ano_GMM_Fn_CWS_DR_exclCeph <- vegan::anosim(distance_GMM_Fn_CWS_DR_exclCeph, metadata_GMM_Fn_dist_CWS_DR_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_DNAReplication_exclCephalothin"] <- ano_GMM_Fn_CWS_DR_exclCeph$signif

# CellWallSynthesis vs DNATranscription
distance_GMM_Fn_CWS_DT_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_DT), c(ind_CWS, ind_DT)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_DT_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_DT), ]
ano_GMM_Fn_CWS_DT_exclCeph <- vegan::anosim(distance_GMM_Fn_CWS_DT_exclCeph, metadata_GMM_Fn_dist_CWS_DT_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_DNATranscription_exclCephalothin"] <- ano_GMM_Fn_CWS_DT_exclCeph$signif

# CellWallSynthesis vs FolicAcidMetabolism
distance_GMM_Fn_CWS_FAM_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_FAM), c(ind_CWS, ind_FAM)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_FAM_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_FAM), ]
ano_GMM_Fn_CWS_FAM_exclCeph <- vegan::anosim(distance_GMM_Fn_CWS_FAM_exclCeph, metadata_GMM_Fn_dist_CWS_FAM_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_FolicAcidMetabolism_exclCephalothin"] <- ano_GMM_Fn_CWS_FAM_exclCeph$signif

# CellWallSynthesis vs Heat
distance_GMM_Fn_CWS_H_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_H), c(ind_CWS, ind_H)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_H_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_H), ]
ano_GMM_Fn_CWS_H_exclCeph <- vegan::anosim(distance_GMM_Fn_CWS_H_exclCeph, metadata_GMM_Fn_dist_CWS_H_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_Heat_exclCephalothin"] <- ano_GMM_Fn_CWS_H_exclCeph$signif

# CellWallSynthesis vs MembraneDisruption
distance_GMM_Fn_CWS_MD_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_MD), c(ind_CWS, ind_MD)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_MD_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_MD), ]
ano_GMM_Fn_CWS_MD_exclCeph <- vegan::anosim(distance_GMM_Fn_CWS_MD_exclCeph, metadata_GMM_Fn_dist_CWS_MD_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_MembraneDisruption_exclCephalothin"] <- ano_GMM_Fn_CWS_MD_exclCeph$signif

# CellWallSynthesis vs ProteinSynthesis30SInhibition
distance_GMM_Fn_CWS_PS30SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PS30SI), c(ind_CWS, ind_PS30SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_PS30SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_PS30SI), ]
ano_GMM_Fn_CWS_PS30SI_exclCeph <- vegan::anosim(distance_GMM_Fn_CWS_PS30SI_exclCeph, metadata_GMM_Fn_dist_CWS_PS30SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_ProteinSynthesis30SInhibition_exclCephalothin"] <- ano_GMM_Fn_CWS_PS30SI_exclCeph$signif

# CellWallSynthesis vs ProteinSynthesis50SInhibition
distance_GMM_Fn_CWS_PS50SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PS50SI), c(ind_CWS, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_PS50SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_PS50SI), ]
ano_GMM_Fn_CWS_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Fn_CWS_PS50SI_exclCeph, metadata_GMM_Fn_dist_CWS_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Fn_CWS_PS50SI_exclCeph$signif

# CellWallSynthesis vs ProteinSynthesistRNAInterference
distance_GMM_Fn_CWS_PStRI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_CWS, ind_PStRI), c(ind_CWS, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_CWS_PStRI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_CWS, ind_PStRI), ]
ano_GMM_Fn_CWS_PStRI_exclCeph <- vegan::anosim(distance_GMM_Fn_CWS_PStRI_exclCeph, metadata_GMM_Fn_dist_CWS_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["CellWallSynthesis_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Fn_CWS_PStRI_exclCeph$signif

# Control vs DNAReplication
distance_GMM_Fn_PC_DR_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_DR), c(ind_PC, ind_DR)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_DR_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_DR), ]
ano_GMM_Fn_PC_DR_exclCeph <- vegan::anosim(distance_GMM_Fn_PC_DR_exclCeph, metadata_GMM_Fn_dist_PC_DR_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_DNAReplication_exclCephalothin"] <- ano_GMM_Fn_PC_DR_exclCeph$signif

# Control vs DNATranscription
distance_GMM_Fn_PC_DT_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_DT), c(ind_PC, ind_DT)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_DT_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_DT), ]
ano_GMM_Fn_PC_DT_exclCeph <- vegan::anosim(distance_GMM_Fn_PC_DT_exclCeph, metadata_GMM_Fn_dist_PC_DT_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_DNATranscription_exclCephalothin"] <- ano_GMM_Fn_PC_DT_exclCeph$signif

# Control vs FolicAcidMetabolism
distance_GMM_Fn_PC_FAM_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_FAM), c(ind_PC, ind_FAM)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_FAM_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_FAM), ]
ano_GMM_Fn_PC_FAM_exclCeph <- vegan::anosim(distance_GMM_Fn_PC_FAM_exclCeph, metadata_GMM_Fn_dist_PC_FAM_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_FolicAcidMetabolism_exclCephalothin"] <- ano_GMM_Fn_PC_FAM_exclCeph$signif

# Control vs Heat
distance_GMM_Fn_PC_H_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_H), c(ind_PC, ind_H)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_H_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_H), ]
ano_GMM_Fn_PC_H_exclCeph <- vegan::anosim(distance_GMM_Fn_PC_H_exclCeph, metadata_GMM_Fn_dist_PC_H_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_Heat_exclCephalothin"] <- ano_GMM_Fn_PC_H_exclCeph$signif

# Control vs MembraneDisruption
distance_GMM_Fn_PC_MD_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_MD), c(ind_PC, ind_MD)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_MD_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_MD), ]
ano_GMM_Fn_PC_MD_exclCeph <- vegan::anosim(distance_GMM_Fn_PC_MD_exclCeph, metadata_GMM_Fn_dist_PC_MD_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_MembraneDisruption_exclCephalothin"] <- ano_GMM_Fn_PC_MD_exclCeph$signif

# Control vs ProteinSynthesis30SInhibition
distance_GMM_Fn_PC_PS30SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_PS30SI), c(ind_PC, ind_PS30SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_PS30SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_PS30SI), ]
ano_GMM_Fn_PC_PS30SI_exclCeph <- vegan::anosim(distance_GMM_Fn_PC_PS30SI_exclCeph, metadata_GMM_Fn_dist_PC_PS30SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_ProteinSynthesis30SInhibition_exclCephalothin"] <- ano_GMM_Fn_PC_PS30SI_exclCeph$signif

# Control vs ProteinSynthesis50SInhibition
distance_GMM_Fn_PC_PS50SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_PS50SI), c(ind_PC, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_PS50SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_PS50SI), ]
ano_GMM_Fn_PC_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Fn_PC_PS50SI_exclCeph, metadata_GMM_Fn_dist_PC_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Fn_PC_PS50SI_exclCeph$signif

# Control vs ProteinSynthesistRNAInterference
distance_GMM_Fn_PC_PStRI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PC, ind_PStRI), c(ind_PC, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PC_PStRI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_PC, ind_PStRI), ]
ano_GMM_Fn_PC_PStRI_exclCeph <- vegan::anosim(distance_GMM_Fn_PC_PStRI_exclCeph, metadata_GMM_Fn_dist_PC_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Control_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Fn_PC_PStRI_exclCeph$signif

# DNAReplication vs DNATranscription
distance_GMM_Fn_DR_DT_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_DT), c(ind_DR, ind_DT)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DR_DT_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_DR, ind_DT), ]
ano_GMM_Fn_DR_DT_exclCeph <- vegan::anosim(distance_GMM_Fn_DR_DT_exclCeph, metadata_GMM_Fn_dist_DR_DT_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNAReplication_DNATranscription_exclCephalothin"] <- ano_GMM_Fn_DR_DT_exclCeph$signif

# DNAReplication vs FolicAcidMetabolism
distance_GMM_Fn_DR_FAM_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_FAM), c(ind_DR, ind_FAM)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DR_FAM_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_DR, ind_FAM), ]
ano_GMM_Fn_DR_FAM_exclCeph <- vegan::anosim(distance_GMM_Fn_DR_FAM_exclCeph, metadata_GMM_Fn_dist_DR_FAM_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNAReplication_FolicAcidMetabolism_exclCephalothin"] <- ano_GMM_Fn_DR_FAM_exclCeph$signif

# DNAReplication vs Heat
distance_GMM_Fn_DR_H_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_H), c(ind_DR, ind_H)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DR_H_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_DR, ind_H), ]
ano_GMM_Fn_DR_H_exclCeph <- vegan::anosim(distance_GMM_Fn_DR_H_exclCeph, metadata_GMM_Fn_dist_DR_H_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNAReplication_Heat_exclCephalothin"] <- ano_GMM_Fn_DR_H_exclCeph$signif

# DNAReplication vs MembraneDisruption
distance_GMM_Fn_DR_MD_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_MD), c(ind_DR, ind_MD)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DR_MD_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_DR, ind_MD), ]
ano_GMM_Fn_DR_MD_exclCeph <- vegan::anosim(distance_GMM_Fn_DR_MD_exclCeph, metadata_GMM_Fn_dist_DR_MD_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNAReplication_MembraneDisruption_exclCephalothin"] <- ano_GMM_Fn_DR_MD_exclCeph$signif

# DNAReplication vs ProteinSynthesis30SInhibition
distance_GMM_Fn_DR_PS30SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_PS30SI), c(ind_DR, ind_PS30SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DR_PS30SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_DR, ind_PS30SI), ]
ano_GMM_Fn_DR_PS30SI_exclCeph <- vegan::anosim(distance_GMM_Fn_DR_PS30SI_exclCeph, metadata_GMM_Fn_dist_DR_PS30SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNAReplication_ProteinSynthesis30SInhibition_exclCephalothin"] <- ano_GMM_Fn_DR_PS30SI_exclCeph$signif

# DNAReplication vs ProteinSynthesis50SInhibition
distance_GMM_Fn_DR_PS50SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_PS50SI), c(ind_DR, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DR_PS50SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_DR, ind_PS50SI), ]
ano_GMM_Fn_DR_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Fn_DR_PS50SI_exclCeph, metadata_GMM_Fn_dist_DR_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNAReplication_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Fn_DR_PS50SI_exclCeph$signif

# DNAReplication vs ProteinSynthesistRNAInterference
distance_GMM_Fn_DR_PStRI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DR, ind_PStRI), c(ind_DR, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DR_PStRI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_DR, ind_PStRI), ]
ano_GMM_Fn_DR_PStRI_exclCeph <- vegan::anosim(distance_GMM_Fn_DR_PStRI_exclCeph, metadata_GMM_Fn_dist_DR_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNAReplication_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Fn_DR_PStRI_exclCeph$signif

# DNATranscription vs FolicAcidMetabolism
distance_GMM_Fn_DT_FAM_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_FAM), c(ind_DT, ind_FAM)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DT_FAM_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_DT, ind_FAM), ]
ano_GMM_Fn_DT_FAM_exclCeph <- vegan::anosim(distance_GMM_Fn_DT_FAM_exclCeph, metadata_GMM_Fn_dist_DT_FAM_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNATranscription_FolicAcidMetabolism_exclCephalothin"] <- ano_GMM_Fn_DT_FAM_exclCeph$signif

# DNATranscription vs Heat
distance_GMM_Fn_DT_H_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_H), c(ind_DT, ind_H)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DT_H_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_DT, ind_H), ]
ano_GMM_Fn_DT_H_exclCeph <- vegan::anosim(distance_GMM_Fn_DT_H_exclCeph, metadata_GMM_Fn_dist_DT_H_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNATranscription_Heat_exclCephalothin"] <- ano_GMM_Fn_DT_H_exclCeph$signif

# DNATranscription vs MembraneDisruption
distance_GMM_Fn_DT_MD_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_MD), c(ind_DT, ind_MD)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DT_MD_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_DT, ind_MD), ]
ano_GMM_Fn_DT_MD_exclCeph <- vegan::anosim(distance_GMM_Fn_DT_MD_exclCeph, metadata_GMM_Fn_dist_DT_MD_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNATranscription_MembraneDisruption_exclCephalothin"] <- ano_GMM_Fn_DT_MD_exclCeph$signif

# DNATranscription vs ProteinSynthesis30SInhibition
distance_GMM_Fn_DT_PS30SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_PS30SI), c(ind_DT, ind_PS30SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DT_PS30SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_DT, ind_PS30SI), ]
ano_GMM_Fn_DT_PS30SI_exclCeph <- vegan::anosim(distance_GMM_Fn_DT_PS30SI_exclCeph, metadata_GMM_Fn_dist_DT_PS30SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNATranscription_ProteinSynthesis30SInhibition_exclCephalothin"] <- ano_GMM_Fn_DT_PS30SI_exclCeph$signif

# DNATranscription vs ProteinSynthesis50SInhibition
distance_GMM_Fn_DT_PS50SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_PS50SI), c(ind_DT, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DT_PS50SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_DT, ind_PS50SI), ]
ano_GMM_Fn_DT_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Fn_DT_PS50SI_exclCeph, metadata_GMM_Fn_dist_DT_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNATranscription_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Fn_DT_PS50SI_exclCeph$signif

# DNATranscription vs ProteinSynthesistRNAInterference
distance_GMM_Fn_DT_PStRI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_DT, ind_PStRI), c(ind_DT, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_DT_PStRI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_DT, ind_PStRI), ]
ano_GMM_Fn_DT_PStRI_exclCeph <- vegan::anosim(distance_GMM_Fn_DT_PStRI_exclCeph, metadata_GMM_Fn_dist_DT_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["DNATranscription_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Fn_DT_PStRI_exclCeph$signif

# FolicAcidMetabolism vs Heat
distance_GMM_Fn_FAM_H_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_H), c(ind_FAM, ind_H)] %>% 
  as.dist()

metadata_GMM_Fn_dist_FAM_H_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_FAM, ind_H), ]
ano_GMM_Fn_FAM_H_exclCeph <- vegan::anosim(distance_GMM_Fn_FAM_H_exclCeph, metadata_GMM_Fn_dist_FAM_H_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["FolicAcidMetabolism_Heat_exclCephalothin"] <- ano_GMM_Fn_FAM_H_exclCeph$signif

# FolicAcidMetabolism vs MembraneDisruption
distance_GMM_Fn_FAM_MD_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_MD), c(ind_FAM, ind_MD)] %>% 
  as.dist()

metadata_GMM_Fn_dist_FAM_MD_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_FAM, ind_MD), ]
ano_GMM_Fn_FAM_MD_exclCeph <- vegan::anosim(distance_GMM_Fn_FAM_MD_exclCeph, metadata_GMM_Fn_dist_FAM_MD_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["FolicAcidMetabolism_MembraneDisruption_exclCephalothin"] <- ano_GMM_Fn_FAM_MD_exclCeph$signif

# FolicAcidMetabolism vs ProteinSynthesis30SInhibition
distance_GMM_Fn_FAM_PS30SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_PS30SI), c(ind_FAM, ind_PS30SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_FAM_PS30SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_FAM, ind_PS30SI), ]
ano_GMM_Fn_FAM_PS30SI_exclCeph <- vegan::anosim(distance_GMM_Fn_FAM_PS30SI_exclCeph, metadata_GMM_Fn_dist_FAM_PS30SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["FolicAcidMetabolism_ProteinSynthesis30SInhibition_exclCephalothin"] <- ano_GMM_Fn_FAM_PS30SI_exclCeph$signif

# FolicAcidMetabolism vs ProteinSynthesis50SInhibition
distance_GMM_Fn_FAM_PS50SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_PS50SI), c(ind_FAM, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_FAM_PS50SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_FAM, ind_PS50SI), ]
ano_GMM_Fn_FAM_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Fn_FAM_PS50SI_exclCeph, metadata_GMM_Fn_dist_FAM_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["FolicAcidMetabolism_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Fn_FAM_PS50SI_exclCeph$signif

# FolicAcidMetabolism vs ProteinSynthesistRNAInterference
distance_GMM_Fn_FAM_PStRI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_FAM, ind_PStRI), c(ind_FAM, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_FAM_PStRI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_FAM, ind_PStRI), ]
ano_GMM_Fn_FAM_PStRI_exclCeph <- vegan::anosim(distance_GMM_Fn_FAM_PStRI_exclCeph, metadata_GMM_Fn_dist_FAM_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["FolicAcidMetabolism_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Fn_FAM_PStRI_exclCeph$signif

# Heat vs MembraneDisruption
distance_GMM_Fn_H_MD_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_MD), c(ind_H, ind_MD)] %>% 
  as.dist()

metadata_GMM_Fn_dist_H_MD_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_H, ind_MD), ]
ano_GMM_Fn_H_MD_exclCeph <- vegan::anosim(distance_GMM_Fn_H_MD_exclCeph, metadata_GMM_Fn_dist_H_MD_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Heat_MembraneDisruption_exclCephalothin"] <- ano_GMM_Fn_H_MD_exclCeph$signif

# Heat vs ProteinSynthesis30SInhibition
distance_GMM_Fn_H_PS30SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_PS30SI), c(ind_H, ind_PS30SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_H_PS30SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_H, ind_PS30SI), ]
ano_GMM_Fn_H_PS30SI_exclCeph <- vegan::anosim(distance_GMM_Fn_H_PS30SI_exclCeph, metadata_GMM_Fn_dist_H_PS30SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Heat_ProteinSynthesis30SInhibition_exclCephalothin"] <- ano_GMM_Fn_H_PS30SI_exclCeph$signif

# Heat vs ProteinSynthesis50SInhibition
distance_GMM_Fn_H_PS50SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_PS50SI), c(ind_H, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_H_PS50SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_H, ind_PS50SI), ]
ano_GMM_Fn_H_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Fn_H_PS50SI_exclCeph, metadata_GMM_Fn_dist_H_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Heat_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Fn_H_PS50SI_exclCeph$signif

# Heat vs ProteinSynthesistRNAInterference
distance_GMM_Fn_H_PStRI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_H, ind_PStRI), c(ind_H, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_H_PStRI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_H, ind_PStRI), ]
ano_GMM_Fn_H_PStRI_exclCeph <- vegan::anosim(distance_GMM_Fn_H_PStRI_exclCeph, metadata_GMM_Fn_dist_H_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["Heat_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Fn_H_PStRI_exclCeph$signif

# MembraneDisruption vs ProteinSynthesis30SInhibition
distance_GMM_Fn_MD_PS30SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_MD, ind_PS30SI), c(ind_MD, ind_PS30SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_MD_PS30SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_MD, ind_PS30SI), ]
ano_GMM_Fn_MD_PS30SI_exclCeph <- vegan::anosim(distance_GMM_Fn_MD_PS30SI_exclCeph, metadata_GMM_Fn_dist_MD_PS30SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["MembraneDisruption_ProteinSynthesis30SInhibition_exclCephalothin"] <- ano_GMM_Fn_MD_PS30SI_exclCeph$signif

# MembraneDisruption vs ProteinSynthesis50SInhibition
distance_GMM_Fn_MD_PS50SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_MD, ind_PS50SI), c(ind_MD, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_MD_PS50SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_MD, ind_PS50SI), ]
ano_GMM_Fn_MD_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Fn_MD_PS50SI_exclCeph, metadata_GMM_Fn_dist_MD_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["MembraneDisruption_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Fn_MD_PS50SI_exclCeph$signif

# MembraneDisruption vs ProteinSynthesistRNAInterference
distance_GMM_Fn_MD_PStRI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_MD, ind_PStRI), c(ind_MD, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_MD_PStRI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_MD, ind_PStRI), ]
ano_GMM_Fn_MD_PStRI_exclCeph <- vegan::anosim(distance_GMM_Fn_MD_PStRI_exclCeph, metadata_GMM_Fn_dist_MD_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["MembraneDisruption_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Fn_MD_PStRI_exclCeph$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesis50SInhibition
distance_GMM_Fn_PS30SI_PS50SI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PS30SI, ind_PS50SI), c(ind_PS30SI, ind_PS50SI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PS30SI_PS50SI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_PS30SI, ind_PS50SI), ]
ano_GMM_Fn_PS30SI_PS50SI_exclCeph <- vegan::anosim(distance_GMM_Fn_PS30SI_PS50SI_exclCeph, metadata_GMM_Fn_dist_PS30SI_PS50SI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["ProteinSynthesis30SInhibition_ProteinSynthesis50SInhibition_exclCephalothin"] <- ano_GMM_Fn_PS30SI_PS50SI_exclCeph$signif

# ProteinSynthesis30SInhibition vs ProteinSynthesistRNAInterference
distance_GMM_Fn_PS30SI_PStRI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PS30SI, ind_PStRI), c(ind_PS30SI, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PS30SI_PStRI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_PS30SI, ind_PStRI), ]
ano_GMM_Fn_PS30SI_PStRI_exclCeph <- vegan::anosim(distance_GMM_Fn_PS30SI_PStRI_exclCeph, metadata_GMM_Fn_dist_PS30SI_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["ProteinSynthesis30SInhibition_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Fn_PS30SI_PStRI_exclCeph$signif

# ProteinSynthesis50SInhibition vs ProteinSynthesistRNAInterference
distance_GMM_Fn_PS50SI_PStRI_exclCeph <- distance_GMM_Fn_novel_compounds_noceph %>% 
  as.matrix() %>% 
  .[c(ind_PS50SI, ind_PStRI), c(ind_PS50SI, ind_PStRI)] %>% 
  as.dist()

metadata_GMM_Fn_dist_PS50SI_PStRI_exclCeph <- metadata_Fn_dist_novel_compounds[c(ind_PS50SI, ind_PStRI), ]
ano_GMM_Fn_PS50SI_PStRI_exclCeph <- vegan::anosim(distance_GMM_Fn_PS50SI_PStRI_exclCeph, metadata_GMM_Fn_dist_PS50SI_PStRI_exclCeph$MOA, distance = "bray", permutations = permutations)
pairwise_anosim_GMM_Fn["ProteinSynthesis50SInhibition_ProteinSynthesistRNAInterference_exclCephalothin"] <- ano_GMM_Fn_PS50SI_PStRI_exclCeph$signif

# Correction for pairwise comparisons
pairwise_anosim_GMM_Fn_bonferroni <- p.adjust(pairwise_anosim_GMM_Fn, method = "bonferroni")
pairwise_anosim_GMM_Fn_BH <- p.adjust(pairwise_anosim_GMM_Fn, method = "BH")

pairwise_anosim_GMM_Fn_df <- as.data.frame(pairwise_anosim_GMM_Fn)
pairwise_anosim_GMM_Fn_bonferroni_df <- as.data.frame(pairwise_anosim_GMM_Fn_bonferroni)
pairwise_anosim_GMM_Fn_BH_df <- as.data.frame(pairwise_anosim_GMM_Fn_BH)

pairwise_anosim_GMM_Fn_export <- merge(pairwise_anosim_GMM_Fn_df, pairwise_anosim_GMM_Fn_bonferroni_df, by = "row.names") %>% 
  column_to_rownames(var = "Row.names") %>% 
  merge(pairwise_anosim_GMM_Fn_BH_df, by = "row.names") %>% 
  rename("classes" = "Row.names", "p" = "pairwise_anosim_GMM_Fn", "p.adjust.bonferroni" = "pairwise_anosim_GMM_Fn_bonferroni", "p.adjust.BH" = "pairwise_anosim_GMM_Fn_BH")
write.csv2(file = "ANOSIM_GMM_Fn_pairwise_novel_compounds.csv", pairwise_anosim_GMM_Fn_export) # Export results as csv file


## 6.3. Both strains combined ----
# Diversity analysis for both strains combined (indicate whether strain or treatment is dominant)
flowData_transformed_GMM_novel_compounds <- flowCore::rbind2(flowData_transformed_Av_GMM_novel_compounds, flowData_transformed_Fn_GMM_novel_compounds)
metadata_GMM_novel_compounds <- rbind(metadata_Av_dist_noceph, metadata_Fn_dist_noceph)
row.names(metadata_GMM_novel_compounds) <- NULL

### 6.3.1. Optimization model ----
# Optimization of amount of phenotypes (Bayesian Information Criterion, BIC)
# Retain usefull information only and pool all samples by MOA and strain
fcs_PGMM_novel_compounds <- Phenoflow::FCS_pool(flowData_transformed_GMM_novel_compounds,
                                                stub = c("20210902_Fabian_MOA_Av_Experiment_CellWallSynthesis",
                                                         "20210902_Fabian_MOA_Av_Experiment_Control",
                                                         "20210902_Fabian_MOA_Av_Experiment_DNAReplication",
                                                         "20210902_Fabian_MOA_Av_Experiment_DNATranscription",
                                                         "20210902_Fabian_MOA_Av_Experiment_FolicAcidMetabolism",
                                                         "20210902_Fabian_MOA_Av_Experiment_Heat",
                                                         "20210902_Fabian_MOA_Av_Experiment_MembraneDisruption",
                                                         "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis30SInhibition",
                                                         "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis50SInhibition",
                                                         "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesistRNAInterference",
                                                         "20210902_Fabian_MOA_Fn_Experiment_CellWallSynthesis",
                                                         "20210902_Fabian_MOA_Fn_Experiment_Control",
                                                         "20210902_Fabian_MOA_Fn_Experiment_DNAReplication",
                                                         "20210902_Fabian_MOA_Fn_Experiment_DNATranscription",
                                                         "20210902_Fabian_MOA_Fn_Experiment_FolicAcidMetabolism",
                                                         "20210902_Fabian_MOA_Fn_Experiment_Heat",
                                                         "20210902_Fabian_MOA_Fn_Experiment_MembraneDisruption",
                                                         "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis30SInhibition",
                                                         "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis50SInhibition",
                                                         "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesistRNAInterference"))
fcs_PGMM_novel_compounds <- FCS_resample(fcs_PGMM_novel_compounds, replace = TRUE, sample = 50000)
fcs_PGMM_novel_compounds <- fcs_PGMM_novel_compounds[, paramGMM]
fcs_PGMM_novel_compounds <- Phenoflow::FCS_pool(fcs_PGMM_novel_compounds, stub = "*")
fcs_PGMM_novel_compounds <- fcs_PGMM_novel_compounds[, paramGMM]
PhenoGMM_noneg_novel_compounds <- Phenoflow::PhenoGMM(fcs_x = fcs_PGMM_novel_compounds, param = paramGMM, downsample = FALSE, nG = 80, auto_nG = TRUE, nG_interval = 10, fcs_scale = FALSE, diagnostic_plot = TRUE)
saveRDS(object = PhenoGMM_noneg_novel_compounds, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/BothStrains_24h_PhenoGMM_8param_10to80per10_novel_compounds.rds")

# Visualization of BIC values
#PhenoGMM_noneg_novel_compounds <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/BothStrains_24h_PhenoGMM_8param_10to80per10_novel_compounds.rds") # Call for results that are to be plotted
NoClusters_PGMM <- dimnames(PhenoGMM_noneg_novel_compounds[[2]]$BIC)[[1]]
BICValues_PGMM <- data.frame(as.matrix(PhenoGMM_noneg_novel_compounds[[2]]$BIC)[1:length(NoClusters_PGMM), ])
BICValues_PGMM$NoClusters <- rownames(BICValues_PGMM)
BICValues_PGMM <- reshape2::melt(BICValues_PGMM, id.vars = "NoClusters")
colnames(BICValues_PGMM) <- c("NoClusters", "ModelType", "BIC")
BICValues_PGMM$NoClusters <- as.numeric(BICValues_PGMM$NoClusters)
BICValues_PGMM <- BICValues_PGMM[!is.na(BICValues_PGMM$BIC), ] # Remove NA values
BICValues_PGMM$ModelType <- droplevels(BICValues_PGMM$ModelType, except = unique(BICValues_PGMM$ModelType)) # Remove levels that are not being used

p_BIC_PGMM_novel_compounds <- BICValues_PGMM %>% 
  ggplot(data = ., aes(x = NoClusters, y = BIC))+
  geom_line(alpha = 1, aes (color = ModelType), show.legend = FALSE)+
  geom_point(shape = 21, size = 3, alpha = 1, aes(fill = ModelType))+
  labs(title = NULL, x = "Number of clusters", y = "BIC", fill = "Model type")+
  theme_bw()+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 34, face = "bold"),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 28),
        legend.text.align = 0)
print(p_BIC_PGMM_novel_compounds)

### 6.3.2. Allocation data to model ----
## Make model for optimum clusters determined in previous step (excluding Cephalothin)
NoClusters_PGMM_novel_compounds <- 60
PhenoGMM_fixedclust_novel_compounds <- Phenoflow::PhenoGMM(fcs_x = fcs_PGMM_novel_compounds, param = paramGMM, downsample = FALSE, nG = NoClusters_PGMM_novel_compounds, fcs_scale = FALSE, diagnostic_plot = TRUE)
saveRDS(object = PhenoGMM_fixedclust_novel_compounds, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/BothStrains_24h_PhenoGMM_8param_60clust_novel_compounds.rds")
#PhenoGMM_fixedclust_novel_compounds <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/GMMFits/BothStrains_24h_PhenoGMM_8param_60clust_novel_compounds.rds") # Call for results that are to be plotted

# Applying mask to data
testPred_novel_compounds <- PhenoMaskGMM(fcs_x = flowData_transformed_GMM_novel_compounds, gmm = PhenoGMM_fixedclust_novel_compounds, fcs_scale = FALSE)

results_novel_compounds <- testPred_novel_compounds[[1]]
rownames(results_novel_compounds) <- sampleNames(flowData_transformed_GMM_novel_compounds)
results_novel_compounds <- select(results_novel_compounds, -c(Sample_names))
results_novel_compounds[is.na(results_novel_compounds)] <- 0 # Replace NA values by 0
results_novel_compounds[1:ncol(results_novel_compounds)] <- lapply(results_novel_compounds[1:ncol(results_novel_compounds)], as.numeric)
# Normalize abundances to sum = 1
results_rel_novel_compounds <- sweep(results_novel_compounds, MARGIN = 1, rowSums(results_novel_compounds), `/`)

### 6.3.3. Ordination of GMM output ----
rownames_GMM_novel_compounds <- row.names(results_novel_compounds)
new_rownames_GMM_novel_compounds <- gsub("^.*?MOA_", "", rownames_GMM_novel_compounds)
new_rownames_GMM_novel_compounds <- gsub("_Experiment", "", new_rownames_GMM_novel_compounds)
new_rownames_GMM_novel_compounds <- gsub("_100_SGPI", "", new_rownames_GMM_novel_compounds)
new_rownames_GMM_novel_compounds <- gsub(".fcs", "", new_rownames_GMM_novel_compounds)
results_ordination_novel_compounds <- results_novel_compounds
rownames(results_ordination_novel_compounds) <- new_rownames_GMM_novel_compounds
results_ordination_novel_compounds <- as.matrix(results_ordination_novel_compounds)

distance_GMM_novel_compounds <- vegan::vegdist(x = results_ordination_novel_compounds, method = "bray", binary = FALSE)
mds_GMM_novel_compounds <- stats::cmdscale(distance_GMM_novel_compounds, k = 2, eig = TRUE, add = TRUE)
NMDS_GMM_novel_compounds <- vegan::metaMDS(distance_GMM_novel_compounds, autotransform = FALSE, k = 2, trymax = 100)

# Plot ordination
plot_PCoA_GMM_novel_compounds <- plot_beta_fcm_custom(mds_GMM_novel_compounds, color = as.factor(metadata_GMM_novel_compounds$MOA), shape = as.factor(metadata_GMM_novel_compounds$Compound), size = as.factor(metadata_GMM_novel_compounds$Strain))+ 
  theme_bw() +
  geom_point(alpha = 0.7)+
  labs(title = NULL, color = "Mode of Action", shape = "Compound", size = "Strain")+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_size_manual(values = c(3, 8),
                    labels = c(expression(italic("Actinomyces viscosus")), expression(italic("Fusobacterium nucleatum"))))+
  scale_shape_manual(values = c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        legend.text.align = 0,
        legend.position = "bottom")
print(plot_PCoA_GMM_novel_compounds)

plot_NMDS_GMM_novel_compounds <- plot_beta_fcm_custom(NMDS_GMM_novel_compounds, color = as.factor(metadata_GMM_novel_compounds$MOA), shape = as.factor(metadata_GMM_novel_compounds$Compound), size = as.factor(metadata_GMM_novel_compounds$Strain))+ 
  theme_bw() +
  geom_point(alpha = 0.7)+
  labs(title = NULL, color = "Mode of Action", shape = "Compound", size = "Strain", x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("purple", "steelblue1", "blue", "seagreen3", "gray0", "red", "chartreuse", "darkgoldenrod1", "yellow4", "yellow1"),
                     labels = c("Cell Wall Synthesis", "DNA Replication", "DNA Transcription", "Folic Acid Metabolism", "Heat", "Membrane Disruption", "Positive Control", "Protein Synthesis: 30S Inhibition", "Protein Synthesis: 50S Inhibition", "Protein Synthesis: tRNA Interference"))+
  scale_size_manual(values = c(3, 8),
                    labels = c(expression(italic("Actinomyces viscosus")), expression(italic("Fusobacterium nucleatum"))))+
  scale_shape_manual(values = c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        legend.text.align = 0,
        legend.position = "bottom")
print(plot_NMDS_GMM_novel_compounds)

plot_PCoA_GMM_novel_compounds_color <- plot_beta_fcm(mds_GMM_novel_compounds, color = as.factor(metadata_GMM_novel_compounds$Strain), shape = as.factor(metadata_GMM_novel_compounds$Compound))+ 
  theme_bw() +
  geom_point(alpha = 0.7)+
  labs(title = NULL, color = "Strain", shape = "Compound")+
  scale_color_manual(values = c("blue", "red"),
                     labels = c(expression(italic("Actinomyces viscosus")), expression(italic("Fusobacterium nucleatum"))))+
  scale_shape_manual(values = c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        legend.text.align = 0,
        legend.position = "bottom")
print(plot_PCoA_GMM_novel_compounds_color)

plot_NMDS_GMM_novel_compounds_color <- plot_beta_fcm(NMDS_GMM_novel_compounds, color = as.factor(metadata_GMM_novel_compounds$Strain), shape = as.factor(metadata_GMM_novel_compounds$Compound))+ 
  theme_bw() +
  geom_point(alpha = 0.7)+
  labs(title = NULL, color = "Strain", shape = "Compound", x = "NMDS1", y = "NMDS2")+
  scale_color_manual(values = c("blue", "red"),
                     labels = c(expression(italic("Actinomyces viscosus")), expression(italic("Fusobacterium nucleatum"))))+
  scale_shape_manual(values = c(1:15, 35, 17, 18, 42, 16, 43, 60, 62, 94))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        legend.text.align = 0,
        legend.position = "bottom")
print(plot_NMDS_GMM_novel_compounds_color)


# 7. Random Forest Classification ----
# Random forest classifiers will be trained based on the output of PhenoGMM.

## 7.1. Matthews correlation coeficient function ----
# Call for function for using MCC as optimization metric in 'caret::train()'
source("/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/MatthewsCC.R")

## 7.2. Av ----

### 7.2.1 Formatting data----
# Add mode of action to results PhenoGMM data frame
results_Av_RF_novel_compounds <- results_Av_rel_novel_compounds %>% 
  add_column(MOA = metadata_Av_noneg$MOA)
results_Av_RF_novel_compounds <- select(results_Av_RF_novel_compounds, -c(Sample_name))

# Change class of MOA
results_Av_RF_novel_compounds$MOA <- as.factor(results_Av_RF_novel_compounds$MOA)

# Change column names for different phenotypes as radomForest function does not take integers as column names
colnames(results_Av_RF_novel_compounds)[1:NoClusters_PGMM_Av_novel_compounds] <- str_c("Cluster", colnames((results_Av_RF_novel_compounds)[1:NoClusters_PGMM_Av_novel_compounds]))

# Remove data for cephalothin (the unseen compound)
results_Av_RF_novel_compounds_no_ceph <- results_Av_RF_novel_compounds[c(1, 5:72), ]


### 7.2.2. Training RF classifier ----

#### 7.2.2.1. Repeated cross-validation ----
# 5-fold repeated cross-validation with 3 repeats
## Train RF model
train_control_RCV <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 3, summaryFunction = MatthewsCC)
mtry_RCV <- c(1:50)
tunegrid_RCV <- base::expand.grid(.mtry = mtry_RCV)
Model_RF_PGMM_Av_caret_RCV <- caret::train(MOA ~ ., data = results_Av_RF_novel_compounds_no_ceph, method = "rf", tuneGrid = tunegrid_RCV, trControl = train_control_RCV, metric = "MCC", ntree = 500)
saveRDS(Model_RF_PGMM_Av_caret_RCV, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/models_rds/Av_24h_8param_PhenoGMM_caret_RCV.rds")

plot_RF_PGMM_Av_caret_RCV <- ggplot(Model_RF_PGMM_Av_caret_RCV)+
  labs(x = "mtry", y = "Matthews Correlation Coefficient (Repeated Cross-Validation)", title = "Model performance (Random Forest caret RCV)")+
  theme_bw()
print(plot_RF_PGMM_Av_caret_RCV)

# Variable importance
variable_importance_Av_caret_RCV <- caret::varImp(Model_RF_PGMM_Av_caret_RCV)
variable_importance_Av_caret_RCV <- as.data.frame(variable_importance_Av_caret_RCV$importance)
variable_importance_Av_caret_RCV$Variable <- rownames(variable_importance_Av_caret_RCV)
# Plot importances
variable_importance_Av_caret_RCV <- variable_importance_Av_caret_RCV[order(variable_importance_Av_caret_RCV$Overall, decreasing = FALSE),]
variable_importance_Av_caret_RCV$Variable <- factor(variable_importance_Av_caret_RCV$Variable , levels = unique(variable_importance_Av_caret_RCV$Variable))
p_VarImp_Av_caret_RCV <- variable_importance_Av_caret_RCV %>% 
  ggplot(data = ., aes(x = Overall, y = Variable)) +
  geom_point(size = 3) +
  labs(title = "Variable importance RF model Av (caret RCV)", x = "Importance (%)", y = "") +
  theme_bw()
print(p_VarImp_Av_caret_RCV)

# Confusion matrix
confusion_RF_PGMM_Av_caret_RCV <- data.frame(Model_RF_PGMM_Av_caret_RCV$finalModel$confusion)
confusion_RF_PGMM_Av_caret_RCV <- confusion_RF_PGMM_Av_caret_RCV[, colnames(confusion_RF_PGMM_Av_caret_RCV) != "class.error"]
confusion_RF_PGMM_Av_caret_RCV <- tibble::rownames_to_column(confusion_RF_PGMM_Av_caret_RCV, "Reference")

melted_confusion_RF_PGMM_Av_caret_RCV <- melt(confusion_RF_PGMM_Av_caret_RCV, id.vars = c("Reference"), variable.name = "Prediction", value.name = "Freq")

# Plot confusion matrix
p_confusion_Av_caret_RCV <- melted_confusion_RF_PGMM_Av_caret_RCV %>% 
  ggplot(data = ., aes(x = Prediction, y = Reference, fill = 100 * Freq/sum(Freq))) +
  geom_raster() +
  geom_text(aes(label = round(100 * Freq/sum(Freq), 0)), size = 6) +
  scale_fill_distiller(name = "% of total samples\n classified\n") +
  theme_bw() +
  scale_x_discrete(position = "top")+
  theme(axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 14, angle = 55, hjust = 0), title = element_text(size = 20),
        plot.margin = unit(c(1.1, 1.1, 1.1, 1.1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = 'transparent', colour = NA))

MCC_Av_RCV_finalModel <- Model_RF_PGMM_Av_caret_RCV$results$MCC[Model_RF_PGMM_Av_caret_RCV$finalModel$mtry]
MCCSD_Av_RCV_finalModel <- Model_RF_PGMM_Av_caret_RCV$results$MCCSD[Model_RF_PGMM_Av_caret_RCV$finalModel$mtry]

table_performance_Av_RCV <- base::round(data.frame(Performance = c(MCC_Av_RCV_finalModel, MCCSD_Av_RCV_finalModel)), 2)
row.names(table_performance_Av_RCV) <- c("Matthews Correlation Coefficient", "Matthews Correlation Coefficient SD")

confusion_RF_Av_RCV_prep <- data.frame(lapply(melted_confusion_RF_PGMM_Av_caret_RCV, rep, melted_confusion_RF_PGMM_Av_caret_RCV$Freq))

confusion_RF_Av_caret_RCV_finalModel <- caret::confusionMatrix(data = as.factor(confusion_RF_Av_RCV_prep$Prediction), reference = as.factor(confusion_RF_Av_RCV_prep$Reference))
table_performance_RF_Av_confusion_RCV <- base::round(data.frame(Performance = confusion_RF_Av_caret_RCV_finalModel$overall), 2)
table_performance_Av_RCV <- rbind(table_performance_Av_RCV, table_performance_RF_Av_confusion_RCV)

p_confusion_Av_caret_RCV_table <- data.frame(melted_confusion_RF_PGMM_Av_caret_RCV) %>% 
  ggplot(data = ., aes(x = Prediction, y = Reference, fill =  100 * Freq/sum(Freq)))+
  theme(axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13, angle = 55, hjust = 0), title = element_text(size = 20),
        plot.margin = unit(c(1.1, 1.1, 1.1, 1.1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = 'transparent', colour = NA))+
  ylab("")+
  xlab("")+
  annotation_custom(gridExtra::tableGrob(table_performance_Av_RCV))

print(cowplot::plot_grid(p_confusion_Av_caret_RCV, p_confusion_Av_caret_RCV_table, align = "h", ncol = 2, rel_widths = c(1/2, 1/4)))

#### 7.2.2.2. Nested repeated cross-validation ----
# 5-fold nested repeated cross-validation: 5 times 20% of data set is kept aside for testing in order to get less biased representation of model performance (final model for predictions of unseen data is the _RCV model from above)
predictions_Av_novel_compounds <- NULL
folds_Av <- createFolds(results_Av_RF_novel_compounds_no_ceph$MOA, k = 5, list = TRUE, returnTrain = FALSE)
train_control_nested <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 3, summaryFunction = MatthewsCC)
mtry_nested <- c(1:50)
tunegrid_nested <- base::expand.grid(.mtry = mtry_nested)
Model_RF_list_Av_nested <- list()
confusion_RF_list_Av_nested <- list()

for (fold in 1:length(folds_Av)){
  # Selected fold
  testindices <- folds_Av[[fold]]
  trainingset <- results_Av_RF_novel_compounds_no_ceph[-testindices, ]
  testset <- results_Av_RF_novel_compounds_no_ceph[testindices, ]
  # Train model according to the train.control scheme
  model <- caret::train(MOA ~ ., data = trainingset, method = "rf", tuneGrid = tunegrid_nested, trControl = train_control_nested, metric = "MCC", ntree = 500)
  Model_RF_list_Av_nested[[fold]] <- model
  # Make predictions on the test set and save the results
  predictionfold <- stats::predict(model, newdata = testset)
  confusion_RF_list_Av_nested[[fold]] <- caret::confusionMatrix(data = predictionfold, as.factor(testset$MOA))
  predictionfold <- cbind.data.frame(fold, predictionfold, testset$MOA)
  predictions_Av_novel_compounds <- rbind.data.frame(predictions_Av_novel_compounds, predictionfold)
}

colnames(predictions_Av_novel_compounds) <- c("Fold", "Predicted", "True")

MCC_Av_nested <- mltools::mcc(preds = predictions_Av_novel_compounds$Predicted, actuals = predictions_Av_novel_compounds$True)
confusion_RF_Av_nested <- caret::confusionMatrix(data = predictions_Av_novel_compounds$Predicted, reference = predictions_Av_novel_compounds$True)

p_confusion_Av_caret_nested <- data.frame(confusion_RF_Av_nested$table) %>% 
  ggplot(data = ., aes(x = Prediction, y = Reference, fill = 100 * Freq/sum(Freq))) +
  geom_raster() +
  geom_text(aes(label = round(100 * Freq/sum(Freq), 0)), size = 6) +
  scale_fill_distiller(name = "% of total samples\n classified\n") +
  theme_bw() +
  scale_x_discrete(position = "top")+
  theme(axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13, angle = 55, hjust = 0),
        title = element_text(size = 20),
        plot.margin = unit(c(1.1, 1.1, 1.1, 1.1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))

table_performance_Av_nested <- base::round(data.frame(Performance = confusion_RF_Av_nested$overall), 2)
MCC_Av_nested_table <- base::round(data.frame(Performance = c(MCC_Av_nested)), 2)
row.names(MCC_Av_nested_table) <- c("MatthewsCorrelationCoefficient")
table_performance_Av_nested <- rbind(MCC_Av_nested_table, table_performance_Av_nested)

p_confusion_Av_caret_nested_table <- data.frame(confusion_RF_Av_nested$table) %>% 
  ggplot(data = ., aes(x = Prediction, y = Reference, fill =  100 * Freq/sum(Freq)))+
  theme(axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13, angle = 55, hjust = 0), title = element_text(size = 20),
        plot.margin = unit(c(1.1, 1.1, 1.1, 1.1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  ylab("")+
  xlab("")+
  annotation_custom(gridExtra::tableGrob(table_performance_Av_nested))

print(cowplot::plot_grid(p_confusion_Av_caret_nested, p_confusion_Av_caret_nested_table, align = "h", ncol = 2, rel_widths = c(1/2, 1/4)))

### 7.2.3. Predict unseen data ----
cephalothin_Av <- results_Av_RF_novel_compounds[c(2:4), ]
prediction_novel_compound_Av <- stats::predict(Model_RF_PGMM_Av_caret_RCV, newdata = cephalothin_Av)
prediction_novel_compound_Av_prob <- stats::predict(Model_RF_PGMM_Av_caret_RCV, newdata = cephalothin_Av, type = "prob")

prediction_novel_compound_Av_exp <- cbind.data.frame(Prediction = prediction_novel_compound_Av, Reference = cephalothin_Av$MOA, Probabiltiy = prediction_novel_compound_Av_prob$CellWallSynthesis)
write.csv2(file = "RF_GMM_Av_pred_Ceph.csv", prediction_novel_compound_Av_exp) # Export results as csv file

### 7.2.4. Random forest using raw fcs data ----
# Define parameters on which RF will build its model
paramRF <- c("FSC-H", "SSC-H", "BL1-H", "BL3-H", "FSC-A", "SSC-A", "BL1-A", "BL3-A")
# Load custom RandomF_FCS function with MCC as optimization metric
source("/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/RandomF_FCS_MCC.R")

# Gather necessary information for training of the algorithm
Sample_Info_Av <- metadata_Av %>% rename(name = Filename) # Necessary to make use of the RandomF_FCS function

# Select the fcs files based on which the model will be trained (cephalothin will be excluded)
fcs_names_Av <- c("20210902_Fabian_MOA_Av_Experiment_CellWallSynthesis_100_SGPI_C1.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_CellWallSynthesis_100_SGPI_C2.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_CellWallSynthesis_100_SGPI_C3.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_CellWallSynthesis_100_SGPI_C4.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_CellWallSynthesis_100_SGPI_C5.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_CellWallSynthesis_100_SGPI_C6.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_CellWallSynthesis_100_SGPI_C7.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_CellWallSynthesis_100_SGPI_C8.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_CellWallSynthesis_100_SGPI_C9.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_Control_100_SGPI_G10.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_Control_100_SGPI_G11.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_Control_100_SGPI_G12.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_Control_100_SGPI_G7.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_Control_100_SGPI_G8.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_Control_100_SGPI_G9.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNAReplication_100_SGPI_D1.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNAReplication_100_SGPI_D2.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNAReplication_100_SGPI_D3.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNAReplication_100_SGPI_D4.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNAReplication_100_SGPI_D5.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNAReplication_100_SGPI_D6.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNAReplication_100_SGPI_D7.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNAReplication_100_SGPI_D8.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNAReplication_100_SGPI_D9.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNATranscription_100_SGPI_D10.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNATranscription_100_SGPI_D11.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNATranscription_100_SGPI_D12.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNATranscription_100_SGPI_E1.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNATranscription_100_SGPI_E2.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_DNATranscription_100_SGPI_E3.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_FolicAcidMetabolism_100_SGPI_G1.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_FolicAcidMetabolism_100_SGPI_G2.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_FolicAcidMetabolism_100_SGPI_G3.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_Heat_100_SGPI_G4.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_Heat_100_SGPI_G5.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_Heat_100_SGPI_G6.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_MembraneDisruption_100_SGPI_B1.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_MembraneDisruption_100_SGPI_B10.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_MembraneDisruption_100_SGPI_B11.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_MembraneDisruption_100_SGPI_B12.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_MembraneDisruption_100_SGPI_B2.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_MembraneDisruption_100_SGPI_B3.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_MembraneDisruption_100_SGPI_B4.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_MembraneDisruption_100_SGPI_B5.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_MembraneDisruption_100_SGPI_B6.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_MembraneDisruption_100_SGPI_B7.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_MembraneDisruption_100_SGPI_B8.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_MembraneDisruption_100_SGPI_B9.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F1.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F2.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F3.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F4.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F5.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F6.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F7.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F8.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F9.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E10.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E11.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E12.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E4.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E5.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E6.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E7.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E8.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E9.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesistRNAInterference_100_SGPI_F10.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesistRNAInterference_100_SGPI_F11.fcs",
                  "20210902_Fabian_MOA_Av_Experiment_ProteinSynthesistRNAInterference_100_SGPI_F12.fcs")

Sample_Info_Av1 <- Sample_Info_Av %>% dplyr::filter(name %in% fcs_names_Av)
Model_RF_Av1_MCC <- RandomF_FCS_MCC(flowData_transformed_Av_gated[fcs_names_Av], sample_info = Sample_Info_Av1, target_label = "MOA", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.8, seed = seed, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = FALSE, metric = "MCC")
saveRDS(object = Model_RF_Av1_MCC, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/models_rds/Av_24h_8param_100000cells_MCC_novel_compounds.rds")
#Model_RF_Av1_MCC <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/models_rds/Av_8param_100000cells_MCC_novel_compounds.rds")

mytable_FCS_MCC_Av <- base::round(data.frame(Performance = Model_RF_Av1_MCC[[2]]$overall), 2)

confusion_RF_Av_FCS_MCC <- data.frame(Model_RF_Av1_MCC[[2]]$table)
confusion_RF_Av_FCS_MCC_prep <- data.frame(lapply(confusion_RF_Av_FCS_MCC, rep, confusion_RF_Av_FCS_MCC$Freq))
MCC_FCS_Av <- mltools::mcc(preds = as.factor(confusion_RF_Av_FCS_MCC_prep$Prediction), actuals = as.factor(confusion_RF_Av_FCS_MCC_prep$Reference))
MCC_FCS_Av_table <- base::round(data.frame(Performance = MCC_FCS_Av), 2)
row.names(MCC_FCS_Av_table) <- c("MatthewsCorrelationCoefficient")
mytable_FCS_MCC_Av_2 <- rbind(MCC_FCS_Av_table, mytable_FCS_MCC_Av)

p_conf_Av_FCS_MCC <- ggplot2::ggplot(data.frame(Model_RF_Av1_MCC[[2]]$table), 
                                     ggplot2::aes(x = Prediction, y = Reference, fill = 100 * Freq/sum(Freq))) +
  ggplot2::geom_raster() + 
  ggplot2::geom_text(ggplot2::aes(label = round(100 * Freq/sum(Freq), 0)), size = 6) +
  ggplot2::scale_fill_distiller(name = "% of total samples\n classified\n") + 
  ggplot2::theme_bw() +
  ggplot2::scale_x_discrete(position = "top") +
  ggplot2::theme(axis.title = element_text(size = 16),
                 strip.text.x = element_text(size = 14),
                 legend.title = element_text(size = 14),
                 legend.text = element_text(size = 14),
                 axis.text.y = element_text(size = 13),
                 axis.text.x = element_text(size = 13, angle = 55, hjust = 0),
                 title = element_text(size = 20),
                 plot.margin = unit(c(1.1, 1.1, 1.1, 1.1), "cm"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_rect(fill = "transparent", colour = NA),
                 plot.background = element_rect(fill = "transparent", colour = NA))

p_conf_table_Av_FCS_MCC <- ggplot2::ggplot(data.frame(Model_RF_Av1_MCC[[2]]$table), 
                                           ggplot2::aes(x = Prediction, y = Reference, fill = 100 * Freq/sum(Freq))) +
  theme(axis.title = ggplot2::element_text(size = 16), 
        strip.text.x = ggplot2::element_text(size = 14), 
        legend.title = ggplot2::element_text(size = 14), 
        legend.text = ggplot2::element_text(size = 14), 
        axis.text.y = ggplot2::element_text(size = 13), 
        axis.text.x = ggplot2::element_text(size = 13, angle = 55, hjust = 0),
        title = ggplot2::element_text(size = 20), 
        plot.margin = ggplot2::unit(c(1.1, 1.1, 1.1, 1.1), "cm"),
        panel.grid.major = ggplot2::element_blank(), 
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(), 
        panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
        plot.background = ggplot2::element_rect(fill = "transparent", colour = NA)) +
  ggplot2::ylab("") +
  ggplot2::xlab("") + 
  ggplot2::annotation_custom(gridExtra::tableGrob(mytable_FCS_MCC_Av_2))

print(cowplot::plot_grid(p_conf_Av_FCS_MCC, p_conf_table_Av_FCS_MCC, align = "h", 
                         ncol = 2, rel_widths = c(1/2, 1/4)))

## 7.3. Fn ----

### 7.3.1 Formatting data----
# Add mode of action to results PhenoGMM data frame
results_Fn_RF_novel_compounds <- results_Fn_rel_novel_compounds %>% 
  add_column(MOA = metadata_Fn_noneg$MOA)
results_Fn_RF_novel_compounds <- select(results_Fn_RF_novel_compounds, -c(Sample_name))

# Change class of MOA
results_Fn_RF_novel_compounds$MOA <- as.factor(results_Fn_RF_novel_compounds$MOA)

# Change column names for different phenotypes as radomForest function does not take integers as column names
colnames(results_Fn_RF_novel_compounds)[1:NoClusters_PGMM_Fn_novel_compounds] <- str_c("Cluster", colnames((results_Fn_RF_novel_compounds)[1:NoClusters_PGMM_Fn_novel_compounds]))

# Remove data for cephalothin (the unseen compound)
results_Fn_RF_novel_compounds_no_ceph <- results_Fn_RF_novel_compounds[c(1, 5:72), ]


### 7.3.2. Training RF classifier ----

#### 7.3.2.1. Repeated cross-validation ----
# 5-fold repeated cross-validation with 3 repeats
## Train RF model
train_control_RCV <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 3, summaryFunction = MatthewsCC)
mtry_RCV <- c(1:50)
tunegrid_RCV <- base::expand.grid(.mtry = mtry_RCV)
Model_RF_PGMM_Fn_caret_RCV <- caret::train(MOA ~ ., data = results_Fn_RF_novel_compounds_no_ceph, method = "rf", tuneGrid = tunegrid_RCV, trControl = train_control_RCV, metric = "MCC", ntree = 500)
saveRDS(Model_RF_PGMM_Fn_caret_RCV, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/models_rds/Fn_24h_8param_PhenoGMM_caret_RCV.rds")

plot_RF_PGMM_Fn_caret_RCV <- ggplot(Model_RF_PGMM_Fn_caret_RCV)+
  labs(x = "mtry", y = "Matthews Correlation Coefficient (Repeated Cross-Validation)", title = "Model performance (Random Forest caret RCV)")+
  theme_bw()
print(plot_RF_PGMM_Fn_caret_RCV)

# Variable importance
variable_importance_Fn_caret_RCV <- caret::varImp(Model_RF_PGMM_Fn_caret_RCV)
variable_importance_Fn_caret_RCV <- as.data.frame(variable_importance_Fn_caret_RCV$importance)
variable_importance_Fn_caret_RCV$Variable <- rownames(variable_importance_Fn_caret_RCV)
# Plot importances
variable_importance_Fn_caret_RCV <- variable_importance_Fn_caret_RCV[order(variable_importance_Fn_caret_RCV$Overall, decreasing = FALSE),]
variable_importance_Fn_caret_RCV$Variable <- factor(variable_importance_Fn_caret_RCV$Variable , levels = unique(variable_importance_Fn_caret_RCV$Variable))
p_VarImp_Fn_caret_RCV <- variable_importance_Fn_caret_RCV %>% 
  ggplot(data = ., aes(x = Overall, y = Variable)) +
  geom_point(size = 3) +
  labs(title = "Variable importance RF model Fn (caret RCV)", x = "Importance (%)", y = "") +
  theme_bw()
print(p_VarImp_Fn_caret_RCV)

# Confusion matrix
confusion_RF_PGMM_Fn_caret_RCV <- data.frame(Model_RF_PGMM_Fn_caret_RCV$finalModel$confusion)
confusion_RF_PGMM_Fn_caret_RCV <- confusion_RF_PGMM_Fn_caret_RCV[, colnames(confusion_RF_PGMM_Fn_caret_RCV) != "class.error"]
confusion_RF_PGMM_Fn_caret_RCV <- tibble::rownames_to_column(confusion_RF_PGMM_Fn_caret_RCV, "Reference")

melted_confusion_RF_PGMM_Fn_caret_RCV <- melt(confusion_RF_PGMM_Fn_caret_RCV, id.vars = c("Reference"), variable.name = "Prediction", value.name = "Freq")

# Plot confusion matrix
p_confusion_Fn_caret_RCV <- melted_confusion_RF_PGMM_Fn_caret_RCV %>% 
  ggplot(data = ., aes(x = Prediction, y = Reference, fill = 100 * Freq/sum(Freq))) +
  geom_raster() +
  geom_text(aes(label = round(100 * Freq/sum(Freq), 0)), size = 6) +
  scale_fill_distiller(name = "% of total samples\n classified\n") +
  theme_bw() +
  scale_x_discrete(position = "top")+
  theme(axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 14, angle = 55, hjust = 0), title = element_text(size = 20),
        plot.margin = unit(c(1.1, 1.1, 1.1, 1.1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = 'transparent', colour = NA))

MCC_Fn_RCV_finalModel <- Model_RF_PGMM_Fn_caret_RCV$results$MCC[Model_RF_PGMM_Fn_caret_RCV$finalModel$mtry]
MCCSD_Fn_RCV_finalModel <- Model_RF_PGMM_Fn_caret_RCV$results$MCCSD[Model_RF_PGMM_Fn_caret_RCV$finalModel$mtry]

table_performance_Fn_RCV <- base::round(data.frame(Performance = c(MCC_Fn_RCV_finalModel, MCCSD_Fn_RCV_finalModel)), 2)
row.names(table_performance_Fn_RCV) <- c("Matthews Correlation Coefficient", "Matthews Correlation Coefficient SD")

confusion_RF_Fn_RCV_prep <- data.frame(lapply(melted_confusion_RF_PGMM_Fn_caret_RCV, rep, melted_confusion_RF_PGMM_Fn_caret_RCV$Freq))

confusion_RF_Fn_caret_RCV_finalModel <- caret::confusionMatrix(data = as.factor(confusion_RF_Fn_RCV_prep$Prediction), reference = as.factor(confusion_RF_Fn_RCV_prep$Reference))
table_performance_RF_Fn_confusion_RCV <- base::round(data.frame(Performance = confusion_RF_Fn_caret_RCV_finalModel$overall), 2)
table_performance_Fn_RCV <- rbind(table_performance_Fn_RCV, table_performance_RF_Fn_confusion_RCV)

p_confusion_Fn_caret_RCV_table <- data.frame(melted_confusion_RF_PGMM_Fn_caret_RCV) %>% 
  ggplot(data = ., aes(x = Prediction, y = Reference, fill =  100 * Freq/sum(Freq)))+
  theme(axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13, angle = 55, hjust = 0), title = element_text(size = 20),
        plot.margin = unit(c(1.1, 1.1, 1.1, 1.1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = 'transparent', colour = NA))+
  ylab("")+
  xlab("")+
  annotation_custom(gridExtra::tableGrob(table_performance_Fn_RCV))

print(cowplot::plot_grid(p_confusion_Fn_caret_RCV, p_confusion_Fn_caret_RCV_table, align = "h", ncol = 2, rel_widths = c(1/2, 1/4)))

#### 7.3.2.2. Nested repeated cross-validation ----
# 5-fold nested repeated cross-validation: 5 times 20% of data set is kept aside for testing in order to get less biased representation of model performance (final model for predictions of unseen data is the _RCV model from above)
predictions_Fn_novel_compounds <- NULL
folds_Fn <- createFolds(results_Fn_RF_novel_compounds_no_ceph$MOA, k = 5, list = TRUE, returnTrain = FALSE)
train_control_nested <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 3, summaryFunction = MatthewsCC)
mtry_nested <- c(1:50)
tunegrid_nested <- base::expand.grid(.mtry = mtry_nested)
Model_RF_list_Fn_nested <- list()
confusion_RF_list_Fn_nested <- list()

for (fold in 1:length(folds_Fn)){
  # Selected fold
  testindices <- folds_Fn[[fold]]
  trainingset <- results_Fn_RF_novel_compounds_no_ceph[-testindices, ]
  testset <- results_Fn_RF_novel_compounds_no_ceph[testindices, ]
  # Train model according to the train.control scheme
  model <- caret::train(MOA ~ ., data = trainingset, method = "rf", tuneGrid = tunegrid_nested, trControl = train_control_nested, metric = "MCC", ntree = 500)
  Model_RF_list_Fn_nested[[fold]] <- model
  # Make predictions on the test set and sFne the results
  predictionfold <- stats::predict(model, newdata = testset)
  confusion_RF_list_Fn_nested[[fold]] <- caret::confusionMatrix(data = predictionfold, as.factor(testset$MOA))
  predictionfold <- cbind.data.frame(fold, predictionfold, testset$MOA)
  predictions_Fn_novel_compounds <- rbind.data.frame(predictions_Fn_novel_compounds, predictionfold)
}

colnames(predictions_Fn_novel_compounds) <- c("Fold", "Predicted", "True")

MCC_Fn_nested <- mltools::mcc(preds = predictions_Fn_novel_compounds$Predicted, actuals = predictions_Fn_novel_compounds$True)
confusion_RF_Fn_nested <- caret::confusionMatrix(data = predictions_Fn_novel_compounds$Predicted, reference = predictions_Fn_novel_compounds$True)

p_confusion_Fn_caret_nested <- data.frame(confusion_RF_Fn_nested$table) %>% 
  ggplot(data = ., aes(x = Prediction, y = Reference, fill = 100 * Freq/sum(Freq))) +
  geom_raster() +
  geom_text(aes(label = round(100 * Freq/sum(Freq), 0)), size = 6) +
  scale_fill_distiller(name = "% of total samples\n classified\n") +
  theme_bw() +
  scale_x_discrete(position = "top")+
  theme(axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13, angle = 55, hjust = 0),
        title = element_text(size = 20),
        plot.margin = unit(c(1.1, 1.1, 1.1, 1.1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))

table_performance_Fn_nested <- base::round(data.frame(Performance = confusion_RF_Fn_nested$overall), 2)
MCC_Fn_nested_table <- base::round(data.frame(Performance = c(MCC_Fn_nested)), 2)
row.names(MCC_Fn_nested_table) <- c("MatthewsCorrelationCoefficient")
table_performance_Fn_nested <- rbind(MCC_Fn_nested_table, table_performance_Fn_nested)

p_confusion_Fn_caret_nested_table <- data.frame(confusion_RF_Fn_nested$table) %>% 
  ggplot(data = ., aes(x = Prediction, y = Reference, fill =  100 * Freq/sum(Freq)))+
  theme(axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13, angle = 55, hjust = 0), title = element_text(size = 20),
        plot.margin = unit(c(1.1, 1.1, 1.1, 1.1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))+
  ylab("")+
  xlab("")+
  annotation_custom(gridExtra::tableGrob(table_performance_Fn_nested))

print(cowplot::plot_grid(p_confusion_Fn_caret_nested, p_confusion_Fn_caret_nested_table, align = "h", ncol = 2, rel_widths = c(1/2, 1/4)))

### 7.3.3. Predict unseen data ----
cephalothin_Fn <- results_Fn_RF_novel_compounds[c(2:4), ]
prediction_novel_compound_Fn <- stats::predict(Model_RF_PGMM_Fn_caret_RCV, newdata = cephalothin_Fn)
prediction_novel_compound_Fn_prob <- stats::predict(Model_RF_PGMM_Fn_caret_RCV, newdata = cephalothin_Fn, type = "prob")

prediction_novel_compound_Fn_exp <- cbind.data.frame(Prediction = prediction_novel_compound_Fn, Reference = cephalothin_Fn$MOA, Probabiltiy_CWS = prediction_novel_compound_Fn_prob$CellWallSynthesis, Probabiltiy_MD = prediction_novel_compound_Fn_prob$MembraneDisruption)
write.csv2(file = "RF_GMM_Fn_pred_Ceph.csv", prediction_novel_compound_Fn_exp) # Export results as csv file

### 7.3.4. Random forest using raw fcs data ----
# Gather necessary information for training of the algorithm
Sample_Info_Fn <- metadata_Fn %>% rename(name = Filename) # Necessary to make use of the RandomF_FCS function

# Select the fcs files based on which the model will be trained (cephalothin will be excluded)
fcs_names_Fn <- c("20210902_Fabian_MOA_Fn_Experiment_CellWallSynthesis_100_SGPI_C1.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_CellWallSynthesis_100_SGPI_C2.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_CellWallSynthesis_100_SGPI_C3.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_CellWallSynthesis_100_SGPI_C4.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_CellWallSynthesis_100_SGPI_C5.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_CellWallSynthesis_100_SGPI_C6.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_CellWallSynthesis_100_SGPI_C7.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_CellWallSynthesis_100_SGPI_C8.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_CellWallSynthesis_100_SGPI_C9.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_Control_100_SGPI_G10.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_Control_100_SGPI_G11.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_Control_100_SGPI_G12.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_Control_100_SGPI_G7.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_Control_100_SGPI_G8.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_Control_100_SGPI_G9.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNAReplication_100_SGPI_D1.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNAReplication_100_SGPI_D2.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNAReplication_100_SGPI_D3.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNAReplication_100_SGPI_D4.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNAReplication_100_SGPI_D5.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNAReplication_100_SGPI_D6.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNAReplication_100_SGPI_D7.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNAReplication_100_SGPI_D8.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNAReplication_100_SGPI_D9.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNATranscription_100_SGPI_D10.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNATranscription_100_SGPI_D11.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNATranscription_100_SGPI_D12.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNATranscription_100_SGPI_E1.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNATranscription_100_SGPI_E2.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_DNATranscription_100_SGPI_E3.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_FolicAcidMetabolism_100_SGPI_G1.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_FolicAcidMetabolism_100_SGPI_G2.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_FolicAcidMetabolism_100_SGPI_G3.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_Heat_100_SGPI_G4.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_Heat_100_SGPI_G5.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_Heat_100_SGPI_G6.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_MembraneDisruption_100_SGPI_B1.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_MembraneDisruption_100_SGPI_B10.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_MembraneDisruption_100_SGPI_B11.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_MembraneDisruption_100_SGPI_B12.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_MembraneDisruption_100_SGPI_B2.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_MembraneDisruption_100_SGPI_B3.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_MembraneDisruption_100_SGPI_B4.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_MembraneDisruption_100_SGPI_B5.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_MembraneDisruption_100_SGPI_B6.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_MembraneDisruption_100_SGPI_B7.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_MembraneDisruption_100_SGPI_B8.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_MembraneDisruption_100_SGPI_B9.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F1.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F2.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F3.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F4.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F5.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F6.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F7.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F8.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis30SInhibition_100_SGPI_F9.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E10.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E11.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E12.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E4.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E5.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E6.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E7.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E8.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesis50SInhibition_100_SGPI_E9.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesistRNAInterference_100_SGPI_F10.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesistRNAInterference_100_SGPI_F11.fcs",
                  "20210902_Fabian_MOA_Fn_Experiment_ProteinSynthesistRNAInterference_100_SGPI_F12.fcs")

Sample_Info_Fn1 <- Sample_Info_Fn %>% dplyr::filter(name %in% fcs_names_Fn)
Model_RF_Fn1_MCC <- RandomF_FCS_MCC(flowData_transformed_Fn_gated[fcs_names_Fn], sample_info = Sample_Info_Fn1, target_label = "MOA", downsample = 9775, classification_type = "sample", param = paramRF , p_train = 0.8, seed = seed, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = FALSE, metric = "MCC")
saveRDS(object = Model_RF_Fn1_MCC, file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/models_rds/Fn_24h_8param_9775cells_MCC_novel_compounds.rds")
#Model_RF_Fn1_MCC <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/MOA_Antibiotics/MOA_Antibiotics_Project/models_rds/Fn_24h_8param_9775cells_MCC_novel_compounds.rds")

mytable_FCS_MCC_Fn <- base::round(data.frame(Performance = Model_RF_Fn1_MCC[[2]]$overall), 2)

confusion_RF_Fn_FCS_MCC <- data.frame(Model_RF_Fn1_MCC[[2]]$table)
confusion_RF_Fn_FCS_MCC_prep <- data.frame(lapply(confusion_RF_Fn_FCS_MCC, rep, confusion_RF_Fn_FCS_MCC$Freq))
MCC_FCS_Fn <- mltools::mcc(preds = as.factor(confusion_RF_Fn_FCS_MCC_prep$Prediction), actuals = as.factor(confusion_RF_Fn_FCS_MCC_prep$Reference))
MCC_FCS_Fn_table <- base::round(data.frame(Performance = MCC_FCS_Fn), 2)
row.names(MCC_FCS_Fn_table) <- c("MatthewsCorrelationCoefficient")
mytable_FCS_MCC_Fn_2 <- rbind(MCC_FCS_Fn_table, mytable_FCS_MCC_Fn)

p_conf_Fn_FCS_MCC <- ggplot2::ggplot(data.frame(Model_RF_Fn1_MCC[[2]]$table), 
                                     ggplot2::aes(x = Prediction, y = Reference, fill = 100 * Freq/sum(Freq))) +
  ggplot2::geom_raster() + 
  ggplot2::geom_text(ggplot2::aes(label = round(100 * Freq/sum(Freq), 0)), size = 6) +
  ggplot2::scale_fill_distiller(name = "% of total samples\n classified\n") + 
  ggplot2::theme_bw() +
  ggplot2::scale_x_discrete(position = "top") +
  ggplot2::theme(axis.title = element_text(size = 16),
                 strip.text.x = element_text(size = 14),
                 legend.title = element_text(size = 14),
                 legend.text = element_text(size = 14),
                 axis.text.y = element_text(size = 13),
                 axis.text.x = element_text(size = 13, angle = 55, hjust = 0),
                 title = element_text(size = 20),
                 plot.margin = unit(c(1.1, 1.1, 1.1, 1.1), "cm"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_rect(fill = "transparent", colour = NA),
                 plot.background = element_rect(fill = "transparent", colour = NA))

p_conf_table_Fn_FCS_MCC <- ggplot2::ggplot(data.frame(Model_RF_Fn1_MCC[[2]]$table), 
                                           ggplot2::aes(x = Prediction, y = Reference, fill = 100 * Freq/sum(Freq))) +
  theme(axis.title = ggplot2::element_text(size = 16), 
        strip.text.x = ggplot2::element_text(size = 14), 
        legend.title = ggplot2::element_text(size = 14), 
        legend.text = ggplot2::element_text(size = 14), 
        axis.text.y = ggplot2::element_text(size = 13), 
        axis.text.x = ggplot2::element_text(size = 13, angle = 55, hjust = 0),
        title = ggplot2::element_text(size = 20), 
        plot.margin = ggplot2::unit(c(1.1, 1.1, 1.1, 1.1), "cm"),
        panel.grid.major = ggplot2::element_blank(), 
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(), 
        panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
        plot.background = ggplot2::element_rect(fill = "transparent", colour = NA)) +
  ggplot2::ylab("") +
  ggplot2::xlab("") + 
  ggplot2::annotation_custom(gridExtra::tableGrob(mytable_FCS_MCC_Fn_2))

print(cowplot::plot_grid(p_conf_Fn_FCS_MCC, p_conf_table_Fn_FCS_MCC, align = "h", 
                         ncol = 2, rel_widths = c(1/2, 1/4)))