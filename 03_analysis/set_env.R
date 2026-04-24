# Author: Lionel Butry (@lbutry)
# Purpose: Set up packages and functions

############################
# Load packages
############################

rm(list = ls())

library(ppcor)         # Partial correlations
library(tidyverse)     # Data handling
library(janitor)       # Data handling
library(furrr)         # Parallel mapping functions
library(parallel)      # Parallel computing
library(ggprism)       # Plotting
library(igraph)        # Graph theory analysis
library(neuroCombat)   # Feature-based post-hoc harmonization
library(rstatix)       # Stats
library(emmeans)       # Stats
library(skimr)         # Summary stats
library(ggh4x)         # Plotting: ggplot2 hacks (for axis margins)
library(patchwork)     # Plotting: Combine plots
library(ggsci)         # Plotting: Color palette
library(ggraph)        # Plotting: Graphs
library(brainconn)     # Plotting: Graphs in glass brain
library(broom)

############################
# Source functions
############################

source("02_scripts/03_analysis/functions_utils.R")
