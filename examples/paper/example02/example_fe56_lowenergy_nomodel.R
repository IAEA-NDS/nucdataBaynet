#
#  Example evaluation of neutron-induced cross sections of Fe56
#  between 1 and 2 MeV using a Bayesian network
#

library(data.table)
library(Matrix)
library(ggplot2)
library(plotly)
library(nucdataBaynet)

source("data_preparation.R")
source("nodes_definition.R")
source("mappings_definition.R")
source("inference.R")

# call plots.R interactively to produce the visualizations
