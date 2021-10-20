## nucdataBaynet: Bayesian networks for nuclear data evaluation - v0.2.0

This R package enables the creation of and inference in Bayesian networks with multivariate normal conditional distributions and linear and non-linear functional relationships between different nodes.
It was designed with the specifics of nuclear data evaluation in mind but can be used in other application scenarios that are compatible with the modeling assumptions mentioned in the previous sentence. 

*Note: The development of this package is in an early stage and the interface of the functions cannot be expected to be stable yet and documentation needs to be extended as well.*

### Installation

The packages *Matrix*, *data.table*, *numDeriv* and *mathjaxr* are prerequisites and the packages *igraph* and *ggplot2* are recommended auxiliary packages. These packages are available on the CRAN network and can be installed by

    install.packages(c("Matrix", "data.table", "numDeriv", "mathjaxr"))
    install.packages(c("igraph", "ggplot2"))

The *nucdataBaynet* package can be installed from the command line, e.g., by

    git clone https://github.com/IAEA-NDS/nucdataBaynet.git
    R CMD INSTALL nucdataBaynet

### General workflow

The workflow with this package can be divided into the following steps:

1.  Define a data table with the information about the nodes including prior estimates and uncertainties (or a prior covariance matrix) and the values of observed nodes
2.  Define the individual mappings between the nodes and combine them in a so-called compound mapping
3.  Use the customized Levenberg-Marquardt algorithm which takes as argument the compound mapping object to obtain a *Maximum a posterior* (MAP) estimate of the values of the nodes
4.  After having obtained the posterior estimate, use additional functions to compute blocks of the approximate posterior covariance matrix or get a sample from the approximate posterior distribution

### Getting started

The `examples/` folder contains tutorials and Bayesian network examples that were given in the paper *Nuclear data evaluation with Bayesian networks* on *arxiv*. The tutorial implementing a simple linear regression with Bayesian networks is a good starting point to learn how the general workflow is implemented in practice.
