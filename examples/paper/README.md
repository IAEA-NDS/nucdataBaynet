# Examples given in the paper

The folders contain the examples given in the paper
*Nuclear data evaluation with Bayesian networks* available on *arxiv*.

The workflow to create the Bayesian networks and do inference in them always
follows the same steps:

1.  Define a datatable with the information about the nodes including prior estimates and uncertainties

2.  Define the mappings, i.e., functional relationships, between the nodes

3.  Perform inference in the Bayesian network to obtain posterior estimates, uncertainties and covariances of variables of interest

These steps are split into separate files: Step one is implemented in `nodes_definition.R`,
step two in `mappings_definition.R` and step three in `inference.R`.

To run the complete example, one needs to change into the directory of the example, e.g., using
the instruction `setwd(<example_dir>)` and then execute the script prefixed by `example_`.

Afterwards, the Bayesian network and results can be visualized by running the code in `plots.R`.
