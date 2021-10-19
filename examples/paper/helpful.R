library(igraph)

# function to normalize coordinates
rescale_layout <- function(layout) {
  apply(layout, 2, function(x) {
    len <- max(x)-min(x)
    mid <- (max(x) + min(x))/2
    (x - mid) / (len/2)
  })
}

# get the Bayesian network structure
#grph <- get_network_structure(compmap$getMaps(), node_dt$NODE, obs)
#layout <- layout.graphopt(grph)

# plot and manually edit the Bayesian network
# the manual editing can be stopped by clicking in the
# bottom-left of the plot in the 'move node' step
edit_network <- function(grph, layout=NULL, layout.fun=layout.graphopt) {
  if (is.null(layout)) {
    layout <- layout.fun(grph)
  }
  layout <- rescale_layout(layout)
  plot.igraph(grph, layout = layout)
  while(TRUE) {
    cat('select node\n')
    node_idx <- identify(layout, n=1, plot=FALSE)
    cat('move node\n')
    newpos <- locator(n=1)
    if (newpos[1] < -0.9 & newpos[2] < -0.9) break
    layout[node_idx, 1:2] <- c(newpos$x, newpos$y)
    layout <- rescale_layout(layout)
    plot.igraph(grph, layout = layout)
  }
}
