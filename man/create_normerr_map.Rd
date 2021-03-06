% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_normerr.R
\name{create_normerr_map}
\alias{create_normerr_map}
\title{Create a mapping for absolute energy-independent normalization errors}
\usage{
create_normerr_map()
}
\value{
Returns a list of functions to operate with the mapping, see \code{\link{create_maptype_map}}.
}
\description{
Creates a map to distribute the absolute normalization errors given
at the source indices to the experimental data at the target indices.
}
\details{
The following fields are required in the parameter list to initialize the mapping:
\tabular{ll}{
\code{mapname} \tab Name of the mapping \cr
\code{maptype} \tab Must be \code{"normerr_map"} \cr
\code{src_idx} \tab Vector of source indices \cr
\code{tar_idx} \tab Vector of target indices \cr
\code{src_feat} \tab Vector with the names of the normalization errors at the
source indices. This vector must not contain duplicated strings. \cr
\code{tar_feat} \tab Vector with the names of the normalization error associated
with the variables at the target indices.
}
}
\note{
This mapping may be used in a product mapping
to define \emph{relative} normalization errors,
see \code{\link{create_product_mapping}}
}
\examples{
params <- list(
  mapname = "mymap",
  maptype = "normerr_map",
  src_idx = 1:3,
  tar_idx = 4:9,
  src_feat = c("expnorm1", "expnorm2", "expnorm3"),
  tar_feat = c(rep("expnorm1", 2), rep("expnorm2", 3),
               rep("expnorm3", 1))
)
mymap <- create_normerr_map()
mymap$setup(params)
x <- c(1, 5, 7, rep(0, 6))
mymap$propagate(x)
mymap$jacobian(x)

}
\seealso{
Other mappings: 
\code{\link{create_compound_map}()},
\code{\link{create_convolution_with_xtrafo_map}()},
\code{\link{create_derivative2nd_map}()},
\code{\link{create_derivative_map}()},
\code{\link{create_linearinterpol_map}()},
\code{\link{create_linearinterpol_with_xtrafo_map}()},
\code{\link{create_linearmap_map}()},
\code{\link{create_maptype_map}()},
\code{\link{create_map}()},
\code{\link{create_nonlinear_map}()},
\code{\link{create_product_map}()}
}
\concept{mappings}
