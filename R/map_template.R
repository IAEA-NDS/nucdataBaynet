create_maptype_map <- function() {

  setup <- function(params) {
    # TODO
  }


  getType <- function() {
    # TODO
  }


  getName <- function() {
    # TODO
  }


  getDescription <- function() {
    # TODO
  }


  get_src_idx <- function() {
    # TODO
  }


  get_tar_idx <- function() {
    # TODO
  }


  propagate <- function(x, with.id=TRUE) {
    # TODO
  }


  jacobian <- function(x, with.id=TRUE) {
    # TODO
  }


  return(list(
    setup = setup,
    getType = getType,
    getName = getName,
    getDescription = getDescription,
    get_src_idx = get_src_idx,
    get_tar_idx = get_tar_idx,
    propagate = propagate,
    jacobian = jacobian
  ))
}
