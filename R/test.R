# Delimiters used in PED files
delims <- "[ \t]"

initialize <- function(.Object, path, n = NULL, p = NULL) {
  path <- path.expand(path)
  if (!file.exists(path)) {
    # Try to add extension (common in PLINK)
    path <- paste0(path, ".bed")
    if (!file.exists(path)) {
      stop("File not found.")
    }
  }
  dir <- substr(path, 1L, nchar(path) - 4L)
  # Create Rcpp object
  .Object@xptr <- .Call("pcaMatrix__new", path, n, p)
  .Object@path <- path
  .Object@dims <- c(n, p)
  return(.Object)
}


#' @aliases pcaMatrix-class
#' @export pcaMatrix
#' @exportClass pcaMatrix
pcaMatrix <- setClass("pcaMatrix", slots = c(xptr = "externalptr", dims = "integer", path = "character"))

#' @export
setMethod("initialize", signature(.Object = "pcaMatrix"), initialize)

#' #' @export
#' get_genotype <- function(x, i, j) {
#'   .Call("pcaMatrix__get_genotype", x@xptr, i, j)
#' }

#' @export
extract_vector <- function(x, i) {
  .Call("pcaMatrix__extract_vector", x@xptr, i)
}

#' @export
extract_matrix <- function(x, i, j) {
  subset <- .Call("pcaMatrix__extract_matrix", x@xptr, i, j)
  return(subset)
}

#' @export
prodvec <- function(x, vec, nIND, nSNP) {
  out <- .Call("pcaMatrix__prodvect", x@xptr, vec, nIND, nSNP)
  return(out)
}

#' @export
crossprodvec <- function(x, vec, nIND, nSNP) {
  out <- .Call("pcaMatrix__crossprodvect", x@xptr, vec, nIND, nSNP)
  return(out)
}

#' @export
# single core implementation
ma_rsvd <- function(X, k, nIND, nSNP) {
  it <- 0
  A <- function(x, args) {
    prodvec(X, x, nIND, nSNP)
  }
  Atrans <- function(x, args) {
    crossprodvec(X, x, nIND, nSNP)
  }
  res <- RSpectra::svds(A, k, nu = k, nv = k, Atrans = Atrans, dim = c(nIND, nSNP))
  
  return(res)
}

