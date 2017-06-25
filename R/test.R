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