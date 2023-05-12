#' Random tree
#'
#' set.seed(666)
#' tree <- ape::rtree(10)
#' tree$tip.label <- paste0(LETTERS[1:10], c(1:10))
#'
#' @format ## `tree`
#' An object of class "phylo"
#' \describe{
#'   \item{edge}{integer matrix where column one holds parent node and column
#'   two holds child node}
#'   \item{tip.label}{character vector of tip labels}
#'   \item{Nnode}{integer number of internal nodes}
#'   \item{edge.length}{numeric vector of edge lengths, ordered according to
#'   edge matrix}
#' }
"tree"
