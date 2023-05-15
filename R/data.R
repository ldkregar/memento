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

#' Simulated object
#'
#' Obtained using Nick's run_sim_infer_example on KX001 (tree downsampled to 30 tips)
#'
#' @format ## `sim_object`
#' A list of items
#' \describe{
#'   \item{details}{data.frame}
#'   \item{G}{matrix}
#'   \item{tree}{ape "phylo" object}
#'   \item{meta}{data.frame with sample meta data for gamma and delta}
#'   \item{MF}{matrix of methylated reads on the forward strand}
#'   \item{MR}{matrix of methylated reads on the reverse strand}
#'   \item{UR}{matrix of unmethylated reads on the reverse strand}
#'   \item{UF}{matrix of unmethylated reads on the forward strand}
#' }
"sim_object"
