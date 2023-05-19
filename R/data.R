#' Simulated object
#'
#' Obtained using Nick's run_sim_infer_example on KX001 (tree downsampled to 30 tips, rate = 0.001)
#'
#' @format ## `object`
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
"object"
