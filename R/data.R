#' Simulated object
#'
#' Obtained using Nick's run_sim_infer_example on KX001 (tree downsampled to 30 tips, rate = 0.001, depth = 20)
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

#' Simulated phased object
#'
#' Obtained using Nick's run_sim_infer_example on KX001 (tree downsampled to 30 tips, rate = 0.001, depth = 20, numstates=4)
#'
#' @format ## `object_phased`
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
#'   \item{MT_M}{matrix of methylated reads on the maternal strand}
#'   \item{MT_U}{matrix of unmethylated reads on the maternal strand}
#'   \item{PT_M}{matrix of methylated reads on the paternal strand}
#'   \item{PT_U}{matrix of unmethylated reads on the paternal strand}
#'   \item{U_M}{matrix of methylated reads of the unphased reads}
#'   \item{U_U}{matrix of unmethylated reads of the unphased reads}
#'   \item{EXTRAFIELDS}{character vector of additional matrix names}
#' }
"object_phased"
