the <- new.env(parent = emptyenv())
the$number_of_states <- 3
the$snv_rate_per_year <- 20
the$num_bases_in_genome <- 3e9
the$human.snv.rate.per.base<- the$snv_rate_per_year/the$num_bases_in_genome
the$epi_fwd = 1e5
the$epi_bck = 1e5
the$all_possible_states <- c(1:the$number_of_states)
the$state_transition_mtx <- matrix(1:the$number_of_states**2,
                                   the$number_of_states,
                                   the$number_of_states,
                                   byrow = TRUE)
colnames(the$state_transition_mtx) <- rownames(the$state_transition_mtx) <- as.character(1:the$number_of_states)

the$state_transition_vec_out <- rep(c(1:the$number_of_states), the$number_of_states)
the$tree_branch_is_snv <- FALSE
the$full_tree <- NULL
the$global_edges_vector_sort <- NULL
the$global_tips_vector <- NULL
the$subtree_sub_nodes <- NULL
the$global_subtree_priors_affected_tips <- NULL
the$top_per_site_name <- NULL
the$list_of_subtree_sol_name_per_ss <- NULL

rebase <- function(b, base = 10){
  base <- as.integer(base)

  znj <- lapply(b, function(x){
    n   <- ceiling(log(x, base))
    vec <- numeric()
    val <- x

    while(n >= 0){
      rem <- val %/% base^n
      val <- val - rem * base^n
      vec <- c(vec, rem)
      n <- n - 1
    }

    while(vec[1] == 0 & length(vec) > 1) vec <- vec[-1]
    return(vec)
  })
  return(unlist(znj))
}

add_padding <- function(v, new_length){
  while(length(v) < new_length){
    v <- c(0, v)
  }
  return(v)
}

get_change_combinations_4s <- function(n_changes_, ns = NULL){
  if(is.null(ns)){
    ns <- get("number_of_states", envir = the)
  }
  if(ns == 2 & n_changes_ != 0){
    change_combination <- matrix(rep(1, n_changes_), ncol = n_changes_)
    names_of_changes <- paste0(rep("A", n_changes_), collapse = "")
  }else{
    possible_perm <- 1:((ns-1)**n_changes_)-1
    change_combination <- sapply(possible_perm, function(x){ rebase(x, (ns -1))} )
    last_item <- change_combination[[length(change_combination)]]
    change_combination2 <- lapply(change_combination, add_padding, new_length = length(last_item))
    change_combination <- do.call("rbind", change_combination2)
    names_of_changes <- apply(change_combination, 1, function(x){
      y <- as.character(x)
      y[y == "0"] <- "A"
      y[y == "1"] <- "B"
      y[y == "2"] <- "C"
      y[y == "3"] <- "D"
      return(paste0(y, collapse = ""))
    })
    change_combination <- change_combination + 1 # this way you get indices
  }
  return(list("change_combinations" = change_combination,
              "names_of_changes" = names_of_changes))
}

the$change_combs_list <- lapply(c(0:10), function(nc){
  cc <- get_change_combinations_4s(nc)
  return(cc)
})
names(the$change_combs_list) <- as.character(0:(length(the$change_combs_list)-1))

getRateDf <- function(numStates){
  do.call("rbind",
          lapply(0:(numStates-1),
                 function(s) data.frame(s.parent=s,s.child=setdiff(0:(numStates-1),s))
          ))
}

the$df_rate <- getRateDf(the$number_of_states)
the$df_rate$rate <- rep(0.001, nrow(the$df_rate))

getQMatrix <- function(df){
  N=max(df$s.parent)+1
  out=matrix(NA,ncol=N,nrow=N)
  out[N*df$s.child+(df$s.parent+1)]=df$rate
  diag(out)=-rowSums(out,na.rm = TRUE)
  out
}

the$Q <- getQMatrix(the$df_rate)
the$useQ <- FALSE
the$use_ff <- FALSE
the$use_Rfast <- FALSE
the$use_rgpumult <- TRUE

#' Report my default settings
#'
#' @export
memento_settings <- function() {
  return(list("number_of_states" = the$number_of_states,
              "epi_fwd" = the$epi_fwd,
              "epi_bck" = the$epi_bck,
              "all_possible_states" = the$all_possible_states,
              "state_transition_mtx" = the$state_transition_mtx,
              "state_transition_vec_out" = the$state_transition_vec_out,
              "length_change_comb_list" = length(the$change_combs_list),
              "df_rate" = the$df_rate,
              "Q" = the$Q,
              "useQ" = the$useQ,
              "use_ff" = the$use_ff,
              "tree_branch_is_snv" = the$tree_branch_is_snv,
              "snv_rate_per_year" = the$snv_rate_per_year,
              "num_bases_in_genome" = the$num_bases_in_genome,
              "human.snv.rate.per.base" = the$human.snv.rate.per.base,
              "use_Rfast" = the$use_Rfast,
              "use_rgpumult" = the$use_rgpumult,
              "full_tree_is_null" = is.null(the$full_tree),
              "global_edges_vector_sort_is_null" = is.null(the$global_edges_vector_sort),
              "subtree_sub_nodes_is_null" = is.null(the$global_tips_vector),
              "subtree_sub_nodes_is_null" = is.null(the$subtree_sub_nodes),
              "global_subtree_priors_affected_tips_is_null" = is.null(the$global_subtree_priors_affected_tips),
              "top_per_site_name_is_null" = is.null(the$top_per_site_name),
              "list_of_subtree_sol_name_per_ss_is_null" = is.null(the$list_of_subtree_sol_name_per_ss)))
}

## @param human.snv.rate.per.base human snv rate per base (usually 20/3e9)

#' Set memento settings
#'
#' @param new_number_of_states integer value for number of states: 3 or 4 (requires phased data)
#' @param new_epi_fwd numeric value for forward rate (gain of methylation): scales the human SNV rate, e.g. 1e5 times greater than SNV rate per base
#' @param new_epi_bck numeric value for backward rate (loss of methylation): scales the human SNV rate
#' @param new_rates data.frame of new rates
#' @param new_Q new Q matrix
#' @param new_useQ logical whether to use Q matrix
#' @param new_use_ff logical whether to use ff for matrix multiplication instead of cpp
#' @param new_combs_list supply a new combinations list manually
#' @param new_tree_branch_is_snv logical whether tree branch lengths are in SNVs
#' @param new_snv_rate_per_year numeric number of SNVs per year (usually 20)
#' @param new_num_bases_in_genome numeric number of bases in genome (3e9 for human)
#' @param new_use_Rfast logical whether to use Rfast for matrix multiplication
#'
#' @export
set_memento_settings <- function(new_number_of_states = NULL,
                                 #human.snv.rate.per.base = NULL,
                                 new_epi_fwd = NULL,
                                 new_epi_bck = NULL,
                                 new_rates = NULL,
                                 new_Q = NULL,
                                 new_useQ = NULL,
                                 new_use_ff = NULL,
                                 new_combs_list = NULL,
                                 new_tree_branch_is_snv = NULL,
                                 new_snv_rate_per_year = NULL,
                                 new_num_bases_in_genome = NULL,
                                 new_use_Rfast = NULL,
                                 new_use_rgpumult = NULL){
  old_num_states <- the$number_of_states
  if(!is.null(new_number_of_states)){
    the$number_of_states <- new_number_of_states
    the$all_possible_states <- c(1:the$number_of_states)
    the$state_transition_mtx <- matrix(1:the$number_of_states**2,
                                       the$number_of_states,
                                       the$number_of_states,
                                       byrow = TRUE)
    colnames(the$state_transition_mtx) <- rownames(the$state_transition_mtx) <- as.character(1:the$number_of_states)

    the$state_transition_vec_out <- rep(c(1:the$number_of_states), the$number_of_states)

    the$change_combs_list <- lapply(c(0:10), function(nc){
      cc <- get_change_combinations_4s(nc, ns = new_number_of_states)
      return(cc)
    })
    names(the$change_combs_list) <- as.character(0:(length(the$change_combs_list)-1))

    the$df_rate <- getRateDf(the$number_of_states)
    the$df_rate$rate <- rep(0.001, nrow(the$df_rate))
    the$Q <- getQMatrix(the$df_rate)
  }

  if(!is.null(new_epi_fwd)){
    the$epi_fwd <- new_epi_fwd
  }
  if(!is.null(new_epi_bck)){
    the$epi_bck <- new_epi_bck
  }
  if(!is.null(new_rates)){
    the$df_rate <- new_rates
    the$Q <- getQMatrix(the$df_rate)
  }
  if(!is.null(new_useQ)){
    the$useQ <- new_useQ
  }
  if(!is.null(new_Q)){
    the$df_rate <- NULL
    the$Q <- new_Q
    cat("New Q matrix, but df_rate is NULL")
  }
  if(!is.null(new_use_ff)){
    if(new_use_ff){ # if TRUE
      if (rlang::is_installed("HDF5Array") &
          rlang::is_installed("S4Arrays") &
          rlang::is_installed("DelayedArray")) {
        the$use_ff <- new_use_ff
      } else {
        message("To use ff functionalities, you must install HDF5Array, S4Arrays and DelayedArray")
      }
    }else{ # if FALSE - do not require additional libraries
      the$use_ff <- new_use_ff
    }
  }
  if(!is.null(new_combs_list)){
    the$change_combs_list <- new_combs_list
  }
  if(!is.null(new_tree_branch_is_snv)){
    the$tree_branch_is_snv <- new_tree_branch_is_snv
  }
  if(!is.null(new_snv_rate_per_year)){
    the$snv_rate_per_year <- new_snv_rate_per_year
    the$human.snv.rate.per.base<- the$snv_rate_per_year/the$num_bases_in_genome
  }
  if(!is.null(new_num_bases_in_genome)){
    the$num_bases_in_genome <- new_num_bases_in_genome
    the$human.snv.rate.per.base<- the$snv_rate_per_year/the$num_bases_in_genome
  }
  if(!is.null(new_use_Rfast)){
    the$use_Rfast <- new_use_Rfast
  }
  if(!is.null(new_use_rgpumult)){
    the$use_rgpumult <- new_use_rgpumult
  }
  invisible(old_num_states)
}

#' Set memento settings to default
#'
#' @export
#'
#' @examples
#' set_memento_settings_to_default()
#' memento_settings()
set_memento_settings_to_default <- function(){
  the$number_of_states <- 3
  the$epi_fwd = 1e5
  the$epi_bck = 1e5
  the$snv_rate_per_year <- 20
  the$num_bases_in_genome <- 3e9
  the$human.snv.rate.per.base<- the$snv_rate_per_year/the$num_bases_in_genome
  the$all_possible_states <- c(1:the$number_of_states)
  the$state_transition_mtx <- matrix(1:the$number_of_states**2,
                                     the$number_of_states,
                                     the$number_of_states,
                                     byrow = TRUE)
  colnames(the$state_transition_mtx) <- rownames(the$state_transition_mtx) <- as.character(1:the$number_of_states)

  the$state_transition_vec_out <- rep(c(1:the$number_of_states), the$number_of_states)

  the$change_combs_list <- lapply(c(0:10), function(nc){
    cc <- get_change_combinations_4s(nc)
    return(cc)
  })
  names(the$change_combs_list) <- as.character(0:(length(the$change_combs_list)-1))

  the$df_rate <- getRateDf(the$number_of_states)
  the$df_rate$rate <- rep(0.001, nrow(the$df_rate))

  the$Q <- getQMatrix(the$df_rate)
  the$useQ <- FALSE
  the$use_ff <- FALSE
  the$tree_branch_is_snv <- FALSE
  the$use_Rfast <- FALSE

  the$full_tree <- NULL
  the$global_edges_vector_sort <- NULL
  the$global_tips_vector <- NULL
  the$subtree_sub_nodes <- NULL
  the$global_subtree_priors_affected_tips <- NULL
  the$top_per_site_name <- NULL
  the$list_of_subtree_sol_name_per_ss <- NULL
}

ignore_unused_imports <- function() {
  HDF5Array::H5ADMatrix
}

#' Wipe global variables set by memento
#'
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' hi <- wipe_before_EM()
#' }
wipe_before_EM <- function(){
  the$full_tree <- NULL
  the$global_edges_vector_sort <- NULL
  the$global_tips_vector <- NULL
  the$subtree_sub_nodes <- NULL
  the$global_subtree_priors_affected_tips <- NULL
  the$top_per_site_name <- NULL
  the$list_of_subtree_sol_name_per_ss <- NULL
  invisible()
}
