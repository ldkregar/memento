P00 = function(t, epi_mut_rate_bck, epi_mut_rate_fwd){
  epi_mut_rate_tot = epi_mut_rate_bck + epi_mut_rate_fwd
  (epi_mut_rate_bck + epi_mut_rate_fwd*exp(-epi_mut_rate_tot*t))/epi_mut_rate_tot
}

P10 = function(t, epi_mut_rate_bck, epi_mut_rate_fwd){
  epi_mut_rate_tot=epi_mut_rate_bck+epi_mut_rate_fwd
  p11=(epi_mut_rate_fwd+epi_mut_rate_bck*exp(-epi_mut_rate_tot*t))/epi_mut_rate_tot
  1-p11
}

## 4 State independent alleles model.
#Gets the transition matrix where row i and column j represent
# the probability of being in state j having started at state i a time
# t units earlier.
getP4 = function(t, epi_mut_rate_bck, epi_mut_rate_fwd){
  P = getP2(t, epi_mut_rate_bck, epi_mut_rate_fwd)
  out = matrixcalc::direct.prod(P, P)
  colnames(out) = c("00", "01", "10", "11")
  rownames(out) = colnames(out)
  out
}

# getP4_new = function(t, Q){
#   P = getP(t, Q)
#   out = matrixcalc::direct.prod(P, P)
#   colnames(out) = c("00", "01", "10", "11")
#   rownames(out) = colnames(out)
#   out
# }

getP3=function(t, epi_mut_rate_bck, epi_mut_rate_fwd){
  out = matrix(NA, ncol = 3, nrow = 3)
  p00 = P00(t, epi_mut_rate_bck, epi_mut_rate_fwd)
  p01 = 1 - p00
  p10 = P10(t, epi_mut_rate_bck, epi_mut_rate_fwd)
  p11 = 1 - p10

  out[1, 1:3] = c(p00*p00, 2*p01*p00, p01*p01) # 0_0, 0_0.5, 0_1
  out[2, 1:3] = c(p10*p00, # 0.5 to 0
                p11*p00 + p10*p01, # 0.5_0.5
                p11*p01) # 0.5 to 1
  out[3, 1:3] = c(p10*p10, 2*p10*p11, p11*p11)
  out
}

getP2 = function(t, epi_mut_rate_bck, epi_mut_rate_fwd){
  out = matrix(NA, ncol=2, nrow=2)
  p00 = P00(t, epi_mut_rate_bck, epi_mut_rate_fwd)
  p01 = 1 - p00
  p10 = P10(t, epi_mut_rate_bck, epi_mut_rate_fwd)
  p11 = 1 - p10
  matrix(c(p00, p10, p01, p11), ncol = 2, nrow = 2)
}

getP <- function(t,Q) as.matrix(expm::expm(t*Q))

even <- function(x) base::subset(x, x %% 2 == 0)

odd <- function(x) base::subset(x, x %% 2 == 1)

find_ancestor_edge <- function(tree, edge){
  rt = ape::Ntip(tree)+1
  anc = tree$edge[edge,1]
  if(anc == rt){
    return(-1)
  }else{
    return(which(tree$edge[,2] == anc))
  }
}

add_mtx_to_change_comb_file <- function(n_changes_){
  additional_matrix <- get_change_combinations_4s(n_changes_)
  the$change_combs_list[[as.character(n_changes_)]] <- additional_matrix
  # save(change_combs_list, file = change_comb_file)
  return(additional_matrix)
}

get_outcomes_per_change <- function(change_combination_list_,
                                    tree_,
                                    previous_edges,
                                    changes_on_branch,
                                    mtx_of_children_local){
  # with this function, you get indices of states per branch
  change_combination <- change_combination_list_$change_combinations
  names_of_changes <- change_combination_list_$names_of_changes
  rownames(change_combination) <- names_of_changes

  change_results_new <- vector("list",
                               length = nrow(change_combination)*length(the$all_possible_states))
  names(change_results_new) <- paste0(names_of_changes, "_", rep(the$all_possible_states, each = nrow(change_combination)))

  starting_states <- rep(the$all_possible_states, each = nrow(change_combination))
  temp_combs <- expand.grid(names_of_changes, the$all_possible_states)

  for(j in 1:nrow(temp_combs)){
    # go through all changes on branches
    current_comb <- temp_combs[j,]
    branch_moves <- change_combination[current_comb[,1],]
    last_state <- current_comb[, 2]
    idx_new <- rep(the$state_transition_mtx[last_state, last_state], nrow(tree_$edge))

    for(i in changes_on_branch){
      # find previous edges, get last state - either inherit from root of previous edge
      previous_edge <- previous_edges[i]
      if(previous_edge == -1){ # means you are at root
        last_state <- current_comb[, 2]
      }else{
        last_state <- which(the$state_transition_mtx == idx_new[previous_edge], arr.ind = TRUE)[,2]
      }
      change_type_of_branch <- which(i == changes_on_branch)
      possible_new <- the$all_possible_states[the$all_possible_states != last_state]
      chosen_new <- possible_new[branch_moves[change_type_of_branch]]
      idx_new[i] <- the$state_transition_mtx[last_state, chosen_new]
      descending_branches <- mtx_of_children_local$branch_children_only[i,]
      idx_new[as.logical(descending_branches)] <- the$state_transition_mtx[chosen_new, chosen_new]
    }
    change_results_new[[j]] <- idx_new
  }
  names(change_results_new) <- paste0(temp_combs[,2], ":",
                                      sapply(change_results_new,
                                             find_state_of_affected_branches,
                                             changed_branches = changes_on_branch))
  return(change_results_new)
}

get_idx_and_sum_for_both_directions <- function(cr,
                                                priors_p_branch,
                                                return_indices = FALSE){

  elements <- cbind(c(1:nrow(priors_p_branch)), cr)
  summ <- priors_p_branch[elements]

  results <- sum(summ)

  if(return_indices == FALSE){
    return(results)
  }else{
    r_and_ind <- list("prior" = results,
                      "idx" = cr)
    return(r_and_ind)
  }

}

find_state_of_affected_branches <- function(branch_state, changed_branches){
  if(length(changed_branches) < 1){
    r <- paste0(unique(branch_state), collapse = "_")
  }else{
    r <- paste0(branch_state[changed_branches], collapse = "_")
  }
  return(r)
}

generate_priors_hemimeth <- function(changes_on_branch,
                                     tree_,
                                     # tree_ = tree,
                                     # priors_p_branch = priors_per_branch,
                                     priors_p_branch,
                                     return_idx = FALSE){
  # When a change happens, some edges will experience a fwd or back change,
  # but the rest will remain the same. Depending on change in an ancestor,
  # nodes will change their state.
  # The change_combination matrix will have as many rows as there are possible
  # permutations of half and full changes across the changed branches. This
  # matrix will be processed by row, where the columns will specify the change
  # affecting the column idx of changed_branches

  changes_on_branch <- sort(changes_on_branch)

  n_changes <- length(changes_on_branch)

  previous_edges <- sapply(1:nrow(tree_$edge), find_ancestor_edge, tree = tree_)

  if(!is.null(the$change_combs_list[[as.character(n_changes)]])){
    change_combination_list <- the$change_combs_list[[as.character(n_changes)]]
  }else{
    change_combination_list <- add_mtx_to_change_comb_file(n_changes)
  }

  mtx_of_children_local_one <- generate_tree_children_matrix(tree_)
  change_results_new <- get_outcomes_per_change(change_combination_list,
                                                tree_,
                                                previous_edges,
                                                changes_on_branch,
                                                mtx_of_children_local = mtx_of_children_local_one)

  # submit state changes on branches and get sum of priors for both directions
  get_idx_from_both_starts <- lapply(change_results_new,
                                     get_idx_and_sum_for_both_directions,
                                     priors_p_branch,
                                     return_indices = return_idx)
  return(get_idx_from_both_starts)
}

map_branch_changes_to_affected_tips_hemi <- function(pr, tree_){
  # takes a vector of branch exit states (generated when priors are produced) and prepares a vector of tip exit states
  exit_states <- vector("list", length = length(pr))
  # find terminal edges
  terminal_edges <- which(tree_$edge[,2] < tree_$edge[1,1])
  names(terminal_edges) <- tree_$tip.label

  for(j in 1:length(exit_states)){
    exit <- pr[[j]]$idx[terminal_edges]
    exit2 <- the$state_transition_vec_out[exit]
    names(exit2) <- names(terminal_edges)
    exit_states[[j]] <- exit2
  }
  names(exit_states) <- names(pr)
  return(exit_states)
}

generate_genotypes <- function(a_t,
                               #original_tree = full_tree,
                               original_tree,
                               # root_tree = one_tree
                               root_tree){
  a_t_names <- names(a_t)
  starting_states <- as.numeric(sapply(strsplit(a_t_names, ":"), "[[", 1))

  desc_per_changes <- mapply(a_t, starting_states, FUN = function(tip_or_node, ss){
    nonzero_tips <- which(tip_or_node != 0)

    if(!is.null(root_tree)){
      nodes_of_tree <- root_tree$tip.label[nonzero_tips]
    }else{
      nodes_of_tree <- original_tree$tip.label[nonzero_tips]
    }
    nonzero_vals <- tip_or_node[nonzero_tips]

    desc <- as.list(match(nodes_of_tree, original_tree$tip.label))
    desc[is.na(desc)] <- lapply(nodes_of_tree[is.na(desc)], function(am_na){
      return(phangorn::Descendants(original_tree, node = as.numeric(am_na), "tips")[[1]])
    })
    u <- unlist(desc)
    full_tip_vector <- rep(0, length(original_tree$tip.label))
    if(length(u) == length(nonzero_vals)){
      full_tip_vector[u] <- nonzero_vals
    }else{
      repeated_vals <- rep(nonzero_vals, lengths(desc))
      full_tip_vector[u] <- repeated_vals
    }
    return(full_tip_vector)
  })

  return(desc_per_changes)
}

## @param tree_branches_are_snv logical: whether tree branches in SNVs or years

#' Generate priors for each tree
#'
#' @param all_trees list of "phylo" objects
#' @param epi_fwd numeric value for forward rate (gain of methylation)
#' @param epi_bck numeric value for backward rate (loss of methylation)
#' @param m integer value for number of spontaneous changes on a tree
#' @param full_tree "phylo" object - entire tree
#' @param verbose logical
#' @param Q_mtx Q matrix of transitions
#'
#' @return List with elements priors, affected_tips, affected_branches,
#' tree_for_priors and genotype_vectors
#' @export
#'
#' @examples
#' \dontrun{
#'  number_of_states <- 3
#'  generate_priors(all_trees = list("root" = tree), m = 2)
#' }
generate_priors <- function(full_tree,
                            all_trees = NULL,
                            epi_fwd = NULL,
                            epi_bck = NULL,
                            Q_mtx = NULL,
                            m = 2,
                            #tree_branches_are_snv = TRUE,
                            verbose = TRUE){
  if(is.null(all_trees)){
    all_trees <- list("root" = full_tree)
  }

  if(is.null(full_tree)) stop("Must provide full tree")
  if(!all(grepl("[[:alpha:]]", full_tree$tip.label))){
    stop("Not all tip labels contain a character, rename tip labels so that they contain a character.")
  }
  priors_per_subtree <- vector("list", length = length(all_trees))

  names(priors_per_subtree) <- names(all_trees)

  for(i in 1:length(all_trees)){
    if(verbose) cat(sprintf("Computing priors for tree %s (%i of %i)\n", names(all_trees)[i], i, length(all_trees)))
    one_tree <- all_trees[[i]]
    n_branches <- nrow(one_tree$edge)
    current_m <- min(m, n_branches)
    permuts <- generate_permutations(n_branches, m = current_m) # possible permutations of changes given tree
    stopifnot(all(unname(table(rowSums(permuts))) == sapply(0:current_m, choose, n = n_branches)))
    totsnvs=sum(one_tree$edge.length)
    if(the$tree_branch_is_snv){
      tottime=totsnvs/the$snv_rate_per_year
    }else{
      tottime=totsnvs
    }
    weighted_branch_lenghts <- tottime* one_tree$edge.length/sum(one_tree$edge.length)

    ##  rates expressed as per cpg per year
    ##  Mutation rates are human mutation rate = 20/3e9 per base
    if(is.null(epi_fwd) & is.null(epi_bck)){
      epi_mut_rate_fwd_ <- the$epi_fwd*the$human.snv.rate.per.base #human.epi.rate.per.cpg
      epi_mut_rate_bck_ <- the$epi_bck*the$human.snv.rate.per.base #human.epi.rate.per.cpg
    }else{
      cat("Using custom forward and backward transition rates.\n")
      epi_mut_rate_fwd_ <- epi_fwd*the$human.snv.rate.per.base #human.epi.rate.per.cpg
      epi_mut_rate_bck_ <- epi_bck*the$human.snv.rate.per.base #human.epi.rate.per.cpg
    }

    if(is.null(Q_mtx)){
      Q_mtx_ <- the$Q
    }else{
      cat("Using custom Q matrix.\n")
      Q_mtx_ <- Q_mtx
    }

    if(the$useQ){
      p_list <- lapply(weighted_branch_lenghts,
                       getP,
                       Q = Q_mtx_)
    }else{
      if(the$number_of_states == 2){
        p_list <- lapply(weighted_branch_lenghts,
                         getP2,
                         epi_mut_rate_bck = epi_mut_rate_bck_,
                         epi_mut_rate_fwd = epi_mut_rate_fwd_)
      }else if(the$number_of_states == 3){
        p_list <- lapply(weighted_branch_lenghts,
                         getP3,
                         epi_mut_rate_bck = epi_mut_rate_bck_,
                         epi_mut_rate_fwd = epi_mut_rate_fwd_)
      }else if(the$number_of_states == 4){
        p_list <- lapply(weighted_branch_lenghts,
                         getP4,
                         epi_mut_rate_bck = epi_mut_rate_bck_,
                         epi_mut_rate_fwd = epi_mut_rate_fwd_)
      }
    }
    p_vectors <- lapply(p_list, function(x){
      y <- as.vector(t(x))
      return(y)
    })
    priors_per_branch <- log(do.call("rbind", p_vectors))

    affected_branches <- apply(permuts, 1, function(x) which(x == 1))

    affected_branches_strings <- sapply(affected_branches, function(x){
      if(length(x) > 0){
        return(paste0(x, collapse = "_"))
      }else{
        return("n")
      }
    })

    names(affected_branches) <- affected_branches_strings

    priors <- tryCatch(lapply(affected_branches,
                              generate_priors_hemimeth,
                              tree_ = one_tree,
                              priors_p_branch = priors_per_branch,
                              return_idx = TRUE),
                       warning = function(w) {
                         cat(paste0("warning!! : ", i))
                       })

    one_tree_aslist <- lapply(1:length(priors), function(x) return(one_tree))

    affected_tips <- mapply(map_branch_changes_to_affected_tips_hemi,
                            pr = priors,
                            tree_ = one_tree_aslist, SIMPLIFY = FALSE)
    if(is.list(affected_tips)){
      genotype_vectors <- lapply(affected_tips,
                                 generate_genotypes,
                                 original_tree = full_tree,
                                 root_tree = one_tree)
    }else if(is.matrix(affected_tips)){
      stop("Check why affected_tips is not a list?")
    }else{
      stop("Unexpected affected_tips object type")
    }

    tree_for_priors <- one_tree

    one_prior_result <- list("priors" = priors,
                             "affected_tips" = affected_tips,
                             "affected_branches" = affected_branches,
                             "tree_for_priors" = tree_for_priors,
                             "genotype_vectors" = genotype_vectors)
    priors_per_subtree[[i]] <- one_prior_result
  }
  return(priors_per_subtree)
}
