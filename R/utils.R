return_edge_descendants <- function(tree){
  lapply(tree$edge[,2], function(nn) paste0(tree$tip.label[phangorn::Descendants(tree, node = nn, type = "tips")[[1]]], collapse = ":"))
}

remap_edges <- function(one_subtree_x,
                        full_tree){
  full_tree_children_by_edge <- return_edge_descendants(full_tree)
  # find children of subtree
  subtree_edge_desc <- return_edge_descendants(one_subtree_x)
  # replace the purely numeric node (nested tree) with colony names
  # find all nodes
  i_have_numeric <- sapply(subtree_edge_desc, function(x){
    splitted <- strsplit(x, ":")[[1]]
    any(grepl("[A-Za-z]", splitted) == FALSE)
  })
  while(any(i_have_numeric)){
    problematic <- which(i_have_numeric)
    problematic_node <- subtree_edge_desc[problematic]
    clean_problematic <- lapply(problematic_node, function(p){
      s <- strsplit(p, ":")[[1]]
      only_num <- grep("[A-Za-z]", s, invert = TRUE, value = TRUE)
      qwe <- lapply(only_num, function(on){
        full_tree$tip.label[phangorn::Descendants(full_tree, node = as.numeric(on), "tips")[[1]]]
      })
      pretty <- s
      pretty[match(only_num, pretty)] <- qwe
      return(paste0(unlist(pretty), collapse = ":"))
    })

    subtree_edge_desc[problematic] <- clean_problematic

    i_have_numeric <- sapply(subtree_edge_desc, function(x){
      splitted <- strsplit(x, ":")[[1]]
      any(grepl("[A-Za-z]", splitted) == FALSE)
    })
  }
  # check intersect with full tree descendants
  orig_edge <- sapply(subtree_edge_desc, function(oneset){
    which(oneset == full_tree_children_by_edge)
  })
  return(orig_edge)
}

#' Input dataframes of assignment per subtree, along with subtrees and full tree
#' to obtain candidate assignments on full tree: returns all solutions
#'
#' @param assignment_dts list of assignment results per subtree
#' @param all_trees list of all subtrees
#' @param full_tree_ full tree
#' @param verbose1 logical whether to print progress
#'
#' @return data table of candidate results on full tree
#' @export
#'
#' @examples
#' \dontrun{
#' dt_combined <- convert_branches(assignment_dts = assignments,
#' all_trees = subtrees, full_tree_ = full_tree, verbose1 = FALSE)
#' }
#'
convert_branches <- function(assignment_dts,
                             all_trees,
                             full_tree_,
                             verbose1){
  chr_pos = sub_tree_id = original_branch = subtree_branch = NULL
  tmp_holder <- vector("list", length = length(all_trees))
  # convert assignment to branch to original tree branch number
  for(i in 1:length(all_trees)){
    if(verbose1){
      cat(sprintf("Converting branches for tree %s: %i of %i\n", names(all_trees)[i], i, length(all_trees)))
    }
    one_subtree <- all_trees[[i]]
    subtree_root <- names(all_trees)[i]

    if(!grepl("[a-zA-Z]", subtree_root)){ # subtree starts at node
      first_edge_in_orig <- which(full_tree_$edge[,1] == as.numeric(subtree_root))[1]
      sub_root_preceding_edge <- which(full_tree_$edge[,2] == as.numeric(subtree_root))
    }else{ # root case
      first_edge_in_orig <- 1
      sub_root_preceding_edge <- "n"
    }

    remapped_e <- remap_edges(one_subtree_x = one_subtree, full_tree = full_tree_)

    all_edges_for_subtree <- data.table::data.table(subtree_edges = c("n", 1:nrow(one_subtree$edge)),
                                                    original_edges = c(paste0("n_", sub_root_preceding_edge), remapped_e))

    if(is.null(assignment_dts)){
      sub_tree_assignments <- as.character(1:nrow(one_subtree$edge))
    }else{
      sub_tree_assignments <- assignment_dts[[i]][["only_branch"]]
    }

    original_tree_assignment <- sapply(sub_tree_assignments, function(one_comb){
      if(is.na(one_comb)){
        return(NA)
      }else{
        split1 <- strsplit(one_comb, "_")[[1]]
        split2 <- gsub("[A-Z]|\\.", "", split1)
        split <- split2
        new <- all_edges_for_subtree[match(split, all_edges_for_subtree[["subtree_edges"]]),][["original_edges"]]
        new_as_one_item <- paste0(new, collapse="_")
        return(new_as_one_item)
      }
    })

    if(is.null(assignment_dts)){
      tmp_holder[[i]] <- data.table::data.table("subtree_branch" = sub_tree_assignments,
                                    "original_branch" = original_tree_assignment,
                                    "sub_tree_id" = subtree_root)
    }else{
      assignment_dts[[i]][,original_branch := original_tree_assignment]
      assignment_dts[[i]][,sub_tree_id := subtree_root]
    }
  }

  if(is.null(assignment_dts)){
    dt_combined <- do.call("rbind", tmp_holder)
    dt_combined[, subtree_branch := as.numeric(subtree_branch)]
    dt_combined[, original_branch := as.numeric(original_branch)]
    data.table::setkey(dt_combined, original_branch)
  }else{
    dt_combined <- do.call("rbind", assignment_dts)
    data.table::setkey(dt_combined, chr_pos)
  }
  return(dt_combined)
}

#' Input candidate assignments across full tree and obtain tip states of full
#' tree
#'
#' @param dt_grouped datatable of all subtree solutions per site
#' @param subtrees list of subtrees
#' @param full ape phylo object of full tree
#'
#' @return tips states of full tree given subtree assignments
#' @export
#'
#' @examples
#' \dontrun{
#' preliminary_solutions <- lapply(combined_by_chr_pos,
#' get_preliminary_solution, subtrees = subtrees, full = full_tree)
#' }
get_preliminary_solution <- function(dt_grouped, subtrees, full){
  subtree_list <- subtrees[dt_grouped$sub_tree_id]
  all_expected_tip_states <- mapply(generate_expected_tip_states,
                                    ss = dt_grouped$starting_state,
                                    changed_branches = dt_grouped$only_branch,
                                    change_types = dt_grouped$branch_state,
                                    tree = subtree_list)
  names(all_expected_tip_states) <- names(subtree_list)
  # find terminal tips
  expected <- unlist(all_expected_tip_states)
  names(expected) <- unlist(lapply(subtree_list, function(x) x$tip.label))

  expected <- expected[match(full$tip.label, names(expected))]
  return(expected)
}

extract_branches <- function(dt_grouped, full){
  sub_tree_id = starting_state = combined_branch = ss_combined = NULL
  original_branch <- dt_grouped[["original_branch"]]

  starting <- dt_grouped[sub_tree_id == "root",][["starting_state"]]

  if(all(is.na(original_branch))){ # this must be first condition
    new_assignments <- NA
  }else if(all(grepl("n_", original_branch)) | all(grepl("n_", original_branch)|is.na(original_branch))){
    new_assignments <- "n"
  }else{
    not_n <- grep("n_", original_branch, value = TRUE, invert = TRUE)
    not_n <- not_n[!is.na(not_n)]
    not_n2 <- lapply(not_n, function(nn){
      strsplit(nn, "_")[[1]]
    })
    not_n3 <- as.numeric(unlist(not_n2))
    new_assignments <- paste0(sort(not_n3), collapse = "_")
  }

  if(length(unique(dt_grouped[["starting_state"]])) > 1){
    # find subtree with different starting state
    dt_mismatch <- dt_grouped[starting_state != starting]
    parent_node_of_odd <- as.numeric(dt_mismatch[["sub_tree_id"]])
    additional_edges <- match(parent_node_of_odd, full$edge[,2])
  }else{
    additional_edges <- c()
  }

  if(length(additional_edges) > 0){
    if(new_assignments != "n"){
      current_assignments <- as.numeric(strsplit(new_assignments, "_")[[1]])
      updated_assignments <- sort(unique(c(current_assignments, additional_edges)))
    }else{
      updated_assignments <- sort(unique(additional_edges))
    }
    updated2 <- paste0(updated_assignments, collapse = "_")
  }else{
    updated2 <- new_assignments
  }

  dt_grouped[, combined_branch := updated2]
  dt_grouped[, ss_combined := starting]
  return(dt_grouped)
}

#' Get proportion of reads in state S given read and depth matrices
#'
#' @param tips_mtx logical matrix of tip states
#' @param reads_mtx matrix of methylated reads
#' @param depth_mtx matrix of depth
#' @param min_in_s minimum number of tips in state s required
#'
#' @return vector of proportion of tips in state S
#' @export
#'
#' @examples
#' \dontrun{
#' # for 3 state model
#' gamma_site <- get_proportion_of_reads(tips_mtx = prelim_s_mtx == 3,
#' reads_mtx = unmeth_matrix, depth_mtx = depth_matrix,
#' min_in_s = min_tips_in_s)
#'
#' delta_site <- get_proportion_of_reads(tips_mtx = prelim_s_mtx == 1,
#' reads_mtx = meth_matrix, depth_mtx = depth_matrix, min_in_s = min_tips_in_s)
#' }
get_proportion_of_reads <- function(tips_mtx, reads_mtx, depth_mtx, min_in_s = 5){
  # Estimate gamma as the proportion of unmethylated reads across sites with genotype=2
  # Estimate delta as the proportion of methylated reads across sites with genotype=0
  if(!all(rownames(reads_mtx) == rownames(tips_mtx))){
    reads_mtx <- reads_mtx[match(rownames(tips_mtx), rownames(reads_mtx)),]
  }
  if(!all(rownames(reads_mtx) == rownames(depth_mtx))){
    depth_mtx <- depth_mtx[match(rownames(reads_mtx), rownames(depth_mtx)),]
  }
  stopifnot(all(rownames(tips_mtx) == rownames(reads_mtx)) & all(rownames(tips_mtx) == rownames(depth_mtx)))
  proportion <- sapply(1:nrow(tips_mtx), function(idx){
    relevant_in_state <- reads_mtx[idx, tips_mtx[idx,]]
    if(length(relevant_in_state) < min_in_s){
      prop <- NA
    }else{
      relevant_total <- depth_mtx[idx, tips_mtx[idx,]]
      prop <- sum(relevant_in_state)/sum(relevant_total)
      if(is.nan(prop)){
        prop <- NA
      }
    }
    return(prop)
  })
  names(proportion) <- rownames(tips_mtx)
  return(proportion)
}


#' Impute gamma and delta vector with its rolling mean*
#'
#' This function will impute values of gamma and delta with a rolling average.
#' The range can be specified. If adjacent sites have too many NAs, mean of
#' entire vector will be used to fill. NB that the input vector must be sorted!
#'
#' @param gd_vector sorted vector of gammas or deltas per site
#' @param range range of sites to consider, centered
#' @param fillwith value to fill NAs with if rolling average doesn't have enough info; if NULL then mean of entire vector is used to fill
#'
#' @return vector of gamma and delta with imputed values
#' @export
#'
#' @examples
#' \dontrun{
#' gd_imputed <- lapply(gd_per_site_sorted, impute_gamma_delta)
#' }
impute_gamma_delta <- function(gd_vector, range = 10, fillwith = NULL){
  running = error = error_imp = NULL
  gd_dt <- data.table::data.table(chr_pos = names(gd_vector),
                                  error = gd_vector)
  if(is.null(fillwith)){
    fillwith <- mean(gd_vector, na.rm = TRUE)
  }
  gd_dt[, running := data.table::frollmean(x = error, n = range,
                                           align = "center", na.rm = TRUE,
                                           hasNA = TRUE, fill = fillwith)]
  gd_dt[, error_imp := ifelse(is.na(error), running, error)]
  gd_new <- gd_dt[["error_imp"]]
  names(gd_new) <- gd_dt[["chr_pos"]]
  return(gd_new)
}


#' Obtain gamma and delta values per site
#'
#' @param assignments either list of preliminary assignments on each subtree or
#' data.frame of post-MLE assignments
#' @param subtrees either list of subtrees (each list element is an ape phylo
#' object) or NULL if computing gamma & delta on post-MLE assignments
#' @param full_tree ape phylo object of full tree
#' @param meth_matrix matrix of methylated reads
#' @param unmeth_matrix matrix of unmethylated reads
#' @param matrix_list list of matrices (when running on phased data)
#' @param min_tips_in_s minimum number of tips that must be in a state
#'
#' @return list of two vectors, gamma and delta
#' @export
#'
#' @examples
#' \dontrun{
#' gd_per_site <- get_gamma_delta(assignments = assignments,
#' subtrees = subtrees, full_tree = tree, meth_matrix = meth_mtx,
#' unmeth_matrix = unmeth_mtx)
#' }
get_gamma_delta <- function(assignments,
                            subtrees = NULL,
                            full_tree,
                            meth_matrix = NULL,
                            unmeth_matrix = NULL,
                            matrix_list = NULL,
                            min_tips_in_s = 5){
  if(is.null(subtrees)){ # gamma & delta on full tree
    # project assignment onto full tree
    list_trees <- lapply(1:nrow(assignments), function(x) return(full_tree))
    children_matrices_tmp <- generate_tree_children_matrix(tree_ = full_tree)
    list_ch_mtx <- lapply(1:nrow(assignments), function(x) return(children_matrices_tmp))
    ts_list <- mapply(generate_expected_tip_states,
                      ss = assignments$starting_state,
                      changed_branches = assignments$only_branch,
                      change_types = assignments$branch_state,
                      tree = list_trees,
                      children_matrices = list_ch_mtx,
                      SIMPLIFY = FALSE)
    prelim_s_mtx <- do.call("rbind", ts_list)
    rownames(prelim_s_mtx) <- assignments$chr_pos
    colnames(prelim_s_mtx) <- full_tree$tip.label
    # prepare prelim_s_mtx
    # TODO: rbind list ts_list ?
  }else{
    dt_combined <- convert_branches(assignment_dts = assignments,
                                    all_trees = subtrees,
                                    full_tree_ = full_tree,
                                    verbose1 = FALSE)
    # filter to rank = 1
    if("rank" %in% colnames(dt_combined)){
      if(length(unique(dt_combined$rank)) > 1){
        cat(sprintf("Rank column suggests there are several top solutions, filtering to rank = 1 only\n"))
        dt_combined <- dt_combined[rank == 1,]
      }
    }

    # combine across sites
    combined_by_chr_pos <- split(dt_combined, by = "chr_pos")

    preliminary_solutions <- lapply(combined_by_chr_pos,
                                    get_preliminary_solution,
                                    subtrees = subtrees,
                                    full = full_tree)
    prelim_s_mtx <- do.call("rbind", preliminary_solutions)
  }

  if(the$number_of_states < 4){
    depth_matrix <- meth_matrix + unmeth_matrix

    if(the$number_of_states == 2){
      gamma_site <- get_proportion_of_reads(tips_mtx = prelim_s_mtx == 2,
                                            reads_mtx = unmeth_matrix,
                                            depth_mtx = depth_matrix,
                                            min_in_s = min_tips_in_s)

      delta_site <- get_proportion_of_reads(tips_mtx = prelim_s_mtx == 1,
                                            reads_mtx = meth_matrix,
                                            depth_mtx = depth_matrix,
                                            min_in_s = min_tips_in_s)
    }else if(the$number_of_states == 3){
      gamma_site <- get_proportion_of_reads(tips_mtx = prelim_s_mtx == 3,
                                            reads_mtx = unmeth_matrix,
                                            depth_mtx = depth_matrix,
                                            min_in_s = min_tips_in_s)

      delta_site <- get_proportion_of_reads(tips_mtx = prelim_s_mtx == 1,
                                            reads_mtx = meth_matrix,
                                            depth_mtx = depth_matrix,
                                            min_in_s = min_tips_in_s)
    }
  }else{
    unmeth_mtx <- Reduce('+', matrix_list[c("MT_U", "PT_U", "U_U")])
    meth_mtx <- Reduce('+', matrix_list[c("MT_M", "PT_M", "U_M")])
    depth_mtx1 <- unmeth_mtx + meth_mtx
    gamma_site <- get_proportion_of_reads(tips_mtx = prelim_s_mtx == 4,
                                          reads_mtx = unmeth_mtx,
                                          depth_mtx = depth_mtx1,
                                          min_in_s = min_tips_in_s)

    delta_site <- get_proportion_of_reads(tips_mtx = prelim_s_mtx == 1,
                                          reads_mtx = meth_mtx,
                                          depth_mtx = depth_mtx1,
                                          min_in_s = min_tips_in_s)
  }
  return(list("gamma" = gamma_site,
              "delta" = delta_site))
}

#' Convert sub tree solution to full tree
#'
#' @param assignments List of assignments per subtree
#' @param subtrees List of subtrees (each is "phylo" object)
#' @param full_tree "phylo" object of full tree
#' @param verbose Logical
#'
#' @return data.table of assignments for the full tree
#' @export
#'
#' @examples
#' \dontrun{
#' converted_dt <- convert_sub_to_full_tree_branches(assignments = assignments,
#' subtrees = subtrees,
#' full_tree = full_tree)
#' }
convert_sub_to_full_tree_branches <- function(assignments,
                                              subtrees,
                                              full_tree,
                                              verbose = TRUE){
  chr_pos = rank = n_assignments = sub_tree_id = combined_branch = NULL
  dt_combined <- convert_branches(assignment_dts = assignments,
                                  all_trees = subtrees,
                                  full_tree_ = full_tree,
                                  verbose1 = verbose)
  # filter to rank = 1
  if("rank" %in% colnames(dt_combined)){
    if(length(unique(dt_combined$rank)) > 1){
      cat(sprintf("Rank column suggests there are several top solutions, filtering to rank = 1 only\n"))
      dt_combined <- dt_combined[rank == 1,]
    }
  }

  # combine across sites
  combined_by_chr_pos <- split(dt_combined, by = "chr_pos")

  renamed_branch_dts <- lapply(combined_by_chr_pos, extract_branches, full = full_tree)

  if(any(is.null(renamed_branch_dts))){
    im_issue <- which(is.null(renamed_branch_dts))
    cat(sprintf("Failed to convert branches from sub tree to original tree for %s\n", paste0(names(combined_by_chr_pos)[im_issue], collapse = ", ")))
  }

  dt_combined2 <- data.table::rbindlist(renamed_branch_dts)

  dt_combined_unq <- dt_combined2[sub_tree_id == "root", c("chr_pos", "rank", "combined_branch", "ss_combined", "probability")]
  data.table::setnames(dt_combined_unq, c("ss_combined", "probability"), c("ss_root_tree", "probability_root"))
  branch_dt <- dt_combined_unq
  data.table::setkey(branch_dt, chr_pos)

  branch_dt[, n_assignments := sapply(combined_branch, function(x){
    y <- strsplit(x, "_")[[1]]
    if(all(y == "n")){
      return(0)
    }else{
      return(length(y))
    }
  })]
  return(branch_dt)
}

convert_to_matrix <- function(matrix, chr_pos){
  if(length(chr_pos) == 1){
    r <- t(as.matrix(matrix[rownames(matrix) %in% chr_pos,]))
    rownames(r) <- chr_pos
  }else{
    r <- matrix[rownames(matrix) %in% chr_pos,]
  }
  return(r)
}

generate_permutations <- function(n, m = m, only_one = FALSE){
  # prepare all permutations across n positions which have up to m changes
  z <- rep(0, n)
  if(only_one){
    lists <- DescTools::CombSet(x = 1:n, m = m, as.list = TRUE)
    b <- lapply(lists, FUN=function(k) {z[k]=1; z})
    c <- do.call("rbind", b)
  }else{
    if(m %% 2 == 0){
      selection <- even(0:m)
    }else{
      selection <- odd(0:m)
    }
    lists <- lapply(selection, function(x){
      return(DescTools::CombSet(0:n, x, as.list = TRUE))
    })
    a <- do.call("c", lists)
    b <- lapply(a, FUN=function(k) {z[k]=1; z})
    c <- do.call("rbind", b)
  }
  return(c)
}

#' Get rates
#'
#' @param branch solution
#' @param priors_w_idx priors
#' @param full_tree ape phylo object
#'
#' @return rates
#' @export
#'
#' @examples
#' \dontrun{
#'    states_rates <- lapply(br_probs_new$branch, get_rate,
#'    priors_w_idx = br_priors_w_idx, full_tree = full_tree)
#' }
get_rate <- function(branch,
                     priors_w_idx,
                     full_tree){
  branch_split <- tryCatch(strsplit(branch, "\\.")[[1]], error=function(e) NULL)
  if(is.null(branch_split)){
    return(NA)
  }

  only_branch <- branch_split[1]
  start_brstates <- branch_split[2]

  branch_state_indices <- priors_w_idx[[only_branch]][[start_brstates]][[2]]

  time_in_state <- sapply(the$all_possible_states, function(s){
    trans_w_start <- as.integer(the$state_transition_mtx[s,])
    return(sum(full_tree$edge.length[branch_state_indices %in% trans_w_start]))
  })
  names(time_in_state) <- paste0("ss_", the$all_possible_states)

  n_trans <- sapply(c(1:length(the$all_possible_states)**2), function(s){
    return(sum(branch_state_indices == s))
  })
  names(n_trans) <- paste0("nt_", c(1:length(the$all_possible_states)**2))

  rate_result <- c(time_in_state, n_trans)
  return(rate_result)
}

#' Produce matrix of affected tip and branches for the tree
#'
#' @param tree_ ape phylo object
#'
#' @return list of lists, tips: affected tips per branch (rows are branches,
#' columns are affected tips), branch_children_only: all descendant branches of
#' affected branch (branch hit is row, its descendants have 1 in their columns)
#' @export
#'
#' @examples
#' \dontrun{
#'  children_matrices <- generate_tree_children_matrix(tree_ = tree)
#' }
generate_tree_children_matrix <- function(tree_){
  # tips
  children <- lapply(tree_$edge[,2], function(node1){
    desc <- phangorn::Descendants(tree_, node1, type = c("tips"))[[1]]
    vec <- rep(0, length(tree_$tip.label))
    vec[desc] <- 1
    return(vec)
  })
  chld_mtx <- do.call("rbind", children)

  # branches
  branch_children <- lapply(tree_$edge[,2], function(node1){
    desc2 <- phangorn::Descendants(tree_, node1, type = c("all"))
    desc_br <- match(desc2, tree_$edge[,2])
    self <- which(tree_$edge[,2] == node1)
    all_affected <- unique(c(desc_br, self))
    vec2 <- rep(0, nrow(tree_$edge))
    vec2[all_affected] <- 1
    return(vec2)
  })
  branch_mtx <- do.call("rbind", branch_children)

  branch_children_noself <- lapply(tree_$edge[,2], function(node1){
    desc2 <- phangorn::Descendants(tree_, node1, type = c("all"))
    desc_br <- match(desc2, tree_$edge[,2])
    vec2 <- rep(0, nrow(tree_$edge))
    vec2[desc_br] <- 1
    return(vec2)
  })
  branch_mtx_noself <- do.call("rbind", branch_children_noself)

  branch_children_only <- lapply(tree_$edge[,2], function(node1){
    desc2 <- phangorn::Descendants(tree_, node1, type = c("all"))
    desc_br <- match(desc2, tree_$edge[,2])
    vec2 <- rep(0, nrow(tree_$edge))
    self <- which(tree_$edge[,2] == node1)
    only_kids <- unique(setdiff(desc_br, self))
    vec2[only_kids] <- 1
    return(vec2)
  })
  branch_children_only_mtx <- do.call("rbind", branch_children_only)

  return(list("tips"=chld_mtx,
              "branches" = branch_mtx,
              "branches_noself" = branch_mtx_noself,
              "branch_children_only" = branch_children_only_mtx))
}

#' Prepare expected tip states given changes on branch
#'
#' @param ss integer starting state
#' @param changed_branches vector of affected branches, either character (branches separated by _) or integer
#' @param change_types vector of branch states, either character (branches separated by _) or integer
#' @param tree ape phylo object
#' @param children_matrices list of matrices of affected tips/branches
#'
#' @return vector of tip states
#' @export
#'
#' @examples
#' \dontrun{
#'  tips <- generate_expected_tip_states(ss = 3, changed_branches = '29_56',
#'  change_types = '8_8', tree = tree)
#' }
generate_expected_tip_states <- function(ss, changed_branches, change_types,
                                         tree,
                                         children_matrices = NULL){
  if(is.null(children_matrices)){
    children_matrices <- generate_tree_children_matrix(tree_ = tree)
  }
  tip_child_mtx <- children_matrices$tips

  result <- rep(ss, length(tree$tip.label))

  if(changed_branches != "n"){
    if(is.character(changed_branches)){
      changed_branches_split <- as.numeric(strsplit(changed_branches, "_")[[1]])
    }else{
      changed_branches_split <- changed_branches
    }
    if(is.character(change_types)){
      change_types_split <- as.numeric(strsplit(change_types, "_")[[1]])
    }else{
      change_types_split <- change_types
    }

    for(i in 1:length(changed_branches_split)){
      cur_branch <- changed_branches_split[i]
      cur_change <- change_types_split[i]
      new_state <- which(the$state_transition_mtx == cur_change, arr.ind = TRUE)[2]
      child_vector <- tip_child_mtx[cur_branch,]
      result[as.logical(child_vector)] <- new_state
    }
  }
  return(result)
}

#' Generate expected branch states given starting state, affected branches and state changes
#'
#' @param ss integer starting state
#' @param changed_branches vector of affected branches, either character (branches separated by _) or integer
#' @param change_types vector of branch states, either character (branches separated by _) or integer
#' @param tree ape phylo object
#' @param children_matrices list of matrices of affected tips/branches
#'
#' @return vector of branch states
#' @export
#'
#' @examples
#' \dontrun{
#'  branch_states <- generate_expected_branch_states(ss = 3,
#'  changed_branches = "29_56", change_types = "8_8", tree = tree)
#' }
generate_expected_branch_states <- function(ss, changed_branches, change_types,
                                            tree,
                                            children_matrices = NULL){
  if(is.null(children_matrices)){
    children_matrices <- generate_tree_children_matrix(tree_ = tree)
  }
  tip_child_mtx <- children_matrices$branch_children_only

  result <- rep(the$state_transition_mtx[ss,ss], nrow(tree$edge))

  if(changed_branches != "n"){
    if(is.character(changed_branches)){
      changed_branches_split <- as.numeric(strsplit(changed_branches, "_")[[1]])
    }else{
      changed_branches_split <- changed_branches
    }
    if(is.character(change_types)){
      change_types_split <- as.numeric(strsplit(change_types, "_")[[1]])
    }else{
      change_types_split <- change_types
    }
    if(is.unsorted(changed_branches_split)){
      new_order <- order(changed_branches_split)
      changed_branches_split <- changed_branches_split[new_order]
      change_types_split <- change_types_split[new_order]
    }
    for(i in 1:length(changed_branches_split)){
      cur_branch <- changed_branches_split[i]
      cur_change <- change_types_split[i]
      result[cur_branch] <- cur_change
      child_vector <- tip_child_mtx[cur_branch,]
      new_state_p <- which(the$state_transition_mtx == cur_change, arr.ind = TRUE)[2]
      new_state <- the$state_transition_mtx[new_state_p,new_state_p]
      result[as.logical(child_vector)] <- new_state
    }
  }
  return(result)
}

#' Get 3 state rates matrix
#' @details
#' Input the methylation change rate per snv change rate
#' per base pair. Here states are 0, 1, 2 and correspond to 1, 2, 3 respectively. Written by NW.
#'
#' @param r01 rate for transition 0 to 1
#' @param r02 rate for transition 0 to 2
#' @param r10 rate for transition 1 to 0
#' @param r21 rate for transition 2 to 1
#' @param r20 rate for transition 2 to 0
#' @param r12 rate for transition 1 to 2
#'
#' @return matrix of rates
#' @export
#'
#' @examples
#' \dontrun{
#' epi_mut_rate_bck=1e5*20/3e9
#' epi_mut_rate_fwd=1e4*20/3e9
#' Q_mtx <- get3StateRates(r01=2*epi_mut_rate_fwd, r02=0, r10=epi_mut_rate_bck,
#' r12=epi_mut_rate_fwd, r21=2*epi_mut_rate_bck, r20=0)
#' }
get3StateRates <- function(r01,r02,r10,r21=r01,r20=r02,r12=r10){
  matrix(c(-(r01+r02),r01,r02,
           r10,-(r10+r12),r12,
           r20,r21,-(r20+r21)),ncol=3,byrow = TRUE)
}


#' Get 4 state rates matrix
#' @details
#' Input the forward and backward rates (cpg rate times snv rate per basepair)
#' and get Q matrix in return. Written by NW.
#'
#' @param fwd numeric forward rate
#' @param bck numeric backward rate
#'
#' @return Q matrix
#' @export
#'
#' @examples
#' \dontrun{
#' epi_mut_rate_bck=1e5*20/3e9
#' epi_mut_rate_fwd=1e4*20/3e9
#' Q2 <- get4StateRates(fwd = epi_mut_rate_fwd,
#' bck = epi_mut_rate_bck)
#' }
get4StateRates <- function(fwd, bck){
  out <- matrix(c(0, fwd, fwd, 0, bck, 0, 0, fwd, bck, 0, 0, fwd, 0, bck, bck, 0), byrow = TRUE, ncol = 4)
  diag(out) <- -rowSums(out)
  return(out)
}
