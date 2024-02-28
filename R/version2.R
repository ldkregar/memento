#' Function to generate genotype vectors of tips only
#'
#' @param full_tree ape object of full tree
#' @param subtrees list of ape objects - subtrees
#'
#' @return list of logical vectors per sub tree
#' @export
#'
#' @examples
#' \dontrun{
#' real_tip_vectors <- get_true_tip_vector_on_full(full_tree = full_tree,
#' subtrees = subtrees)
#' }
get_true_tip_vector_on_full <- function(full_tree, subtrees){
  tip_vectors <- lapply(subtrees, function(sub1){
    return(full_tree$tip.label %in% sub1$tip.label)
  })
  return(tip_vectors)
}


get_tip_only_genotypes <- function(current_genotypes, node_is_real_tip){
  # for each genotype vector, change tips belonging to subtree clade to 0, return
  new_genotypes <- current_genotypes
  if(nrow(new_genotypes) == length(node_is_real_tip)){
    new_genotypes[!node_is_real_tip,] <- 0
  }else{
    stop("Number of rows of the genotype matrix and lenght of tip vector on full tree do not match")
  }
  return(new_genotypes)
}

# store subtree solutions
get_max_lik_per_ss <- function(lik_mtx){
  cnames <- colnames(lik_mtx)
  test_pattern <- "(?<=\\.).*(?=:)"
  genotype_vec_ss <- regmatches(cnames,
                                gregexpr(test_pattern, cnames, perl=TRUE))
  genotype_vec_ss <- as.numeric(genotype_vec_ss)
  unq_ss <- sort(unique(genotype_vec_ss))

  # subset matrix by ss
  idx_max <- lapply(unq_ss, function(oness){
    tf_vec <- genotype_vec_ss == oness
    tmp <- lik_mtx[,tf_vec,drop=FALSE]
    cn <- colnames(tmp)
    res1 <- apply(X = tmp, MARGIN = 1, FUN = which.max)
    res_name <- sapply(res1, function(x){
      cn[x]
    })
    res_value <- tmp[cbind(c(1:nrow(tmp)), res1)]
    return(list("name" = res_name,
                "value" = res_value))
  })
  names(idx_max) <- unq_ss
  return(idx_max)
}

add_values_to_end <- function(current, new, item, subtree_name = NULL){
  new2 <- new[[item]]
  new3 <- cbind(current, new2)
  if(!is.null(subtree_name)){
    colnames(new3)[ncol(new3)] <- subtree_name
  }
  return(new3)
}

take_br_return_new <- function(br_string, reorder_v){
    if(br_string == "n"){
      new_string <- br_string
    }else{
      new_string <- paste0(reorder_v[as.numeric(strsplit(br_string, "_")[[1]])], collapse = "_")
    }
  return(new_string)
}

reconstruct_recursive <- function(node, row_idx, node_ss){
  if(node == the$full_tree$edge[1,1]){ # root
    node_chr <- "root"
  }else{
    node_chr <- as.character(node)
  }
  node_tips <- the$global_tips_vector[names(the$global_tips_vector) == node_chr]
  node_tips_subtrees <- the$subtree_sub_nodes[names(the$subtree_sub_nodes) == node_chr]
  if(node_chr == "root"){
    assignment <- the$top_per_site_name[row_idx]
  }else{
    assignment <- the$list_of_subtree_sol_name_per_ss[[node_ss]][row_idx, node_chr]
  }
  assignment_split <- strsplit(assignment, split = "\\.|:")[[1]]
  node_ss <- as.numeric(assignment_split[2])
  node_branches <- strsplit(assignment_split[1], "_")[[1]]
  if(all(node_branches == "n")){
    new_assignment <- NULL
  }else{
    node_branches_num <- as.numeric(node_branches)
    edges_of_this_node <- which(names(the$global_edges_vector_sort) == node_chr)
    real_edges <- edges_of_this_node[node_branches_num]
    new_assignment <- sprintf("%s:%s", paste0(real_edges, collapse = "_"), assignment_split[3])
  }
  if(length(node_tips_subtrees) == 0){ # tippy tree
    return(new_assignment)
  }else{ # go down subtree
    node_tip_out_states <- the$global_subtree_priors_affected_tips[[node_chr]][names(the$global_subtree_priors_affected_tips[[node_chr]]) == assignment][[1]]
    next_nodes <- node_tip_out_states[grep("[[:alpha:]]", names(node_tip_out_states), invert = TRUE)]
    return(c(new_assignment,
             mapply(FUN = reconstruct_recursive,
                    node = names(next_nodes),
                    row_idx = rep(row_idx, length(next_nodes)),
                    node_ss = next_nodes, SIMPLIFY = FALSE)))
  }
}

clean_up <- function(items, root_ss){
  unlisted <- unlist(items)
  if(is.null(unlisted)){
    final_assign <- sprintf("n.%i:%i", root_ss, the$state_transition_mtx[root_ss, root_ss])
  }else{
    split_prelim <- strsplit(unlisted, ":")
    if(length(unlisted) == 1){
      final_assign <- sprintf("%s.%i:%s", split_prelim[[1]][1], root_ss, split_prelim[[1]][2])
    }else{
      all_branches <- sapply(split_prelim, "[[", 1)
      all_states <- sapply(split_prelim, "[[", 2)
      all_n_num <- as.numeric(unlist(strsplit(all_branches, "_")))
      all_s_num <- as.numeric(unlist(strsplit(all_states, "_")))
      b_ord <- order(all_n_num)
      final_assign <- sprintf("%s.%i:%s", paste0(all_n_num[b_ord], collapse = "_"),
                              root_ss, paste0(all_s_num[b_ord], collapse = "_"))
    }
  }
  return(final_assign)
}


#' Infer set of changes with expectation-maximisation: this function will infer
#' assignments per subtree and then combine them. It will also infer gamma and
#' delta for each assignment and then run expectation maximisation to infer the
#' most likely assignment.
#'
#' @param full_tree ape phylo object of the full tree
#' @param subtrees list of ape phylo objects
#' @param priors_p_subtree list of priors per subtree
#' @param meth_mtx for 2 or 3 state model: matrix of methylated reads (rows are CpG sites, columns are samples)
#' @param unmeth_mtx for 2 or 3 state model: matrix of unmethylated reads (rows are CpG sites, columns are samples)
#' @param matrix_list for 4 state model: list of phased reads matrices
#' @param ro error rate
#' @param gammas error rate: proportion of unmethylated reads when site is methylated (state 1)
#' @param deltas error rate: proportion of methylated reads when site is unmethylated (state 0)
#' @param em_rounds integer: number of EM rounds
#' @param return_intermediate_assignments logical: whether to return solutions from each EM round
#' @param gamma_prior numeric: prior for gamma (if not enough tips in state to estimate gamma, this will be used)
#' @param delta_prior numeric: prior for delta (if not enough tips in state to estimate delta, this will be used)
#' @param min_gamma numeric: minimum value for gamma
#' @param min_delta numeric: minimum value for delta
#' @param max_gamma numeric: maximum value for gamma
#' @param max_delta numeric: maximum value for delta
#' @param min_tips_in_state integer: minimum number of tips in state 0/1 to infer gamma and delta
#'
#' @return data.table of final assignments or list of data.table of assignments and intermediate solutions
#' @export
#'
#' @examples
#' \dontrun{
#' priors_p_subtree <- generate_priors(all_trees = subtrees, full_tree = tree)
#' gv <- rep(0.05, nrow(meth_mtx))
#' names(gv) <- rownames(meth_mtx)
#' dv <- gv
#' assignments <- run_EM_on_subtrees(full_tree = tree, subtrees = subtrees,
#' priors_p_subtree = priors_p_subtree, meth_mtx = meth_mtx,
#' unmeth_mtx = unmeth_mtx, gammas = gv, deltas = dv)
#' }
run_EM_on_subtrees <- function(full_tree,
                               subtrees,
                               priors_p_subtree = NULL,
                               meth_mtx = NULL,
                               unmeth_mtx = NULL,
                               matrix_list = NULL,
                               ro = NA,
                               gammas = NA,
                               deltas = NA,
                               em_rounds = 5,
                               return_intermediate_assignments = FALSE,
                               gamma_prior = 0.05,
                               delta_prior = 0.05,
                               min_gamma = 0.01,
                               min_delta = 0.01,
                               max_gamma=0.4,
                               max_delta=0.4,
                               min_tips_in_state=5){
  chr_pos = branch = only_branch = starting_state = n_assignments = Chrom = NULL
  Pos = probability = em_round = delta = gamma = NULL

  if(is.null(priors_p_subtree)){
    priors_p_subtree <- generate_priors(all_trees = subtrees, full_tree = full_tree)
  }

  invisible <- wipe_before_EM()

  the$full_tree <- full_tree

  global_edges <- lapply(subtrees, FUN = remap_edges, full_tree)
  global_edges_vector <- unlist(global_edges)
  names(global_edges_vector) <- unlist(mapply(FUN = function(node_name, len){
    return(rep(node_name, len))
  }, node_name = names(global_edges), len = lengths(global_edges)))
  global_edges_vector_sort <- global_edges_vector[order(global_edges_vector)]
  the$global_edges_vector_sort <- global_edges_vector_sort

  global_tips <- lapply(subtrees, function(x) x$tip.label)
  global_tips_vector <- unlist(global_tips)
  names(global_tips_vector) <- unlist(mapply(FUN = function(node_name, len){
    return(rep(node_name, len))
  }, node_name = names(global_tips), len = lengths(global_tips)))
  the$global_tips_vector <- global_tips_vector

  the$subtree_sub_nodes <- global_tips_vector[grep("[[:alpha:]]", global_tips_vector, invert = TRUE)]

  global_subtree_priors_affected_tips <- lapply(priors_p_subtree, function(x){
    node_priors <- unlist(x$affected_tips, recursive = FALSE)
    return(node_priors)
  })
  the$global_subtree_priors_affected_tips <- global_subtree_priors_affected_tips

  real_tip_vectors <- get_true_tip_vector_on_full(full_tree = full_tree,
                                                  subtrees = subtrees)
  while_counter <- 1
  save_data <- vector("list", length = em_rounds)
  final_list <- vector("list", length = em_rounds)
  if(is.null(matrix_list)){
    site_names <- rownames(meth_mtx)
  }else{
    site_names <- rownames(matrix_list[[1]])
  }
  original_sites <- data.table::data.table(chr_pos = site_names)

  while(while_counter <= em_rounds){

    if(is.null(matrix_list)){ # 2 or 3 states
      meth_list <- lapply(c(1:length(subtrees)), function(x){
        return(meth_mtx)
      })
      unmeth_list <- lapply(c(1:length(subtrees)), function(x){
        return(unmeth_mtx)
      })
      nsites <- nrow(meth_mtx)
      site_names <- rownames(meth_mtx)
    }else{
      phased_per_subtree <- lapply(c(1:length(subtrees)), function(x){
        return(matrix_list)
      })
      nsites <- nrow(matrix_list[[1]])
      site_names <- rownames(matrix_list[[1]])
    }
    dt_chunk_sites <- data.table::data.table(chr_pos = site_names)

    if(!is.na(ro)){
      ro_list <- rep(ro, length(priors_p_subtree))
    }else{
      ro_list <- lapply(priors_p_subtree, function(x){
        return(NA)
      })
      if(length(gammas) == 1 & length(deltas) == 1){
        gamma_list <- rep(gammas, length(priors_p_subtree))
        delta_list <- rep(deltas, length(priors_p_subtree))
      }else{
        gamma_list <- lapply(priors_p_subtree, function(x){
          return(gammas)
        })
        delta_list <- lapply(priors_p_subtree, function(x){
          return(deltas)
        })
      }
    }

    if(is.null(matrix_list)){
      likelihoods <- mapply(FUN = compute_probabilities3,
                            priors_per_tree = priors_p_subtree,
                            m_reads_matrix = meth_list,
                            u_reads_matrix = unmeth_list,
                            ro=ro_list,
                            gamma = gamma_list,
                            delta = delta_list,
                            true_tip_vector = real_tip_vectors,
                            only_return_likelihoods = rep(TRUE, length(priors_p_subtree)),
                            add_priors_to_likelihoods = rep(TRUE, length(priors_p_subtree)),
                            SIMPLIFY = FALSE)
    }else{ # 4 state model
      likelihoods <- mapply(FUN = compute_probabilities3,
                            priors_per_tree = priors_p_subtree,
                            matrix_list = phased_per_subtree,
                            ro=ro_list,
                            gamma = gamma_list,
                            delta = delta_list,
                            true_tip_vector = real_tip_vectors,
                            only_return_likelihoods = rep(TRUE, length(priors_p_subtree)),
                            add_priors_to_likelihoods = rep(TRUE, length(priors_p_subtree)),
                            SIMPLIFY = FALSE)
    }

    nodes_to_visit <- names(subtrees)
    nodes_to_visit[nodes_to_visit == "root"] <- full_tree$edge[1,1]
    nodes_to_visit <- as.numeric(nodes_to_visit)
    names(nodes_to_visit) <- names(subtrees)
    nodes_to_visit <- nodes_to_visit[order(nodes_to_visit, decreasing = TRUE)]

    list_of_subtree_sol_name_per_ss <- vector("list", length = length(the$number_of_states))
    list_of_subtree_sol_loglik_per_ss <- vector("list", length = length(the$number_of_states))

    while(length(nodes_to_visit) > 0){
      last_node <- nodes_to_visit[1]
      nodes_to_visit <- nodes_to_visit[nodes_to_visit != last_node]
      subtree_tips <- subtrees[[names(last_node)]]$tip.label
      max_lik_ss <- get_max_lik_per_ss(lik_mtx = likelihoods[[names(last_node)]])
      if(all(grepl("[[:alpha:]]", subtree_tips))){ # subtree has tips only
        list_of_subtree_sol_name_per_ss <- mapply(FUN = add_values_to_end,
                                                  current = list_of_subtree_sol_name_per_ss,
                                                  new = max_lik_ss,
                                                  item = rep(1, length(the$number_of_states)),
                                                  subtree_name = rep(last_node, length(the$number_of_states)),
                                                  SIMPLIFY = FALSE)
        list_of_subtree_sol_loglik_per_ss <- mapply(FUN = add_values_to_end,
                                                    current = list_of_subtree_sol_loglik_per_ss,
                                                    new = max_lik_ss,
                                                    item = rep(2, length(the$number_of_states)),
                                                    subtree_name = rep(last_node, length(the$number_of_states)),
                                                    SIMPLIFY = FALSE)
      }else{ # subtree contains a subtree
        # for each assignment, find subtree states
        node_priors <- priors_p_subtree[[names(last_node)]]
        node_affected_tips <- unlist(node_priors$affected_tips, recursive = FALSE)
        # which tips are subtrees
        subtree_parents <- grep("[[:alpha:]]", subtree_tips, invert = TRUE, value = TRUE)
        subtree_parents_idx <- match(subtree_parents, names(node_affected_tips[[1]]))
        list_subtree_genotypes <- lapply(node_affected_tips, function(x){
          return(x[subtree_parents_idx])
        })
        matrix_subtree_gnt <- do.call("rbind", list_subtree_genotypes)
        matrix_subtree_gnt_ss_split <- lapply(c(1:the$number_of_states), function(y){
          znj <- matrix_subtree_gnt == y
          mode(znj) <- "numeric"
          return(znj)
        })

        subtree_likelihoods <- matrix(0,
                                      nrow = nrow(likelihoods[[names(last_node)]]),
                                      ncol = ncol(likelihoods[[names(last_node)]]))
        for(i in the$all_possible_states){
          subtree_likelihoods <- subtree_likelihoods + list_of_subtree_sol_loglik_per_ss[[i]][,subtree_parents,drop=FALSE] %*% t(matrix_subtree_gnt_ss_split[[i]][,subtree_parents,drop=FALSE])
        }
        full_likelihood <- subtree_likelihoods + likelihoods[[names(last_node)]]

        # add solution to list of solutions
        max_lik_ss <- get_max_lik_per_ss(lik_mtx = full_likelihood)
        list_of_subtree_sol_name_per_ss <- mapply(FUN = add_values_to_end,
                                                  current = list_of_subtree_sol_name_per_ss,
                                                  new = max_lik_ss,
                                                  item = rep(1, length(the$number_of_states)),
                                                  subtree_name = rep(last_node, length(the$number_of_states)),
                                                  SIMPLIFY = FALSE)
        list_of_subtree_sol_loglik_per_ss <- mapply(FUN = add_values_to_end,
                                                    current = list_of_subtree_sol_loglik_per_ss,
                                                    new = max_lik_ss,
                                                    item = rep(2, length(the$number_of_states)),
                                                    subtree_name = rep(last_node, length(the$number_of_states)),
                                                    SIMPLIFY = FALSE)
      }
    }

    top_per_site_val <- apply(full_likelihood, MARGIN = 1, FUN = max)
    raw_rowsum <- rowSums(full_likelihood)
    max_subtracted <- sweep(full_likelihood, 1, top_per_site_val)
    max_prob <- exp(max_subtracted) / apply(exp(max_subtracted), 1, sum)

    top_per_site_idx <- apply(full_likelihood, MARGIN = 1, FUN = which.max)
    the$top_per_site_name <- colnames(full_likelihood)[top_per_site_idx]

    top_per_site_val_pr <- apply(max_prob, 1, max)

    the$list_of_subtree_sol_name_per_ss <- list_of_subtree_sol_name_per_ss

    split_pat <- "(?<=\\.).*(?=:)"
    winning_ss_root <- as.numeric(regmatches(the$top_per_site_name,
                                             gregexpr(split_pat, the$top_per_site_name, perl=TRUE)))

    reconstructed_assignments <- mapply(reconstruct_recursive,
                                        node = rep(full_tree$edge[1,1], nsites),
                                        row_idx = c(1:nsites),
                                        node_ss = winning_ss_root,
                                        SIMPLIFY = FALSE)

    cleaned_up <- mapply(clean_up,
                         items = reconstructed_assignments,
                         root_ss = winning_ss_root,
                         SIMPLIFY = FALSE)

    mle_tmp <- data.table::data.table(chr_pos = site_names)
    mle_tmp[, c("Chrom", "Pos") := tstrsplit(chr_pos, "_")]
    mle_tmp[, Pos := as.integer(Pos)]
    mle_tmp[, branch := unlist(cleaned_up)]
    mle_tmp[, c("only_branch", "starting_state", "branch_state") := tstrsplit(branch, "\\.|:")]
    mle_tmp[, starting_state := as.numeric(starting_state)]
    mle_tmp[, probability := top_per_site_val_pr]
    mle_tmp[, em_round := while_counter]

    mle_tmp[, gamma := gammas]
    mle_tmp[, delta := deltas]

    save_data[[while_counter]] <- mle_tmp
    # re-estimate the gamma and delta from current best solution
    if(is.null(matrix_list)){
      gd_val <- get_gamma_delta(assignments = mle_tmp,
                                full_tree = full_tree,
                                meth_matrix = meth_mtx,
                                unmeth_matrix = unmeth_mtx,
                                min_tips_in_s = min_tips_in_state)
    }else{
      gd_val <- get_gamma_delta(assignments = mle_tmp,
                                full_tree = full_tree,
                                matrix_list = matrix_list,
                                min_tips_in_s = min_tips_in_state)
    }

    g_tmp <- gd_val$gamma
    g_tmp[is.na(g_tmp)] <- gamma_prior
    g_tmp[g_tmp < min_gamma] <- min_gamma
    g_tmp[g_tmp > max_gamma] <- max_gamma

    d_tmp <- gd_val$delta
    d_tmp[is.na(d_tmp)] <- delta_prior
    d_tmp[d_tmp < min_delta] <- min_delta
    d_tmp[d_tmp > max_delta] <- max_delta

    gammas <- g_tmp
    deltas <- d_tmp

    not_null_items <- Filter(Negate(is.null), save_data)

    if(length(not_null_items) > 1 & em_rounds > 1){
      # if stable - add to final
      extract_stable <- lapply(mle_tmp$chr_pos, function(asite){
        assign_all <- lapply(not_null_items, function(onerun){
          return(onerun[chr_pos == asite,])
        })
        dt_assign_all <- data.table::rbindlist(assign_all)
        last <- dt_assign_all[which.max(em_round),][["branch"]]
        if(is.list(last)){
          last <- unlist(last)
        }
        intersecting <- last %in% dt_assign_all[!which.max(em_round),][["branch"]]
        max_item <- NULL
        if(intersecting){
          # add to final set
          dt_assign_all_int <- dt_assign_all[branch == last,]
          max_item <- dt_assign_all_int[which.max(gamma+delta),]
        }
        return(max_item)
      })
      final_sites <- Filter(Negate(is.null), extract_stable)
      final_sites_dt <- data.table::rbindlist(final_sites)
      final_list[[while_counter]] <- final_sites_dt
      full_final <- data.table::rbindlist(final_list)
      # if unstable - redo probabilities
      subchunk <- dt_chunk_sites[!(chr_pos %in% full_final$chr_pos),]
      if(nrow(subchunk) == 0){
        while_counter <- em_rounds + 1 # so that you leave loop
      }else{
        # Fix inputs
        if(the$number_of_states < 4){
          meth_mtx <- meth_mtx[rownames(meth_mtx) %in% subchunk$chr_pos,,drop=FALSE]
          unmeth_mtx <- unmeth_mtx[rownames(unmeth_mtx) %in% subchunk$chr_pos,,drop=FALSE]
        }else{
          phased_reads_list2 <- lapply(matrix_list, function(oi){
            return(oi[rownames(oi) %in% subchunk$chr_pos,,drop=FALSE])
          })
          matrix_list <- phased_reads_list2
        }
        gammas <- gammas[names(gammas) %in% subchunk$chr_pos]
        deltas <- deltas[names(deltas) %in% subchunk$chr_pos]
      }
    }else if(em_rounds == 1){
      final_list[[1]] <- mle_tmp
    }

    while_counter <- while_counter + 1
  }
  from_final <- Filter(Negate(is.null), final_list)
  if(length(from_final) > 0){
    full_final <- data.table::rbindlist(final_list)
  }

  check_done <- match(original_sites$chr_pos, full_final$chr_pos)
  if(any(is.na(check_done))){
    # find missing sites
    dt_missed <- original_sites[!(chr_pos %in% full_final$chr_pos),]
    # extract from save_data
    not_null_items <- Filter(Negate(is.null), save_data)
    dt_not_null <- data.table::rbindlist(not_null_items)
    # take most conservative
    dt_missed_prelim <- dt_not_null[chr_pos %in% dt_missed$chr_pos,]
    dt_missed_prelim1 <- dt_missed_prelim[dt_missed_prelim[, .I[which.max(gamma+delta)], by=chr_pos]$V1]
    full_final2 <- rbind(full_final, dt_missed_prelim1)
    full_final <- full_final2
  }
  if(nrow(full_final) > nrow(original_sites)){
    stop("Extra rows???")
  }
  converted_dt <- data.table::merge.data.table(full_final, original_sites, by = "chr_pos")
  converted_dt <- converted_dt[order(Chrom, Pos),]
  converted_dt[, n_assignments := sapply(only_branch, function(x){
    y <- strsplit(x, "_")[[1]]
    if(all(y == "n")){
      return(0)
    }else{
      return(length(y))
    }
  })]
  if(return_intermediate_assignments){
    return(list("converted_dt" = converted_dt,
                "intermediate" = Filter(Negate(is.null), save_data)))
  }else{
    return(converted_dt)
  }
}
