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


#' Title
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

  # dt_combined_all_trees <<- dt_combined

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


