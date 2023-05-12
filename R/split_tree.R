bfs_node_order <- function(tree){
  root <- tree$edge[1,1]

  explored <- vector("numeric", length = max(tree$edge))
  counter <- 1
  explored[root] <- counter
  counter <- counter + 1

  level <- rep(NA, max(tree$edge))
  level[root] <- 0

  first_nodes <- tree$edge[tree$edge[,1] == root, 2]
  q <- first_nodes
  while(length(q) > 0){
    v <- q[1] # dequeue vertex
    q <- q[-1]
    explored[v] <- counter
    counter <- counter + 1
    v_parent <- tree$edge[,1][tree$edge[,2] == v] # find parent of vertex
    if(is.na(level[v])){
      level[v] <- level[v_parent] + 1
    }
    v_children <- tree$edge[,2][tree$edge[,1] == v]
    if(length(v_children) > 0){
      q[((length(q) + 1):(length(q)+length(v_children)))] <- v_children
    }
  }
  explored <- order(explored)
  level <- level[explored]
  return(list("explored" = explored,
              "levels" = level))
}

get_num_descendants <- function(node_, tree = tree){
  length(phangorn::Descendants(tree, node = node_, type = "tips")[[1]])
}

generate_tree_partition <- function(tree,
                                    cutoff_upper = 50,
                                    cutoff_lower = 5,
                                    original_tree = NULL,
                                    custom_nodes = NULL){
  tree_root <- tree$edge[1,1]

  # breadth first search to get list of node order
  r1 <- bfs_node_order(tree)
  node_order <- r1$explored
  # node_level <- r1$levels
  # node_level_ordered <- node_level[node_order]

  num_descendants <- sapply(node_order, get_num_descendants, tree = tree)
  num_desc_idx <- sapply(1:max(tree$edge), get_num_descendants, tree = tree)

  # visit each node and accept if number of descendents lower than cutoff_upper
  used_tips <- length(tree$tip.label)
  picked_nodes <- c()
  nodes_considered <- node_order
  while(used_tips > 0 & length(nodes_considered) > 0){
    node1 <- nodes_considered[1]
    nodes_considered <- nodes_considered[-1]
    n_des <- num_desc_idx[node1]
    if(n_des > cutoff_upper){ # if node too big, try next

    }else if(n_des <= cutoff_upper){ # node below cutoff_upper, can add to list of picked
      picked_nodes <- c(picked_nodes, node1)
      used_tips <- used_tips - n_des # update the num of nodes that are in use
      node_desc <- phytools::getDescendants(tree, node = node1) # remove descendants of picked node from the list of nodes to consider
      nodes_considered <- nodes_considered[!nodes_considered %in% node_desc]
    }
  }

  # extract subtrees
  if(is.null(custom_nodes)){
    internal_nodes <- picked_nodes[picked_nodes >= tree_root & num_desc_idx[picked_nodes] > cutoff_lower]
    if(length(internal_nodes) == 0){
      tree_unit <- tree
      tree_unit$edge.length <- rep(1, length(tree_unit$edge.length))
      plot.phylo(tree_unit)
      stop(sprintf("Cannot create tree - modify cutoffs. Upper: %i, lower: %i\n", cutoff_upper, cutoff_lower ))
    }
  }else{
    internal_nodes <- as.integer(custom_nodes)
  }

  subtrees <- lapply(internal_nodes, function(node_){
    sub_tree <- ape::extract.clade(tree, node = node_)
    return(sub_tree)
  })
  names(subtrees) <- as.character(internal_nodes)

  # prepare root structure
  terminal_nodes <- picked_nodes[picked_nodes < tree_root]

  tree_df_orig <- tidytree::as_tibble(tree)
  tree_df <- tree_df_orig
  tree_df$label <- ifelse(grepl("PD|[[:alpha:]]", tree_df$label), tree_df$label, tree_df$node)
  for(i in 1:length(internal_nodes)){
    cur_node <- internal_nodes[i]
    all_desc <- phangorn::Descendants(tree, node = cur_node, type = "all")
    tree_df <- tree_df[!(tree_df$node %in% all_desc),] # remove all descendants of subtrees from tree_df
  }
  # renumber nodes
  # find terminal nodes which are numeric - they contain a subtree that already exists
  nested_root <- tree$tip.label[as.numeric(grep("[[:alpha:]]", tree_df$label, value = TRUE, invert = TRUE))]
  nested_root <- nested_root[!is.na(nested_root)]
  new_root <- length(internal_nodes) + sum(grepl("^[[:alpha:]]", tree_df$label)) + length(nested_root) + 1

  # first find internal nodes of root structure - renumber them from root # onwards
  root_structure_internal <- tree_df$node[tree_df$node >= tree_root]
  root_structure_internal <- root_structure_internal[!root_structure_internal %in% internal_nodes]
  names(root_structure_internal) <- new_root:(new_root + length(root_structure_internal) - 1)

  # find new terminal nodes of root structure
  root_structure_terminal <- tree_df$node[tree_df$node < tree_root | tree_df$node %in% internal_nodes]
  old_dfs_edge_order <- match(root_structure_terminal, tree$edge[,2])
  names(root_structure_terminal) <- old_dfs_edge_order
  root_structure_terminal <- root_structure_terminal[order(as.numeric(names(root_structure_terminal)))]
  names(root_structure_terminal) <- 1:length(root_structure_terminal)

  pairs_new <- paste0(tree_df$parent, "_", tree_df$node)

  full_new_nodes <- c(root_structure_terminal, root_structure_internal)
  tree_df$parent_new <- 0
  tree_df$node_new <- 0
  for(i in 1:nrow(tree_df)){
    tree_df$parent_new[i] <- as.numeric(names(which(full_new_nodes == tree_df$parent[i])))
    tree_df$node_new[i] <- as.numeric(names(which(full_new_nodes == tree_df$node[i])))
  }

  tree_df$pairs_old <- paste0(tree_df$parent, "_", tree_df$node)
  tree_df$pairs_new <- paste0(tree_df$parent_new, "_", tree_df$node_new)

  pairs_old <- paste0(tree$edge[,1], "_", tree$edge[,2])

  ordering <- pairs_new[match(pairs_old, pairs_new)] # find pairs of root tree in the original tree
  ordering <- ordering[!is.na(ordering)]
  #ordering <- pairs_old

  best_yet <- match(ordering, tree_df$pairs_old)
  root_row <- which(is.na(tree_df$branch.length))
  tree_df <- tree_df[c(best_yet, root_row),]

  tree_df_copy <- tree_df
  tree_df_copy$parent <- tree_df_copy$parent_new
  tree_df_copy$node <- tree_df_copy$node_new
  tree_df_copy <- tree_df_copy[order(tree_df_copy$node),]

  new_root_phylo <- as.phylo(tree_df_copy)

  pairs_extra_new <- paste0(new_root_phylo$edge[,1], "_", new_root_phylo$edge[,2])
  idx <- match(tree_df$pairs_new, pairs_extra_new)
  idx <- idx[!is.na(idx)]
  new_root_phylo$edge <- new_root_phylo$edge[idx,]
  new_root_phylo$edge.length <- new_root_phylo$edge.length[idx]

  if(!is.null(original_tree)){ #if original tree is provided - map to original numbering
    nonPD_tips <- as.numeric(grep("[[:alpha:]]", new_root_phylo$tip.label, value = TRUE, invert = TRUE))
    # tree already had subtree at that node
    already_subtree <- tree$tip.label[nonPD_tips]
    tmp_idx <- which(!is.na(already_subtree))
    wrong_node_name <- nonPD_tips[tmp_idx]
    correct_node_name <- already_subtree[tmp_idx]
    new_root_phylo$tip.label[match(as.numeric(wrong_node_name), new_root_phylo$tip.label)] <- correct_node_name

    # tip was chosen for subtree within this call to function
    # go through all trees, take all descendants (and also descendants of nested trees)
    # then check which node in original tree has node which has identical descendants

    # the following might be problematic - different way to remap?
    nodes_to_rename <- c(names(subtrees), new_root_phylo$node.label)
    is_subtree <- c(rep(TRUE, length(subtrees)), rep(FALSE, length(new_root_phylo$node.label)))
    new_subtrees <- subtrees
    for(j in 1:length(nodes_to_rename)){
      if(is_subtree[j]){
        current_name <- names(subtrees)[j]
        all_cur_tips <- subtrees[[j]]$tip.label
      }else{
        current_name <- nodes_to_rename[j]
        number_of_node_rel <- which(current_name == new_root_phylo$node.label)
        number_of_node_abs <- new_root_phylo$edge[1,1] + number_of_node_rel - 1
        all_cur_tips_nums <- phangorn::Descendants(new_root_phylo, node = number_of_node_abs, type = "tips")[[1]]
        all_cur_tips <- new_root_phylo$tip.label[all_cur_tips_nums]
      }
      all_descendants <- lapply(all_cur_tips, function(node){
        if(grepl("[[:alpha:]]", node)){
          return(node)
        }else{
          children <- phangorn::Descendants(x = original_tree,
                                  node = as.numeric(node),
                                  type = "tips")[[1]]
          PD_children <- original_tree$tip.label[children]
          return(PD_children)
        }
      })
      unlisted_desc <- unlist(all_descendants)
      # find one nodes path
      path_to_root <- phangorn:: Ancestors(original_tree,
                                           node = which(original_tree$tip.label == unlisted_desc[1]),
                                           type = "all")
      # concatenate all descendants
      local_desc <- phangorn::Descendants(original_tree,
                                          path_to_root,
                                          type = "tips")

      # find node in original tree where the match is 100%
      corresponding_node <- sapply(local_desc, function(ld){
        if(length(unlisted_desc) == length(ld)){
          all(sort(unlisted_desc) == sort(original_tree$tip.label[ld]))
        }else{
          return(FALSE)
        }
      })
      stopifnot(sum(corresponding_node) == 1)

      if(is_subtree[j]){
        names(new_subtrees)[j] <- path_to_root[corresponding_node]
        new_root_phylo$tip.label[new_root_phylo$tip.label == current_name] <- path_to_root[corresponding_node]
      }else{
        new_root_phylo$node.label[new_root_phylo$node.label == current_name] <- path_to_root[corresponding_node]
      }
    }
  }else{
    new_subtrees <- subtrees
  }

  # return all trees
  tmp_all_trees <- c(new_subtrees, "root" = list(new_root_phylo))
  return(tmp_all_trees)
}

get_tree_ntips <- function(x){
  length(x$tip.label)
}

find_orig_node_in_current <- function(full_tree, original_node, current_tree){
  descendants <- phangorn::Descendants(full_tree, node = original_node, type = "tips")[[1]]
  named_desc <- full_tree$tip.label[descendants]
  named_desc_col <- paste0(named_desc, collapse = ":")

  internal <- current_tree$edge[,2]
  internal <- internal[current_tree$edge[1,1] <= internal]
  int_desc <- lapply(internal, function(x){
    dd <- phangorn::Descendants(current_tree, node = x, type = "tips")[[1]]
    nd <- as.list(current_tree$tip.label[dd])
    to_replace_idx <- grep("[[:alpha:]]", nd, invert = TRUE)
    to_replace_num <- nd[to_replace_idx]
    new_ones <- lapply(to_replace_num, function(y){
      dorig <- phangorn::Descendants(full_tree, node = as.integer(y), type = "tips")[[1]]
      dorig_n <- full_tree$tip.label[dorig]
      return(dorig_n)
    })
    nd[to_replace_idx] <- new_ones
    return(paste0(unlist(nd), collapse = ":"))
  })
  new_node <- match(named_desc_col, int_desc)
  return(internal[new_node])
}

divide_tree_into_subtrees <- function(full_tree,
                                      min_subtree_ = 5,
                                      max_subtree_ = 20,
                                      custom_nodes = NULL){
  tmp_list <- list("root"=full_tree)
  nest_level <- list()

  if(is.null(custom_nodes)){
    subtree_sizes <- sapply(tmp_list, get_tree_ntips)
    tree_correct_size <- subtree_sizes >= min_subtree_ & subtree_sizes <= max_subtree_

    while_counter <- 1
    while(!all(tree_correct_size)){
      problematic_tree <- tmp_list[[which(tree_correct_size == FALSE)]]
      sub_root_trees <- generate_tree_partition(tree = problematic_tree,
                                                cutoff_upper = max_subtree_,
                                                cutoff_lower = min_subtree_,
                                                original_tree = full_tree)
      tmp_list[[which(tree_correct_size == FALSE)]] <- NULL
      tmp_list <- c(tmp_list, sub_root_trees)
      nest_level[[while_counter]] <- names(sub_root_trees)
      tmp_subtree_nodes_id <- paste0(names(tmp_list), collapse = "_")

      sub_root_trees_sizes <- sapply(tmp_list, get_tree_ntips)
      tree_correct_size <- sub_root_trees_sizes >= min_subtree_ & sub_root_trees_sizes <= max_subtree_
      while_counter <- while_counter + 1
      if(while_counter > 20){
        stop("Stuck in while loop?")
      }
    }
  }else{
    custom_order <- rev(sort(custom_nodes))
    atree <- full_tree
    for(i in 1:length(custom_order)){
      sub_root_trees2 <- generate_tree_partition(tree = atree,
                                                 cutoff_upper = max_subtree_,
                                                 cutoff_lower = min_subtree_,
                                                 original_tree = full_tree,
                                                 custom_nodes = custom_order[i])
      atree <- sub_root_trees2$root
      if(!is.na(custom_order[i+1])){ # next element exists
        # find original node in current root tree
        custom_order[i+1] <- find_orig_node_in_current(full_tree = full_tree,
                                                       original_node = custom_order[i+1],
                                                       current_tree = atree)
      }
      tmp_list[["root"]] <- NULL
      tmp_list <- c(tmp_list, sub_root_trees2)
      nest_level[[i]] <- names(sub_root_trees2)
    }
    tmp_subtree_nodes_id <- paste0(names(tmp_list), collapse = "_")
  }

  all_trees <- tmp_list

  return(list("all_trees" = all_trees,
              "nest_level" = nest_level))
}
