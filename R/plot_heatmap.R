get_df_depth_mm <- function(mtx_list, chr_pos){
  depth_mtx <- mtx_list$meth_mtx + mtx_list$unmeth_mtx
  mm_mtx <- mtx_list$meth_mtx

  depth_df <- depth_mtx[chr_pos,,drop=FALSE]
  meth_df <- mm_mtx[chr_pos,,drop=FALSE]

  stopifnot(all(colnames(depth_df) == colnames(meth_df)))

  hello <- data.frame(label = colnames(depth_df),
                      id = colnames(depth_df),
                      meth_r = meth_df[1,],
                      depth = depth_df[1,])
  return(hello)
}



#' Plot heatmap of assignments on sites
#'
#' @param assignment character: an assignment, branch.start_state:state_changes
#' @param full_dt data.table of all assignments
#' @param tree "phylo" object of tree
#' @param meth_rate_mtx matrix of methylation rates across sites and tips
#' @param include_tip_labels logical whether to include tip labels
#' @param max_pos integer: maximum number of positions to plot
#' @param specific_positions character vector of "chr_pos" sites to show
#' @param df_col_meta data frame of colony metadata
#' @param custom_title character
#' @param custom_subtitle character
#' @param highlight_additional_nodes numeric vector of nodes you wish to highlight in addition, will be coloured grey
#' @param show_bar logical: if you would like to plot a bar chart for just 1 site, requires meth_reads_list
#' @param meth_reads_list list of methlyated and unmethylated read counts
#' @param meth_reads_df data frame of methlyated and unmethylated read counts
#' @param fill_by_exp_tip_state logical - when plotting single sites - fill the bar by expected tip state (alternative is filling by methylation rate)
#' @param offset_by numeric distance between phylogeny and heatmap
#' @param sort_by_gd character ("min", "max" or "no") to sort sites by sum of gamm and delta
#' @param change_names named vector of all transitions, where the names are text that will replace corresponding transition
#' @param highlight_extend how much to extend the node highlight by (default 20)
#' @param heatmap_width width of heatmap relative to tree, default is 1
#'
#' @return ggplot plot object
#' @export
#'
#' @examples
#' \dontrun{
#'   hm_plot <- plot_heatmap(assignment = "n.3:9",
#'   full_dt = assignments_full[rank == 1,],
#'   meth_rate_mtx = rate_mtx,
#'   tree = tree)
#' }
plot_heatmap <- function(assignment,
                         full_dt,
                         tree,
                         meth_rate_mtx,
                         include_tip_labels = TRUE,
                         max_pos = 30,
                         specific_positions = NULL,
                         df_col_meta = NULL,
                         custom_title = NULL,
                         custom_subtitle = NULL,
                         highlight_additional_nodes = NULL,
                         show_bar = FALSE,
                         meth_reads_list = NULL,
                         meth_reads_df = NULL,
                         fill_by_exp_tip_state = TRUE,
                         offset_by = 30,
                         highlight_extend = 20,
                         sort_by_gd = "no",
                         change_names = NULL,
                         heatmap_width = 1){
  branch = chr_pos = chrom_fct = Chrom = Affected = color = edge_text = NULL
  text1 = label = meth_r = depth = rate = gamma = delta = inferred_s =Pos = NULL

  full_dt_subset <- full_dt[branch == assignment,]

  only_branch <- unique(full_dt_subset[["only_branch"]])
  only_starting_state <- unique(full_dt_subset[["starting_state"]])
  only_states <- unique(full_dt_subset[["branch_state"]])

  if(is.null(specific_positions)){
    specific_positions <-  full_dt_subset[["chr_pos"]]
  }

  if(the$number_of_states == 2){
    state_explain <- "1: 0->0, 2: 0->1, 3: 1->0, 4: 1->1"
  }else if(the$number_of_states == 3){
    state_explain <- "1: 0->0, 2: 0->0.5, 3: 0->1, 4: 0.5-> 0, 5: 0.5->0.5, 6: 0.5->1, 7: 1->0, 8: 1->0.5, 9: 1->1"
  }else if(the$number_of_states == 4){
    state_explain <- "Stable states: MP: 1:00, 6:01, 11:10, 16:11\nChange types: MP: 2: 00->01, 3: 00->10, 4: 00->11, 5:01->00, 7:01->10, 8:01->11,\n9:10->00, 10:10->01, 12:10->11, 13:11->00, 14:11->01, 15:11->10"
  }

  if(is.null(custom_subtitle)){
    if(the$useQ){
      custom_subtitle <- sprintf("Rates: used Q matrix. Transitions:\n%s",
                                 state_explain)
    }else{
      custom_subtitle <- sprintf("Rates: forward: %1.2e, backward: %1.2e. Transitions:\n%s",
                                 the$epi_fwd, the$epi_bck,
                                 state_explain)
    }
  }
  if(sort_by_gd == "no"){
    if(length(specific_positions) > max_pos){
      specific_positions <- sample(specific_positions,
                                   size = min(max_pos,
                                              length(specific_positions)),
                                   replace = FALSE)
      full_dt_subset <- full_dt_subset[chr_pos %in% specific_positions]
    }
  }else{
    if(!all(c("gamma", "delta") %in% colnames(full_dt_subset))){
      stop("To plot by min/max sum of gamma & delta, both columns gamma and delta must be present in the input data frame")
    }
    if(sort_by_gd == "min"){
      full_dt_subset_min <- full_dt_subset[order(gamma+delta),][1:min(max_pos, length(specific_positions)),]
      full_dt_subset <- full_dt_subset_min
    }else if(sort_by_gd == "max"){
      full_dt_subset_max <- full_dt_subset[order(gamma+delta, decreasing = TRUE),][1:min(max_pos, length(specific_positions)),]
      full_dt_subset <- full_dt_subset_max
    }
    specific_positions <- full_dt_subset$chr_pos
  }

  # reorder positions by chromosome
  chrom_levels <- c(as.character(1:22), "X", "Y")

  full_dt_subset <- full_dt_subset[, chrom_fct := factor(Chrom, levels = chrom_levels, ordered = TRUE)]

  full_dt_subset <- full_dt_subset[order(chrom_fct, Pos),]
  if(is.matrix(meth_rate_mtx)){
    pos_rates <- convert_to_matrix(meth_rate_mtx, full_dt_subset[["chr_pos"]])
    pos_rates <- pos_rates[,match(tree$tip.label, colnames(pos_rates)), drop=FALSE]
    stopifnot(all(tree$tip.label == colnames(pos_rates)))
    pos_rates_flip <- t(pos_rates)
  }else if(is.list(meth_rate_mtx)){
    pos_rates_list <- lapply(meth_rate_mtx, convert_to_matrix, chr_pos = full_dt_subset[["chr_pos"]])
    pos_rates_list <- lapply(pos_rates_list, function(one){
      two <- one[, match(tree$tip.label, colnames(one)), drop = FALSE]
      stopifnot(all(tree$tip.label == colnames(two)))
      return(t(two))
    })
    pos_rates_mtx <- do.call("cbind", pos_rates_list)
    eg <- expand.grid(full_dt_subset[["chr_pos"]], names(pos_rates_list))
    new_colnames <- apply(eg, 1, function(n){
      n <- rev(n)
      paste0(n, collapse = ":")
    })
    colnames(pos_rates_mtx) <- new_colnames
    # new_order <- lapply(specific_positions, function(p){
    #   grep(p, colnames(pos_rates_mtx), value = TRUE)
    # })
    # pos_rates_mtx <- pos_rates_mtx[,unlist(new_order)]
    pos_rates_flip <- pos_rates_mtx
  }
  pos_rates_flip[is.nan(pos_rates_flip)] <- NA
  if(sort_by_gd %in% c("min", "max")){
    pos_rates_flip <- pos_rates_flip[,specific_positions,drop=FALSE]
  }
  if(the$number_of_states == 2){
    conversion <- c(0, 1)
    names(conversion) <- as.character(c(1:2))
  }else if(the$number_of_states == 3){
    conversion <- c(0, 0.5, 1)
    names(conversion) <- as.character(c(1:3))
  }else if(the$number_of_states == 4){
    conversion <- c("00", "01", "10", "11")
    names(conversion) <- as.character(c(1:4))
  }
  ss <- conversion[as.character(only_starting_state) == names(conversion)]
  plot_assignment_vector_char <- generate_expected_tip_states(ss = only_starting_state,
                                                          changed_branches = only_branch,
                                                          change_types = only_states,
                                                          tree = tree)

  plot_assignment_vector <- conversion[plot_assignment_vector_char]
  names(plot_assignment_vector) <- tree$tip.label
  # plot_assignment_vector <- plot_assignment_vector[names(plot_assignment_vector) %in% tree$tip.label]
  plot_assignment_df <- data.frame(row.names = names(plot_assignment_vector),
                                   label = names(plot_assignment_vector),
                                   Affected = factor(plot_assignment_vector, levels = conversion, ordered = TRUE))
  tb <- tidytree::as_tibble(tree)
  tb2 <- full_join(tb, plot_assignment_df, by = "label")

  new_tree <- tidytree::as.treedata(tb2)

  if(the$number_of_states == 2){
    affected_tip_colors <- c("#313695", "#a50026")
    names(affected_tip_colors) <- as.factor(c(0, 1))

    edge_change_colours <- c("#d73027", "#4575b4")
    names(edge_change_colours) <- as.factor(c(2, 3))
  }else if(the$number_of_states == 3){
    affected_tip_colors <- c("#313695","#fdae61", "#a50026")
    names(affected_tip_colors) <- as.factor(c(0, 0.5, 1))

    edge_change_colours <- c("#fdae61","#d73027", "#abd9e9", "#fdae61",
                             "#4575b4","#abd9e9", "##808080")
    names(edge_change_colours) <- as.factor(c(2, 3, 4, 6, 7, 8, 9))

  }else if(the$number_of_states == 4){
    affected_tip_colors <- c("#313695","#fdae61", "#fb9a99", "#a50026")
    names(affected_tip_colors) <- as.factor(c("00", "01", "10", "11"))

    edge_change_colours <- c("#fdae61", "#fb9a99", "#d73027",
                             "#cab2d6", "#bdbdbd", "#fb9a99",
                             "#abd9e9", "#bdbdbd", "#fdae61",
                             "#4575b4", "#abd9e9", "#cab2d6")
    names(edge_change_colours) <- as.factor(c(2, 3, 4,
                                              5, 7, 8,
                                              9, 10, 12,
                                              13, 14, 15))
  }
  #max_x <- ifelse(donor_id == "KX004", ifelse(include_tip_labels, 1600, 1500), 80)
  p1 <- ggtree::ggtree(new_tree, ladderize = FALSE) +
    ggtree::geom_tippoint(aes(color=Affected)) +
    ggtree::scale_color_manual(values = affected_tip_colors, drop=FALSE) +
    ggplot2::guides(color=ggplot2::guide_legend(title=sprintf("Expected\ntip state%s",
                                                              ifelse(the$number_of_states == 4,
                                                                     " (MP)",
                                                                     ""))))
  #p1 <- p1 + ggtree::theme_tree2()
  if(include_tip_labels){
    p1 <- p1 +
      ggtree::geom_tiplab(size = 2, hjust = -0.2)
  }
  edges <- strsplit(only_branch, "_")[[1]]
  edge_changes <- strsplit(only_states, "_")[[1]]
  if(the$number_of_states == 3){
    edge_changes <- edge_changes[!(edge_changes %in% c("1", "5", "9"))]
  }else if(the$number_of_states == 4){
    edge_changes <- edge_changes[!(edge_changes %in% c("1", "6", "11", "16"))]
  }

  if(only_branch != "n"){
    nodes <- tree$edge[as.numeric(edges),2]
    p1 <- p1 +
      ggtree::geom_hilight(node = nodes,
                   fill = edge_change_colours[edge_changes],
                   extend = highlight_extend,
                   alpha = 0.4)
  }
  if(!is.null(highlight_additional_nodes)){
    p1 <- p1 +
      ggtree::geom_hilight(node = highlight_additional_nodes,
                           fill = "#808080",
                           extend = highlight_extend,
                           alpha = 0.4)
  }
  if(is.null(custom_title)){
    custom_title <-sprintf("Assigned to %s;\n%s transition(s): %s, respectively;\nstarting state%s: %s",
                           paste0(edges, collapse = ", "),
                           paste0(ifelse(only_branch == "n", 0, length(edges))),
                           paste0(strsplit(only_states, "_")[[1]], collapse = ", "),
                           ifelse(the$number_of_states == 4, " (MP)", ""),
                           as.character(ss))

    if(!is.null(change_names)){
      just_trans_txt <- strsplit(only_states, "_")[[1]]
      if("n" %in% just_trans_txt){
        warning("Cannot rename transitions because there's a character present, must be all numeric")
      }else{
        just_trans_num <- as.numeric(just_trans_txt)
        renamed_trans <- names(change_names)[match(just_trans_num, change_names)]
        custom_title <-sprintf("Assigned to %s;\n%s transition(s): %s, respectively;\nstarting state%s: %s",
                               paste0(edges, collapse = ", "),
                               paste0(ifelse(only_branch == "n", 0, length(edges))),
                               paste0(renamed_trans, collapse = ", "),
                               ifelse(the$number_of_states == 4, " (MP)", ""),
                               as.character(ss))
      }
    }
  }

  if(only_branch != "n" | !is.null(highlight_additional_nodes)){
    p1_tmp <- p1
    if(!exists("nodes")){
      nodes <- c()
    }
    if(!exists("highlight_additional_nodes")){
      highlight_additional_nodes <- c()
    }
    all_nodes <- c(highlight_additional_nodes, nodes)
    associated_edges <- match(all_nodes, tree$edge[,2])
    names(all_nodes) <- as.character(associated_edges)
    #p1_tmp$data$text1 <- names(all_nodes)[match(p1_tmp$data$node, all_nodes)]
    #p1 <- p1_tmp + ggtree::geom_nodelab(geom='label', aes(label = text1), node = "all", size = 2)
    p1_tmp$data$edge_text <- names(all_nodes)[match(p1_tmp$data$node, all_nodes)]
    p1 <- p1_tmp  + ggtree::geom_label2(aes(x=branch, label=edge_text), size = 2) #, subset = !is.na(p1_tmp$data$edge_text)) #geom='label', aes(label = text1), node = "all", size = 2)
  }

  if(show_bar){
    if(nrow(full_dt_subset) > 1){
      cat("Will sample random site to plot\n")
      full_dt_subset <- full_dt_subset[sample(1:nrow(full_dt_subset), size = 1)]
    }

    if(is.null(meth_reads_list) & is.null(meth_reads_df)){
      stop("Must provide dataframe or list of methylated and unmethylated read counts")
    }else{
      if(!is.null(meth_reads_df)){
        df_bar <- meth_reads_df
        # meth_reads_list = NULL,
        # meth_reads_df = NULL
      }else{
        onesite <- full_dt_subset$chr_pos
        df_bar <- get_df_depth_mm(mtx_list = meth_reads_list,
                                  chr_pos = onesite)
      }
      df_bar <- df_bar[match(tree$tip.label, df_bar$label),]
      df_bar$label <- factor(df_bar$label,
                             levels = tree$tip.label,
                             ordered = TRUE)
      df_bar$rate <- df_bar$meth_r/df_bar$depth

      df_bar$inferred_s <- as.character(plot_assignment_vector[match(df_bar$id, names(plot_assignment_vector))])

      p2 <- ggplot2::ggplot(df_bar, aes(label, meth_r)) +
        ggplot2::geom_col() +
        ggplot2::coord_flip()
      p2_depth <- ggplot2::ggplot(df_bar, aes(label, depth)) +
        ggplot2::geom_col() +
        ggplot2::coord_flip()
      if(fill_by_exp_tip_state){
        p2_combined <- ggplot2::ggplot(df_bar, aes(x = label)) +
          ggplot2::geom_col(aes(y = depth), position = "identity", fill = "#bdbdbd") +
          ggplot2::geom_col(aes(y = meth_r, fill = inferred_s), position = "identity") +
          ggplot2::scale_fill_manual(values = affected_tip_colors) +
          ggplot2::labs(y = sprintf("%s: Read count (methylated overlaid on depth)", onesite)) +
          ggplot2::coord_flip() +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                         axis.ticks.y = ggplot2::element_blank(),
                         axis.title.y = ggplot2::element_blank(),
                         legend.position = "none")
      }else{
        p2_combined <- ggplot2::ggplot(df_bar, aes(x = label)) +
          ggplot2::geom_col(aes(y = depth), position = "identity", fill = "#bdbdbd") +
          ggplot2::geom_col(aes(y = meth_r, fill = rate), position = "identity") +
          ggplot2::scale_fill_gradient2(low = "#313695",
                                        mid = "#ffffbf",
                                        high = "#a50026",
                                        midpoint = 0.5,
                                        name = "Methylation\nratio",
                                        limits = c(0,1)) +
          ggplot2::labs(y = sprintf("%s: Read count (methylated overlaid on depth)", onesite)) +
          ggplot2::coord_flip() +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                         axis.ticks.y = ggplot2::element_blank(),
                         axis.title.y = ggplot2::element_blank())
      }


      p3_int <- p2_combined %>% aplot::insert_left(p1)

      p3 <- aplot::plot_list(p3_int) + patchwork::plot_annotation(title = custom_title,
                                                subtitle = custom_subtitle,
                                                theme = ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 16),
                                                              plot.subtitle = element_text(hjust = 0.5, size = 12)))
    }
  }else{
    p3 <- suppressMessages(ggtree::gheatmap(p1,
                                            data=pos_rates_flip,
                                            offset = offset_by,
                                            colnames_angle=45,
                                            hjust = 1,
                                            color=NA,
                                            font.size = 2,
                                            width = heatmap_width) +
                             ggtree::scale_x_ggtree() +
                             ggtree::vexpand(.1, -1) +
                             ggplot2::scale_fill_gradient2(low = "#313695",
                                                           mid = "#ffffbf",
                                                           high = "#a50026",
                                                           midpoint = 0.5,
                                                           name = "Methylation\nratio",
                                                           limits = c(0,1)) +
                             ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 16),
                                            plot.subtitle = element_text(size = 12, hjust = 0.5)) +
                             ggplot2::ggtitle(label = custom_title,
                                              subtitle = custom_subtitle))
  }


  return(p3)
}
