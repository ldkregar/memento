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
                         custom_subtitle = NULL){
  branch = chr_pos = chrom_fct = Chrom = Pos = Affected = NULL

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
  }else{
    stop("You must call plot_heatmap_phased if you wish to plot the 4 state")
  }

  if(is.null(custom_subtitle)){
    custom_subtitle <- sprintf("Rates: forward: %1.2e, backward: %1.2e. Transitions:\n%s",
                               the$epi_fwd, the$epi_bck,
                               state_explain)
  }

  if(length(specific_positions) > max_pos){
    specific_positions <- sample(specific_positions, size = max_pos, replace = FALSE)
    full_dt_subset <- full_dt_subset[chr_pos %in% specific_positions]
  }
  # reorder positions by chromosome
  chrom_levels <- c(as.character(1:22), "X", "Y")

  full_dt_subset <- full_dt_subset[, chrom_fct := factor(Chrom, levels = chrom_levels, ordered = TRUE)]

  full_dt_subset <- full_dt_subset[order(chrom_fct, Pos),]

  pos_rates <- convert_to_matrix(meth_rate_mtx, full_dt_subset[["chr_pos"]])
  pos_rates <- pos_rates[,match(colnames(pos_rates), tree$tip.label), drop=FALSE]
  pos_rates_flip <- t(pos_rates)
  if(the$number_of_states == 3){
    conversion <- c(0, 0.5, 1)
    names(conversion) <- as.character(c(1:3))
  }else if(the$number_of_states == 2){
    conversion <- c(0, 1)
    names(conversion) <- as.character(c(1:2))
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
                                   Affected = as.factor(plot_assignment_vector))
  tb <- tidytree::as_tibble(tree)
  tb2 <- full_join(tb, plot_assignment_df, by = "label")

  new_tree <- tidytree::as.treedata(tb2)

  if(the$number_of_states == 2){
    affected_tip_colors <- c("#313695", "#a50026")
    names(affected_tip_colors) <- as.factor(c(0, 1))
  }else if(the$number_of_states == 3){
    affected_tip_colors <- c("#313695","#fdae61", "#a50026")
    names(affected_tip_colors) <- as.factor(c(0, 0.5, 1))
  }

  edge_change_colours <- c("#fdae61","#d73027", "#abd9e9", "#fdae61", "#4575b4","#abd9e9", "##808080")
  names(edge_change_colours) <- as.factor(c(2, 3, 4, 6, 7, 8, 9))

  #max_x <- ifelse(donor_id == "KX004", ifelse(include_tip_labels, 1600, 1500), 80)
  color = NULL
  p1 <- ggtree::ggtree(new_tree, ladderize = FALSE) +
    ggtree::geom_tippoint(aes(color=Affected)) +
    ggtree::scale_color_manual(values = affected_tip_colors) +
    ggplot2::guides(color=ggplot2::guide_legend(title="Expected\ntip state"))

  if(include_tip_labels){
    p1 <- p1 +
      ggtree::geom_tiplab(size = 2, hjust = -0.2)
  }
  edges <- strsplit(only_branch, "_")[[1]]
  edge_changes <- strsplit(only_states, "_")[[1]]
  edge_changes <- edge_changes[edge_changes != "9" & edge_changes != "5" & edge_changes!= "1"]
  if(only_branch != "n"){
    nodes <- tree$edge[as.numeric(edges),2]
    p1 <- p1 +
      ggtree::geom_hilight(node = nodes,
                   fill = edge_change_colours[edge_changes],
                   extend = 20,
                   alpha = 0.4)
  }
  if(is.null(custom_title)){
    custom_title <-sprintf("Assigned to %s;\ntransition(s): %s, respectively;\nstarting state: %s",
                           paste0(edges, collapse = ", "),
                           paste0(strsplit(only_states, "_")[[1]], collapse = ", "),
                           as.character(ss))
  }
  p3 <- suppressMessages(ggtree::gheatmap(p1,
                     data=pos_rates_flip,
                     offset = 30,
                     # offset = ifelse(donor_id == "KX004", ifelse(include_tip_labels, 180, 30), 30)
                     colnames_angle=45,
                     hjust = 1,
                     color=NA,
                     font.size = 2) +
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
  return(p3)
}
