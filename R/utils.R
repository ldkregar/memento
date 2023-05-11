return_edge_descendants <- function(tree){
  lapply(tree$edge[,2], function(nn) paste0(tree$tip.label[phangorn::Descendants(tree, node = nn, type = "tips")[[1]]], collapse = ":"))
}
