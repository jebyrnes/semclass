library(DiagrammeR)


plot.psem <- function(mod, return=FALSE,
                      node_attrs = data.frame(shape = "rectangle", color = "black",
                                              fillcolor = "white"),
                      edge_attrs = data.frame(style = "solid", color="black"),
                      ns_dashed = T, alpha=0.05,
                      show = "std", digits = 3, 
                      add_edge_label_spaces = TRUE, ...
                      ){
  #get the coefficients table
  ctab <- coefs(mod)
  ctab$Response <- as.character(ctab$Response)
  ctab$Predictor <- as.character(ctab$Predictor)
  
  #make a nodes DF
  unique_nodes <- unique(c(ctab$Response, ctab$Predictor))
  nodes <- create_node_df(n = length(unique_nodes),
                          nodes = unique_nodes,
                          type = "lower",
                          label = unique_nodes)
  nodes <- cbind(nodes, node_attrs)
  nodes[] <- lapply(nodes, as.character)
  nodes$id <- as.numeric(nodes$id)
  
  
  #make an edges DF
  edges <- create_edge_df(
                          from = match(ctab$Predictor, unique_nodes),
                          to = match(ctab$Response, unique_nodes))
  
  edges <- data.frame(edges, edge_attrs)
  edges[] <- lapply(edges, as.character)
  edges$id <- as.numeric(edges$id)
  edges$from <- as.numeric(edges$from)
  edges$to <- as.numeric(edges$to)
  if(ns_dashed) edges$style[which(ctab$P.Value>alpha)] <- "dashed"
  if(show == "std") edges$label = round(ctab$`Std.Estimate`, digits)
  if(show == "unstd") edges$label = round(ctab$Estimate, digits)
  if(add_edge_label_spaces) edges$label = paste0(" ", edges$label, " ")
  
  
  #turn into a graph
  sem_graph <- create_graph(nodes, edges, directed=TRUE)
  
  if(return) return(sem_graph)
  
  render_graph(sem_graph, ...)
  
  
}
