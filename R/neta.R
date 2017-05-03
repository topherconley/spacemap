


adj2ig <- function(YY, XY = NULL, weighted = F, dropnull = TRUE) { 
  
  if (!is.null(XY)) { 
    adjnet <- rbind(cbind(matrix(0, nrow(XY), nrow(XY)), XY), 
                      cbind(t(XY), YY))
  } else { 
    adjnet <- YY 
  } 
  
  if (weighted) { 
    ig <- igraph::graph_from_adjacency_matrix(adjnet, mode = "max", weighted = weighted, diag = FALSE, add.colnames = NA) 
  } else { 
    ig <- igraph::graph_from_adjacency_matrix(adjnet, mode = "max", diag = FALSE, add.colnames = NA)  
  }
  
  #remove vertices with no edges
  if(dropnull) { 
    ig <- igraph::delete_vertices(graph  = ig, v = V(ig)[igraph::degree(ig) == 0])
  }
  ig
}







netaHub <- function(net, xinfo, yinfo, bdeg) { 
  
}