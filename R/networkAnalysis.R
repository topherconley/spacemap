#' Convert adjacency matrix to igraph with annotations.
#' 
#' @description Convert adjacency matrix output from conditional graphical model output (xy, yy) into 
#' into an igraph object and optionally layer attributes onto the vertices. 
#' @param yy Adjacency matrix encoding y--y edges of conditional graphical model.
#' Format can be either a sparse matrix (i.e. from Matrix package) or a regular matrix. 
#' @param xy Adjacency matrix encoding x--y edges of conditional graphical model. 
#' Format can be either a sparse matrix (i.e. from Matrix package) or a regular matrix.
#' If NULL, then a network encoding just y--y edges is made.
#' @param yinfo Data.frame encoding the attributes of the y vertices. If non-null, 
#' the column 'id' must exist for labelling vertices with a unique identifier. 
#' The order of the rows must match the order of the \code{yy} columns/rows. 
#' Other attributes in the form of data.frame columns are optional and 
#' can have missing data for some of the vertices. 
#' See 'Details' for further explanation.
#' @param xinfo Data.frame encoding the attributes of the x vertices. If non-null, 
#' the column 'id' must exist for labelling vertices with a unique identifier. 
#' The order of the rows must match the order of the \code{xy} columns/rows. 
#' Other attributes in the form of data.frame columns are optional and 
#' can have missing data for some of the vertices. 
#' See 'Details' for further explanation.
#' @param weighted Same argument passed to \code{igraph::graph_from_adjacency_matrix}. 
#' Defaults to NULL. 
#' 
#' @details The annotation parameters \code{xinfo} and \code{yinfo}
#'  are of type `data.frame`, where each row is expected to be ordered 
#'  according to the rows of \code{xy} and \code{yy}, respectively. 
#'  The first column `id` is always required and must be unique for each node.
#'  The second column `alias` is an optional column that reserves more 
#'  human-readable labels to be applied to nodes, and can be duplicated. 
#'  For example, two protein isoforms may have different ID's,
#'  but the same gene symbol alias. If gene coordinates are supplied for 
#'  cis/trans identification, they require at least three variables: 
#'  'chr' a character for the chromosome location, 
#'  'start' an integer for the beginning of the genomic feature, 
#'  and 'end' for the end location of the feature. 
#'  The 'strand' argument is optional. 
#'  Additional annotations such as `description` or
#'  whatever other annotations are deemed important 
#'  can be added as node attributes. Note that if any vertex
#'  is missing some or all annotations, an `NA` should be used in its place. 
#' 
#' @return An igraph network object.
#' @seealso \code{\link{rankHub}}, \code{\link{cisTrans}}, 
#' \code{\link{reportHubs}}, \code{\link{xHubEnrich}}, \code{\link{modEnrich}}
#' @examples 
#' 
#' #Load BCPLS CNA-protein network output
#' library(spacemap)
#' data("bcpls")
#' #Attach "bcpls" data objects to the R Global Environment
#' attach(bcpls)
#' 
#' suppressPackageStartupMessages(library(igraph))
#' #Label y and x nodes with attributes
#' ig <- spacemap::adj2igraph(yy = net$yy, xy = net$xy, 
#'                            yinfo = yinfo, xinfo = xinfo)
#' #Label y with attributes
#' igy <- spacemap::adj2igraph(yy = net$yy, 
#'                            yinfo = yinfo)
#' @importFrom igraph graph_from_adjacency_matrix V E set_vertex_attr %--%
#' @importFrom Matrix Matrix t                    
#' @export
adj2igraph <- function(yy, xy = NULL, yinfo = NULL, xinfo = NULL, weighted = NULL) { 
  
  requireNamespace("igraph")
  
  if(!inherits(yy, "Matrix")) { 
    if(!is.matrix(yy)) { 
      stop("yy must inherit class matrix or Matrix.")
    }
  }
  
  givenX <- !is.null(xy)
  givenY <- !is.null(yy)
  
  if(givenX & !inherits(xy, "Matrix")) { 
    if(givenX & !is.matrix(xy)) { 
      stop("xy must inherit class matrix or Matrix.")
    }
  }
  if (givenX & !givenY) { 
    stop("yy is NULL, but xy is non-NULL. Input not supported. ")  
  }
  
  if (givenX) { 
    Xindex <- seq_len(nrow(xy))
    Yindex <- (length(Xindex) + 1):(nrow(xy) + nrow(yy))
    requireNamespace("Matrix")
    xy <- Matrix(xy)
    yy <- Matrix(yy)
    adjnet <- rbind(cbind(Matrix(0, nrow(xy), nrow(xy)), xy), 
                    cbind(Matrix::t(xy), yy))
    
  } else { 
    Yindex <- seq_len(nrow(yy))
    adjnet <- Matrix(yy)
  } 
   
  ig <- igraph::graph_from_adjacency_matrix(adjnet, mode = "max", weighted = weighted, diag = FALSE, add.colnames = NA) 


  if (!is.null(xinfo)) { 
    if(!givenX) { 
        stop("xy cannot be null and xinfo non-null.")
    }
    if (nrow(xinfo) != length(Xindex)) { 
      stop("xinfo and xy are not of the same dimension.")  
    }
    if (is.null(xinfo$id)) { 
      stop("xinfo$id is null.")
    } 
    ig <- layerInfo(ig = ig, info = xinfo, vs =  V(ig)[Xindex])
    
  } 
  
  if (!is.null(yinfo)) { 
    if (nrow(yinfo) != length(Yindex)) { 
      stop("yinfo and yy are not of the same dimension.")  
    }
    if (is.null(yinfo$id)) { 
      stop("yinfo$id is null.")
    }
    ig <- layerInfo(ig = ig, info = yinfo, vs =  V(ig)[Yindex])
  }
  
  if (givenX) { 
    ig <- set_edge_attr(graph = ig, name = "levels",
                        index = E(ig)[Xindex %--% Yindex], value = "x-y")  
    ig <- set_vertex_attr(graph = ig, name = "levels", index = V(ig)[Xindex], value = "x")
  }
  ig <- set_edge_attr(graph = ig, name = "levels",
                      index = E(ig)[Yindex %--% Yindex], value = "y-y")  
  ig <- set_vertex_attr(graph = ig, name = "levels", index = V(ig)[Yindex], value = "y")
  
  #remove vertices with no edges
  #if(dropnull) { 
  #  ig <- igraph::delete.vertices(graph  = ig, v = V(ig)[igraph::degree(ig) == 0])
  #}
  ig
}

layerInfo <- function(ig, info, vs) { 
  requireNamespace("igraph")
  for(var in names(info)) { 
    #require factors to be characters
    if(is.factor(info[,var])) { 
      info[,var] <- as.character(info[,var])
    }
    if(var == "id") { 
      ig <- igraph::set_vertex_attr(graph = ig, name = "name", 
                                    index = vs,
                                    value = info[,var])  
    } else { 
      ig <- igraph::set_vertex_attr(graph = ig, name = var, 
                                    index = vs,
                                    value = info[,var])  
    }
  }
  ig
}

#' Prioritize networks hubs by degree.
#' 
#' @description Rank hubs by degree either with the final network or 
#' according to how consistently the degree is high across bootstrap replicates. 
#' @param ig An igraph network object output from \code{\link{adj2igraph}}
#' @param bdeg If NULL, the ranking of hubs nodes is done by their 
#' highest degree in the \code{ig} network within their respective node 
#' \code{level = "y"} or \code{level = "x"}. 
#' 
#' If non-null, List output from list element of \code{\link{bootVote}}. 
#' Element \code{bdeg$yy[b,]} is an integer vector representing 
#' the degree distribution for the \eqn{b}th bootstrap replicate 
#' across the y nodes.  
#' Similarly, element \code{bdeg$xy[b,]} is an integer vector 
#' representing the degree distribution for the \eqn{b}th bootstrap 
#' replicate across the x nodes. 
#' If `bdeg` is available, the hubs will be prioritized according to 
#' their mean rank of degree across the \eqn{B} bootstrap replicates. 
#' Highly ranked hubs consistently have a larger degree than other 
#' nodes across the bootstrap replicates.
#' @param level Character value either 'y' or 'x' (defaults to 'x') 
#' specifying to rank hubs that are x or y nodes. Ranking levels 
#' x and y together is not supported.
#' @return The igraph object \code{ig} updated with a 
#' 'rank_hub' attribute containing each nodes' respective rank. 
#' If \code{bdeg} is non-null the mean degree rank across
#' bootstrap replicates is labeled with attribute 'mean_rank_hub',
#' and its corresponding standard deviation is labeled
#' as attribute 'sd_rank_hub'.
#' @importFrom igraph V E set_vertex_attr
#' @importFrom stats sd 
#' @export
#' @seealso \code{\link{adj2igraph}}, \code{\link{cisTrans}}, 
#' \code{\link{reportHubs}}, \code{\link{xHubEnrich}}, 
rankHub <- function(ig, bdeg = NULL, level = c("x", "y")) { 
  
  requireNamespace("igraph")
  
  level <- match.arg(level)
  vs <- V(ig)[levels %in% level]
  
  if (is.null(bdeg)) { 
    deg <- igraph::degree(graph = ig, v = vs)
    ig <- set_vertex_attr(graph = ig, name = "rank_hub", index = vs, value = rank(-deg))
  } else {
    brank <- apply(-bdeg, 1, rank)
    #average the ranks across the bootstrap replicates
    avg_bdeg_ranks <- rowMeans(brank)
    sd_bdeg_ranks <- apply(brank, 1, sd)
    ig <- set_vertex_attr(graph = ig, name = "rank_hub", 
                          index = vs, value = rank(avg_bdeg_ranks))
    ig <- set_vertex_attr(graph = ig, name = "mean_rank_hub", 
                          index = vs, value = avg_bdeg_ranks)
    ig <- set_vertex_attr(graph = ig, name = "sd_rank_hub", 
                          index = vs, value = sd_bdeg_ranks)
    #bdeg_pres <-  paste(round(avg_bdeg_ranks,2), 
    #" (", round(sd_bdeg_ranks, 2), ")")
  }
  ig
}


splitEdgeVector <- function(evec, ig) { 
  
  requireNamespace("igraph")
  lev <- strsplit(evec, "|", fixed  = TRUE)
  tab <- data.frame(left = sapply(lev, function(x) x[1]),
                    right = sapply(lev, function(x) x[2]), stringsAsFactors = F)
  
  #no nodes of levels x should be in the right column, 
  #flip to the left column so that it is always (x, y) or (y, y)
  rev_idx <- which(tab$right %in% as_ids(V(ig)[levels %in% "x"]))
  if(length(rev_idx) > 0) { 
    tab_rev <- tab[rev_idx,]
    tab[rev_idx,] <- tab[rev_idx,c("right", "left")]
  }
  tab
}

#' Identify x-y cis and trans edges.
#'
#' @description Label igraph object with cis/trans 
#' status for x-y edges. 
#' 
#' @param ig An igraph network object output from \code{\link{rankHub}}
#' @param level Character vector currently only supporting 'x-y' edge 
#' labeling of cis/trans status. 
#' @param cw Numeric value denoting the half length of the genomic interval 
#' used for calling a response node (i.e. of levels y) as cis regulated by a predictor 
#' node (i.e. of levels x). Defaults to 2Mb. 
#' @param ignoreStrand Logical defaults to TRUE specifying that strand specificity 
#' is not required to call cis/trans edges.
#' 
#' @details This function requires the GenomicRanges package, 
#' specifically 'findOverlaps' and 'countOverlaps'. The \code{ig} 
#' network parameter requires vertex attributes for genome coordinates 
#' specified as \code{chr,start,end} 
#' (see details section of \code{\link{adj2igraph}}).
#' 
#' @return The igraph object \code{ig} with x-y edges 
#' updated with attribute 'cis_trans' indicating cis/trans status. 
#' 
#' Also x node vertices are updated with attributes 
#' 'ntrans' reporting the number of trans regulations by a node
#'  (similarly for attribute 'ncis') and 'regulates_in_cis'
#'  lists specific y nodes regulated in cis by x nodes.
#' @importFrom igraph as_data_frame set_vertex_attr V E edge_attr as_ids set_edge_attr
#' @importFrom stats na.omit
#' @importFrom GenomicRanges GRanges findOverlaps countOverlaps 
#' @importClassesFrom S4Vectors Rle
#' @export
#' @seealso \code{\link{rankHub}}, \code{\link{adj2igraph}}, 
#' \code{\link{reportHubs}}, \code{\link{xHubEnrich}}
cisTrans <- function(ig, level = c("x-y"), cw = 2e6, ignoreStrand = TRUE) { 
  
  requireNamespace("igraph")
  #get edges and vertex attributes for specified levels
  de <- as_data_frame(x = ig, what = "edges")
  sublevel <- de[,"levels"] %in% level
  if (all(!sublevel)) { 
    stop("level parameter did not have any edges")  
  }
  de <- de[sublevel,]
  #de <- splitEdgeVector(evec, ig)
  
  dv <- as_data_frame(x = ig, what = "vertices")
  dvn <- names(dv)
  
  #define features for which missing values is not tolerated 
  requiredAttr <-   c("name","chr", "start", "end")
  if(!ignoreStrand) {
    requiredAttr <- c(requiredAttr, "strand")
  } 
  for(ra in requiredAttr)
    if (!(ra %in% dvn)) { 
      stop("Vertex attribute", ra, "is required, but is missing for all vertices.")  
    }

  #remove features that have missing values in genome coordinates that 
  #does not permite cis/trans identification
  rmNA <- na.omit(dv[,requiredAttr])
  dgc <- dv[-attr(rmNA, "na.action"),requiredAttr]
  
  if(!is.null(dgc$strand)) { 
    dgc$strand <- as.character(dgc$strand)
    dgc$strand <- ifelse(is.na(dgc$strand), "*", dgc$strand)
  } else { 
    dgc$strand <- "*"  
  }
  
  if(nrow(dgc) == 0) { 
    stop("Too many missing values in gene coordinates to compute cis/trans activity.")
  }

  requireNamespace("GenomicRanges")
  gr <- GenomicRanges::GRanges(seqnames = Rle(dgc$chr),
                ranges = IRanges(start = pmax(dgc$start - cw, 1), 
                                 end = dgc$end + cw, 
                                 names = dgc$name),
                strand = Rle(strand(dgc$strand)))

  ## cis/trans edge identification
  
  #identify those vertices with edges and gene coordinates
  fromgc <- unique(de$from[de$from %in% names(gr)])
  togc <- unique(de$to[de$to %in% names(gr)])
  
  #find overlapping
  ov <- GenomicRanges::findOverlaps(query = gr[fromgc,], 
                                    subject = gr[togc,], 
                                    ignore.strand = ignoreStrand)
  candidateCisEdges <- paste0(fromgc[from(ov)], "|", togc[to(ov)])
  igEdges <- as_ids(E(ig))
  cis <- igEdges[igEdges %in%  candidateCisEdges]
  trans <- setdiff(igEdges, cis)
  ig <- set_edge_attr(ig, name = "cis_trans", 
                     index = E(ig)[ igEdges %in% cis], 
                     value = "cis")
  ig <- set_edge_attr(ig, name = "cis_trans", 
                      index = E(ig)[ igEdges %in% trans], 
                      value = "trans")
  
  #cis/trans distribution
  vseqx <- V(ig)[match(fromgc, V(ig)$name)]	
  ig <- setCisTransDistr(ig = ig, vseq = vseqx)	
  
  #potential cis
  allYNodes <- as_ids(V(ig)[levels %in% "y"])
  allYNodesGC <- allYNodes[allYNodes %in% names(gr)]
  potcis <- GenomicRanges::countOverlaps(query = gr[fromgc,], 
                                         subject = gr[allYNodesGC,], 
                                         ignore.strand = ignoreStrand)
  ig <- set_vertex_attr(graph = ig, name = "potential_cis", 
                      index = V(ig)[match(names(potcis), V(ig)$name)], 
                      value = potcis)
  
  ig
}

setCisTransDistr <- function(ig, vseq) {
  requireNamespace("igraph")
  vids <- as_ids(vseq)
  dct <- t(sapply(vseq, function(v) {
    #
    #vid <- igraph::as_ids(v)
    eseq <- as_ids(E(ig)[ inc(v) ])
    vct <- edge_attr(ig, name = "cis_trans", index = eseq)  
    de <- splitEdgeVector(eseq, ig)
    cidx <- vct %in% "cis"
    ncis <- sum(cidx)
    if (ncis > 0) { 
      dec <- de[cidx,]
      chits <- setdiff(c(dec[,"right"], dec[,"left"]),vids)
      cisgenes <- paste(V(ig)[match(chits, V(ig)$name)]$alias, 
                        collapse = ", ")
    } else { 
      cisgenes <- "--"
    }
    c(ncis = ncis, ntrans = length(vct) - ncis, cisgenes = cisgenes)
  }))
  dct <- as.data.frame(dct, stringsAsFactors = F)
  dct$ncis <- as.integer(dct$ncis)
  dct$ntrans <- as.integer(dct$ntrans)
  #idx <- match(rownames(dct), V(ig)$name)
  #V(ig)[idx]
  ig <- set_vertex_attr(graph = ig, name = "ncis", index = vseq, value = dct$ncis)
  ig <- set_vertex_attr(graph = ig, name = "ntrans", index = vseq, value = dct$ntrans)
  ig <- set_vertex_attr(graph = ig, name = "regulates_in_cis", index = vseq, value = dct$cisgenes)
  ig
}

#' Report top hub nodes with attributes.
#' 
#' @description Organize x and y hub analysis into a publish-ready format.
#' @param ig An igraph network object returned from either 
#' \code{\link{rankHub}} or \code{\link{cisTrans}}.
#' @param top An integer specifying the number of top-degree hubs to report.
#' @param level Character specifying nodes of level x or y, but not both. 
#' Defaults to x.
#' @param unit If gene coordinates for x-hubs exist, then this will 
#' label the chosen genome location scale in the table. Defaults to mega-base (Mb).
#' @param unit_base Numeric divisor of the genome coordinates and should 
#' correspond to \code{unit} label. Defaults to 1e6 for mega-base.
#' @param sig Integer for controlling precision of genome coordinates.
#' 
#' @return A table of top ranked node hubs with attributes deemed 
#' relevant to report. See vignette on network analysis for examples. 
#' @importFrom igraph as_data_frame degree
#' @export
#' @seealso \code{\link{rankHub}}, \code{\link{cisTrans}}, 
#' \code{\link{adj2igraph}}, \code{\link{xHubEnrich}}, \code{\link{modEnrich}}
reportHubs <- function(ig, top = 10, level = c("x", "y"), 
                       unit = "Mb", unit_base = 1e6, sig = 2) {
  
  requireNamespace("igraph")

  level <- match.arg(level)
  hubidx <- which(degree(ig) != 0 & V(ig)$levels == level)
  rawtab <- as_data_frame(x = ig, what = "vertices")[hubidx,]
  #remove those columns that are all NA's
  rawtab <- rawtab[, colSums(is.na(rawtab)) != nrow(rawtab)]
  
  #Name Column
  if(is.null(rawtab$alias)) { 
    label <- rawtab[,"name"]  
  } else  { 
    label <- rawtab[,"alias"]  
  }
  tab <- data.frame(hub = label, row.names = rawtab[,"name"])
  
  if (level == "x") { 
    if(!is.null(rawtab$start) & !is.null(rawtab$end) & !is.null(rawtab$chr)) { 
      ub_start <- signif(as.integer(rawtab$start) / unit_base,sig)
      ub_end <- signif(as.integer(rawtab$end)/ unit_base, sig)
      tab$hub <- paste(label, 
                       paste0("(", paste(ub_start, ub_end, sep = "-"), 
                              paste0(" ", unit, ")")))
    }
    
    if (!is.null(rawtab$ncis) & !is.null(rawtab$ntrans))
      tab$`# cis/ # trans` =  paste0(rawtab$ncis, " / ", rawtab$ntrans) 
    
    
    if(!is.null(rawtab$potential_cis)) 
      tab$`Potential # cis` <- rawtab$potential_cis
    
    if(!is.null(rawtab$regulates_in_cis)) 
      tab$`cis genes` <- rawtab$regulates_in_cis
    
  } else if (level == "y") { 
    
    tab$degree <- igraph::degree(graph = ig, v = V(ig)[name %in% rawtab$name])
    
  }
  
  if(!is.null(rawtab$description)) { 
    tab$description <- rawtab$description
  }
  
  tab <- tab[order(rawtab$rank_hub),]
  tab[seq_len(min(top,nrow(tab))),]
}

#' Functional enrichment of x-hub neighborhood.
#' 
#' @description A means to assess functional 
#' impact of x-hub perturbations
#' 
#' @param ig An igraph network object output from either 
#' \code{\link{adj2igraph}, \link{rankHub}, \link{cisTrans}}
#' @param go2eg  Named list where the names denote a
#'  biological process (e.g. Gene Ontology ID) and
#' the elements of the list is a vector of members 
#' belonging to the biological process. 
#' The list ought to be non-redundant in names. 
#' 
#' @details A x-hub is defined to be any x node with at least one edge to a y node. 
#' A x neighborhood comprises all y nodes that are directly connected to a x-hub by an edge. 
#' x neighborhoods are of interest because they represent direct perturbations to y nodes.
#' To quantitatively assess how much those perturbations are functionally meaningful,
#' we compute a score called the GO-neighbor percentage. Two y nodes are called GO-neighbors 
#' if they share a common Gene Ontology (GO) term in the same x neighborhood. We postulate 
#' that a high percentage of GO-neighbors within an x neighborhood associates 
#' the x-hub with more functional meaning. These scores, 
#' as presented in figure 5 of the spaceMap publication, 
#' can be generated with a GO mapping
#' 
#' @return Data.frame with each x-hub, out-degree, and the 
#' neighborhood percentage. 
#' @importFrom igraph as_ids V E degree as_adj
#' @importFrom utils combn
#' @export
#' @seealso \code{\link{rankHub}}, \code{\link{cisTrans}}, 
#' \code{\link{reportHubs}}, \code{\link{adj2igraph}}
xHubEnrich <- function(ig, go2eg) { 
  
  requireNamespace("igraph")
  
  xytab <- splitEdgeVector(igraph::as_ids(E(ig)[levels %in% "x-y"]), ig)
  hubHood <- split(x = xytab$right, f = as.factor(xytab$left))
  hubHoodSize2 <- hubHood[sapply(hubHood, length) > 1]
  prop <- sapply(hubHoodSize2, function(gt1) { 
    
    #make a go adjacency matrix
    ngt <- length(gt1)
    mgo <- matrix(0, nrow = ngt, ncol = ngt)
    pw <- combn(1:ngt,2)
    for(j in 1:ncol(pw)) { 
      ep <- gt1[pw[,j]]
      mgo[pw[1,j],pw[2,j]] <- mgo[pw[2,j],pw[1,j]] <- any(sapply(go2eg, function(egs) all(ep %in% egs)))
    }
    m3 <- sum(colSums(mgo) > 0) / ngt
    m3
  })
  vidx <- match(names(hubHoodSize2), V(ig)$name)
  vseq <- V(ig)[vidx]
  
  if(!is.null(V(ig)$alias)) { 
    hub <- V(ig)$alias[vidx]
  } else { 
    hub <- as_ids(vseq)
  }
  data.frame(hub = hub,
             degree = igraph::degree(graph = ig, v = vseq), 
             neighbor_percentage = prop*100)
}

# rndCGGraph <- function(xy, yy, info) { 
#   
#   requireNamespace("igraph")
#   
#   nz_xy <- which(rowSums(xy) > 0)
#   nz_yy <- which(rowSums(yy) > 0)
#   rxy <- xy
#   Q <- ncol(xy)
#   rxy[nz_xy,] <- t(sapply(nz_xy, function(i) sample(x = xy[i,], size = Q, replace = F)))
#   igy <- adj2igraph(yy = yy, xy = NULL, info = info[info$levels == "y",], dropnull = FALSE)
#   rigy <- rewire(graph = igy, with = keeping_degseq(niter = vcount(igy) * 10))
#   ryy <- as_adj(graph = rigy)
#   rig <- adj2igraph(yy = ryy, xy = rxy, info = info, weighted = F)
#   rig <- setEdgeTypeAttr(rig, "x", "y")
#   rig <- setEdgeTypeAttr(rig, "y", "y")
#   #rig <- xyCisTrans(ig = rig, et = c("x-y", "y-y"))
#   rig
# }
# 
# rndHubGoProportion <- function(xy, yy, info, go2eg) { 
#   rig <- rndCGGraph(xy = xy, yy = yy, info = info)
#   rhgp <- hubGOproportion(rig, go2eg)
#   #weighted.mean(x = rhgp$`Proportion of GO Pairs`, w = rhgp$`Out Degree`)
#   colMeans(rhgp[,c("m1","m2","m3")])
# }

#' Module enrichment analysis.
#' 
#' @description Identify biological processes 
#' which are significantly-enriched in network modules. 
#' 
#' @param ig An igraph network object output from either 
#' \code{\link{adj2igraph}, \link{rankHub}, \link{cisTrans}}.
#' @param mods modules with three accepted formats. 
#' \itemize{ 
#'   \item An igraph object of class 'communities', 
#'   which is typically the output of the cluster_* functions of igraph.
#'   \item A list where each element is a character vector that only contains 
#'   the node identifiers from either \code{yinfo$id} or  \code{xinfo$id} input to 
#'   function \code{\link{adj2igraph}}.
#'   \item An integer vector with names as node identifiers and values as an integer. 
#'   For example \code{c(id1 = 1, id2 = 2, id3 = 2, id4 = 1)}.
#' }
#' @param levels Any given module can contain x nodes or y nodes.
#'   If both predictors and responses have a functional mapping 
#'   in the \code{go2eg} argument, then specify \code{levels = c("x","y")}. 
#'   Otherwise, specify only those nodes that have a functional mapping. 
#'   See details for more discussion.
#' @param go2eg  Named list where the names denote a
#'  biological process (e.g. Gene Ontology ID) and
#' the elements of the list is a vector of members 
#' belonging to the biological process. 
#' The list ought to be non-redundant in names. For example, 
#' \code{list(bio_proc_1 = c("gene1", "gene2", "gene3"),
#'            bio_proc_2 = c("gene4", "gene5", "gene6")
#'            )}
#' 
#' @param glb Integer defining the smallest possible size of a module
#' in order for the module to be tested for enrichment.
#' @param minGO  Integer defining the smallest possible number of 
#' nodes represented in a biological process to be called a 
#' significant enrichment of that biological process. 
#' @param thresh Numeric between 0 and 1 indicating the threshold at which 
#' adjusted P-values should be considered significant.
#' @param adjust Character of type \code{stats::p.adjust.methods} for 
#' specifying the type of multiple comparison adjustment desired.
#' @param prefix Character to prefix module identifiers.
#' @param process_alias Vector mapping biological process identifiers in 
#' \code{go2eg} with biologically meaningful descriptions. The vector 
#' \code{process_alias} must have names
#' as the same names in \code{go2eg} and the elements are the  
#' biologically meaningful descriptions.
#' 
#' @details The hyper-geometric test is used to test for 
#' over-representation of a biological process. In the 
#' \code{phyper} R function, parameter \code{q} is the 
#' overlap between the biological process group and the module, 
#' where the module is reduced to only its y node members if \code{level = "y"}.
#' Parameter \code{m} is the size of the biological process.
#' Parameter \code{n} is the number of nodes in the network not 
#' in the biological process. This excluding node levels that do not have a functional mapping.
#' In other words, if no x nodes do not appear in the mapping of \code{go2eg}, 
#' corresponding to \code{level = "y"}, then x nodes are not counted,
#'  but "y" nodes without a mapping are counted, because most y nodes do have a mapping.  
#'  
#' @return A three-item list: 
#' 
#' \itemize{
#' \item The element "ig" is the input igraph network object 
#' updated with a "process_id" attribute 
#' for nodes affiliated with a significant GO-term.
#' The "process_id" and "module" attributes together can be especially
#' useful for visualizing which nodes of a module 
#' are enriched for a specific biological function. 
#' \item The element "etab" is the polished module enrichment table
#'  conveniently organized to report significant GO terms in modules,
#'   the representation of the GO term in the module relative to the 
#'   size of the GO term, and what x-hubs may belong to the module. 
#' \item The element "eraw" contains details for each (module, GO-term) pair
#'  that was subjected to the hyper-geometric test. 
#'  This output gives the user more control (if desired) over enrichment 
#'  by reporting all tests, the relative over-representation 
#'  of a GO-term in that module, the raw P-value, and the adjusted P-value.
#' }
#' @importFrom igraph set_vertex_attr V membership as_ids vertex_attr
#' @importFrom foreach %:% %do%
#' @importFrom stats p.adjust p.adjust.methods phyper
#' @export
#' @seealso \code{\link{adj2igraph}}, \code{\link{rankHub}}, \code{\link{cisTrans}}, 
#' \code{\link{reportHubs}}, \code{\link{xHubEnrich}}
modEnrich <- function(ig, mods, levels = c("x", "y"), go2eg, glb = 15, minGO = 5, 
                      thresh = 0.05, adjust = "BH",  
                      prefix = "M", process_alias = NULL) { 

  requireNamespace("igraph")
  
  levels <- match.arg(levels, choices = c("x", "y"), several.ok = TRUE)
  adjust <- match.arg(adjust, choices = stats::p.adjust.methods)
  
  #list is easier to work with than a vector of members
  ll <- mods2list(ig = ig, mods = mods, glb = glb, prefix = prefix)
  lmods <- ll$all
  ylmods <- ll$yy

  #annotate nodes with module membership
  for(i in seq_along(lmods)) {
    ig <- set_vertex_attr(graph = ig, name= "module", 
                          index = V(ig)[ name %in% lmods[[i]] ], 
                          value = names(lmods)[i])
  }
  
  if(sum(levels %in% c("x", "y")) == 2) { 
    inmods <- lmods
  } else { 
    inmods <- ylmods
  }

  
  #hypothesis testing
  emods <- hypergeoTest(ig = ig, levels = levels, lmods = inmods, go2eg = go2eg,
                            thresh = thresh, adjust = adjust,
                            process_alias = process_alias)
  #make pretty result
  etab <- renderModuleTable(ig = ig, lmods = lmods, 
                            emods = emods, minGO = minGO)
  #store result in igraph
  ig <- setGOTags(lmods = lmods, emods = emods, ig = ig, go2eg = go2eg)
  
  #organize raw results 
  if (!is.null(process_alias)) { 
    proc_alias_name <- "process_alias"
  } else { 
    proc_alias_name <- NULL
  }
  eraw <- emods[,c("modname", "process_id", proc_alias_name, "representation", "pvalue", "fdr")]
  names(eraw)  <- c("Module", "process_id", proc_alias_name, "GO Obs./Total", "P_value", "Adjusted_P_value")
  
  list(ig = ig, etab = etab, eraw = eraw)
}

mods2list <- function(ig, mods, glb, prefix) { 
  
  requireNamespace("igraph")
  
  if (class(mods) == "communities") { 
    modmems <- membership(mods)
    modtab <- sizes(mods)
    #Group the modules and report the modules that are of at least size glb.
    modsBigEnough <- modtab[modtab >= glb]
    #form list of modules
    lmods <- lapply(names(modsBigEnough), function(mbe) { 
      mbe <- as.integer(mbe)
      names(modmems)[modmems == mbe]
    })
    #name list
    names(lmods) <- sapply(seq_along(lmods), function(i) paste0(prefix, i))
  } else if (is.list(mods)) { 
    lmods <- mods
    #missing all names OR missing some names 
    if (is.null(names(lmods)) | any(nchar(names(lmods)) == 0)) { 
      warning("Renaming module names because list of modules not correctly specified.")
      names(lmods) <- sapply(seq_along(lmods), function(i) paste0(prefix, i))
    }
  } else if (is.vector(mods)) { 
    if(!is.null(names(mods)) | any(nchar(names(mods)) == 0)) { 
      stop("A name for each element of vector 'mods' not specified.")  
    }
    mods <- as.integer(mods)
    #all integer valued or all numeric valued
    if(!all(sapply(mods, is.integer))) { 
      stop("Input 'mods' ought to be an integer vector.")
    }
    modtab <- table(mods)
    modsBigEnough <- modtab[modtab >= glb]
    #form list of modules
    lmods <- lapply(names(modsBigEnough), function(mbe) { 
      mbe <- as.integer(mbe)
      names(mods)[mods == mbe]
    })
  } else { 
    stop("Module input 'mods' is not currently supported.")  
  }
  
  if(any(sapply(lmods, function(m) any(is.na(m))))) { 
    stop("No missing values allowed  in 'mods' input.")
  }
  
  #assure all members of module belong to node ids
  valid_mods <- sapply(lmods, function(modmems) {
    modmems <- as.character(modmems)
    all(modmems %in% V(ig)$name)
  })
  if (!all(valid_mods)) { 
    stop("There exists members of the modules that are not contained in graph 'ig'.")
  }
  
  #no duplicates allowed in modules
  #but currently allow for overlapping fuzzy modules 
  if (!all(sapply(lmods, anyDuplicated) == 0)) { 
    stop("No duplicates allowed within a modules.")  
  }  
  
  ylmods <- lapply(lmods, function(mod) { 
    as_ids(V(ig)[(levels %in% "y") & name %in% mod])
  })
  
  list(all = lmods, yy = ylmods)
}

hypergeoTest <- function(ig, levels, lmods, go2eg, thresh, adjust, process_alias) { 
  
  #number of y nodes 
  requireNamespace("igraph")
  ng <- sum(V(ig)$levels %in% levels)
  requireNamespace("foreach")
  #for R CMD check NOTE passing
  i <- NULL; j <- NULL;
  modgo <- foreach::foreach(i = seq_along(lmods)) %:%
    foreach::foreach(j = seq_along(go2eg), .combine = 'rbind') %do% { 
      mod <- lmods[[i]]
      gt <- go2eg[[j]]
      #overlap
      ov <- length(intersect(mod,gt))
      #minus -1 because of the upper tail
      pval <- phyper(q = ov - 1L, 
                     m = length(gt), 
                     n = ng - length(gt), 
                     k = length(mod),lower.tail = FALSE)
      data.frame(process_id = names(go2eg)[j], pvalue = pval, ov = ov,
                 representation = paste(ov, length(gt), sep = "/"),
                 modname = names(lmods)[i], stringsAsFactors = F)
    }
  amodgo <- do.call(rbind, modgo)
  #multiple-testing correction with Benjamini Hochberg
  amodgo$fdr <- p.adjust(p = amodgo$pvalue, method = adjust)
  #identify significant terms
  amodgo$sig <- amodgo$fdr < thresh
  #order by module name and then by significance
  amodgo <- amodgo[order(amodgo$modname, amodgo$fdr),]
  if (!is.null(process_alias)) { 
    amodgo$process_alias <- process_alias[match(amodgo$process_id, names(process_alias))]
    amodgo <- amodgo[,c("modname", "process_id", "process_alias", "representation", "fdr", "pvalue", "sig", "ov")]
  } else { 
    amodgo <- amodgo[,c("modname", "process_id", "representation", "fdr", "pvalue", "sig", "ov")]
  }
  amodgo
}

setGOTags <- function(lmods, emods, ig, go2eg) { 
  
  requireNamespace("igraph")
  
  for (mod_name in names(lmods)) {
    mod_refseq <- lmods[[mod_name]]
    mod_attr_id <- V(ig)$name[match(mod_refseq, V(ig)$name)]
    #mod_attr_id <- info[match(mod_refseq, info$id),"id"]
    
    eidx <- (emods$modname %in% mod_name) & emods$sig
    if (sum(eidx) == 0) next; 
    enriched_goids <- emods[,"process_id"]
    #record enriched process_id's into the igraph object
    
    goid_enriched <-  sapply(mod_attr_id, function(prot) 
      paste0(enriched_goids[sapply(go2eg[enriched_goids], 
                                   function(g) prot %in% g)], collapse = "; "))
    goid_enriched_clean <- goid_enriched[nchar(goid_enriched) != 0]
    if(length(goid_enriched_clean) == 0) next;
    ig <- set_vertex_attr(graph = ig, name = "process_id", 
                          index = V(ig)[name %in% names(goid_enriched_clean)],
                          value = goid_enriched_clean)
    
  }
  ig
}

renderModuleTable <- function(ig, lmods, emods, minGO = 5) { 
  
  requireNamespace("igraph")
  
  #only apply this to significant terms
  emods <- emods[emods$sig,]
  
  if(nrow(emods)  == 0) { 
    return((message("No significant enrichment of modules.")))
  }
  
  emods$modname <- as.character(emods$modname)
  #restrict mininum number of GO overlap to be reported. 
  emods <- emods[emods$ov >= minGO,]
  
  #find associated major CNA hubs
  modids <- emods$modname
  moduids <- unique(modids)
  emods$xhubs <- ""
  
  for(moduid in moduids) { 
    xmembers <- V(ig)[name %in% lmods[[moduid]] & levels %in% "x"]
    char_xmembers <- as_ids(xmembers)
    if(length(char_xmembers) == 0) { 
      next;
    } else { 
      ymembers <- V(ig)[name %in% lmods[[moduid]] & levels %in% "y"]
      xyedges <- as_ids(E(ig)[xmembers %--% ymembers])
      xhits <- sapply(char_xmembers, function(xm) sum(grepl(pattern = xm, x = xyedges)))
      xhits_be <- xhits[xhits >= minGO]
      if (length(xhits_be) == 0) next;
      #only report first two
      if(length(xhits_be) > 1) { 
        xhits_be <- sort(xhits_be, decreasing = TRUE)[1:2]  
      }
      
      if(!is.null(V(ig)$alias)) { 
        label <- "alias"  
      } else { 
        label <- "name"  
      }
      
      xhubs <- paste0(paste0(vertex_attr(ig, name = label, 
                                            index = V(ig)[name %in% names(xhits_be)]), 
                                " (", xhits_be, ")"), collapse = "; ")
      emods$xhubs[modids %in% moduid] <- xhubs
    }
  }
  modsize <- sapply(lmods,length)
  emods$modname_size <- sapply(emods$modname, function(mn) paste0(mn, " (", modsize[mn], ")"))
  #remove duplicates
  emods$modname_size[duplicated(emods$modname_size)] <- ""
  emods$xhubs[duplicated(emods$xhubs)] <- ""
  
  if (!is.null(emods$process_alias)) { 
    col2 <- "process_alias"  
  } else { 
    col2 <- "process_id"
  }
  
  colorder <- c("modname_size", col2, "representation", "xhubs")
  emods <- emods[,names(emods)[match(colorder, names(emods))]]
  names(emods) <- c("Module (size)" , "GO Category", "GO Obs./ Total ","X-hubs (hits)")
  emods
}

getModuleInfo <- function(mod_name, ig, go2eg, lmods, emods, info, cross_edge_tab, minGO = 5)  { 
  
  mod_refseq <- lmods[[mod_name]]
  #mod_name <- paste0(mod_name, " (", length(mod_refseq), ")")
  
  mod_attr <- info[match(mod_refseq, info$id),]
  
  eidx <- (emods$modname %in% mod_name) & (emods$sig) & (emods$ov >= minGO)
  if (sum(eidx) == 0)  { 
    mod_attr$go_enriched <- ""
  }  else { 
    
    stop("Fix me here.")
    if (!is.null(emods$process_alias)) { 
      col2 <- "process_alias"  
    } else { 
      col2 <- "process_id"
    }
    
    
    enriched_go_terms <- emods[eidx,"process_alias"]
    enriched_goids <- emods[eidx,"process_id"] 
    mod_attr$go_enriched <- sapply(mod_attr$id, function(prot) 
      paste0(enriched_go_terms[sapply(go2eg[enriched_goids], 
                                      function(g) prot %in% g)], collapse = "; "))
  }
  
  mod_attr_tab <- mod_attr[,c("cytoscape", "id", "chr", 
                              "start", "end", "strand", "go_enriched")]
  #make sure degree matches up right
  mod_attr_tab$degree <- igraph::degree(ig, v = V(ig)[name %in% mod_attr_tab$id ])
  stopifnot(identical(names(V(ig)[name %in% mod_attr_tab$id ]), mod_attr_tab$id))
  
  #hub perturbation
  hub_perturbers_id <- mod_attr$id[mod_attr$levels %in% "x"]
  hub_perturbers_cytoscape <- mod_attr$cytoscape[mod_attr$levels %in% "x"]
  prot_all_edges <- splitEdgeVector(as_ids(E(ig)), ig)
  map_perturbers <- prot_all_edges[prot_all_edges$left %in% hub_perturbers_id,]
  mod_attr_tab$CNA_perturbation <- sapply(mod_attr_tab$id, function(id) { 
    query <- map_perturbers$left[map_perturbers$right %in% id]
    if(length(query)) { 
      cytoscape_query <- mod_attr_tab$cytoscape[mod_attr_tab$id %in% query]
      paste(cytoscape_query, collapse = ";")
    } else { 
      "--"
    }
  })
  
  #E(ig)[ V(ig)[name %in% hub_perturbers_id] %--% map_perturbers$right ]$cis_trans
  xy_cis_edges <- as_ids(E(ig)[ E(ig)$cis_trans == "cis" & E(ig)$levels == "x-y"])
  mod_attr_tab$cis_perturbation <- sapply(mod_attr_tab$id, function(id) { 
    query <- map_perturbers$left[map_perturbers$right %in% id]
    if(length(query)) { 
      ct_query <- ifelse(paste0(query, "|", id) %in% xy_cis_edges, "c", "t")
      paste(ct_query, collapse = ";")
    } else { 
      "--"
    }
  })
  
  #Cross edge detection
  mod_attr_tab$cross_edge <- sapply(mod_attr$id, function(prot) { 
    paste0(cross_edge_tab$left[cross_edge_tab$right %in% prot], collapse = ";")
  })
  ce_exists <- nchar(mod_attr_tab$cross_edge) > 0
  mod_attr_tab$cross_edge_symbol <- "--"
  mod_attr_tab$cross_edge_module <- "--"
  if (sum(ce_exists) > 0) { 
    ce_match_index <- match(mod_attr_tab$cross_edge[ce_exists], V(ig)$name)
    mod_attr_tab$cross_edge_symbol[ce_exists] <- V(ig)$cytoscape[ce_match_index]  
    mod_attr_tab$cross_edge_module[ce_exists] <- V(ig)$module[ce_match_index]
  }
  
  names(mod_attr_tab) <- c("Node Alias", "Unique_ID", "Chromosome", "Start", "End", 
                           "Strand", "GO Enriched", "Degree", "CNA Perturbation", "Cis/Trans", 
                           "Cross Edge Node ID", "Cross Edge Node Alias", "Cross Edge Module")
  tab <- mod_attr_tab[order(mod_attr_tab$Degree, decreasing=T),
                      c(1,2,8,9,10,7,11,12,13,3,4,5,6)]
  tab
}