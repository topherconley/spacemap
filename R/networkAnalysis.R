#' Network analysis toolkit
#' 
#' Downstream analysis of spaceMap networks
#'
#' @param net list containing two adjacency matrices, 
#' named \code{yy} and \code{xy}, encoding the 
#' response--response interactions and the predictor--response 
#' interactions, respectively. 
#' @param bdeg list containing two matrices, \code{xy} and \code{yy}, 
#' where each row represents a bootstrap relicate of the degree distribution and 
#' the columns contain the degree across all x nodes (in the case of \code{xy} )
#'  or all y nodes (in the case of \code{yy}). 
#' @param info data.frame containing several columns that describe the 
#' attributes of the nodes in the network. The columns are currently required to 
#' the following variables, but they can be missing values for all variables but the ID
#' and levels of node. The variables are:
#' id (a unique ID for the gene node); chr (chromosome location); 
#' start (starting coordinate on the chromosome); end (ending coordinate); 
#' cyto (cyotband abbreviation, e.g. 5q34); cytoscape (a gene symbol); 
#' strand (the +/- strand); levels (either a "x" or "y" corresponding to predictor or response
#' nodes).  
#' @param go2eg named list where the names denote a biological process (e.g. GO ID) and
#' the elements of the list is a vector of members belonging to the biological process. 
#' The list ought to be non-redundant in names. 
#' @param cis_window numeric denoting the half length of the genomic interval 
#' used for calling a response node (i.e. of levels y) as cis regulated by a predictor 
#' node (i.e. of levels x)..
#' @param lmods named list of modules in the network. When \code{lmdods=NULL}, the 
#' Girvan-Newman edge betweenness algorithm is used to find network modules. 
#' @param min_module_size integer defining the minimum size a module must be to 
#' be tested for enrichment. 
#' @param mod_prefix character to prefix module ID's. 
#' @param thresh numeric between 0 and 1 to denote the thresh threshold in the 
#' Bonferroni-Hochberg correction. 
#' @param min_go integer defining the minimum number of nodes required to be 
#' to call a significant enrichment of a biological process based on the 
#' hypergeometric test.  
#' 
#' @return list
#' \itemize{
#'  \item \code{ig} The igraph object encoding the network and its attributes. 
#'  \item \code{xhubtab} data.frame reporting prioritized predictor hubs and 
#'  their attributes.
#'  \item \code{emodstab} data.frame reporting which modules were enriched 
#'  with biological processes and what they were. If the module contains 
#'  predictor nodes (e.g. CNA) with at least \code{min_go} edges to other 
#'  module members it is also reported for the module.
#'  \item \code{lmodtab} list of data.frames corresponding to each module of size at least 
#'  \code{min_module_size}. Each data.frame reports the nodes of the module and 
#'  other information such as degree and which nodes participate in specific enrichment
#'  of biological processes. 
#'  \item \code{goneighbors} data.frame where each row corresponds to a predictor hub. 
#'  The predictor hub induces a neighborhood of response nodes with which it has an 
#'  edge. A certain proportion of response node neighbors share the same 
#'  biological processes and this is reported.  
#'  }
neta <- function(net, bdeg, info, go2eg, cis_window = 2e6, lmods = NULL, 
                 min_module_size = 15, mod_prefix = "P", thresh = 0.05, 
                 min_go = 5) { 
  
  #make igraph object
  message("Importing network attributes ...")
  ig <- adj2igraph(yy = net$yy, 	
                   xy = net$xy, 	
                   info = info)	
  #subset just the yy info for GO enrichment and cis/trans 
  #functions later
  infoy <- info[info$levels == "y",]
  
  message("Identify cis and trans regulation ...") 	
  ig <- xyCisTrans(ig = ig, et = c("x-y", "y-y"))	
  #E(ig)$cis_trans	
  
  #obtain the distribution of cis/trans edges for 
  #reporting in X hub table
  vseqx <- V(ig)[levels %in% "x"]	
  xdct <- setCisTransDistr(ig = ig, vseq = vseqx)	
  ig <- xdct$ig	
  dct <- xdct$dct	
  
  #Which CNA are most stabl?e	
  message("X hub prioritization")
  pc <- potentialCis(pact = as_ids(V(ig)[levels %in% "x"]),	
                     preg = infoy$id, 	
                     info = info,	
                     yy = net$yy, xy = net$xy, 	
                     cw = cis_window, ignoreStrand = T)	
  #REPORT TABLE 
  xhubtab <- renderXhubTable(ig = ig, dct = dct, pc = pc, bdeg = bdeg$xy)	
  
  message("Identifying modules ... ")
  #Obtain network modules using the `Girvan-Newman` algorithm for clust	
  if(is.null(lmods)) { 
    cebout <- cebModules(ig, glb = min_module_size, prefix = mod_prefix)	
    ig <- cebout$ig	
    lmods <- cebout$module_list	
    yy_lmods <- cebout$yy_module_list	
    ceb <- cebout$community_obj	
    cross_edges <- crossing(communities = ceb, graph = ig)	
    cross_edges2 <- names(cross_edges)[cross_edges]	
    cross_edge_tab <- splitEdgeVector(cross_edges2, ig)	
  } else { 
    message("User supplied module list not yet supported.")
    stop()
  }
  
  message("GO enrichment analysis")
  #GO neighbor enrichment:
  #Identify the GO pair proportion for each CNA hub by considering any two response nodes with a common CNA hub as adjacent.	
  hgp <- hubGOproportion(ig, go2eg)	
  
  #GO module enrichment:
  #need to also consider nodes with no edges so use a different network 
  #object
  igy <- adj2igraph(yy = net$yy, info = infoy, dropnull = F)	
  emods <- moduleEnrichment(lmods = yy_lmods, go2eg = go2eg, igy = igy, thresh = thresh)	
  #head(emods)
  #library(ggplot2)	
  #qplot(emods$pvalue) + xlab("P-value of Module Enrichment Analysis")	
  
  #REPORT TABLE
  #format the module enrichment results for displaying in a table.	 	
  emodstab <- renderModuleTable(ig, lmods, emods, minGO = min_go) 	

    #embed GO enrichment in igraph object 	
  ig <- setGOTags(lmods = lmods, emods = emods, 	
                  ig = ig, go2eg = go2eg, info = info)	

  #REPORT LIST OF TABLES
  #gather attributes of members each module in a data.frame 	
  lmod_attr <- lapply(names(lmods), getModuleInfo,  	
                      ig = ig,	
                      go2eg = go2eg,	
                      lmods = lmods,	
                      emods = emods,	
                      info = info, 	
                      cross_edge_tab = cross_edge_tab, 	
                      minGO = min_go)	
  names(lmod_attr) <- names(lmods)	
  
  #return: xhubtab 
  ret <- list(ig  = ig, 
              xhubtab = xhubtab, 
              emodstab = emodstab,
              lmodtab = lmod_attr, 
              goneighbors = hgp)
  return(ret)
}

#' Initialize network's vertex attributes
#' 
#' @description Label vertex attributes of a  conditional graphical model's
#'  output. 
#' @param yy adjacency matrix encoding y--y edges of conditional graphical model. Format can be 
#' @param xy adjacency matrix encoding x--y edges of conditional graphical model. If NULL,
#' then a network encoding just y--y edges is made and \code{info} should just contain attributes for (y(1), ..., y(Q)).  
#' @param info data.frame encoding the attributes of the vertices 
#' in the order of (x(1), ..., x(p), y(1), ..., y(q)) if both yy and xy are specified.
#' If just yy is specified, then \code{info} should just contain attributes for (y(1), ..., y(Q)). 
#' See details for required attributes to include and optional attributes. 
#' @param weighted Same argument passed to \code{igraph::graph_from_adjacency_matrix}. Defaults to NULL.
#' weights of the vertices. Defaults to FALSE.
#' @export
adj2igraph <- function(yy, xy = NULL, yinfo = NULL, xinfo = NULL, weighted = NULL) { 
  
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
    library(Matrix)
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
    ig <- layerInfo(ig = ig, info = xinfo, vs =  V(ig)[Xindex])
    
  }
  if (!is.null(yinfo)) { 
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
  for(var in names(info)) { 
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

rankHub <- function(ig, bdeg = NULL, level = c("y", "x")) { 
  
  level <- match.arg(level)
  vs <- V(ig)[levels %in% level]
  
  if (is.null(bdeg)) { 
    deg <- igraph::degree(graph = ig, v = vs)
    ig <- set_vertex_attr(graph = ig, name = "mean_rank_hub", index = vs, value = rank(-deg))
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

cisTrans <- function(ig, level = c("x-y"), cw = 2e6, ignoreStrand = TRUE) { 
  
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

  suppressMessages(library(GenomicRanges))
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
  library(igraph)
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


reportHubs <- function(top = 10, level = "x", unit = "Mb", unit_base = 1e6, sig = 2) {
  
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

xHubEnrich <- function(ig, go2eg) { 
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

rndCGGraph <- function(xy, yy, info) { 
  nz_xy <- which(rowSums(xy) > 0)
  nz_yy <- which(rowSums(yy) > 0)
  rxy <- xy
  Q <- ncol(xy)
  rxy[nz_xy,] <- t(sapply(nz_xy, function(i) sample(x = xy[i,], size = Q, replace = F)))
  igy <- adj2igraph(yy = yy, xy = NULL, info = info[info$levels == "y",], dropnull = FALSE)
  rigy <- rewire(graph = igy, with = keeping_degseq(niter = vcount(igy) * 10))
  ryy <- as_adj(graph = rigy)
  rig <- adj2igraph(yy = ryy, xy = rxy, info = info, weighted = F)
  rig <- setEdgeTypeAttr(rig, "x", "y")
  rig <- setEdgeTypeAttr(rig, "y", "y")
  #rig <- xyCisTrans(ig = rig, et = c("x-y", "y-y"))
  rig
}

rndHubGoProportion <- function(xy, yy, info, go2eg) { 
  rig <- rndCGGraph(xy = xy, yy = yy, info = info)
  rhgp <- hubGOproportion(rig, go2eg)
  #weighted.mean(x = rhgp$`Proportion of GO Pairs`, w = rhgp$`Out Degree`)
  colMeans(rhgp[,c("m1","m2","m3")])
}

modAnalysis <- function(ig, mods, go2eg, glb = 15, prefix = "M", 
                        thresh = 0.05, adjust = "BH", minGO = 5, 
                        process_alias = NULL) { 

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
  
  #hypothesis testing
  emods <- moduleEnrichment(ig = ig, lmods = ylmods, go2eg = go2eg,
                            thresh = thresh, adjust = adjust,
                            process_alias = process_alias)
  #make pretty result
  etab <- renderModuleTable(ig = ig, lmods = lmods, 
                            emods = emods, minGO = minGO)
  #store result in igraph
  ig <- setGOTags(lmods = lmods, emods = emods, ig = ig, go2eg = go2eg)
  
  list(ig = ig, etab = etab, eraw = emods)
}

mods2list <- function(ig, mods, glb, prefix) { 
  
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

moduleEnrichment <- function(ig, lmods, go2eg, thresh, adjust, process_alias) { 
  
  #number of y nodes 
  ny <- sum(V(ig)$levels == "y")
  library(foreach)
  modgo <- foreach(i = seq_along(lmods)) %:%
    foreach(j = seq_along(go2eg), .combine = 'rbind') %do% { 
      mod <- lmods[[i]]
      gt <- go2eg[[j]]
      #overlap
      ov <- length(intersect(mod,gt))
      #minus -1 because of the upper tail
      pval <- phyper(q = ov - 1L, 
                     m = length(gt), 
                     n = ny - length(gt), 
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
  names(emods) <- c("Module (size)" , "GO Category", "GO Obs./ Total ","X hubs (hits)")
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