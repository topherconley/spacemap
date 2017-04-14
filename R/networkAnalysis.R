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
#' and type of node. The variables are:
#' id (a unique ID for the gene node); chr (chromosome location); 
#' start (starting coordinate on the chromosome); end (ending coordinate); 
#' cyto (cyotband abbreviation, e.g. 5q34); cytoscape (a gene symbol); 
#' strand (the +/- strand); type (either a "x" or "y" corresponding to predictor or response
#' nodes).  
#' @param go2eg named list where the names denote a biological process (e.g. GO ID) and
#' the elements of the list is a vector of members belonging to the biological process. 
#' The list ought to be non-redundant in names. 
#' @param cis_window numeric denoting the half length of the genomic interval 
#' used for calling a response node (i.e. of type y) as cis regulated by a predictor 
#' node (i.e. of type x)..
#' @param lmods named list of modules in the network. When \code{lmdods=NULL}, the 
#' Girvan-Newman edge betweenness algorithm is used to find network modules. 
#' @param min_module_size integer defining the minimum size a module must be to 
#' be tested for enrichment. 
#' @param mod_prefix character to prefix module ID's. 
#' @param fdr_thresh numeric between 0 and 1 to denote the FDR threshold in the 
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
#'  \item \code{modenrichtab} data.frame reporting which modules were enriched 
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
                 min_module_size = 15, mod_prefix = "P", fdr_thresh = 0.05, 
                 min_go = 5) { 
  
  #make igraph object
  message("Importing network attributes ...")
  ig <- adj2igraph(YY = net$yy, 	
                   XY = net$xy, 	
                   info = info)	
  #subset just the yy info for GO enrichment and cis/trans 
  #functions later
  infoy <- info[info$type == "y",]
  
  message("Distinguishing edges by node type ...")
  ig <- setEdgeTypeAttr(ig, "x", "y")
  #E(ig)[edge_type %in% "x--y"]
  ig <- setEdgeTypeAttr(ig, "y", "y")
  #E(ig)[edge_type %in% "y--y"]	
  
  message("Identify cis and trans regulation ...") 	
  ig <- setCisTransAttr(ig = ig, et = c("x--y", "y--y"))	
  #E(ig)$cis_trans	
  
  #obtain the distribution of cis/trans edges for 
  #reporting in X hub table
  vseqx <- V(ig)[type %in% "x"]	
  xdct <- setCisTransDistr(ig = ig, vseq = vseqx)	
  ig <- xdct$ig	
  dct <- xdct$dct	
  
  #Which CNA are most stabl?e	
  message("X hub prioritization")
  pc <- potentialCis(pact = as_ids(V(ig)[type %in% "x"]),	
                     preg = infoy$id, 	
                     info = info,	
                     YY = net$yy, XY = net$xy, 	
                     cw = cis_window, is = T)	
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
  igy <- adj2igraph(YY = net$yy, info = infoy, dropnull = F)	
  modenrich <- moduleEnrichment(lmods = yy_lmods, go2eg = go2eg, igy = igy, fdr_thresh = fdr_thresh)	
  #head(modenrich)
  #library(ggplot2)	
  #qplot(modenrich$pvalue) + xlab("P-value of Module Enrichment Analysis")	
  
  #REPORT TABLE
  #format the module enrichment results for displaying in a table.	 	
  modenrichtab <- renderModuleTable(ig, lmods, modenrich, minGO = min_go) 	

    #embed GO enrichment in igraph object 	
  ig <- setGOTags(lmods = lmods, modenrich = modenrich, 	
                  ig = ig, go2eg = go2eg, info = info)	

  #REPORT LIST OF TABLES
  #gather attributes of members each module in a data.frame 	
  lmod_attr <- lapply(names(lmods), getModuleInfo,  	
                      ig = ig,	
                      go2eg = go2eg,	
                      lmods = lmods,	
                      modenrich = modenrich,	
                      info = info, 	
                      cross_edge_tab = cross_edge_tab, 	
                      minGO = min_go)	
  names(lmod_attr) <- names(lmods)	
  
  #return: xhubtab 
  ret <- list(ig  = ig, 
              xhubtab = xhubtab, 
              modenrichtab = modenrichtab,
              lmodtab = lmod_attr, 
              goneighbors = hgp)
  return(ret)
}

#' Initialize network's vertex attributes
#' 
#' 
#' 
#' @description Label vertex attributes of a  conditional graphical model's
#'  output. 
#' @param YY adjacency matrix encoding y--y edges of conditional graphical model. If a 
#' logical matrix, TRUE denotes the existence of an edge. Format should conform to 
#' \code{igraph::graph_from_adjacency_matrix} input. 
#' @param XY adjacency matrix encoding x--y edges of conditional graphical model. If NULL,
#' then a network encoding just y--y edges is made and \code{info} should just contain attributes for (y(1), ..., y(Q)).  
#' @param info data.frame encoding the attributes of the vertices 
#' in the order of (x(1), ..., x(p), y(1), ..., y(q)) if both YY and XY are specified.
#' If just YY is specified, then \code{info} should just contain attributes for (y(1), ..., y(Q)). 
#' See details for required attributes to include and optional attributes. 
#' @param weighted logical indicating whether the values of YY and XY should indicate 
#' weights of the vertices. Defaults to FALSE.
#' @param dropnull logical indicating to remove vertices with degree zero from the graph. Defaults to TRUE.
#' @export
adj2igraph <- function(YY, XY = NULL, info, weighted = F, dropnull = TRUE) { 

  if (!is.null(XY)) { 
    adj_yyxy <- rbind(cbind(matrix(FALSE, nrow(XY), nrow(XY)), XY), 
                      cbind(t(XY), YY))
  } else { 
    adj_yyxy <- YY 
  } 
    
  if (weighted) { 
    ig <- graph_from_adjacency_matrix(adj_yyxy, mode = "max", weighted = weighted, diag = FALSE, add.colnames = NA) 
  } else { 
    ig <- graph_from_adjacency_matrix(adj_yyxy, mode = "max", diag = FALSE, add.colnames = NA)  
  }
  
  #label with info
  for(var in names(info)) { 
    if (var == "id") { 
      ig <- igraph::set_vertex_attr(graph = ig, name = "name", value = info[,var])
    } else { 
      ig <- igraph::set_vertex_attr(graph = ig, name = var, value = info[,var])
    }
  }
  
  #remove vertices with no edges
  if(dropnull) { 
    ig <- igraph::delete.vertices(graph  = ig, v = V(ig)[igraph::degree(ig) == 0])
  }
  ig
}

iGRanges <- function(ig, vtype = NULL, cw, is) { 
  library(foreach)
  #get gene coordinates
  gcnames <- c("chr", "start", "end", "strand", "name")
  names(gcnames) <- gcnames
  gc <- lapply(gcnames, function(gcn) { 
    if (is.null(vtype)) { 
      vertex_attr(graph = ig, name = gcn) 
    } else { 
      vertex_attr(graph = ig, name = gcn, 
                  index = V(ig)[type %in% vtype])
    }
  })
  dgc <- as.data.frame(gc, stringsAsFactors = F)
  
  
  #define features for which missing values is not tolerated 
  if(is) {
    #allow for missing values in the strand
    mfi <-   c("chr", "start", "end")
  } else { 
    #do not allow for missing values in the strand
    mfi <-   c("chr", "start", "end", "strand")
  }
  #remove features that have missing values in genome coordinates that 
  #does not permite cis/trans identification
  notNA <- rowSums(!is.na(dgc[,gcnames[mfi]])) == length(mfi)  
  dgc <- dgc[notNA,]
  dgc$strand <- as.character(dgc$strand)
  dgc$strand <- ifelse(is.na(dgc$strand), "*", dgc$strand)
  if(nrow(dgc) == 0) { 
    stop("Too many missing values in gene coordinates to compute cis/trans activity.")
  }
  
  suppressMessages(library(GenomicRanges))
  gr <- GRanges(seqnames = Rle(dgc$chr),
                ranges = IRanges(start = pmax(dgc$start - cw, 1), 
                                 end = dgc$end + cw, 
                                 names = dgc$name),
                strand = Rle(strand(dgc$strand)))
  gr
}

#' Label edges with type attribute
#' 
#' 
setEdgeTypeAttr <- function(ig, x="x", y="y") { 
  yidx <- V(ig)[type %in% y]
  xidx <- V(ig)[type %in% x]
  set_edge_attr(graph = ig, name = "edge_type", index = E(ig)[xidx %--% yidx], value = paste0(x,"--",y))
}

splitEdgeVector <- function(evec, ig) { 
  lev <- strsplit(evec, "|", fixed  = TRUE)
  tab <- data.frame(left = sapply(lev, function(x) x[1]),
                    right = sapply(lev, function(x) x[2]), stringsAsFactors = F)
  
  #no nodes of type x should be in the right column, 
  #flip to the left column so that it is always (x, y) or (y, y)
  rev_idx <- which(tab$right %in% as_ids(V(ig)[type %in% "x"]))
  if(length(rev_idx) > 0) { 
    tab_rev <- tab[rev_idx,]
    tab[rev_idx,] <- tab[rev_idx,c("right", "left")]
  }
  tab
}

#as_ids(E(ig)[type %in% "x--y"])
#ldf <- as_long_data_frame(ig)
#head(ldf)

setCisTransAttr <- function(ig, et = NULL, cw = 2e6, is = TRUE) { 
  if(!is.null(et)) { 
    evec <- as_ids(E(ig)[edge_type %in% et])
    vtype <- unique(unlist(strsplit(x = et, split = "--")))
  } else { 
    evec <- as_ids(E(ig))
    vtype <- NULL
  }
  de <- splitEdgeVector(evec, ig)
  gr <- iGRanges(ig=ig, vtype=vtype, cw=cw, is=is)
  left <- de$left[de$left %in% names(gr)]
  right <- de$right[de$right %in% names(gr)]
  ov <- findOverlaps(gr[left,], gr[right,], ignore.strand = is)
  candidateCis <- paste0(left[ov@from], "|", right[ov@to])
  cis <- evec[evec %in%  candidateCis]
  trans <- setdiff(evec, cis)
  ig <- set_edge_attr(ig, name = "cis_trans", 
                     index = E(ig)[ as_ids(E(ig)) %in% cis], 
                     value = "cis")
  ig <- set_edge_attr(ig, name = "cis_trans", 
                      index = E(ig)[ as_ids(E(ig)) %in% trans], 
                      value = "trans")
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
      cisgenes <- paste(V(ig)[match(chits, V(ig)$name)]$cytoscape, 
                        collapse = ", ")
    } else { 
      cisgenes <- "--"
    }
    c(ncis = ncis, ntrans = length(vct) - ncis, cisgenes = cisgenes)
  }))
  dct <- as.data.frame(dct, stringsAsFactors = F)
  dct$ncis <- as.integer(dct$ncis)
  dct$ntrans <- as.integer(dct$ntrans)
  idx <- match(rownames(dct), V(ig)$name)
  ig <- set_vertex_attr(graph = ig, name = "ncis", index = V(ig)[idx], value = dct$ncis)
  ig <- set_vertex_attr(graph = ig, name = "ntrans", index = V(ig)[idx], value = dct$ntrans)
  list(ig = ig, dct = dct)
}


potentialCis <- function(pact, preg, info, XY, YY, cw, is) { 
  
  ig <- adj2igraph(YY = YY, 
                   XY = XY,
                   info = info, 
                   dropnull = F)
  gr <- iGRanges(ig=ig, vtype=NULL, cw=cw, is=is)
  library(GenomicRanges)
  pactgc <- pact[pact %in% names(gr)]
  preggc <- preg[preg %in% names(gr)]
  pregquery <- GenomicRanges::findOverlaps(query = gr[preggc,], subject = gr[pactgc,], ignore.strand = TRUE)
  pactgctab <- table(pactgc[pregquery@to])
  potcis <- vector(mode = "numeric", length = length(pact))
  names(potcis) <- pact
  potcis[match(names(pactgctab), names(potcis))] <- pactgctab
  potcis
}

renderXhubTable <- function(ig, dct, pc, bdeg, sig = 2, unit_base = 1e6, unit = "Mb") { 
  idx <- match(rownames(dct), V(ig)$name)
  tab <- data.frame(Cytoband = V(ig)[idx]$cytoscape, 
                    chr = V(ig)[idx]$chr,
                    start = V(ig)[idx]$start, 
                    end = V(ig)[idx]$end, 
                    `Cis/Trans` =  paste0(dct$ncis, " / ", dct$ntrans), 
                    `Potential Cis`=  pc[match(names(pc), rownames(dct))],
                    `Cis Proteins` = dct$cisgenes)
  
  ###degree-based ranking
  #rank the degree's for each bootstrapped degree distribution
  #smaller rank is larger degree, larger
  brank <- apply(-bdeg, 1, rank)
  #average the ranks across the bootstrap replicates
  avg_bdeg_ranks <- rowMeans(brank)
  sd_bdeg_ranks <- apply(brank, 1, sd)
  bdeg_pres <-  paste(round(avg_bdeg_ranks,2), " (", round(sd_bdeg_ranks, 2), ")")
  names(bdeg_pres) <- names(avg_bdeg_ranks)
  tab$bdeg_pres <- bdeg_pres[rownames(tab)]
  tab$rank_bootstrap_rank <- rank(avg_bdeg_ranks)[rownames(tab)]
  tab <- tab[order(tab$rank_bootstrap_rank, tab$chr, tab$start),]
  
  #add in gene coordinate location
  start <- tab$start
  end <-  tab$end
  ub_start <- signif(start / unit_base,sig)
  ub_end <- signif(end/ unit_base, sig)
  tab$Cytoband <- paste(tab$Cytoband, 
                        paste0("(", paste(ub_start, ub_end, sep = "-"), 
                               paste0(" ", unit, ")")))
  #head(tab, 40)
  tab <- tab[,!(names(tab) %in%"rank_bootstrap_rank")]
  names(tab) <- c("Cytoband", "Chromosome", "Start", "End",
                  "Cis/Trans", "Potential Cis",
                  "Cis Proteins", "Rank of CNA degree")
  tab
}


hubGOproportion <- function(ig, go2eg) { 
  xytab <- splitEdgeVector(igraph::as_ids(E(ig)[edge_type %in% "x--y"]), ig)
  cna_nbrhood <- split(x = xytab$right, f = as.factor(xytab$left))
  gt1_cna_nbrhood <- cna_nbrhood[sapply(cna_nbrhood, length) > 1]
  cna_prop <- lapply(gt1_cna_nbrhood, function(gt1) { 
    
    #make a go adjacency matrix
    ngt <- length(gt1)
    mgo <- matrix(0, nrow = ngt, ncol = ngt)
    pw <- combn(1:ngt,2)
    for(j in 1:ncol(pw)) { 
      ep <- gt1[pw[,j]]
      mgo[pw[1,j],pw[2,j]] <- mgo[pw[2,j],pw[1,j]] <- any(sapply(go2eg, function(egs) all(ep %in% egs)))
    }
    m1 <- sum(mgo[upper.tri(mgo)]) / ncol(pw)
    m2 <- sum(colSums(mgo)) / ngt
    m3 <- sum(colSums(mgo) > 0) / ngt
    
    # pw <- combn(gt1, 2)
    # ngopairs <- sum(sapply(1:ncol(pw), function(i) { 
    #   ep <- pw[,i]
    #   any(sapply(go2eg, function(egs) all(ep %in% egs)))
    # }))
    # data.frame(total_edges = ncol(pw), ngopairs  = ngopairs, prop = ngopairs/ncol(pw))
    data.frame(m1 = m1, m2 = m2, m3 = m3)
  })
  df_cna_prop <- do.call(rbind, cna_prop)
  df_cna_prop$id <- names(gt1_cna_nbrhood)
  vids <- V(ig)[name %in% df_cna_prop$id]
  cytoband <- igraph::vertex_attr(graph = ig, name = "cyto", index = vids)
  deg <- df_cna_prop$degree <- igraph::degree(graph = ig, v = vids)
  char_vids <- as_ids(vids)
  correct_order_idx <- match(df_cna_prop$id, char_vids)
  df_cna_prop$degree <- deg[correct_order_idx]
  df_cna_prop$cytoband <- cytoband[correct_order_idx]
  out <- df_cna_prop[,c("cytoband", "id", "degree", "m3")]
  rownames(out) <- NULL
  names(out) <- c("Cytoband", "ID", "Out_Degree", "Proportion")
  out
}

rndCGGraph <- function(xy, yy, info) { 
  nz_xy <- which(rowSums(xy) > 0)
  nz_yy <- which(rowSums(yy) > 0)
  rxy <- xy
  Q <- ncol(xy)
  rxy[nz_xy,] <- t(sapply(nz_xy, function(i) sample(x = xy[i,], size = Q, replace = F)))
  igy <- adj2igraph(YY = yy, XY = NULL, info = info[info$type == "y",], dropnull = FALSE)
  rigy <- rewire(graph = igy, with = keeping_degseq(niter = vcount(igy) * 10))
  ryy <- as_adj(graph = rigy)
  rig <- adj2igraph(YY = ryy, XY = rxy, info = info, weighted = F)
  rig <- setEdgeTypeAttr(rig, "x", "y")
  rig <- setEdgeTypeAttr(rig, "y", "y")
  #rig <- setCisTransAttr(ig = rig, et = c("x--y", "y--y"))
  rig
}

rndHubGoProportion <- function(xy, yy, info, go2eg) { 
  rig <- rndCGGraph(xy = xy, yy = yy, info = info)
  rhgp <- hubGOproportion(rig, go2eg)
  #weighted.mean(x = rhgp$`Proportion of GO Pairs`, w = rhgp$`Out Degree`)
  colMeans(rhgp[,c("m1","m2","m3")])
}

cebModules <- function(ig, glb = 15, prefix = "M") { 
  ceb <- cluster_edge_betweenness(ig)
  mods <- membership(ceb)
  #Group the modules and report the modules that are of at least size glb.
  modsBigEnough <- sizes(ceb)[sizes(ceb) >= glb]
  
  lmods <- lapply(names(modsBigEnough), function(mbe) { 
    mbe <- as.integer(mbe)
    names(mods)[mods == mbe]
  })
  names(lmods) <- sapply(seq_along(lmods), function(i) paste0(prefix, i))
  yy_lmods <- lapply(lmods, function(mod) { 
    as_ids(V(ig)[(type %in% "y") & name %in% mod])
  })
  
  #each module is a unique list
  stopifnot(all(sapply(lmods, anyDuplicated) == 0))
  
  for(i in seq_along(lmods)) {
    ig <- set_vertex_attr(graph = ig, name= "module", 
                          index = V(ig)[ name %in% lmods[[i]] ], 
                          value = names(lmods)[i])
  }
  
  
  
  
  #?communities
  #?compare.communities
  list(ig = ig, modularity_score = modularity(ceb), 
       module_list = lmods, yy_module_list = yy_lmods, community_obj = ceb, glb = glb)
}

moduleEnrichment <- function(lmods, go2eg, igy, fdr_thresh = 0.05) { 
  modgo <- foreach(i = seq_along(lmods)) %:%
    foreach(j = seq_along(go2eg), .combine = 'rbind') %do% { 
      mod <- lmods[[i]]
      gt <- go2eg[[j]]
      #overlap
      ov <- length(intersect(mod,gt))
      #minus -1 because of the upper tail
      pval <- phyper(q = ov - 1L, 
                     m = length(gt), 
                     n =  vcount(igy) - length(gt), 
                     k = length(mod),lower.tail = FALSE)
      data.frame(goid = names(go2eg)[j], pvalue = pval, ov = ov,
                 representation = paste(ov, length(gt), sep = "/"),
                 modname = names(lmods)[i])
    }
  amodgo <- do.call(rbind, modgo)
  #multiple-testing correction with Benjamini Hochberg
  amodgo$fdr <- p.adjust(p = amodgo$pvalue, method = "BH")
  #identify significant terms
  amodgo$sig <- amodgo$fdr < fdr_thresh
  #order by module name and then by significance
  amodgo <- amodgo[order(amodgo$modname, amodgo$fdr),]
  library(GO.db)
  amodgo$GOTERM <- Term(as.character(amodgo$goid))
  amodgo <- amodgo[,c("modname", "goid", "GOTERM", "representation", "fdr", "pvalue", "sig", "ov")]
  amodgo
}

################

setGOTags <- function(lmods, modenrich, ig, go2eg, info) { 
  for (mod_name in names(lmods)) {
    mod_refseq <- lmods[[mod_name]]
    mod_attr_id <- info[match(mod_refseq, info$id),"id"]
    
    eidx <- (modenrich$modname %in% mod_name) & modenrich$sig
    if (sum(eidx) == 0) next; 
    enriched_goids <- modenrich[,"goid"]
    #record enriched GOID's into the igraph object
    
    goid_enriched <-  sapply(mod_attr_id, function(prot) 
      paste0(enriched_goids[sapply(go2eg[enriched_goids], 
                                   function(g) prot %in% g)], collapse = "; "))
    goid_enriched_clean <- goid_enriched[nchar(goid_enriched) != 0]
    if(length(goid_enriched_clean) == 0) next;
    ig <- set_vertex_attr(graph = ig, name = "GO_ID", 
                          index = V(ig)[name %in% names(goid_enriched_clean)],
                          value = goid_enriched_clean)
    
  }
  ig
}


renderModuleTable <- function(ig, lmods, modenrich, minGO = 5) { 
  
  #only apply this to significant terms
  modenrich <- modenrich[modenrich$sig,]
  
  if(nrow(modenrich)  == 0) { 
    return((message("No significant enrichment of modules.")))
  }
  
  modenrich$modname <- as.character(modenrich$modname)
  #restrict mininum number of GO overlap to be reported. 
  modenrich <- modenrich[modenrich$ov >= minGO,]
  
  #find associated major CNA hubs
  modids <- modenrich$modname
  moduids <- unique(modids)
  modenrich$cna_hubs <- ""
  
  for(moduid in moduids) { 
    xmembers <- V(ig)[name %in% lmods[[moduid]] & type %in% "x"]
    char_xmembers <- as_ids(xmembers)
    if(length(char_xmembers) == 0) { 
      next;
    } else { 
      ymembers <- V(ig)[name %in% lmods[[moduid]] & type %in% "y"]
      xyedges <- as_ids(E(ig)[xmembers %--% ymembers])
      xhits <- sapply(char_xmembers, function(xm) sum(grepl(pattern = xm, x = xyedges)))
      xhits_be <- xhits[xhits >= minGO]
      if (length(xhits_be) == 0) next;
      #only report first two
      if(length(xhits_be) > 1) { 
        xhits_be <- sort(xhits_be, decreasing = TRUE)[1:2]  
      }
      cna_hubs <- paste0(paste0(vertex_attr(ig, name = "cyto", 
                                            index = V(ig)[name %in% names(xhits_be)]), 
                                " (", xhits_be, ")"), collapse = "; ")
      modenrich$cna_hubs[modids %in% moduid] <- cna_hubs
    }
  }
  modsize <- sapply(lmods,length)
  modenrich$modname_size <- sapply(modenrich$modname, function(mn) paste0(mn, " (", modsize[mn], ")"))
  #remove duplicates
  modenrich$modname_size[duplicated(modenrich$modname_size)] <- ""
  modenrich$cna_hubs[duplicated(modenrich$cna_hubs)] <- ""
  
  colorder <- c("modname_size", "GOTERM", "representation", "cna_hubs")
  modenrich <- modenrich[,names(modenrich)[match(colorder, names(modenrich))]]
  names(modenrich) <- c("Module (size)" , "GO Category", "GO Obs./ Total ","CNA hubs (hits)")
  modenrich
}


getModuleInfo <- function(mod_name, ig, go2eg, lmods, modenrich, info, cross_edge_tab, minGO = 5)  { 
  
  mod_refseq <- lmods[[mod_name]]
  #mod_name <- paste0(mod_name, " (", length(mod_refseq), ")")
  
  mod_attr <- info[match(mod_refseq, info$id),]
  
  eidx <- (modenrich$modname %in% mod_name) & (modenrich$sig) & (modenrich$ov >= minGO)
  if (sum(eidx) == 0)  { 
    mod_attr$go_enriched <- ""
  }  else { 
    enriched_go_terms <- modenrich[eidx,"GOTERM"]
    enriched_goids <- modenrich[eidx,"goid"] 
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
  hub_perturbers_id <- mod_attr$id[mod_attr$type %in% "x"]
  hub_perturbers_cytoscape <- mod_attr$cytoscape[mod_attr$type %in% "x"]
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
  xy_cis_edges <- as_ids(E(ig)[ E(ig)$cis_trans == "cis" & E(ig)$edge_type == "x--y"])
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