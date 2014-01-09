# OK
getLargestComp = function(g) {
	dg = decompose.graph(g)
	dg[order(sapply(dg, vcount), decreasing=TRUE)][[1]]
}

.reduceGraph = function (g, target, query, max.dist = 3)
{
  target = target[ target %in% V(g)$name ]
  query = query[ query %in% V(g)$name ]
  
  ## TODO: check the nodes exists...? (they should)
  
  res = list()
  for (t in target) {
    for (q in query) {
      if (t != q) {
        sp = get.all.shortest.paths(g, from = t, to = q)
        sp$res = sp$res[sapply(sp$res, function(r) length(r)) <= max.dist]
        res = c(res, sp$res)
      }
    }
  }
    
  v = unique(unlist(res))
  message("removing: ", vcount(g)-length(v), " nodes out of ", vcount(g), " [keeping ", length(v), " nodes]")
  if (length(v) > 0) 
    induced.subgraph(g, v)
  else NULL
}

.findNeigh = function(g, n) {
  i = V(g)[ n ]
	V(g)[ c(i, neighbors(g, i)) ]$name
}


##
# extended: allows to include the query nodes into the restrictions to .reduceGraph(), returning larger TRMs.
# strict: returns a single connected component by restricting it to the one containing the target node. If none it returns NULL
findTRM = function (g, target, query, method = "nsa", max.bridge = 1, extended = FALSE, strict = FALSE, type = "igraph") 
{
  type = match.arg(type, c("igraph", "graphNEL"))
  
  if (!any(target %in% V(g)$name))
    stop("none of the target nodes can be found in the network!")
  
  if (!any(query %in% V(g)$name))
    stop("none of the query nodes can be found in the network!")
  
  sel = target %in% V(g)$name
  target = target[sel]
  if (length(which(!sel)) > 0) 
      message(length(which(!sel)), " target nodes NOT FOUND in network-- removed")
  
  sel = query %in% V(g)$name
  query = query[sel]
  if (length(which(!sel)) > 0) 
    message(length(which(!sel)), " query nodes NOT FOUND in network-- removed")
  
  switch(method, all = {
      m = g
  }, first = {
      v = .findNeigh(g, target)
      m = induced.subgraph(g, v)
  }, second = {
      v = .findNeigh(g, target)
      all_v = v
      for (cv in all_v) {
          v = c(v, .findNeigh(g, cv))
      }
      m = induced.subgraph(g, unique(v))
  }, nsa = {
      v = unique(unlist(lapply(unique(c(target, query)), 
          function(x) .findNeigh(g, x))))
      m = induced.subgraph(g, v)
      if(vcount(m)>1) {
        if(extended)
          m = .reduceGraph(m, unique(c(target, query)), unique(c(target, query)), max.dist = max.bridge + 2)
        else
          m = .reduceGraph(m, target, query, max.dist = max.bridge + 2)
      }
  })
  if(!is.null(m)) {
    if(strict) {
      dm = decompose.graph(m)
      m = dm[sapply(dm, function(d) target %in% V(d)$name)][[1]]
    }
    m = annotateTRM(m, target = target) 
  }
  if(type == "graphNEL") m = igraph.to.graphNEL(m)
  m
}


removeVertices = function(g, filter, keep.hanging = FALSE) {
	g = delete.vertices(g, filter)
	if(keep.hanging) g
	else {
		getLargestComp(g)
	}
}