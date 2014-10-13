## Package: rTRM (c)
##
## File: plot.R
## Description: support plotting functions for rTRM
##
## Author: Diego Diez
## Contact: diego10ruiz@gmail.com
##


## plots the degree distribution of the graph object.
plotDegree = function(g) {
	d = table(igraph::degree(g))
	x = as.numeric(names(d))
	y = as.numeric(d)
		
	x0 = x[! x == 0 ]
	y0 = y[! x == 0 ]
	
	plot(x0, y0, log = "xy", xlab = "Degree", ylab = "# Nodes", pch = 21, bg = "gray")
	
	dl = lm(log10(y0) ~ log10(x0))
	abline(dl)
	l = .getLegend(dl)
	legend("topright", l, bg = "white", lty = "solid", cex = 0.7)
}

.getLegend <- function(m) {
	s <- summary(m)
	r2 = format(s$r.squared, nsmall = 2, digits = 2)
	a = format(m$coef[1], nsmall = 2, digits = 2)
	if (m$coef[2] > 0) {
		b = format(m$coef[2], nsmall = 2, digits = 2)
		l <- eval(substitute(expression(paste(R^2 == r2, 
			", ", y == a + b * x)), list(r2 = r2, a = a, 
			b = b)))
	}
	else {
		b = format(abs(m$coef[2]), nsmall = 2, digits = 2)
		l <- eval(substitute(expression(paste(r^2 == r2, 
			", ", y == a - b * x)), list(r2 = r2, a = a, 
			b = b)))
	}
}


.lum = function(col) {
    tmp = col2rgb(col)
    apply(tmp, 2, function(x) {
    (0.2126*x[1]) + (0.7152*x[2]) + (0.0722*x[3])
    })
}


.lum2 = function (col) 
{
    tmp = col2rgb(col)
    apply(tmp, 2, function(x) {
        (0.2126 * x[1]^2) + (0.7152 * x[2]^2) + (0.0722 * x[3]^2)
    })
}


.getNiceColors = function (n, eq.spaced = TRUE, order.by = TRUE, min.lum, max.lum) 
{
    cs = grep("^grey|^gray", colors(), value = TRUE, invert = TRUE)
    cs = grep("white", cs, value = TRUE, invert = TRUE)
    cols = c(cs, grey.colors(10))
  
    cols.lum = .lum2(cols)
      
    if (order.by) {
    	#cols.lum = lum2(cols)
      sel = order(cols.lum, decreasing = TRUE)
      cols = cols[sel]
      cols.lum = cols.lum[sel]
    }
    
    if(!missing(min.lum))
    	cols = cols[cols.lum >= min.lum]
    
    if(!missing(max.lum))
    	cols = cols[cols.lum <= max.lum]
    
    if(!missing(n)) {
    	if(eq.spaced) {
	    	ns = floor(length(cols)/n)
    		sel = seq(1, ns * n, ns)
    		cols = cols[sel]
    	} else cols = cols[1:n]
    }
    cols
}

.getTFclassColors = function(subset = "Class") {
  terms = getTFterms(subset)
  col = .getNiceColors(length(unique(terms)))
  names(col) = unique(terms)
  col
}

.annotateTFclass = function(g) {
  family = getTFclassFromEntrezgene(V(g)$name)
  family[sapply(family, is.null)] = "unclassified"
  V(g)$family = family
  
  col = c("unclassified" = "white", .getTFclassColors())
  V(g)$pie.color = lapply(V(g)$family, function(f) col[f])
  V(g)$pie = lapply(V(g)$pie.color, function(f) rep(1,length(f)))
  g
}

.annotateTarget = function(g, target) {
  target = target[ target %in% V(g)$name ]
  V(g)$frame.color = "grey"
  V(g)[ target ]$frame.color = "black"
  V(g)$frame.width = 1
  V(g)[ target ]$frame.width = 2
  g
}

annotateTRM = function(g, target) {
  g = .annotateTarget(g, target)
  .annotateTFclass(g)
}

.checkParam = function(p1, p2, default, multi) {
  if(!missing(p1))
    if(!missing(multi) & length(p1) == 1)
      return(rep(p1, multi))
  else
    return(p1)
  if(!is.null(p2))
    return(p2)
  if(!missing(multi))
    rep(default, multi)
  else
    default
}

plotTRM = function(g, layout = layout.fruchterman.reingold, mar = .5, vertex.col, vertex.cex, vertex.lwd, edge.col, edge.lwd, edge.lty, label = TRUE, label.cex, label.col, label.pos = NULL, label.offset = 1.5, adjust.label.col = FALSE, normalize.layout=TRUE) {
  
  if(class(layout) == "function")
    l = layout(g)
  else
    l = layout

  # normalize layout.
  if(normalize.layout)
    l = layout.norm(l, -1, 1, -1, 1)
  
  vertex.col = .checkParam(vertex.col, V(g)$frame.color, "black", multi = vcount(g))
  vertex.cex = .checkParam(vertex.cex, V(g)$size, 10, multi = vcount(g))
  vertex.lwd = .checkParam(vertex.lwd, V(g)$frame.width, 1, multi = vcount(g))
  
  edge.col = .checkParam(edge.col, E(g)$color, "grey")
  edge.lwd = .checkParam(edge.lwd, E(g)$width, 1, multi = ecount(g))
  edge.lty = .checkParam(edge.lty, E(g)$lty, "solid", multi = ecount(g))
  
  label.col = .checkParam(label.col, V(g)$label.color, "black")
  label.cex = .checkParam(label.cex, V(g)$label.cex, 1)
  if (adjust.label.col) {
    col.range = range(.lum2(colors()))
    col.cut = round(diff(col.range)/2)
    lum.mean = sapply(V(g)$pie.color, function(x) mean(.lum2(x)))
    label.col = ifelse(lum.mean < col.cut, "snow", "gray20")
  }
  
  mat = get.edgelist(g)
  
  op = par(mar = rep(mar,4), xpd = TRUE)
  plot(0, xlim = range(l[,1]), ylim = range(l[,2]), cex = 0, axes = FALSE, xlab = "", ylab = "", asp = 1)
  for(i in 1:nrow(mat)){
    ni = as.numeric(V(g)[ mat[i,] ])
    x1 = l[ni[1] ,]
    x2 = l[ni[2] ,]
    lines(c(x1[1], x2[1]), c(x1[2], x2[2]), col = edge.col, lwd = edge.lwd[i], lty = edge.lty[i])
  }
  for(i in 1:nrow(l)) {
    np=V(g)$pie[[i]]
    col=V(g)$pie.color[[i]]
    .floating.pie(l[i,1], l[i,2], x=np, col=col, radius=vertex.cex[i]/100, frame.width=vertex.lwd[i], frame.color=vertex.col[i])
  }
  if(label) {
    ll = as.character(V(g))
    if(!is.null(V(g)$name))
      ll = V(g)$name
    if(!is.null(V(g)$label))
      ll = V(g)$label
    text(l, labels = ll, cex = label.cex, col = label.col, pos = label.pos, offset = label.offset) 
  }
  par(op)
}

plotTRMlegend = function (x, title = NULL, cex = 1) 
{
	
	if(class(x) == "igraph")
  	family = sort(unique(unlist(V(x)$family)))
	else
		family = x
  
  col = .getTFclassColors()
  col = c(col, unclassified = "white")
  
  op = par(mar = rep(0, 4))
  plot(0, col = "transparent", axes = FALSE)
  legend("center", legend = family, fill = col[family], bty = "n", cex = cex, title = title)
  par(op)
}

.floating.pie = function (xpos, ypos, x, edges = 200, radius = 0.8, col = NULL, border = TRUE, lty = NULL, frame.color = "black", frame.width = 1)
{
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("'x' values must be positive.")
  
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  
  if (!is.null(col)) 
    col <- rep(col, length.out = nx)
    
  if (!is.null(lty)) 
    lty <- rep(lty, length.out = nx)
  
  
  init.angle = 90
  
  t2xy = function(t) {
    t2p <- -2 * pi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x+xpos, xpos), c(P$y+ypos, ypos), border = NA, col = col[i], lty = lty[i])
  }
  if (border) 
    symbols(xpos, ypos, circles = radius, inches = FALSE, add = TRUE, fg = frame.color, lwd = frame.width)
}

plotGraph = function(g, layout = layout.fruchterman.reingold, mar = .5, vertex.pch = 21, vertex.cex, vertex.col, vertex.bg, vertex.lwd, edge.col, edge.lwd, edge.lty, label = TRUE, label.col, label.cex, label.pos = NULL, label.offset = 1.5, adjust.label.col = FALSE, normalize.layout=TRUE) { 
  if(class(layout) == "function")
    l = layout(g)
  else
    l = layout
  
  # normalize layout.
  if(normalize.layout)
    l = layout.norm(l, -1, 1, -1, 1)
  
  vertex.cex = .checkParam(vertex.cex, V(g)$size, 5)
  vertex.col = .checkParam(vertex.col, V(g)$frame.color, "grey")
  vertex.bg = .checkParam(vertex.bg, V(g)$color, "white")
  vertex.lwd = .checkParam(vertex.lwd, V(g)$frame.width, 1)
  
  edge.col = .checkParam(edge.col, E(g)$color, "grey", multi = ecount(g))
  edge.lwd = .checkParam(edge.lwd, E(g)$width, 1, multi = ecount(g))
  edge.lty = .checkParam(edge.lty, E(g)$lty, "solid", multi = ecount(g))
  
  label.col = .checkParam(label.col, V(g)$label.color, "black")
  label.cex = .checkParam(label.cex, V(g)$label.cex, 1)
  if (adjust.label.col) {
    col.range = range(.lum2(colors()))
    col.cut = round(diff(col.range)/2)
    lum.mean = sapply(V(g)$pie.color, function(x) mean(.lum2(x)))
    label.col = ifelse(lum.mean < col.cut, "snow", "gray20")
  }
  
  mat = get.edgelist(g)
  op = par(mar = rep(mar,4), xpd = TRUE)
  plot(0, xlim = range(l[,1]), ylim = range(l[,2]), cex = 0, axes = FALSE, xlab = "", ylab = "", asp = 1)
  for(i in 1:nrow(mat)){
    ni = as.numeric(V(g)[ mat[i,] ])
    x1 = l[ni[1] ,]
    x2 = l[ni[2] ,]
    lines(c(x1[1], x2[1]), c(x1[2], x2[2]), col = edge.col[i], lwd = edge.lwd[i], lty = edge.lty[i])
  }
  points(l, pch = vertex.pch, cex = vertex.cex, bg = vertex.bg, col = vertex.col, lwd = vertex.lwd)
  if(label) {
    ll = as.character(V(g))
    if(!is.null(V(g)$name))
      ll = V(g)$name
    if(!is.null(V(g)$label))
      ll = V(g)$label
    text(l, labels = ll, cex = label.cex, pos = label.pos, offset=label.offset) 
  }
  par(op)
}


## annotation functions to compare modules.
# annotate a graph based on the information in other module, as well as expression and targets, etc.
annotateModule = function(g, enrich, trm, targets, ppi, exprs, tfs) {
  trm_genes = .checkInNetwork(g, V(trm)$name)
  enrich_genes = .checkInNetwork(g, enrich)
  targets_found = targets[targets %in% trm]
  targets_found = .checkInNetwork(g, targets_found)
  exprs = .checkInNetwork(g, exprs)
                             
	V(g)$color = "white"
	V(g)[ trm_genes ]$color = "white"
	V(g)[ enrich_genes ]$color = "steelblue2"
	
	V(g)[ targets_found ]$color = "steelblue4"
	V(g)$size = 1
	V(g)[ exprs ]$size = 15
	
	V(g)$frame.color = "gray"
	V(g)$frame.width = 1
	V(g)[ trm_genes ]$frame.color = "black"
	V(g)[ trm_genes ]$frame.width = 3
		
	if(!missing(tfs)) {
	  tfs = .checkInNetwork(g, tfs)
		V(g)$shape = "circle"
		V(g)[ tfs ]$shape = "square"
	}
	
	E(g)$lty = "dotted"
	el = get.edgelist(graph.intersection(ppi, g))
	for(j in 1:nrow(el)) {
		E(g, P = which(V(g)$name %in% el[j,]), directed = FALSE)$lty = "solid"
	}
	
	E(g)$color = "grey"
	E(g)$width = 1
	el = get.edgelist(graph.intersection(trm, g))
	for(j in 1:nrow(el)) {
		E(g, P = which(V(g)$name %in% el[j,]), directed = FALSE)$color = "black"
		#		E(g, P = which(V(g)$name %in% el[j,]), directed = FALSE)$width = 3
	}
	
	g
}


# annotate a graph based on the frequency the nodes and edges appear in 
# a set of other graphs, provided as graph_list.
annotateFreq = function(g, graph_list) {
	s = sapply(V(g)$name, function(x) {
		10*length(which(sapply(graph_list, function(s) if(!is.null(s)) any(V(s)$name %in% x) else FALSE)))/length(graph_list[!sapply(graph_list, is.null)])
	})
	
	ew = apply(get.edgelist(g), 1, function(e) {
		sum(sapply(names(graph_list), function(n) {

			s = graph_list[[n]]
			if(!is.null(s)) {
        e = .checkInNetwork(s, e)
				if(length(V(s)[ e ]) == 2)
					if(are.connected(s, e[1], e[2])) return(1)
			}
			return(0)
		}))
	})
	
	et = rep("solid", length(ew))
	et[ew == 0] = "dotted"
	ew[ew == 0] = 1
	
	V(g)$size = s
	E(g)$width = ew
	E(g)$lty = et
	g
}

.checkInNetwork = function(g, x) {
  x[ x %in% V(g)$name ]
}
