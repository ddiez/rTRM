.getCoordinates = function(x, r) {
	l = length(x)
	d = 360/l
	c1 = seq(0, 360, d)
	c1 = c1[1:(length(c1)-1)]
	tmp = t(sapply(c1, function(cc) c(cos(cc*pi/180)*r, sin(cc*pi/180)*r)))
	rownames(tmp) = x
	tmp
}

.checkValid = function(x) {
	if(any(table(x) > 1)) FALSE else TRUE
}

# TODO: make layout independent the $name $label (i.e. based on indexes.)
layout.concentric = function (g, concentric = NULL, radius = NULL, order.by) 
{
	if(is.null(concentric))
    concentric = list(V(g)$name)

	all_c = unlist(concentric, use.names = FALSE)
	
  if (!.checkValid(all_c))
	  stop("Duplicated nodes in layers!")
	
	if (!.checkValid(radius)) 
	  stop("Duplicated radius in layers!")
	
	all_n = V(g)$name
	sel_other = all_n[ ! all_n %in% all_c ]
	
	if(length(sel_other) > 0)
	  concentric[[length(concentric)+1]] = sel_other
    
	if(is.null(radius)) {
		radius = seq(0, 1, 1/(length(concentric)))
		if(length(concentric[[1]]) == 1)
			radius = radius[-length(radius)]
		else
			radius = radius[-1]
	}
	
	if( ! missing(order.by) )
		order.values = lapply(order.by, function(b) get.vertex.attribute(g, b))
		
	res = matrix(NA, nrow = length(all_n), ncol = 2)
	for(k in 1:length(concentric)) {
		r = radius[k]
		l = concentric[[k]]
    
		i = which(V(g)$name %in% l) - 1
		i_o = i
		if (!missing(order.by)) {
			ob = lapply(order.values, function(v) v[i + 1])
			ord = do.call(order, ob)
			i_o = i_o[ord]
		}
		res[i_o+1, ] = .getCoordinates(i_o, r)

	}
	res
}

getConcentricList = function(g, t, e, max.size = 60, order.by = "label") {
	sel.all = V(g)$name
	
	# filter out not in graph.
	t = t[t %in% sel.all]
	e = e[e %in% sel.all]
	
	sel.e = V(g)[ e ]$name
	sel.t = V(g)[ t ]$name
	sel.t = sel.t[sel.t %in% sel.e ] # choose only target that are enriched.
	sel.e = sel.e[! sel.e %in% sel.t ]
	
	sel.b = sel.all[! sel.all %in% c(sel.t, sel.e) ]
	
	tmp = list(sel.t, sel.b, sel.e)

  if(!is.null(order.by)) {
    tmp = lapply(tmp, function(l) {
      l[order(get.vertex.attribute(g, order.by, V(g)[ l ]))]
    })
  }

	res = list()
	for(k in 1:length(tmp)) {
		r = tmp[[k]]
    if(length(r)>max.size) {
      s = ceiling(length(r)/max.size)
      #r1 = split(r,rep(1:s, s,length.out = length(r)))
      v = rep(1:s, each = ceiling(length(r)/s))
      v = v[1:length(r)]
      r1 = split(r,v)
      for(kk in r1) {
        res[[length(res)+1]] = kk
      } 
    }	else res[[length(res)+1]] = r
	}
	res
}

layout.arc = function (g, target, query)
{
  n = vcount(g)
  if(! all(target %in% V(g)$name)) {
    warning("some targets not in graph, removing them.")
    target=target[target %in% V(g)$name]
  }
  
  #
  target = target[target %in% query]
  
  V(g)$type = "bridge"
  all_name=V(g)$name # could be V(g)[ name %in% query ] but want to avoid "note's" in R CMD check.
  V(g)[all_name %in% query]$type = "query"
  V(g)[all_name %in% target]$type = "target"
  
  g_con = g
  n_left = character()
  
  if(! is.connected(g)) {
    g_con = getLargestComp(g)
    n_left = setdiff(V(g)$name,V(g_con)$name)
  }
  
  all_type=V(g_con)$type
  all_name=V(g_con)$name
  set = list(target = target, bridge=V(g_con)[all_type == "bridge"]$name, query1 = character(), query2 = character(), query3 = character(), left=n_left)
  for(q in V(g_con)[all_type == "query"]$name) {
    sp = get.all.shortest.paths(g_con,from=V(g_con)[q],to=V(g_con)[all_name %in% target])
    #print(sp)
    sp_min = min(sapply(sp$res,length))
    #print(sp_min)
    if(sp_min == 2) {
      set$query1 = c(set$query1, V(g_con)[q]$name)
    }
    else set$query2 = c(set$query2, V(g_con)[q]$name)
  }
  
  set = lapply(set,function(s) {
    ns = V(g)[s]$label
    s[order(ns,decreasing=TRUE)]
  })
  
  x0 = c(left=-2,query1=-1,target=0,bridge=1,query2=2)
  y0 = sapply(set,function(x){
    -1 * floor(length(x)/2)
  })
  
  res = matrix(NA, nrow = n, ncol = 2)
  all_n = unlist(set)
  for(my_n in all_n) {
    k = which(V(g)$name == my_n)
    my_type = names(set)[sapply(set,function(x) my_n %in% x)]
    x1 = x0[my_type]
    y1 = y0[my_type]
    y0[my_type] = y0[my_type] + 1
    res[k,1] = x1
    res[k,2] = y1
  }
  res
}