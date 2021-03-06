# write tab delimited report with information about a TRM nodes.
writeTRMreport = function(graph, file, organism, target, query, sort.by = "symbol") {
	x = V(graph)$name
  
  map=.getMapFromOrg(organism)
  res=select(map,keys=x,columns=c("SYMBOL","GENENAME"))
  
  S=res$SYMBOL
	S[is.na(S)] = ""
  
  D=res$GENENAME
  D[is.na(D)] = ""
	
	grole = rep("bridge", length(x))
	grole[x %in% query] = "enriched"
	
	gtype = rep("query", length(x))
	gtype[x %in% target] = "target"
	
	family = sapply(getTFclassFromEntrezgene(x), function(z) if(length(z) > 0) paste(z, sep = " | ") else "")
	
	d = data.frame("entrezgene" = x, symbol = S, role = grole, type = gtype, description = D, family = family, check.names = FALSE)
	
	d = d[order(d[, sort.by], decreasing = TRUE),]
	
	if(!missing(file))
		write.table(d, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
	invisible(d)
}
# compare a list of graph all-against-all.
getSimilarityMatrix = function(g_list, type = "edges") {
	type = match.arg(type, c("nodes", "edges"))
	m = matrix(NA, nrow = length(g_list), ncol = length(g_list), dimnames = list(names(g_list), names(g_list)))
	for(i in 1:length(g_list)) {
		for(j in i:length(g_list)) {
			#message(i, " vs ", j)
			i_n = names(g_list)[i]
			j_n = names(g_list)[j]
			#message(i_n, " vs ", j_n)
			i_g = g_list[[i_n]]
			j_g = g_list[[j_n]]
			if(!is.null(i_g) & !is.null(j_g)) {
				gi = graph.intersection(i_g, j_g)
				gu = graph.union(i_g, j_g)
				switch(type, 
							 nodes = {
							 	t_n = vcount(gu)
							 	c_n = vcount(gi)
							 	p = 100*c_n/t_n
							 	#if(length(p) == 0) p = 0
							 	m[i_n, j_n] = p
							 	m[j_n, i_n] = p
							 },
							 edges = {
							 	t_e = ecount(gu)
							 	c_e = ecount(gi)
							 	p = 100*c_e/t_e
#							 	if(length(p) == 0) p = 0
							 	m[i_n, j_n] = p
							 	m[j_n, i_n] = p
							 }
				)
			}
		}
	}
	m
}

.dcor = function (x, use = "pairwise") 
{
	as.dist(1 - cor(t(x), use = use))
}

getSequencesFromGenome = function(BED, genome, append.id) {
  requireNamespace("Biostrings")
  s = Biostrings::getSeq(genome, names = BED$chr, start = BED$start, end = BED$end)
  sn = paste(ifelse(!missing(append.id), append.id, ""), BED$chr, BED$start, BED$end, sep = "_")
  names(s) = sn
  s
}