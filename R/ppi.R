.getTaxidFromOrg = function(org) {
	switch(org,
		human = "9606",
		mouse = "10090"
	)
}

getBiogridData = function(release) {
	tmp = tempfile()
	file = paste("BIOGRID-ALL-", release, ".tab2", sep = "")
	url = paste("http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-", release, "/", file, ".zip", sep = "")
	download.file(url, destfile = tmp)
	db = read.delim(unz(tmp, paste(file, ".txt", sep = "")), check.names = FALSE, colClasses = "character")
	unlink(tmp)
	list(db = db, release = release, date = Sys.Date())
}

processBiogrid = function(dblist, org = "human", simplify = TRUE, type = "physical", mimic.old = FALSE) {
  
  db = dblist$db
  
  # get taxid
  txid = .getTaxidFromOrg(org)
  # restrict to interactions *only* in the target organism.
  # TODO: add more flexibility to this? depends on the problem to address...
  if(mimic.old)
    db = db[(db$"Organism Interactor A" == txid) | (db$"Organism Interactor B" == txid),]
  else
    db = db[(db$"Organism Interactor A" == txid) & (db$"Organism Interactor B" == txid),]
  
  # filter by interaction type [default: physical]
  # TODO: add flexibility to this? depends on the problem to address...
  db = db[db[, "Experimental System Type"] %in% type, ]
  
  dbtmp = as.matrix(db[, c("Entrez Gene Interactor A", "Entrez Gene Interactor B")])
  #dbtmp

  biogrid = graph.edgelist(dbtmp, directed = FALSE)
   
  # add gene annotations
  map = .getMapFromOrg(org)
  sym = unlist(AnnotationDbi::mget(V(biogrid)$name, map, ifnotfound = NA))
  sym[is.na(sym)] = paste("eg:", names(sym[is.na(sym)]), sep = "")
  
  V(biogrid)$label = sym[V(biogrid)$name]
  
  # add other annotations.
  biogrid = set.edge.attribute(biogrid, "exp_sys", value = db[, "Experimental System"])
  biogrid = set.edge.attribute(biogrid, "exp_sys_type", value = db[, "Experimental System Type"])
  biogrid = set.edge.attribute(biogrid, "pubmed_id", value = db[, "Pubmed ID"])
  
  # create simplified graph.
  if (simplify) {
    E(biogrid)$biogrid_count <- count.multiple(biogrid)
    biogrid = simplify(biogrid)
  }
  
  # print stats.
  print(paste("# nodes: ", vcount(biogrid), sep = ""))
  print(paste("# links: ", ecount(biogrid), sep = ""))
  
  # add graph level annotations.
  biogrid$info$release = dblist$release
  biogrid$info$date = dblist$date
  biogrid$info$organism = org
  
  biogrid
}