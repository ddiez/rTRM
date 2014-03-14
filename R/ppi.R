.getTaxidFromOrg = function(org) {
	switch(org,
		human = "9606",
		mouse = "10090"
	)
}

# NOTE: new implementation to allow automatic downloading of latest release.
getBiogridData = function(release) {
  tmp = tempfile()
  
  if(missing(release)) {
    release = "LATEST"
    dir = "Latest%20Release"
  } else {
    dir = paste("Release%20Archive/BIOGRID-", release, sep ="")
  }
  
  file = paste("BIOGRID-ALL-", release, ".tab2", sep = "")
  url = paste("http://thebiogrid.org/downloads/archives/", dir, "/", file, ".zip", sep = "")
  
  download.file(url, destfile = tmp)
  if(release=="LATEST") {
    f=unzip(tmp,list=TRUE)
    release=sub(".tab2.txt", "", sub("BIOGRID-ALL-", "", f$Name[1]))
    file = paste("BIOGRID-ALL-", release, ".tab2", sep = "")
  }
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
  res = select(map, keys=V(biogrid)$name, columns="SYMBOL")
  sym = res$SYMBOL
  names(sym)=res$ENTREZID
  sym[is.na(sym)] = paste("eg:", res$ENTREZID[is.na(sym)], sep="")
  
  V(biogrid)$label = sym[V(biogrid)$name]
  
  # add other annotations.
  biogrid = set.edge.attribute(biogrid, "exp_sys", value = db[, "Experimental System"])
  biogrid = set.edge.attribute(biogrid, "exp_sys_type", value = db[, "Experimental System Type"])
  biogrid = set.edge.attribute(biogrid, "pubmed_id", value = db[, "Pubmed ID"])
  
  # create simplified graph.
  if (simplify) {
    biogrid = igraph::simplify(biogrid,edge.attr.comb="concat")
    E(biogrid)$biogrid_count=sapply(E(biogrid)$pubmed_id,function(x) length(x))
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
