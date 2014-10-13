.getNode = function(x, subset, tfclass, dbname = NULL) {
	if(missing(tfclass)) {
		sql = "SELECT * FROM tfclass"
		tfclass = .getSQLresult(sql, dbname = dbname)
	}
  
  if(! x %in% tfclass$id) return(NULL)
  
  new_ids = tfclass$is_a[tfclass$id == x]
  unique(sapply(new_ids, function(id) {
    if(tfclass$subset[tfclass$id == id] == subset)
      id
    else
      .getNode(id, subset, tfclass)
  }))
}

# TODO: rename getTFclassFromEntrezgene() and integrate with the one at database.R.
getTFclassFromEntrezgene = function(x, subset = "Class", tfclass, dbname = NULL) {
  if(missing(tfclass)) {
    sql = paste("SELECT * FROM tfclass")
    tfclass = .getSQLresult(sql, dbname = dbname) 
  }
  
  if(length(x) > 1)
    lapply(x, function(xx) unique(getTFclassFromEntrezgene(xx, subset, tfclass)))
  else {
    # check if there is a node
    if(! any(tfclass$id == x)) {
      o_o = getOrthologs()
      e_map = o_o$entrezgene[o_o$map_entrezgene %in% x]
      x = getOrthologs(e_map, "human")$map_entrezgene
    }
    
    # remove those not in d:
    x = x[x %in% tfclass$id]
    
    if(length(x)>0)
      c(unique(sapply(x, function(xx) {
        tmp = .getNode(xx, subset, tfclass)
        tfclass$name[ tfclass$id %in% tmp ]
      })))
    else NULL
  }
}

getTFterms = function (subset = "Class", dbname = NULL) {
	sql = paste("SELECT * FROM tfclass")
	d = .getSQLresult(sql, dbname = dbname)
	d$name[d$subset == subset]
}

### tests:
# getEntrezgeneByClass = function(subset = "Class") {
#   # load tfclass_g
#   if(! exists("tfclass", envir = .datacache))
#     load(system.file("db/tfclass.rda", package = "rTRM"), envir = .datacache)
#   g = get("tfclass_g", envir = .datacache)
#   
#   genes = V(g)[ subset == "Map"]$name
#   lapply(genes, function(gene) {
#     message(gene)
#     .getNode(gene, g, subset)
#   })
# }
