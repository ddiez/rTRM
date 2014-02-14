### basica functions.

# getMatrices
# returns a list with the PWM matrices.
# filter: name of PWM to obtain.
# matrices are returned as 'matrix' object (before was 'pwm' object (MotIV package) but now I am 
# moving to standard bioconductor classes and methods.
getMatrices = function (filter, dbname = NULL) 
{
  if (missing(filter)) 
    sql = "SELECT * FROM matrix"
  else {
    filter = paste("('", paste(filter, collapse = "','"), 
                   "')", sep = "")
    sql = paste("SELECT * FROM matrix WHERE pwm_id IN", filter)
  }
  tmp = .getSQLresult(sql, dbname = dbname)
  ids = unique(tmp$pwm_id)
  res = lapply(ids, function(id) {
    sel = tmp$pwm_id == id
    tmp2 = t(tmp[sel, c("A", "C", "G", "T")])
    colnames(tmp2) = 1:ncol(tmp2)
    tmp2
  })
  names(res) = ids
  res
}


getAnnotations = function(filter, dbname = NULL) {	
	if(missing(filter))
		sql = "SELECT * FROM pwm"
	else {
		filter = paste("('", paste(filter, collapse = "','"), "')", sep = "")
		sql = paste("SELECT * FROM pwm WHERE pwm_id IN", filter)
	}

	.getSQLresult(sql, dbname = dbname)
}

getMaps = function(filter, dbname = NULL) {	
	if(missing(filter))
		sql = "SELECT * FROM map"
	else {
		filter = paste("('", paste(filter, collapse = "','"), "')", sep = "")
		sql = paste("SELECT * FROM map WHERE pwm_id IN", filter)
	}

	.getSQLresult(sql, dbname = dbname)
}

getOrthologs = function (filter, organism, dbname = NULL) 
{
  sql = paste("SELECT * FROM ortholog", sep = "")
  if (!missing(organism))
    sql = paste(sql, " WHERE map_organism='", organism, "'", sep = "")
  if (!missing(filter)) {
    filter = paste("('", paste(filter, collapse = "','"), "')", sep = "")
    sql = paste(sql, "AND entrezgene IN", filter)
  }
  .getSQLresult(sql, dbname = dbname)
}

getTFclass = function(dbname = NULL) {
  sql = "SELECT * FROM tfclass"
  .getSQLresult(sql, dbname = dbname)
}

.getDefaultDb = function() {
	system.file("db/database.db", package = "rTRM")
}

.getSQLresult = function(sql, drv = "SQLite", dbname = NULL) {
	if(is.null(dbname))
		dbname = .getDefaultDb()
	
	drv <- dbDriver(drv)
	con = dbConnect(drv, dbname = dbname)
	res = dbSendQuery(con, sql)
	r = fetch(res, n = -1)
	dbClearResult(res)
	dbDisconnect(con)
	r
}


######## utility functions.

# getOrthologFromMatrix
# returns a vector with with corresponding entrezgene ids for the
# target organism from a list of motif ids.
getOrthologFromMatrix = function(filter, organism = "human", dbname = NULL) {
  m = getMaps(filter, dbname = dbname)
  o = unique(getOrthologs(m$entrezgene, organism = organism)$map_entrezgene)
  o[o != ""]
}

# get corresponding motifs for entrezgene ids.
getMotifsFromEntrezgene = function(e, organism) {
  o = getOrthologs(organism = organism)
  e_map = unique(o$entrezgene[o$map_entrezgene == e])
  m = getMaps()
  m$pwm_id[m$entrezgene %in% e_map]
}

# get corresponding motifs from entrezgene symbol
getMotifsFromSymbol = function(s, organism) {
  map=.getMapFromOrg(organism)
  res=select(map,keys=s,columns="ENTREZID",keytype="ALIAS")
  e=unique(na.omit(res$ENTREZID))
  
  o = getOrthologs(organism = organism)
  e_map = unique(o$entrezgene[o$map_entrezgene == e])
  m = getMaps()
  m$pwm_id[m$entrezgene %in% e_map]
}