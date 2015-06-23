# 
initBiomart = function(filter, biomart = "ensembl", host) {
	requireNamespace("biomaRt")
	if(missing(host)) host = "www.biomart.org"
	
	ds = list(
		"human" = "hsapiens_gene_ensembl",
		"mouse" = "mmusculus_gene_ensembl"
	)
	
	if(! missing(filter))
		ds = ds[filter]
		
	lapply(ds, function(d) useMart(biomart = biomart, dataset = d, host = host))
}


.getMapFromOrg = function(org) {#, map = "SYMBOL") {
	switch(org,
    human = {
      if(requireNamespace("org.Hs.eg.db")) get("org.Hs.eg.db")
    },
    mouse = {
      if(requireNamespace("org.Mm.eg.db")) get("org.Mm.eg.db")
    }
	)
}


getOrthologsFromBiomart = function(eg, target_org, mart) {
	org = list(
		"human" = "hsapiens_homolog",
		"mouse" = "mmusculus_homolog"
	)
	query = "ensembl_gene_id"
	query = c(query, sapply(org[[target_org]], function(x) paste(x, c("ensembl_gene"), sep = "_")))
	res = getBM(query, filters = "entrezgene", values = eg, mart = mart)
	res = res[,2]
	res = res[!is.na(res)]
	res = unique(res[res != ""])
	if(length(res) > 0) {
    map=.getMapFromOrg(target_org)
    res=select(map,keys=res,columns="ENTREZID",keytype="ENSEMBL")
		unique(na.omit(res$ENTREZID))
	}
}
