biogrid_mm = function() {
	data(biogrid_mm, envir = .datacache)
	message("Mouse PPI network data of class igraph")
	message("Number of nodes: ", vcount(get("biogrid_mm", envir = .datacache)))
	message("Number of edges: ", ecount(get("biogrid_mm", envir = .datacache)))
	message("Source: The BioGRID (http://www.thebiogrid.org)")
	message("Release: ", get("biogrid_mm", envir = .datacache)$info$release)
	message("Downloaded: ", get("biogrid_mm", envir = .datacache)$info$date)
	message("Use data(biogrid_mm) to load it")
	rm(biogrid_mm, envir = .datacache)
}

biogrid_hs = function() {
	data(biogrid_hs, envir = .datacache)
	message("Human PPI network data of class igraph")
	message("Number of nodes: ", vcount(get("biogrid_hs", envir = .datacache)))
	message("Number of edges: ", ecount(get("biogrid_hs", envir = .datacache)))
	message("Source: The BioGRID (http://www.thebiogrid.org)")
	message("Release: ", get("biogrid_hs", envir = .datacache)$info$release)
	message("Downloaded: ", get("biogrid_hs", envir = .datacache)$info$date)
	message("Use data(biogrid_hs) to load it")
	rm(biogrid_hs, envir = .datacache)
}
