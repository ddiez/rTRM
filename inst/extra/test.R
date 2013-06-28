.testTRM = function(extended = FALSE, strict = FALSE) {
  require(org.Mm.eg.db)
  # get motifs.
  motif_file = system.file("extra/sox2_motif_list.rda", package = "rTRM")
  load(motif_file)
  
  # get tfs.
  tfs_list = getOrthologFromMatrix(motif_list, organism = "mouse")
  tfs_list = unique(unlist(tfs_list, use.names = FALSE))
  
  # get ppi.
  message("PPI")
  data(biogrid_mm)
  ppi = biogrid_mm
  
  # remove outliers.
  f = c("Ubc", "Sumo1", "Sumo2", "Sumo3")
  f = unlist(AnnotationDbi::mget(f, revmap(org.Mm.egSYMBOL)))
  ppi = removeVertices(ppi, f)
  
  # get expressed genes.
  message("EXPRS")
  eg_esc_file = system.file("extra/ESC-expressed.txt", package = "rTRM")
  eg_esc = scan(eg_esc_file, what = "")
  eg_esc = c(eg_esc, "71950") # fix Nanog...

  # filter data by expression.
  tfs_list_esc = tfs_list[tfs_list %in% eg_esc]
  
  ppi_esc = induced.subgraph(ppi, eg_esc[ eg_esc %in% V(ppi)$name ])
  ppi_esc = getLargestComp(ppi_esc)
  
  # define target.
  target = unlist(AnnotationDbi::mget("Sox2", org.Mm.egSYMBOL2EG))
  
  # find TRM
  message("TRM")
  s = findTRM(ppi_esc, target, tfs_list_esc, method = "nsa", max.bridge = 1, extended = extended, strict = strict)
  
  list(trm = s, exprs = eg_esc, ppi = ppi_esc, motif = motif_list, target = target, query = tfs_list_esc)
}