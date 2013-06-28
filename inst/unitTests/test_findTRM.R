test_findTRM = function() {
  load(system.file(package = "rTRM", "extra/example.rda"))
  target = "N6"
  query = c("N7", "N12", "N28")
  s = findTRM(g, target = target, query = query, method = "nsa", max.bridge = 1)
  checkEquals(vcount(s), 4)
  checkTrue(identical(V(s)$name, c("N6", "N7", "N9", "N12")))
}