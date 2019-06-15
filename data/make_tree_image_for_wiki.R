set.seed(1)
temp <- rtree(4)
temp$tip.label <- c("sample_1", "sample_2", "sample_3", "sample_4")
plot(temp, use.edge.length = FALSE, edge.width = 5, cex = 2, no.margin = TRUE)

jpeg("./data/tree_for_wiki.jpg", width = 400, height = 400)
plot(temp, use.edge.length = FALSE, edge.width = 5, cex = 2, no.margin = TRUE)
dev.off()
