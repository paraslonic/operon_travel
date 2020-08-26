pdu_tree <- read.tree("/data/ps/operon_travel/tree_region/pdu/pdu.treefile")
plot(pdu_tree)

core_tree <- read.tree("//data/ps/operon_travel/orthosnake/pdu/Results/coreogs_nucleotide.treefile")
plot(core_tree)


library(ggtree)
install.packages("phytools")
install.packages("phyext")
library("dendextend")
library("ape")

# create random trees
# firstTree <- (pdu_tree)
# secondTree <- (core_tree)
firstTree <- midpoint(pdu_tree)
length(firstTree$tip.label)
# firstTree <- di2multi(firstTree, tol = 1e-2)
# plot(firstTree)
# firstTree <- multi2di(firstTree,random = TRUE)
length(firstTree$tip.label)

which(firstTree$edge.length == 0)
firstTree$edge[49, ]
MRCA(firstTree, c(209,73))

firstTree$edge.length[which(firstTree$edge.length == 0)] <- 1e-20
firstTree <- chronopl(firstTree, lambda = 0.01, tol = 1e-19, iter.max = 2000, eval.max = 2000)
firstTree <- as.dendrogram(firstTree)

secondTree <- midpoint(core_tree)
secondTree$edge.length[which(secondTree$edge.length == 0)] <- 1e-20
secondTree <- chronopl(secondTree, lambda = 0.01, tol = 1e-19, iter.max = 2000, eval.max = 2000)
#secondTree <- untangle_step_rotate_1side(as.dendrogram(secondTree), firstTree)[[1]]
#secondTree <- untangle_step_rotate_2side(as.dendrogram(secondTree), firstTree)[[1]]
secondTree <- untangle(as.dendrogram(secondTree), firstTree, method = "random", R = 100)[[1]]

dends <- dendlist(firstTree, secondTree)
#plot(dends)

x <- dends %>% untangle(method = "step2side")
x %>% plot(main = paste("entanglement =", round(entanglement(x), 2)), common_subtrees_color_branches = TRUE)
