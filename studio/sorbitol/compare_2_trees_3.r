library(ggtree)
library("dendextend")
library("ape")


setwd("/data12/bio/runs-manolov/operon_travel/studio/sorbitol/")

region_tree <- read.tree("sorbitol_aligned.fasta.raxml.bestTree")
core_tree <- read.tree("coreogaligned.fasta.raxml.bestTree")

firstTree <- region_tree

firstTree <- midpoint(firstTree)
secondTree <- midpoint(secondTree)

firstTree$edge.length[which(firstTree$edge.length == 0)] <- 1e-20
firstTree <- chronopl(firstTree, lambda = 0.01, tol = 1e-19, iter.max = 2000, eval.max = 2000)
firstTree <- as.dendrogram(firstTree)

secondTree$edge.length[which(secondTree$edge.length == 0)] <- 1e-20
secondTree <- chronopl(secondTree, lambda = 0.01, tol = 1e-19, iter.max = 2000, eval.max = 2000)
secondTree <- as.dendrogram(secondTree)

dends <- dendlist(firstTree, secondTree)
#plot(dends)

x <- dends %>% untangle(method = "step2side")
x %>% plot(main = paste("entanglement =", round(entanglement(x), 2)), common_subtrees_color_branches = TRUE)
