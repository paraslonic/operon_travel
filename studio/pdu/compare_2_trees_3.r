library(ggtree)
library("dendextend")
library("ape")


setwd("/data12/bio/runs-manolov/operon_travel/studio/sorbitol/")

region_tree <- read.tree("sorbitol_aligned.fasta.raxml.bestTree")
core_tree <- read.tree("coreogaligned.fasta.raxml.bestTree")

makeDendros <- function(firstTree, secondTree)
{
  firstTree <- midpoint(firstTree)
  secondTree <- midpoint(secondTree)
  firstTree$edge.length[which(firstTree$edge.length == 0)] <- 1e-20
  firstTree <- chronopl(firstTree, lambda = 0.01, tol = 1e-19, iter.max = 2000, eval.max = 2000)
  firstTree <- as.dendrogram(firstTree)
  secondTree$edge.length[which(secondTree$edge.length == 0)] <- 1e-20
  secondTree <- chronopl(secondTree, lambda = 0.01, tol = 1e-19, iter.max = 2000, eval.max = 2000)
  secondTree <- as.dendrogram(secondTree)
  dends <- dendlist(firstTree, secondTree)
  dends <- dends %>% untangle(method = "step2side")
  return(dends)
}

plot_tree_compare <- function(dends)
{
  dends %>% plot(main = paste("entanglement =", round(entanglement(x), 2)), common_subtrees_color_branches = TRUE)
}

dends <- makeDendros(core_tree, region_tree)

png("compare_sorbitol_full.png", type="cairo", width = 1200, height = 1200)
  plot_tree_compare(dends)
dev.off()

### H I S T O
taxas <- region_tree$tip.label
N <- 50
e_values <- rep(NA, N)
for(i in 1:N)
{
  res = try({
    the_taxas <- sample(taxas, 20)
    to_del_taxas <- setdiff(taxas, the_taxas)
    region_subtree <- drop.tip(region_tree, to_del_taxas)
    core_subtree <- drop.tip(core_tree, to_del_taxas)
    dends <- makeDendros(core_subtree, region_subtree)
    e_values[i] <- entanglement(dends)
    #plot_tree_compare(dends)
  })
  if(isTRUE(class(res)=="try-error")) { next } 
}
hist(e_values, breaks = 20)

