# This function creates species clusters based on dissimilarity matrix based on species traits
# Modified from Prima et al. 2024

## This function creates clusters of species based on a distance matrix between those species. 
## Several metrics are computed, the metrics and a figure is produced to select the cluster

speciesclusters <- function(species_matrixdist, no_clust_max) {
library(foreach)
avail.methods = c("complete", "ward.D", "ward.D2", "single",
                  "average", "mcquitty", "median", "centroid")
trait_matlist <- list(species_matrixdist)
clust.choice = foreach(clust.method = avail.methods) %do% {
  ## CALCULATE dendrograms from distance matrices
  clust.dendrograms = lapply(trait_matlist, function(x) {
    hclust(as.dist(x), method = clust.method)
  })
  ## CALCULATE THE DISTANCES corresponding to these dendrograms
  clust.DIST = lapply(clust.dendrograms, cophenetic)
  
  ## CALCULATE Mouchet measure
  clust.choice = sapply(1:length(clust.DIST), function(x){
    return(1 - (cor(as.dist(clust.DIST[[x]]), as.dist(trait_matlist[[x]])) *
                  cor(as.dist(clust.DIST[[x]]), as.dist(trait_matlist[[x]]))))
  })
  
  return(data.frame(clust.method = clust.method
                    , metric = clust.choice
                    , stringsAsFactors = FALSE))
}
clust.choice = do.call(rbind, clust.choice)
## CHOICE OF CLUSTERING METHOD ----------------------------------------------
clust.method = clust.choice$clust.method[which.min(clust.choice$metric)] 

## CALCULATE dendrograms from distance matrices
clust.dendrograms = lapply(trait_matlist, function(x) {
  hclust(as.dist(x), method = clust.method)
})

## COMPUTATION OF SEVERAL INDICES TO EVALUATE THE 'QUALITY' OF CLUSTERING
clust.evaluation = foreach(no.clusters = 2:no_clust_max) %do%
  {
    k1 = no.clusters
    k2 = no.clusters + 1
    c1 = cutree(clust.dendrograms[[1]], k = k1)
    c2 = cutree(clust.dendrograms[[1]], k = k2)
    stats = fpc::cluster.stats(trait_matlist[[1]], c1, c2)
    ## Dunn index : ratio of the smallest distance between observations
    ## not in the same cluster to the largest intra-cluster distance.
    ## Value between zero and infinity, and should be maximized.
    mdunn = clValid::dunn(trait_matlist[[1]], c1)
    ## Meila's VI index (Variation of Information) : measures the amount of 
    ## information lost and gained in changing between 2 clusterings.
    ## Should be minimized (?)
    mVI = stats$vi
    ## Value between zero and one. Should be maximized.
    R2 = stats$average.between / (stats$average.between + stats$average.within)
    ## Calinski and Harabasz index : 
    ## The higher the value, the "better" is the solution.
    ch = stats$ch
    ## Corrected rand index : measure of the similarity between two data clusterings.
    ## Value between 0 and 1, with 0 indicating that the two data clusters do not agree
    ## on any pair of points and 1 indicating that the data clusters are exactly the same.
    Rand = stats$corrected.rand
    ## Average silhouette width :
    ## Observations with a large s(i) (almost 1) are very well clustered,
    ## a small s(i) (around 0) means that the observation lies between two clusters,
    ## and observations with a negative s(i) are probably placed in the wrong cluster.
    ## Should be maximized.
    av.sil = stats$avg.silwidth
    return(data.frame(no.clusters, 
                      mdunn, mVI, R2, ch, Rand, av.sil, 
                      stringsAsFactors = FALSE))
  }
clust.evaluation = do.call(rbind, clust.evaluation)
clust.evaluation = reshape2::melt(clust.evaluation, id.vars = c("no.clusters"))
clust.evaluation.optim = split(clust.evaluation
                               , list(clust.evaluation$variable))
clust.evaluation.optim = foreach(ii = 1:length(clust.evaluation.optim), .combine = "rbind") %do%
  {
    tab = clust.evaluation.optim[[ii]]
    ord = ifelse(length(grep("mVI", names(clust.evaluation.optim)[ii])) > 0, FALSE, TRUE)
    tab$ORDER = NA
    tab$ORDER[order(tab$value, decreasing = ord)] = nrow(tab):1
    return(tab)
  }
### Graphical visualization of the results
pp2 <- ggplot(clust.evaluation.optim) +
  facet_grid("variable ~ .", scales = "free") +
  geom_vline(
    aes(xintercept = no.clusters, color = ORDER, alpha = ORDER), lwd = 4) +
  scale_color_viridis_c(guide = "none") +
  scale_alpha(guide = "none", range = c(0.1, 0.8)) +
  geom_point(aes(x = no.clusters, y = value)) +
  geom_line(aes(x = no.clusters, y = value)) +
  labs(x = "", y = "",
       title = "Choice of cluster number",
       subtitle = paste0(
         "Evolution of clustering evaluation variables with ",
         "the number of clusters.\n",
         "All values except that of mVI must be maximized ",
         "(check script for more details about the measures).\n",
         "Values are highlighted to help finding the number of clusters to keep : ",
         "the brighter (yellow-ish) the better."
       )
  )
plot(pp2)
 
############################################
return(list(
  clust.dendrograms = clust.dendrograms,
  clust.evaluation = clust.evaluation,
  plot.clustNo = pp2
))
}
