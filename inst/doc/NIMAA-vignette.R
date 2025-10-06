## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(NIMAA)

## ----show data,echo=FALSE, results='asis'-------------------------------------
knitr::kable(NIMAA::beatAML[1:10,], caption='The first ten rows of beatAML dataset')

## ----read the data------------------------------------------------------------
# read the data
beatAML_data <- NIMAA::beatAML

## ----plotIncMatrix function, results='asis'-----------------------------------
beatAML_incidence_matrix <- plotIncMatrix(
  x = beatAML_data, # original data with 3 columns
  index_nominal = c(2,1), # the first two columns are nominal data
  index_numeric = 3,  # the third column is numeric data
  print_skim = FALSE, # if you want to check the skim output, set this as TRUE
  plot_weight = TRUE, # when plotting the weighted incidence matrix
  verbose = FALSE # NOT save the figures to local folder
  )

## ----echo=FALSE, out.width="100%", fig.cap="The beatAML dataset as an incidence matrix"----
knitr::include_graphics("patient_id-inhibitor.png")

## ----plotBipartite, fig.height = 7, fig.width = 7, fig.align = "center"-------
bipartGraph <- plotBipartite(inc_mat = beatAML_incidence_matrix, vertex.label.display = T)

# show the igraph object
bipartGraph

## ----plotBipartiteInteractive, results='hide'---------------------------------
plotBipartiteInteractive(inc_mat = beatAML_incidence_matrix)

## ----echo=FALSE, out.width="100%", fig.cap="beatAML dataset as incidence matrix"----
knitr::include_graphics("interactiveplot.png")

## -----------------------------------------------------------------------------
analysis_reuslt <- analyseNetwork(bipartGraph)

# showing the general measures for network topology
analysis_reuslt$general_stats

## ----extractSubMatrix,eval=TRUE, echo=TRUE------------------------------------
sub_matrices <- extractSubMatrix(
  x = beatAML_incidence_matrix,
  shape = c("Square", "Rectangular_element_max"), # the selected shapes of submatrices
  row.vars = "patient_id",
  col.vars = "inhibitor",
  plot_weight = TRUE,
  print_skim = FALSE
  )

## ----echo=FALSE, out.width="30%", fig.cap="Row-wise arrangement"--------------
knitr::include_graphics("Row_wise_arrangement.png")

## ----echo=FALSE, out.width="100%", fig.cap="Column-wise arrangement"----------
knitr::include_graphics("Column_wise_arrangement.png")

## ----eval=TRUE, echo= TRUE, results='hide',error=FALSE, warning=FALSE, message=FALSE, fig.height = 3, fig.width = 7, fig.align = "center"----
cls1 <- findCluster(
  sub_matrices$Rectangular_element_max,
  part = 1,
  method = "all", # all available clustering methods
  normalization = TRUE, # normalize the input matrix
  rm_weak_edges = TRUE, # remove the weak edges in graph
  rm_method = 'delete', # delete the weak edges instead of lowering their weights to 0.
  threshold = 'median', # Use median of edges' weights as threshold
  set_remaining_to_1 = TRUE, # set the weights of remaining edges to 1
  )

## ----generate extra features, eval=TRUE, echo=TRUE, results='asis'------------
external_feature <- data.frame(row.names = cls1$walktrap$names)
external_feature[,'membership'] <- paste('group',
                                         findCluster(sub_matrices$Rectangular_element_max,
                                                     method = c('walktrap'),
                                                     rm_weak_edges = T,
                                                     set_remaining_to_1 = F,
                                                     threshold = 'median',
                                                     comparison = F)$walktrap$membership)
dset1 <- utils::head(external_feature)
knitr::kable(dset1, format = "html")


## ----findcluster 1, eval=TRUE, fig.height = 6, fig.width = 7, fig.align = "center"----
cls <- findCluster(
  sub_matrices$Rectangular_element_max,
  part = 1,
  method = "all", 
  normalization = TRUE, 
  rm_weak_edges = TRUE, 
  rm_method = 'delete', 
  threshold = 'median', 
  set_remaining_to_1 = TRUE, 
  extra_feature = external_feature # add an extra feature reference here
  )

## ----eval=TRUE,fig.height = 3, fig.width = 7, fig.align = "center"------------
cls2 <- findCluster(
  sub_matrices$Rectangular_element_max, # the same submatrix
  part = 2 # set to 2 to use the other part of bipartite network
  )

## ----plotcluster,eval = FALSE-------------------------------------------------
# plotCluster(graph=cls2$graph,cluster = cls2$louvain)

## ----echo=FALSE, out.width="100%", fig.cap="Interactive plot for function plotCluster()"----
knitr::include_graphics("plotCluster.png")

## ----visualClusterInBipartite, eval = FALSE-----------------------------------
# visualClusterInBipartite(
#   data = beatAML_data,
#   community_left = cls1$leading_eigen,
#   community_right = cls2$fast_greedy,
#   name_left = 'patient_id',
#   name_right = 'inhibitor')

## ----echo=FALSE, out.width="100%", fig.cap="Interactive plot for function visualClusterInBipartite()"----
knitr::include_graphics("visualClusterInBipartite.png")

## ----eval=FALSE, echo=TRUE----------------------------------------------------
# scoreCluster(community = cls2$infomap,
#              graph = cls2$graph,
#              dist_mat = cls2$distance_matrix)

## ----eval=TRUE----------------------------------------------------------------
validateCluster(community = cls$leading_eigen,
                extra_feature = external_feature,
                dist_mat = cls$distance_matrix)

## -----------------------------------------------------------------------------
imputations <- predictEdge(inc_mat = beatAML_incidence_matrix,
                                  method = c('svd','median','als','CA'))
# show the result format
summary(imputations)

## ----eval=TRUE, fig.height = 6, fig.width = 7, fig.align = "center"-----------
validateEdgePrediction(imputation = imputations,
                   refer_community = cls1$fast_greedy,
                   clustering_args = cls1$clustering_args)

