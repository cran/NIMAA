## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(NIMAA)

## ----show data,echo=FALSE, results='asis'-------------------------------------
knitr::kable(NIMAA::beatAML[1:10,], caption='beatAML dataset samples')

## ----read the data------------------------------------------------------------
# read the data
beatAML_data <- NIMAA::beatAML

## ----plotInput function,results='asis'----------------------------------------
beatAML_incidence_matrix <- plotInput(
  x = beatAML_data, # original data with 3 columns
  index_nominal = c(2,1), # the first two columns are nominal data
  index_numeric = 3,  # the third column inumeric data
  print_skim = FALSE, # if you want to check the skim output, set this as TRUE(Default)
  plot_weight = TRUE, # when plotting the figure, show the weights
  verbose = FALSE # NOT save the figures to local folder
  )

## ---- echo=FALSE, out.width="100%", fig.cap="beatAML dataset as incidence matrix"----
knitr::include_graphics("patient_id-inhibitor.png")

## ----plotBipartite, fig.height = 7, fig.width = 7, fig.align = "center"-------
graph <- plotBipartite(inc_mat = beatAML_incidence_matrix,vertex.label.display = T)

# show the igraph graph object
graph

## ----plotBipartiteInteractive, results='hide'---------------------------------
plotBipartiteInteractive(inc_mat = beatAML_incidence_matrix)

## -----------------------------------------------------------------------------
analysis_reuslt <- analyseNetwork(graph)

## ----extractSubMatrix,eval=FALSE----------------------------------------------
#  sub_matrices <- extractSubMatrix(
#    x = beatAML_incidence_matrix,
#    shape = c("Square", "Rectangular_element_max"), # the shapes wanted
#    row.vars = "patient_id",
#    col.vars = "inhibitor",
#    plot_weight = TRUE,
#    verbose = FALSE,
#    print_skim = TRUE # just to reduce the length of vignette
#    )

## ---- echo=FALSE, out.width="30%", fig.cap="Row-wise arrangement"-------------
knitr::include_graphics("Row_wise_arrangement.png")

## ---- echo=FALSE, out.width="100%", fig.cap="Column-wise arrangement"---------
knitr::include_graphics("Column_wise_arrangement.png")

## ---- eval=FALSE, fig.height = 3, fig.width = 7, fig.align = "center"---------
#  cls <- findCluster(
#    sub_matrices$Rectangular_element_max,
#    dim = 1,
#    method = "all", # clustering mehod
#    normalization = TRUE, # normalize the input matrix
#    rm_weak_edges = TRUE, # remove the weak edges in graph
#    rm_method = 'delete', # removing method is deleting the edges
#    threshold = 'median', # edges with weights under the median of all edges' weight are weak edges
#    set_remaining_to_1 = TRUE, # set the weights of remaining edges to 1
#    )

## ----generate extra features,eval=FALSE,echo=FALSE, results='asis'------------
#  external_feature <- data.frame(row.names = cls$walktrap$names)
#  external_feature[,'membership'] <- paste('group', findCluster(sub_matrices$Rectangular_element_max,method = c('walktrap'),rm_weak_edges = T,set_remaining_to_1 = F,threshold = 'median',comparison = F)$walktrap$membership)
#  dset1 <- utils::head(external_feature)
#  knitr::kable(dset1, format = "html")

## ----findcluster 1, eval=FALSE,fig.height = 6, fig.width = 7, fig.align = "center"----
#  cls <- findCluster(
#    sub_matrices$Rectangular_element_max,
#    dim = 1,
#    method = "all", # clustering mehod
#    normalization = TRUE, # normalize the input matrix
#    rm_weak_edges = TRUE, # remove the weak edges in graph
#    rm_method = 'delete', # removing method is deleting the edges
#    threshold = 'median', # edges with weights under the median of all edges' weight are weak edges
#    set_remaining_to_1 = TRUE, # set the weights of remaining edges to 1
#    extra_feature =external_feature # ADD A EXTRA FEATURE REFRRENCE HERE!
#    )

## ----eval=FALSE,fig.height = 3, fig.width = 7, fig.align = "center"-----------
#  cls2 <- findCluster(
#    sub_matrices$Rectangular_element_max, # the same sub-matrix
#    dim = 2 # set to 2 to use the other part of graph
#    )

## ----plotcluster,eval = FALSE-------------------------------------------------
#  plotCluster(graph=cls2$graph,cluster = cls2$louvain)

## ---- eval=FALSE, echo=FALSE, out.width="100%", fig.cap="Interactive plot for function plotCluster()"----
#  knitr::include_graphics("plotCluster.png")

## ----visualClusterInBipartite,eval = FALSE------------------------------------
#  visualClusterInBipartite(
#    data = beatAML_data,
#    community_left = cls$leading_eigen,
#    community_right = cls2$fast_greedy,
#    name_left = 'patient_id',
#    name_right = 'inhibitor')

## ---- echo=FALSE, out.width="100%", fig.cap="Interactive plot for function visualClusterInBipartite()"----
knitr::include_graphics("visualClusterInBipartite.png")

## ----eval=FALSE---------------------------------------------------------------
#  scoreCluster(community = cls2$infomap,
#               graph = cls2$graph,
#               distance_matrix = cls2$distance_matrix)

## ----eval=FALSE---------------------------------------------------------------
#  validateCluster(dist_mat = cls$distance_matrix,
#                  extra_feature = external_feature,
#                  community = cls$leading_eigen)

## -----------------------------------------------------------------------------
imputations <- imputeMissingValue(
  inc_mat = beatAML_incidence_matrix,
  method = c('svd','median','als','CA')
  )
# show the result format
names(imputations)

## ---- eval=FALSE, fig.height = 6, fig.width = 7, fig.align = "center"---------
#  validation_of_imputation <- validateImputation(
#    imputation = imputations,
#    refer_community = cls$fast_greedy,
#    clustering_args = cls$clustering_args
#    )

## ----eval=FALSE---------------------------------------------------------------
#  NA %in% imputations$als

## ----eval=FALSE---------------------------------------------------------------
#  dim(imputations$als) == dim(beatAML_incidence_matrix)

