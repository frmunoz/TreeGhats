\name{qian}
\alias{qian}
\docType{data}
\title{
Megaphylogeny of plants
}
\description{
The dataset includes the phylogenetic tree and the table of node names provided by Qian and Jin (2016).}
\usage{data("qian")}
\format{
  The \code{qian} object includes two elements, \code{$phylo.all} and \code{$nodes}:
  
  \code{phylo.all} is a phylogenetic tree in \link{phylo} format,
  
  \code{nodes} is a data frame including basic information on the nodes and the related branch lengths. This information is used by \link{S.PhyloMaker} to include new taxa in the megaphylogeny depending on their taxonomy and to define branch lengths in generated phylogenies.
}
\source{
Qian H. & Jin Y. (2016). An updated megaphylogeny of plants, a tool for generating plant phylogenies, and an analysis of phylogenetic community structure. Journal of Plant Ecology. 9(2):233-239.
}
\seealso{
\link{S.PhyloMaker}
}
\keyword{datasets}