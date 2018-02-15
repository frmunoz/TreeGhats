\name{qian}
\alias{qian}
\docType{data}
\title{
Megaphylogeny of plants
}
\description{
The dataset includes the phylogenetic tree and the table of node names provided by Qian and Jin (2015)}
\usage{data("qian")}
\format{
  The \code{qian} object includes two elements, \code{$phylo.all} and \code{$nodes}:
  
  \code{phylo.all} is a phylogenetic tree in \code{\link{phylo}} format,
  
  \code{nodes} is a data.frame including basic information on the nodes and the related branch lengths. This information is used by \code{TreeGhats:S.PhyloMaker} to include new taxa in the megaphylogeny depending on their taxonomy and to define branch lengths in the generated phylogeny.
}
\source{
Qian H. & Jin Y. (2015). An updated megaphylogeny of plants, a tool for generating plant phylogenies, and an analysis of phylogenetic community structure. Journal of Plant Ecology.
}
\seealso{
\code{\link{S.PhyloMaker}}
}
\keyword{datasets}
