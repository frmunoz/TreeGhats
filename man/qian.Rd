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
  The format is:
List of 2
 $ phylo.all:List of 5
  ..$ edge       : int [1:62775, 1:2] 31390 31391 31391 31392 31392 31393 31393 31394 31394 31395 ...
  ..$ Nnode      : int 31387
  ..$ tip.label  : chr [1:31389] "Blasia_pusilla" "Lunularia_cruciata" "Marchantia_polymorpha" "Riccia_fluitans" ...
  ..$ edge.length: num [1:62775] 73.2 367.1 54 313.1 177.7 ...
  ..$ node.label : chr [1:31387] "" "" "" "" ...
  ..- attr(*, "class")= chr "phylo"
  ..- attr(*, "order")= chr "cladewise"
 $ nodes    :'data.frame':	8493 obs. of  7 variables:
  ..$ level     : chr [1:8493] "G" "G" "G" "G" ...
  ..$ genus     : chr [1:8493] "Hamatocaulis" "Jaegerina" "Actinostrobus" "Helleborus" ...
  ..$ family    : chr [1:8493] "" "" "Cupressaceae" "Ranunculaceae" ...
  ..$ node.label: chr [1:8493] "N100" "N100" "N1000" "N10008" ...
  ..$ node.age  : num [1:8493] 22.7 22.7 18.7 20.8 23.7 ...
  ..$ clade.size: int [1:8493] 1 1 2 13 39 1 2 4 2 32 ...
  ..$ taxa      : chr [1:8493] "" "" "" "" ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Qian H. & Jin Y. (2015). An updated megaphylogeny of plants, a tool for generating plant phylogenies, and an analysis of phylogenetic community structure. Journal of Plant Ecology.
}
\references{

}
\examples{
data(qian)
## maybe str(qian) ; plot(qian) ...
}
\keyword{datasets}
