\name{S.PhyloMaker}
\alias{S.PhyloMaker}
\title{
Function of Qian and Jin (2016) to generate a phylogenetic tree}
\description{
The function builds a phylogenetic tree based on taxonomic information on the taxa to be included.
}
\usage{
S.PhyloMaker(splist, tree, nodes, output.splist = T, scenarios = "S3")
}
\arguments{
  \item{splist}{
a data frame containing taxonomic information for a set of plant species (rows), for which \code{S.PhyloMaker} generates a phylogeny. It should include three columns named "species", "genus", and "family", and providing the species names (e.g. Acacia berlandieri), genus names (e.g. Acacia), and family names (e.g. Fabaceae), respectively.
}
  \item{tree}{
a \code{phylo} object including the reference megaphylogeny of vascular plants (available in \code{\link{qian}}, i.e., \code{qian$phylo.all}). 
}
  \item{nodes}{
a data frame containing the basal node information of every family and genus of the megaphylogeny (available in \code{\link{qian}}, i.e., \code{qian$nodes}).
}
  \item{output.splist}{
if \code{TRUE}, the output species list will include in \code{splist} an extra column showing which species have matched with species in the megaphylogeny.
}
  \item{scenarios}{
defines which scenario(s) will be used to generate phylogenies. See Qian and Jin (2016) for more information on the alternative options \code{S1}, \code{S2} and \code{S3}. By default, scenario \code{S3} is performed.
}
}
\details{
(The function description is derived from the Appendix of Qian and Jin (2016))

Using the PhytoPhylo megaphylogeny as a backbone, this function takes a user-specified species list and matches species in the species list with species in the megaphylogeny of Qian and Jin (2016). If a species is found in the megaphylogeny, the species will be selected for pruning; if the species is not found but its genus or family is found in the megaphylogeny, the species will be added to the backbone at a certain place within the genus or family depending on which scenarios are taken. Finally, the selected and added species are pruned from the megaphylogeny, which results in three phylogenies corresponding to the three scenarios (when default is chosen). 
}
\value{
  \item{Scenario.x}{
a phylogenetic tree in \code{phylo} format, built according to the \code{scenarios} defined in input. \code{x} is replaced with 1, 2 or 3 depending on whether scenarios \code{S1}, \code{S2} or \code{S3} have been defined.
}
  \item{Species.list}{
a species list including original columns of the \code{splist} data frame and, if \code{output.splist = T}, an additional column called \code{status} showing which species have matched with species of the megaphylogeny. Species matching between the species list and the megaphylogeny are indicated as \code{match(prune)}, species that are added to the megaphylogeny before pruning are indicated as \code{match(add)}, while species whose genera or families are not found in the megaphylogeny are indicated as \code{unmatch}.
}
}
\references{
Qian H. & Jin Y. (2016). An updated megaphylogeny of plants, a tool for generating plant phylogenies, and an analysis of phylogenetic community structure. Journal of Plant Ecology. 9(2):233-239.
}
\author{
F. Munoz, adapted from Qian and Jin (2016)
}
\seealso{
\code{\link{create.phylo}}
}
\examples{
example <- data.frame(species = c("Dipterocarpus indicus", "Memecylon umbellatum",
"Michelia champaca"), genus = c("Dipterocarpus", "Memecylon", "Michelia"), 
family = c("Dipterocarpaceae", "Melastomataceae", "Magnoliaceae"))

# run the function S.PhyloMaker
result <- S.PhyloMaker(splist = example, tree = qian$phylo.all, nodes = qian$nodes)
str(result)       # the structure of the ouput of S.PhyloMaker.

plot(result$Scenario.3, cex = 1.1)
}
