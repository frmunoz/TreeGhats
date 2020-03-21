\name{create.phylo}
\alias{create.phylo}
\alias{create.phylo.all}
\title{
Generating a phylogeny of selected Western Ghats trees}
\description{
Creates a phylogeny of tree taxa based on species, genus and APGIII family information.
}
\usage{
create.phylo(names = NULL, scenarios = "S3")
create.phylo.all(scenarios = "S3")
}
\arguments{
  \item{names}{a data frame including taxonomic information for each taxon on row. The columns must be named "Binome", "Genus" and "Family_APGIII", and include latine genus species combination, genus name and APGIII family name, respectively. Genus and species names in "Binome" must be separated with a space.
  Default value is \code{NULL}, the phylogeny is generated for all taxa in the TreeGhatsData database.
}
  \item{scenarios}{one of the scenarios of Jin and Qian (2019) proposed for including taxa that are not already present in the megaphylogeny. See \code{\link[V.PhyloMaker]{phylo.maker}} for more information on the available options.
}
}
\details{
\code{create.phylo} establishes a baseline phylogeny with the taxa in \code{names} that are already present in a reference megaphylogeny, and further includes the missing taxa by branching them to the corresponding genus node or, if the genus is absent from the megaphylogeny, to the corresponding family node.

\code{create.phylo.all} is used to generate the phylogeny of all the taxa of the database.

}
\value{
The function returns an object of class \code{link[ape]{phylo}}, which can be analyzed and plotted with tools available in \code{ape} package.
}
\references{
Jin, Y., Qian, H., 2019. V.PhyloMaker: an R package that can generate very large phylogenies for vascular plants. Ecography 42, 1353–1359. https://doi.org/10.1111/ecog.04434
}
\author{
F. Munoz}
\seealso{
\code{\link[V.PhyloMaker]{phylo.maker}}
}
\examples{
example <- data.frame(Binome = c("Dipterocarpus indicus", "Michelia champaca", 
"Drypetes roxburghii", "Corypha umbraculifera"), Genus = c("Dipterocarpus", "Michelia", 
"Drypetes", "Corypha"), Family_APGIII = c("Dipterocarpaceae", "Magnoliaceae", 
"Euphorbiaceae", "Aceraceae"), stringsAsFactors = FALSE)

phylo <- create.phylo(example)
plot(phylo$scenario.3)
}
