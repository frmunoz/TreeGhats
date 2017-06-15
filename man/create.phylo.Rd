\name{create.phylo}
\alias{create.phylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating a phylogeny of selected Western Ghats trees}
\description{
Creates a phylogenic tree of taxa based on species, genus and APGIII family information}
\usage{
create.phylo(names = NULL, scenarios = "S3")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{names}{Table including taxonomic information for each taxon on row, and columnes denoted "Binome", "Genus" and "Family_APGIII", with latine genus species combination, genus name and APGIII complying family name, respectively. Genus and species names in "Binome" must be separated with a space.
  With default value NULL, the phylogeny is generated for all the taxa of the TreeGhats database.
}
  \item{scenarios}{One of the scenarios of Qian and Jin (2015) proposed for including taxa that are not already present in the worldwide phylogeny.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Qian H. & Jin Y. (2015). An updated megaphylogeny of plants, a tool for generating plant phylogenies, and an analysis of phylogenetic community structure. Journal of Plant Ecology.}
\author{
F. Munoz}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link{qian}}.
}
\examples{
example <- data.frame(Binome=c("Dipterocarpus indicus","Michelia champaca","Drypetes roxburghii","Corypha umbraculifera"),Genus=c("Dipterocarpus","Michelia","Drypetes","Corypha"),Family_APGIII=c("Dipterocarpaceae","Magnoliaceae","Euphorbiaceae","Aceraceae"),stringsAsFactors = F)

phylo <- create.phylo(example)
plot(phylo$Scenario.3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
