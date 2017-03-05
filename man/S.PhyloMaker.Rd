\name{S.PhyloMaker}
\alias{S.PhyloMaker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function of Qian and Jin (2015) to generate a phylogenetic tree}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
S.PhyloMaker(splist, tree, nodes, output.splist = T, scenarios = "S3")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{splist}{
%%     ~~Describe \code{splist} here~~
A user specified data frame which contains a list of seed plant species for which S.PhyloMaker generates phylogenies. This data frame should include three columns with field names as species, genus, and family. The three columns are, respectively, for species name (e.g. Acacia berlandieri), genus name (e.g. Acacia), and family name (e.g. Fabaceae)
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
A phylo object; the megaphylogeny of vascular plants (i.e. PhytoPhylo). 
}
  \item{nodes}{
%%     ~~Describe \code{nodes} here~~
A data.frame; it contains the basal node information of every family and genus on the megaphylogeny (i.e. PhytoPhylo).
}
  \item{output.splist}{
%%     ~~Describe \code{output.splist} here~~
Logical; the output species list includes the columns of the "splist" data frame and an extra column showing which species in "splist" have matched with species in the megaphylogeny. 
}
  \item{scenarios}{
%%     ~~Describe \code{scenarios} here~~
Optional; defining which scenario(s) will be used to generate phylogenies. By default, scenario "S3" is performed. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
(The function description is derived from the Appendix of Qian and Jin (2015))

Using the PhytoPhylo megaphylogeny as a backbone, this function takes a user-specified species list and matches species in the species list with species in the megaphylogeny of Qian and Jin (2015). If a species is found in the megaphylogeny, the species will be selected for pruning; if the species is not found but its genus or family is found in the megaphylogeny, the species will be added to the backbone at a certain place within the genus or family depending on which scenarios are taken. Finally, the selected and added species are pruned from the megaphylogeny, which results in three phylogenies corresponding to the three scenarios (when default is chosen). 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The output contains a phylogenetetic tree in phylo format corresponding to the scenario defined in input, and a species list, which includes the original columns of the "spList" data frame and an additional column (called "status") showing which species in the user’s species list have matched with species in the megaphylogeny. Species matched between the species list and the megaphylogeny are indicated as "match(prune)"; species that are added to the megaphylogeny before the pruning are indicated as "match(add)"; species whose genera or families are not found in the megaphylogeny are indicated as "unmatch". 
}
\references{
Qian H. & Jin Y. (2015). An updated megaphylogeny of plants, a tool for generating plant phylogenies, and an analysis of phylogenetic community structure. Journal of Plant Ecology.
}
\author{
F. Munoz, adapted from Qian and Jin (2015)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library("phytools") # load the "phytools" package.

example <- 

data(qian)
result<-S.PhyloMaker(splist=example, tree=qian$phylo.all, nodes=qian$nodes)      # run the function S.PhyloMaker.
str(result)       # the structure of the ouput of S.PhyloMaker.

plot(result$Scenario.1,cex=1.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
