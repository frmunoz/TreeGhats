\name{TreeGhatsData}
\alias{TreeGhatsData}
\docType{data}
\title{Western Ghats trees database
}
\description{
A comprehensive database of tree species of the Western Ghats with an up-to-date taxonomic reference. It contains all the taxon names that have been commonly used in forest ecology research on Western Ghats over the last decades. It considers infra taxa. The Database is dynamic and will be updated. 
}
\usage{data("TreeGhatsData")}
\value{
  A data frame with the following information.
  
  \item{Full_name}{Plant names without authority.}
  
  \item{Spelling_variants}{Other spelling found for plant names}
  \item{Family_TBGRI}{Family names according to the Tropical Botanic Garden Research Institute (TBGRI).The field is filled when the name is present in TBGRI}
  \item{Family_APGIII}{Family names according to APGIII.}
  
  \item{Genus}{Genus of plant names.}
  
  \item{Species}{Specific epithet of plant names.}
  
  \item{InfrataxonRank}{Intraspecific rank of plant names.}
  
  \item{InfrataxonName}{Intraspecific epithet of plant names.}
  
  \item{Authority_InitialDatabases}{Names of the scientist(s) who first published the plant name.}
  \item{Authority_TPL}{Names of the scientist(s) who first published the plant name, according to TPL.}
  \item{Complet_name}{Plant names with authorities.}
  \item{ID_Tropicos}{Tropicos record 'id'.'Absent' indicates plant names not found in Tropicos.}
  \item{ID_TPL}{The Plant List record 'id'.'Absent' indicates plant names not found in TPL.}
  
  \item{Status_TPL}{Taxonomic status found in TPL, either 'Accepted', 'Synonym', 'Unresolved', or 'Misapplied'. 'Absent' indicates plant names not found in TPL.}
  
  \item{ReferenceName_TPL}{Valid forms of plant names, extrated from TPL. The field is filled with \code{NA} when the name is absent in TPL.}
  
  \item{ReferenceAuthority_TPL}{Names of the scientist(s) who first published valid plant names, extracted from TPL. The field is filled with \code{NA} when the name is absent in TPL.}
  
  
  \item{Status_TBGRI}{Taxonomic status according to TBGRI, either 'Absent', 'Accepted', 'Synonym', 'Unresolved'. 'Absent' designes names not found in TBGRI.}
  
  \item{ReferenceName_TBGRI}{Valid forms of plant names, according to TBGRI.The field is filled with \code{NA} when the name is absent in TBGRI.}
  
  \item{ReferenceAuthorithy_TBGRI}{Names of the scientist(s) who first published plant names, according to TBGRI.The field is filled with \code{NA} when the name is absent in TBGRI.}
  
  \item{Status_proposed}{Taxonomic status proposed, either 'Accepted', 'Synonym', 'Unresolved'.}
  
  \item{ReferenceName_proposed}{Valid forms of refrence names proposed.}
  
  \item{ReferenceAuthorithy_proposed}{Names of the scientist(s) who first published reference names proposed.}
  
  \item{Origin}{The origin of plants according to TBGRI or IFP expertise, either 'exotic', 'endemic of Western Ghats', 'endemic of india' or 'native'.}
  
  \item{Habit}{The habit of plants according to TBGRI or IFP expertise, either 'liana', 'shrub', 'tree','...','unknown' ,or any combination of these terms.}
  
  \item{Phenology}{The leaves phenology of plants according to IFP expertise, either 'evergreen' or 'deciduous'.}
  
  \item{IUCN}{The IUCN status of plants, either 'CR', 'DD', 'EN', 'EX', 'LC', 'LR/cd', 'LR/lc', 'LR/nt', 'NE', 'NT', 'VU'.}
  \item{Ecology_habitat}{Some informations are dispensed in this field about particular habitat: 'Planted/Cultivated', 'Coastal plant', 'Mangrove plant', Riparian}
}
\details{
The database is a compilation of different taxa lists :
Species for which the word "tree" appeared in its habit in the book Flowering plants of the western Ghats edited by the Tropical Botanic Garden Research Institute (TBGRI); Species recently described; Species names used in studies done in different parts of the Western Ghats by the French institut of Pondicherry (IFP) or from 2000 and for which the list of taxa is published; Species names used in the Indian Biodiversity Portail (\url{http://indiabiodiversity.org/}). 
}
\source{
Nayar, T., Rasiya Beegam, A. & Sibi, M. 2015. Flowering plants of the Western Ghats, India, Volume 1 Dicots; Volume 2 Monocots. Jawaharlal Nehru Tropical Botanic Garden and Research Insitute.

The International Plant Names Index. 2017. \url{http://www.ipni.org}.

Tropicos.org. Missouri Botanical Garden. 2017. \url{http://www.tropicos.org}.


}
\author{
V. Raevel & F. Munoz
}
\references{
%%  ~~ possibly secondary sources and usages ~~
To be completed
}
\examples{
data(TreeGhatsData)
str(TreeGhatsData)
}
\keyword{datasets}
\keyword{ Western Ghats }
\keyword{ Tree }
\keyword{ Flowering plants }