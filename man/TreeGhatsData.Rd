\name{TreeGhatsData}
\alias{TreeGhatsData}
\docType{data}
\title{Western Ghats tree database
}
\description{
A comprehensive database of tree taxa of the Western Ghats with up-to-date taxonomic reference. It contains taxon names that have been commonly used in forest ecology research on Western Ghats. It also considers infrataxon level. The database is dynamic and will be updated. 
}
\usage{data("TreeGhatsData")}
\value{
  A data frame with the following information.
  
  \item{Name}{Taxon names without authority.}
  
  \item{Spelling_variants}{Other spelling found for taxon names.}
  \item{Family_TBGRI}{Family names according to the Tropical Botanic Garden Research Institute (TBGRI).}
  \item{Family_APGIII}{Family names according to APGIII.}
  
  \item{Genus}{Genus of taxon names.}
  
  \item{Species}{Specific epithet of taxon names.}
  
  \item{InfrataxonRank}{Intraspecific rank of taxon names.}
  
  \item{InfrataxonName}{Intraspecific epithet of taxon names.}
  
  \item{Authority_InitialDatabases}{Names of the author(s) who first published the taxon name presented in the reference lists to build TreeGhatsData.}
  \item{Authority_TPL}{Names of the author(s) who first published the taxon name, according to TPL.}
  \item{Full_name}{Taxon names with authorities.}
  \item{ID_Tropicos}{Tropicos record 'id'. 'Absent' indicates taxon names not found in Tropicos.}
  \item{ID_TPL}{The Plant List record 'id'. 'Absent' indicates taxon names not found in TPL.}
  
  \item{Status_TPL}{Taxonomic status found in TPL, either 'Accepted', 'Synonym' or 'Unresolved' (see \url{http://www.theplantlist.org/1.1/about/#status}). 'Absent' indicates taxon names not found in TPL.}
  
  \item{ReferenceName_TPL}{Reference taxon names, extracted from TPL. The field is filled with \code{NA} when the name is absent in TPL.}
  
  \item{ReferenceAuthority_TPL}{Names of the author(s) who first published the TPL reference taxon names. The field is filled with \code{NA} when the name is absent in TPL.}
  
  
  \item{Status_TBGRI}{Taxonomic status according to TBGRI, either 'Absent', 'Accepted', 'Synonym' or 'Unresolved'. 'Absent' designs names not found in TBGRI.}
  
  \item{ReferenceName_TBGRI}{Taxa names, according to TBGRI. The field is filled with \code{NA} when the name is absent in TBGRI.}
  
  \item{ReferenceAuthorithy_TBGRI}{Names of the author(s) who first published the TBGRI reference taxon names. The field is filled with \code{NA} when the name is absent in TBGRI.}
  
  \item{Status_proposed}{Taxonomic status proposed, either 'Accepted', 'Synonym' or 'Unresolved'.}
  
  \item{ReferenceName_proposed}{Reference names proposed.}
  
  \item{ReferenceAuthorithy_proposed}{Names of the author(s) who first published the reference names proposed.}
  
  \item{Origin}{Biogeographic origin according to TBGRI or IFP expertise, either 'exotic', 'endemic of Western Ghats', 'endemic of india' or 'native'.}
  
  \item{Habit}{Habit according to TBGRI or IFP expertise, either 'climbing shrub','liana', 'shrub', 'tree','unknown' ,or any combination of these terms.}
  
  \item{Phenology}{Leaf phenology of plants according to IFP expertise, either 'evergreen' or 'deciduous'.}
  
  \item{IUCN}{IUCN status of plants, either 'CR', 'DD', 'EN', 'EX', 'LC', 'LR/cd', 'LR/lc', 'LR/nt', 'NE', 'VU'. \code{NA} correspond to the non evaluated taxa.}
  \item{Ecology_habitat}{particular habitats are precised, either 'Planted/Cultivated', 'Coastal', 'Mangrove'or 'Riparian'.}
}
\details{
The database is a compilation of different taxa lists :
taxa for which the word "tree" appeared in habit in the book Flowering plants of the western Ghats edited by the Tropical Botanic Garden Research Institute (TBGRI) including planted or cultivated taxa; Taxon described after 2014; taxon names used in studies done in different parts of the Western Ghats by the French institut of Pondicherry (IFP) or from 2000 and for which the list of taxa is published; taxon names reported as "tree" growth form used in the Indian Biodiversity Portail (\url{http://indiabiodiversity.org/}). As most of studies did not distinguish the infraspecific level, for each infra species (sub-species, variety or form), we automatically add the species level as an entry in our taxonomic database. The addition of the species level correspond to a degraded information. Indeed, according to the International code of nomenclature, the publication of an infrataxon name automatically created the infrataxon autonym corresponding to the taxon level.
The version of the database available in TreeGhats package has been published on Zenodo on August 21st (https://zenodo.org/record/846291).
}
\source{
Raevel, V., Ayyappan, N., Balachandran, N., Aravajy, S., Barathan, N., Ramesh, B. R. & Munoz, F. (2017) Taxonomic and ecological database of trees of Western Ghats - TreeGhatsData (ed Zenodo), https://zenodo.org/record/846291

Nayar, T., Rasiya Beegam, A. & Sibi, M. 2015. Flowering plants of the Western Ghats, India, Volume 1 Dicots; Volume 2 Monocots. Jawaharlal Nehru Tropical Botanic Garden and Research Insitute.

Ramesh, B.R., and Jean-Pierre Pascal. 1997. Atlas of Endemics of the Western Ghats (India): Distribution of Tree Species in the Evergreen and Semi-Evergreen Forests. French Institute of Pondicherry.

The International taxon names Index. 2017. \url{http://www.ipni.org}.

Tropicos.org. Missouri Botanical Garden. 2017. \url{http://www.tropicos.org}.

India biodiversity Portal.2017.\url{http://indiabiodiversity.org/}

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
\keyword{database}
\keyword{Western Ghats}
\keyword{Tree}
\keyword{ Flowering plants }

\note{
Two limitations of the database have been identified. First, it does not content all synomym names for a taxon. Second, the attribution of 'tree' as a growth form to a taxon depends on the criteria used in the study and/or on the study area. This can be explain missing taxa in the database.  The database will be updated at least every two years or more often depending on the requests send to <treeghats@gmail.com>. }