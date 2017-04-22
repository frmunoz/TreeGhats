\name{TreeGhats}
\alias{TreeGhats}
\docType{data}
\title{Western Ghats trees database
}
\description{
A comprehensive database of the tree species of the Western Ghats with an up-to-date taxonomic reference. It contains all the taxon names that have been commonly used in forest ecology research over the last decades in South India.
}
\usage{data("TreeGhats")}
\value{
  A data frame with n names and the following information.
  
  \item{Full_name}{Plant name without authority.}
  
  \item{Spelling_variants}{Other spelling found for the plant name}
  \item{Family_TBGRI}{Family name according to TBGRI.}
  \item{Family_APGIII}{Family name according to APGIII.}
  
  \item{Genus}{Genus of plant name.}
  
  \item{Species}{Specific epithet of plant name.}
  
  \item{InfrataxonRank}{Intraspecific rank of plant name.}
  
  \item{InfrataxonName}{Intraspecific epithet of plant name.}
  
  \item{Authority_InitialDatabases}{Names of the scientist(s) who first published the plant name}
  \item{Authority_TPL}{Names of the scientist(s) who first published the plant name, according to TPL}
  \item{Complet_name}{Plant name with authority}
  \item{ID_Tropicos}{Tropicos record 'id'.}
  \item{ID_TPL}{The Plant List record 'id'.}
  
  \item{Status_TPL}{Taxonomic status found in TPL, either 'Accepted', 'Synonym', 'Unresolved', or 'Misapplied'. \code{NA} indicates plant names not found in TPL.}
  
  \item{ReferenceName_TPL}{The valid form of the plant name, extrated from TPL.}
  
  \item{ReferenceAuthority_TPL}{Names of the scientist(s) who first published the valid name, extracted from TPL.}
  
  
  \item{Status_TBGRI}{Taxonomic status according to TBGRI, either 'Absent', 'Accepted', 'Synonym', 'Unresolved'. 'Absent' designes names not found in TBGRI.}
  
  \item{ReferenceName_TBGRI}{The valid form of the plant name, according to TBGRI.}
  
  \item{ReferenceAuthorithy_TBGRI}{Names of the scientist(s) who first published the plant name, according to TBGRI.}
  
  \item{Status_proposed}{Taxonomic status proposed, either 'Accepted', 'Synonym', 'Unresolved'.}
  
  \item{ReferenceName_proposed}{The valid form of the refrence name proposed.}
  
  \item{ReferenceAuthorithy_proposed}{Names of the scientist(s) who first published the reference name proposed.}
  
  \item{Origin}{The origin of the plant according to TBGRI or IFP expertise, either 'exotic', 'endemic of Western Ghats', 'endemic of india' or 'native'.}
  
  \item{Habit}{The habit of the plant according to TBGRI or IFP expertise, either 'liana', 'shrub', 'tree','...','unknown' ,or any combination of these terms.}
  
  \item{Phenology}{The leaves phenology of the plant according to IFP expertise, either 'evergreen' or 'deciduous'.}
  
  \item{IUCN}{The IUCN status of the plant names, either 'CR', 'DD', 'EN', 'EX', 'LC', 'LR/cd', 'LR/lc', 'LR/nt', 'NE', 'NT', 'VU'.}
}  
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\authors{
V. Raevel & F. Munoz
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(TreeGhats)
str(TreeGhats)
}
\keyword{datasets}
\keyword{ Western Ghats }
\keyword{ Tree }
\keyword{ Flowering plants }