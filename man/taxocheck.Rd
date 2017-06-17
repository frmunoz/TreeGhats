\name{taxocheck}
\alias{taxocheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Check the status for a list of plant names according to TBGRI and The Plant List 

}
\description{Check the status of a vector of plant names, according to the TreeGhats database. This database relies on the Tropical Botanic Garden and Research Institute of India (TBGRI) and The Plant List (TPL). The function removes orthographical errors in plant names and provides the accepted names according to TBGRI and TPL. A synthese between TPL and TBGRI is done and a status, a reference name and a reference authority are proposed. The family according to APGIII is provided.  The origin, the habit and the phenology of the taxon are provided if \code{othersinfo = TRUE}, if the plant name is present in the TreeGhats database.The IUCN status of the taxon is provided if \code{iucn = TRUE}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
taxocheck(names, othersinfo = T, iucn = F, max.distance = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{names}{A vector of plant names, each element including the genus, the specific epithet and, additionally, the infraspecific epithet. The elements must be separated with a space. The authority of the taxon name should not be included.
}
\item{othersinfo}{Logical. If 'TRUE' (default), then the origin, the habit and the phenology are provided for the plant names found in the TreeGhats database.
}
\item{iucn}{Logical. If 'TRUE' (default), the IUCN status is provided for the reference name proposed.
}
\item{max.distance}{Maximum distance allowed for a match in agrep function when performing corrections of spelling errors.
}
}

\details{
First, the function standardizes the plant names by removing unnecessary spaces, homogenising abbreviations of intraspecific rank and checking for incomplete names, i.e. missing species epithet or names including a number. Then the function searches plant names in the TreeGhats database and provides their taxonomic status according to TPL and TBGRI. Differences in typos can be detected depending on the arguments 'max.distance' and the correct spelling is provided in foundNames.
Some caution is recommended in increasing 'max.distance' as false match can be done. If the taxon name is absent from the database, the function used a modified version of the function \code{TPLck} (see ) in order to provide the taxonomic status of the plant according to TPL.
The family according to APGIII is provided for each plant name using the function \code{apg_lookup} from the package Taxize. 
If \code{othersinfo = T}, the origin, the habit, the phenology of the reference names proposed are extracted from the TreeGhats database. For plant names absent from the TreeGhats database, the function return \code{NA}.
If \code{iucn = T}, the IUCN status is provided using the function \code{lets.iucn} from the package Taxize.
}
\value{
The function return an object of class \code{data.frame} with the following components:

  \item{FoundName}{Plant names as found TreeGhats database or when it is absent in TPL. Orthographic errors are corrected. The plant names provided by the user is kept in the rownames of the return object.}
  
  \item{Typo}{Logical. If 'TRUE' there was a spelling error in the taxon name provided by the user. The correct spelling can be found in the Foundnames.}
  
  \item{ID_TPL}{The Plant List record 'id'.}
  
  \item{Status_TPL}{Taxonomic status found in TPL, either 'Accepted', 'Synonym', 'Unresolved', or 'Misapplied'. \code{NA} indicates plant names not found in TPL. When TPL found several homonyms for the plant name provided the function returns 'TPLfoundSeveralHomonyms'.}
  
  \item{ReferenceName_TPL}{The valid form of the plant name, extrated from TPL.}
  
  \item{ReferenceAuthority_TPL}{Names of the scientist(s) who first published the reference name, extracted from TPL.}
  
  \item{NewID_TPL}{The Plant List record 'id' for the reference name, extracted from TPL.}
  
  \item{Status_TBGRI}{Taxonomic status according to TBGRI and extract from the TreeGhats data base, either 'Absent', 'Accepted', 'Synonym', 'Unresolved'. 'Absent' designes species names not found in TBGRI.}
  
  \item{ReferenceName_TBGRI}{The valid form of the taxon name, according to TBGRI and extract from the TreeGhats data base.}
  
  \item{ReferenceAuthorithy_TBGRI}{Names of the scientist(s) who first published the reference name, according to TBGRI and extract from the TreeGhats database.}
  
  \item{Status_proposed}{Taxonomic status proposed, either 'Accepted', 'Synonym', 'Unresolved'. When a name is absent from TreeGhats database and TPL found several homonyms for the taxon name provided the function returns 'TPLfoundSeveralHomonyms'. When the taxon name is incomplete (missing species epithet or included a number), the function returns 'IncompleteName' as Status_proposed}
  
  \item{ReferenceName_proposed}{The valid form of the plant name proposed.}
  
  \item{ReferenceAuthorithy_proposed}{Names of the scientist(s) who first published the reference name proposed.}
  
  \item{Family_APGIII}{Family name according to APGIII for the plant reference name proposed.}
  
  \item{Genus}{Genus of plant reference name proposed.}
  
  \item{Species}{Specific epithet  of plant reference name proposed.}
  
  \item{InfrataxonRank}{Intraspecific rank of plant reference name proposed.}
  
  \item{InfrataxonName}{Intraspecific epithet of taxon reference name proposed.}
  
  \item{Origin}{if \code{othersinfo = T}, the origin of the plant extracted from TreeGhats database according to TBGRI or IFP expertise, either 'exotic', 'endemic of Western Ghats', 'endemic of india' or 'native'.  When the plant name is not present in TreeGhats database, the function returns \code{NA}.}
  
  \item{Habit}{if \code{othersinfo = T}, the habit of the plant extracted from TreeGhats database according to TBGRI or IFP expertise, either 'liana', 'shrub', 'tree','...','unknown' ,or any combination of these terms. When the plant name is not present in TreeGhats database, the function returns \code{NA}.}
  
  \item{Phenology}{if \code{othersinfo = T}, the leaves phenology of the plant extracted from TreeGhats database or according to IFP expertise, either 'evergreen' or 'deciduous'.}
  
  \item{IUCN}{if \code{iucn = T}, the IUCN status, either 'CR', 'DD', 'EN', 'EX', 'LC', 'LR/cd', 'LR/lc', 'LR/nt', 'NE', 'NT', 'VU' for the reference name proposed. The search is done at the species level. When the plant name is not present in TreeGhats database, the function returns \code{NA}.}
}


\references{
Cayuela, L., Granzow-de la Cerda, I., Albuquerque, F.S. and Golicher, J.D. 2012. Taxonstand: An R package for species names standardisation in vegetation databases. Methods in Ecology and Evolution, 3(6): 1078-1083.

Kalwij, J.M. 2012. Review of 'The Plant List, a working list of all plant species'. Journal of Vegetation Science, 23(5): 998-1002.

Nayar, T., Rasiya Beegam, A. & Sibi, M. 2015. Flowering plants of the Western Ghats, India, Volume 1 Dicots; Volume 2 Monocots. Jawaharlal Nehru Tropical Botanic Garden and Research Insitute.

Nayar, T.S., Beegam, A.R. & Sibi, M. 2014. Flowering Plants of the Western Ghats, India: Dicots. Jawaharlal Nehru Tropical Botanic Garden and Research Institute.
}

\author{
V. Raevel & F. Munoz
}

\seealso{
 See Also  \code{\link{TreeGhatsData}}, which is the database used by the function.
 
 See Also  \code{\link[Taxonstand]{TPLck}}
 
 See Also  \code{\link{TPLck2}}
 
 See Also  \code{\link[taxize]{get_uid}}
 
 See Also  \code{\link[taxize]{apg_lookup}}
 
 See Also  \code{\link[letsR]{lets.iucn}}
 
}

\examples{
example <- c("Dipterocarpus indicus","Michelia champaca","Drypetes roxburghii",
"Corypha umbraculifera")
\dontrun{
result <- taxocheck(example)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Western Ghats }
\keyword{ Tree }

\note{The function only consider subspecies, varieties and forms.
The function relies on the TreeGhats database, which is included in the package. This database will be updated.
The running time of the function depends on the length of the taxon names vector and which arguments are TRUE (if \code{othersinfo = T} the running time is longer).

The function uses internet connection to check names in TPL, get the family according to APG III and the IUCN status.

The function return the warnings relative to the use of the \code{TLPck2} function.}
