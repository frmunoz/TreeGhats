\name{taxocheck}
\alias{taxocheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Provide taxonomic status of a tree names list 

}
\description{Check the status of a vector of taxon names, according to Tropical Botanic Garden and Research Institute of India (TBGRI) and The Plant List (TPL), from the TreeGhats database (\code{TreeGhatsData}). The function corrects orthographical errors in taxon names and provides the status and the reference names according to TBGRI and TPL. A synthese between TPL and TBGRI is done and a status, a reference name and a reference authority are proposed. The family name according to APGIII is provided.  The origin, the habit, the phenology and the IUCN status of the taxon are provided if \code{otherinfo = TRUE}, if the taxon name is present in (\code{TreeGhatsData}).

}
\usage{
taxocheck(names, otherinfo = T,  max.distance = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{names}{A vector of taxon names, each element including the genus, the specific epithet and, additionally, the infraspecific rank and the infraspecific epithet. The elements must be separated with a space. The authority of the taxon name should not be included.
}
\item{otherinfo}{Logical. If 'TRUE' (default), the biogeographic origin, the habit, the phenology and the IUCN status are provided for the taxon names found in \code{TreeGhatsData}.
}
\item{max.distance}{Maximum distance allowed for a match in \code{agrep} function when performing corrections of spelling errors.
}
}

\details{
First, the function standardizes the taxon names by removing unnecessary spaces, homogenising abbreviations of intraspecific rank and checking for incomplete names, i.e. missing species epithet, or names including a number. Then the function searches taxon names in the \code{TreeGhatsData}  and provides their taxonomic status according to TPL and TBGRI. Typos can be detected depending on the arguments \code{max.distance} and the correct spelling is provided in \code{FoundName}.
Some caution is recommended in increasing \code{max.distance} to avoid incorrect matching. If the taxon name is absent from \code{TreeGhatsData}, the function used a modified version of the function \code{TPLck} (see \code{\link{TPLck2}}) in order to provide the taxonomic status of the taxon name according to TPL.
The family according to APGIII is provided for each taxon name using the function \code{apg_lookup} from the package \pkg{Taxize}. 
If \code{otherinfo = T}, the origin, the habit, the phenology and the IUCN status of the reference names proposed are extracted from the \code{TreeGhatsData}. For taxon names absent from the TreeGhats database, the function return \code{NA}.
}
\value{
The function return an object of class \code{data.frame} with the following components:

  \item{FoundName}{taxon names as found TreeGhats database or when it is absent in TPL. Orthographic errors are corrected. The taxon names provided by the user is kept in the rownames of the return object.}
  
  \item{Typo}{Logical. If 'TRUE' there was a spelling error in the taxon name provided by the user. The correct spelling can be found in the \code{Foundname} column.}
  
  \item{ID_TPL}{The Plant List record 'id'.}
  
  \item{Status_TPL}{Taxonomic status found in The Plant List, either 'Accepted', 'Synonym', 'Unresolved', or 'Misapplied' (see \url{http://www.theplantlist.org/1.1/about/#status}). \code{NA} indicates taxon names not found in TPL. When TPL found several homonyms for the taxon name provided the function returns 'TPLfoundSeveralHomonyms'.}
  
  \item{ReferenceName_TPL}{The taxon name extracted from TPL.}
  
  \item{ReferenceAuthority_TPL}{Names of the author(s) who first published the reference name, extracted from TPL.}
  
  \item{NewID_TPL}{The Plant List record 'id' for the reference name, extracted from TPL.}
  
  \item{Status_TBGRI}{Taxonomic status according to TBGRI and extracted from the TreeGhats database, either 'Absent', 'Accepted', 'Synonym', 'Unresolved'. 'Absent' indicates taxon names not found in TBGRI.}
  
  \item{ReferenceName_TBGRI}{The taxon name according to TBGRI and extracted from the TreeGhats data base.}
  
  \item{ReferenceAuthorithy_TBGRI}{Names of the author(s) who first published the reference name, according to TBGRI and extracted from the TreeGhats database.}
  
  \item{Status_proposed}{Taxonomic status proposed, either 'Accepted', 'Synonym', 'Unresolved'. When a name is absent from TreeGhats database and TPL found several homonyms for the taxon name provided the function returns 'TPLfoundSeveralHomonyms'. When the taxon name is incomplete (missing species epithet or included a number), the function returns 'IncompleteName' as \code{Status_proposed}.}
  
  \item{ReferenceName_proposed}{The taxon name proposed.}
  
  \item{ReferenceAuthorithy_proposed}{Names of the author(s) who first published the reference name proposed.}
  
  \item{Family_APGIII}{Family name according to APGIII for the plant reference name proposed.}
  
  \item{Genus}{Genus of plant reference name proposed.}
  
  \item{Species}{Specific epithet  of plant reference name proposed.}
  
  \item{InfrataxonRank}{Intraspecific rank of plant reference name proposed.}
  
  \item{InfrataxonName}{Intraspecific epithet of taxon reference name proposed.}
  
  \item{Origin}{if \code{otherinfo = T}, the origin of the plant extracted from TreeGhats database according to TBGRI or IFP expertise, either 'exotic', 'endemic of Western Ghats', 'endemic of india' or 'native'.  When the taxon name is not present in TreeGhats database, the function returns \code{NA}.}
  
  \item{Habit}{if \code{otherinfo = T}, the habit of the plant extracted from TreeGhats database according to TBGRI or IFP expertise, either 'liana', 'shrub', 'tree','...','unknown' ,or any combination of these terms. When the taxon name is not present in TreeGhats database, the function returns \code{NA}.}
  
  \item{Phenology}{if \code{otherinfo = T}, the leaf phenology of the plant extracted from TreeGhats database or according to IFP expertise, either 'evergreen' or 'deciduous'.}
  
  \item{IUCN}{if \code{otherinfo = T}, the IUCN status, either 'CR', 'DD', 'EN', 'EX', 'LC', 'LR/cd', 'LR/lc', 'LR/nt', 'NE', 'VU' for the reference name proposed. When the taxon name is not present in TreeGhats database or not evaluated by IUCN, the function returns \code{NA}.}
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

\note{The function only considers subspecies, varieties and forms.
The function relies on the TreeGhats database, which is included in the package. This database will be updated.
The running time of the function depends on the length of the taxon names vector and which arguments are TRUE (if \code{otherinfo = T} the running time is longer).

The function uses internet connection to check names in TPL, get the family according to APG III.

The function return the warnings relative to the use of the \code{TLPck2} function.}
