\name{taxocheck}
\alias{taxocheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Provide taxonomic status of a tree names list 

}

\description{Check the taxonomic status of a vector of taxon names, according to the Tropical Botanic Garden and Research Institute of India (TBGRI) and The Plant List (TPL), from the TreeGhats database (\code{TreeGhatsData}). The function corrects orthographical errors in taxon names and provides the taxonomic status and the reference names according to TBGRI and TPL. From these information a status, a reference name and  authority are proposed. For the reference taxon proposed, the family name according to APGIII is provided.  The biogeographic origin, the habit, the phenology and the IUCN status are also provided, if \code{otherinfo = TRUE} and if the taxon name is present in \code{\link{TreeGhatsData}}.
}
\usage{
taxocheck(names, otherinfo = T,  max.distance = 2, phylo = F)
}
\arguments{
\item{names}{Typically a vector of taxon names, each element including a genus, a specific epithet and, additionally, an infraspecific rank and an infraspecific epithet. The elements must be separated with a space. The authority of the taxon name should not be included.
If the user instead provides a data frame, \code{taxocheck} selects the first column named "names" or "binome" (case insensitive). This column should include taxon names following the same format as presented above.
}
\item{otherinfo}{Logical. If 'TRUE' (default), the biogeographic origin, the habit, the phenology and the IUCN status are provided for the taxon names found in \code{TreeGhatsData}.
}
\item{max.distance}{Maximum distance allowed for a match in \code{agrep} function when performing corrections of spelling errors.
}
\item{phylo}{If 'TRUE', the function will return an object of class \code{phylo} with the phylogenetic tree. Default is 'FALSE'.
}
}

\details{
First, the function standardizes the taxon names by removing unnecessary spaces and homogenising abbreviations of intraspecific rank, and indicates the incomplete names, i.e. if the species epithet is missing or the name includes a number. Second the function searches the taxon names in the \code{TreeGhatsData}  and provides their taxonomic status according to TPL and TBGRI. Typos can be detected depending on the arguments \code{max.distance} and the correct spelling is provided in \code{FoundName}.Orthographic errors can be corrected only in specific epithets. Some caution is recommended in increasing \code{max.distance} to avoid incorrect matching. Third and if the taxon name is absent from \code{TreeGhatsData}, the function provides the taxonomic status of the taxon name according to TPL, by using a modified version of the function \code{\link{TPLck}} (see \code{\link{TPLck2}}), which indicates if several homonyms exists in TPL (i.e. names with identical spelling but decribed by differents authority and so corresponding to different taxa).  
From the status found in TPL and TBGRI, a status, a reference name and authority are proposed, according to the following rules: if a taxon name is absent from \code{TreeGhatsData}, we retain the information given by TPL otherwise we retain the taxonomic status poposed in \code{TreeGhatsData}. For taxon names with several homonyms in TPL or incomplete, the function returns respectevely 'SeveralHomonyms' and 'IncompleteName' as the status and \code{NA} for reference name and authority; For a taxon name provided by the user at the species level, if only one infrataxon is present in \code{TreeGhatsData} and so in Western Ghats, the function returns the information concerning this infrataxon otherwise the function returns the information concerning the species level and  'SeveralInfrataxa' in the comment column, as the infrataxon level was not identified.   
The family of the reference taxon proposed according to APGIII is provided using the function \code{apg_lookup} from the package \pkg{Taxize}. If \code{otherinfo = T}, the biogeographic origin, the habit, the phenology and the IUCN status of the reference taxon proposed are extracted from the \code{TreeGhatsData}.
}
\value{
If \code{phylo = F} (default), the function returns an object of class \code{data.frame} with the following components:

  \item{FoundName}{taxon names as found in \code{TreeGhatsData} or when it is absent in TPL. Orthographic errors are corrected. The taxon names provided by the user is kept in the rownames of the return object. 'IncompleteName' is return when the user provide an name without species epithet, or including a number. \code{NA} indicates taxon names found neither in \code{TreeGhatsData} nor in TPL.}
  
  \item{Typo}{Logical. If 'TRUE' there was a spelling error in a taxon name provided by the user. The correct spelling can be found in the \code{Foundname} column, while the spelling provided by the user is kept as row name.}
  
  \item{Status_TPL}{Taxonomic status found in The Plant List, either 'Accepted', 'Synonym', 'Unresolved', or 'Misapplied' (see \url{http://www.theplantlist.org/1.1/about/#status}). \code{NA} indicates taxon names not found in TPL. 'SeveralHomonyms' indicates several homonyms in TPL for the taxon name provided.}
  
  \item{ReferenceName_TPL}{The reference taxon name extracted from TPL.}
  
  \item{ReferenceAuthority_TPL}{Names of the author(s) who first published the reference name extracted from TPL.}
  
  \item{Status_TBGRI}{Taxonomic status according to TBGRI and extracted from \code{TreeGhatsData}, either 'Absent', 'Accepted', 'Synonym', 'Unresolved'. 'Absent' indicates taxon names not found in TBGRI. \code{NA} indicates taxon names not found in \code{TreeGhatsData} i.e. taxon names not checked in TBGRI.}
  
  \item{ReferenceName_TBGRI}{The reference taxon name according to TBGRI and extracted from \code{TreeGhatsData}.}
  
  \item{ReferenceAuthorithy_TBGRI}{Names of the author(s) who first published the reference name, according to TBGRI and extracted from \code{TreeGhatsData}.}
  
  \item{Status_proposed}{Taxonomic status proposed, either 'Accepted', 'Synonym', 'Unresolved','IncompleteName' or 'SeveralHomonyms'. 'IncompleteName' indicates incomplete taxon names provided by the user. 'SeveralHomonyms' indicates the presence of several homonyms in TPL.}
  
  \item{ReferenceName_proposed}{The taxon name proposed.}
  
  \item{ReferenceAuthorithy_proposed}{Names of the author(s) who first published the reference taxon proposed.}
  
  \item{Comment}{Indicate if the taxon name given by the user has several infra taxa present in the WG according to TBGRI. If this field is empty, it is not shown.}
  
  \item{Family_APGIII}{Family name according to APGIII for the reference taxon proposed.}
  
  \item{Genus}{Genus of reference taxon proposed.}
  
  \item{Species}{Specific epithet of reference taxon proposed.}
  
  \item{InfrataxonRank}{Intraspecific rank of reference taxon proposed. If this field is empty, it is not shown.}
  
  \item{InfrataxonName}{Intraspecific epithet of reference taxon proposed. If this field is empty, it is not shown.}
  
  \item{Origin}{if \code{otherinfo = T}, the biogegraphical origin of the reference taxon proposed extracted from \code{TreeGhatsData}, either 'exotic', 'endemic of Western Ghats', 'endemic of india' or 'native'. \code{NA} indicates that the reference taxon proposed is absent from \code{TreeGhatsData}.}
  
  \item{Habit}{if \code{otherinfo = T}, the habit of the reference taxon proposed extracted from from \code{TreeGhatsData}, either 'liana', 'shrub', 'tree', '...', 'unknown' ,or any combination of these terms. \code{NA} indicates that the reference taxon proposed is absent from \code{TreeGhatsData}.}
  
  \item{Phenology}{if \code{otherinfo = T}, the leaf phenology of the reference taxon proposed extracted from \code{TreeGhatsData}, either 'evergreen' or 'deciduous'.}
  
  \item{IUCN}{if \code{otherinfo = T}, the IUCN status, either 'CR', 'DD', 'EN', 'EX', 'LC', 'LR/cd', 'LR/lc', 'LR/nt', 'NE', 'VU' for the reference taxon proposed. When the taxon name is not present in \code{TreeGhatsData} or not evaluated by IUCN, the function returns \code{NA}.}
  
If \code{phylo = T}, the function returns a list with the same data frame as above in \code{$tab}, plus an object of class \code{phylo} in \code{$phylo} providing the phylogenetic tree of the requested taxa.

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
\dontrun{
example <- c("Dipterocarpus indicus","Michelia champaca","Drypetes roxburghii",
"Corypha umbraculifera")
result <- taxocheck(example)

# Tree species from the Uppangala permanent forest plot
test <- taxocheck(UPSP_Species, phylo = TRUE)
plot(test$phylo)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Western Ghats }
\keyword{ Tree }

\note{
The function relies on \code{TreeGhatsData} included in the package, which will be update.
The running time of the function depends on the length of the taxon names vector and the internet connection used to check names in TPL and get the family according to APG III. The function return the warnings relative to the use of the \code{TLPck2} function.
Various limitations have been identified to date in the implementation of this function. 
First, homonyms (i.e. a name for a taxon that is identical in spelling to another such name, but belongs to a different taxon) cannot be identified. Second, if the input infraspecific epithet does not match any of the infraspecific epithets provided in TPL, then the first accepted name with no infraspecific epithet is selected as the best match. If all names provided by TPL have infraspecific epithets and are all synonyms (e.g. Pottia starkeana), then there is no best match and the output will match the original name.}