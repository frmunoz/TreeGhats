\name{taxocheck}
\alias{taxocheck}
\title{Provide taxonomic status of a tree names list 
}

\description{Check the taxonomic status of a vector of taxon names, according to the Tropical Botanic Garden and Research Institute of India (TBGRI) and The Plant List (TPL), from the \code{\link{TreeGhatsData}} database. The function corrects orthographical errors in taxon names and provides the taxonomic status and the reference names according to TBGRI and TPL. From these information a status, a reference name and  authority are proposed. For the reference taxon proposed, the family name according to APGIII is provided.  The biogeographic origin, the habit, the phenology and the IUCN status are also provided, if \code{otherinfo = TRUE} and if the taxon name is present in \code{\link{TreeGhatsData}}.A phylogenetic tree of the data is provided if \code{phylo = TRUE}.
}
\usage{
taxocheck(names, otherinfo = T,  max.distance = 2, phylo = F)
}
\arguments{
\item{names}{typically a vector of taxon names, each element including a genus, a specific epithet and, additionally, an infraspecific rank and an infraspecific epithet. The elements must be separated with a space. The authority of the taxon name should not be included.
If the user instead provides a data frame, \code{taxocheck} selects the first column named "names" or "binome" (case insensitive). This column should include taxon names following the same format as presented above.
}
\item{otherinfo}{logical. If 'TRUE' (default), the biogeographic origin, the habit, the phenology and the IUCN status are provided for the taxon names found in \code{\link{TreeGhatsData}}.
}
\item{max.distance}{maximum distance allowed in spelling correction with \code{\link[base]{agrep}} function.
}
\item{phylo}{if 'TRUE', the function will return an object of class \code{phylo} with the phylogenetic tree. Default is 'FALSE'.
}
}

\details{
First, the function standardizes the taxon names by removing unnecessary spaces and homogenising the abbreviations of intraspecific rank. It identifies incomplete names, i.e. if the species epithet is missing or the name includes a number instead of a species name. Second, the function searches the taxon names in the database \code{\link{TreeGhatsData}} and provides their taxonomic status according to TPL and TBGRI. Typos can be detected depending on the arguments \code{max.distance} and the correct spelling is provided in \code{FoundName}. Orthographic errors can be corrected only in specific epithets. \code{max.distance} should not be too large, to avoid incorrect matching. Third, for taxon names absent from \code{TreeGhatsData}, the function provides taxonomic status according to TPL, by using a modified version of the function \code{\link[Taxonstand]{TPLck}} (see \code{\link{TPLck2}}). The modified version indicates if several homonyms are present in TPL (i.e., names with identical spelling but decribed by different authorities, thus corresponding to different taxa).  

Based on the statuses found in TPL and TBGRI, a status, a reference name and authority are proposed, with the following rules: if a taxon name is absent from \code{\link{TreeGhatsData}}, we retain information given by TPL, and otherwise we retain the taxonomic status poposed in \code{\link{TreeGhatsData}}. For taxon names with several homonyms in TPL or incomplete, the function returns 'SeveralHomonyms' and 'IncompleteName' as the status, respectively, and \code{NA} as reference name and authority. For a taxon name provided by the user at species level, if only one infrataxon is present in \code{\link{TreeGhatsData}}, the function returns the information concerning this infrataxon, otherwise the function returns the information concerning the species level and includes 'SeveralInfrataxa' in the comment column. In the later case, the infrataxon level cannot be identified.   
If \code{otherinfo = T}, the biogeographic origin, the habit, the phenology and the IUCN status of the reference taxon proposed are extracted from the \code{\link{TreeGhatsData}}.
}
\value{
If \code{phylo = F} (default), the function returns an object of class \code{data.frame} with the following components:

  \item{FoundName}{taxon names as found in \code{\link{TreeGhatsData}} or in TPL. Orthographic errors are corrected. The taxon names provided by the user is kept in the rownames of the return object. 'IncompleteName' is return when the user provide a name without species epithet, or including a number. \code{NA} indicates taxon names found neither in \code{\link{TreeGhatsData}} nor in TPL.}
  
  \item{Typo}{logical. If 'TRUE', there was a spelling error in a taxon name provided by the user. The correct spelling can be found in the \code{Foundname} column, while the spelling provided by the user is kept as row name.}
  
  \item{Status_TPL}{taxonomic status found in The Plant List, either 'Accepted', 'Synonym', 'Unresolved', or 'Misapplied' (see \url{http://www.theplantlist.org/1.1/about/#status}). \code{NA} indicates taxon names not found in TPL. 'SeveralHomonyms' indicates that there are several homonyms in TPL for the taxon name provided.}
  
  \item{ReferenceName_TPL}{the reference taxon name extracted from TPL.}
  
  \item{ReferenceAuthority_TPL}{authority of \code{ReferenceName_TPL}.}
  
  \item{Status_TBGRI}{taxonomic status according to TBGRI and extracted from \code{\link{TreeGhatsData}}, either 'Absent', 'Accepted', 'Synonym'. 'Absent' indicates taxon names not found in TBGRI.}
  
  \item{ReferenceName_TBGRI}{the reference taxon name according to TBGRI extracted from \code{\link{TreeGhatsData}}.}
  
  \item{ReferenceAuthorithy_TBGRI}{authority of \code{ReferenceName_TBGRI}.}
  
  \item{Status_proposed}{proposed taxonomic status, either 'Accepted', 'Synonym', 'Unresolved', 'IncompleteName' or 'SeveralHomonyms'. 'IncompleteName' indicates incomplete taxon names provided by the user. 'SeveralHomonyms' indicates that there are several homonyms in TPL.}
  
  \item{ReferenceName_proposed}{the taxon name proposed.}
  
  \item{ReferenceAuthorithy_proposed}{authority of \code{ReferenceName_proposed}.}
  
  \item{Infrataxon}{indicates if a taxon name given by the user includes  one or several infra-taxa on the WG according to TBGRI. If it not the case for any of the taxa, the field is absent.}
  
  \item{Family_APGIII}{family name according to APGIII for the reference taxon proposed.}
  
  \item{Genus}{genus of reference taxon proposed.}
  
  \item{Species}{specific epithet of reference taxon proposed.}
  
  \item{InfrataxonRank}{infraspecific rank of reference taxon proposed. If there is no infrataxon in the result, this field is absent.}
  
  \item{InfrataxonName}{infraspecific epithet of reference taxon proposed. If there is no infrataxon in the result, this field is absent.}
  
  \item{Origin}{if \code{otherinfo = T}, the biogegraphical origin of the reference taxon proposed. The information is extracted from \code{\link{TreeGhatsData}}, and is either 'exotic', 'endemic of Western Ghats', 'endemic of india' or 'native'. \code{NA} indicates that the reference taxon proposed is absent from \code{\link{TreeGhatsData}}.}
  
  \item{Habit}{if \code{otherinfo = T}, the habit of the reference taxon proposed. It is extracted from from \code{\link{TreeGhatsData}}, and is either 'liana', 'shrub', 'tree', 'unknown', or any combination of the terms. \code{NA} indicates that the reference taxon proposed is absent from \code{\link{TreeGhatsData}}.}
  
  \item{Phenology}{if \code{otherinfo = T}, the leaf phenology of the reference taxon proposed. It is extracted from \code{\link{TreeGhatsData}}, and is either 'evergreen' or 'deciduous'.}
  
  \item{IUCN}{if \code{otherinfo = T}, the IUCN status, either 'CR', 'DD', 'EN', 'EX', 'LC', 'LR/cd', 'LR/lc', 'LR/nt', 'NE', 'VU' for the reference taxon proposed. When the taxon name is not present in \code{\link{TreeGhatsData}} or not evaluated by IUCN, the function returns \code{NA}.}
  
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
 \code{\link{TreeGhatsData}}, the database used by the function.
 \code{\link[Taxonstand]{TPLck}}
 \code{\link{TPLck2}}
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

\note{
The function relies on the database \code{\link{TreeGhatsData}} included in the package, which will be updated regularly.
The running time of the function depends on the length of the taxon names vector.
An internet connection is needed to check names in TPL.
\code{\link{TPLck2}} is a modified version of \code{\link[Taxonstand]{TPLck}}, which allows handling homonyms.}