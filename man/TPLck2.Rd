\name{TPLck2}
\alias{TPLck2}
\title{Modified version of \code{TPLck} indicating the number of homonyms of a taxon name 
}
\description{
A modified version of the \code{\link[Taxonstand]{TPLck}} function from package \code{Taxonstand}. The function basically provides the same information than \code{TPLck}, but further keeps track of homonyms. Although a warning message indicates in \code{TPLck} when there are several homonyms, the function returns the first homonym in alphabetic order only. 
\code{TPLck2} returns the number of homonyms for a given taxon name.
}
\usage{
TPLck2(sp, corr = TRUE, diffchar = 2, max.distance = 1, infra = TRUE, abbrev = TRUE,
version = "1.1", encoding = "UTF-8")
}
\arguments{
\item{sp}{a character vector specifying the input taxon, i.e. genus and specific epithet and, potentially, author name and infraspecific abbreviation and epithet.}
\item{corr}{logical. If TRUE (default), spelling errors are corrected (only) in the specific and infraspecific epithets prior to taxonomic standardization.}
\item{diffchar}{a number indicating the maximum difference between the number of characters in corrected and original taxon names. Not used if corr = FALSE.}
\item{max.distance}{a number indicating the maximum distance allowed for a match in agrep when performing corrections of spelling errors in specific epithets. Not used if corr = FALSE.}
\item{infra}{logical. If TRUE (default), infraspecific epithets are used to match taxon names in TPL.}
\item{abbrev}{logical. If 'TRUE' (default), abbreviations (aff., cf., subsp., var.) and their variants are removed prior to taxonomic standardization.}
\item{version}{a character vector indicating whether to connect to the newest version of TPL (1.1) or to the older one (1.0). Default to "1.1".}
\item{encoding}{encoding to be assumed for input strings from TPL website; defaults to "UTF-8" (see \code{\link{read.csv}} and \code{\link{Encoding}}).}
}
\details{
\code{TPLck2} overcomes two issues in \code{TPLck}. First, valid homonyms (i.e., a name for a taxon that is identical in spelling to another such name, but belongs to a different taxon) cannot be identified. 
Second, if the infraspecific epithet does not match any of the infraspecific epithets considered in TPL, then the first accepted name with no infraspecific epithet is selected as the best match. If all names provided by TPL have infraspecific epithets and are all synonyms (e.g. Pottia starkeana), then there is no best match and the output will match the original name.
}
\value{
The function return an object of class \code{data.frame} including the same elements than those provided by \code{\link[Taxonstand]{TPLck}}, plus the element \code{$Homonym} indicating the number of homonyms if any. 
The value of \code{$Homonym} is \code{NA} if there is no homonym.
}
\author{
V. Raevel, adapted from Luis Cayuela & Jari Oksanen
}
\seealso{
 \code{\link[Taxonstand]{TPLck}}
}
\examples{
\dontrun{
Taxonstand::TPLck("Acacia sundra") 
test <- TPLck2("Acacia sundra") 
test
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.