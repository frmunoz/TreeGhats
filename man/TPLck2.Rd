\name{TPLck2}
\alias{TPLck2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Modified version of \code{TPLck} indicating the number of homonyms of a plant name 
}
\description{
A modified version of \code{\link[Taxonstand]{TPLck}} function from the package \code{Taxonstand} in order to connect to The Plant list (TPL),and check for the validity of a plant name. The function provides the same informations than \code{TPLck}. However we add a new value indicating the number of homonyms. Indeed, even the presence of several homonyms is indicating in \code{Warning message}, by  \code{the plant name has more than one valid synonym},  the function \code{TPLck} returns the values for the first homonyms. 
}
\usage{
TPLck2(sp, corr = TRUE, diffchar = 2, max.distance = 1, infra = TRUE, abbrev = TRUE,
version = "1.1", encoding = "UTF-8")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp}{A character specifying the genus and specific epithet and, additionally, the infraspecific epithet.
}
  \item{corr}{Logical. If 'TRUE' (default), then removal of spelling errors is performed (only) on the specific and infraspecific epithets prior to taxonomic standardization.
}
  \item{diffchar}{
%%     ~~Describe \code{diffchar} here~~
}
  \item{max.distance}{
%%     ~~Describe \code{max.distance} here~~
}
  \item{infra}{
%%     ~~Describe \code{infra} here~~
}
  \item{abbrev}{
%%     ~~Describe \code{abbrev} here~~
}
  \item{version}{
%%     ~~Describe \code{version} here~~
}
  \item{encoding}{
%%     ~~Describe \code{encoding} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
To be completed
}
\value{
The function return an object of class \code{data.frame} with the following components:
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
To be completed
}
\author{
V. Raevel adapted from Luis Cayuela & Jari Oksanen
}
\note{
%%  ~~further notes~~
To be completed
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 See Also  \code{\link[Taxonstand]{TPLck}}
}
\examples{
example <- c("Dipterocarpus indicus")
result <- TPLck2(example)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
