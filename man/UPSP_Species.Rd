\name{UPSP_Species}
\alias{UPSP_Species}
\docType{data}
\title{
List of tree species in the Uppangala rainforest plot
}
\description{
The data frame includes a list of tree species names reported in the Uppangala Permanent Sampling Plot (UPSP), Karnataka.
}
\usage{data("UPSP_Species")}
\value{
  A data frame including the following columns:
  \item{SpCode}{
  a four-letter code including the two first letters of genus and two first letter of species names.
  }
  \item{LatinName}{
  genus names (capitalized).
  }
  \item{species}{
  species names.
  }
  \item{family}{
  family names.
  }
  \item{VoucherRef}{
  the reference code of samples collected and conserved in the Herbarium of French Institute of Pondicherry (HFIP).
  }
  \item{PotStrata}{
  the potential vertical stratum of tree species.
  }
  \item{binome}{
  binome name (genus plus species).
  }   
}
\details{
The dataset includes all the tree species reaching a diameter at breadth height above 10cm in the Uppangala Permanent Forest Plot(UPSP). The context of the forest plot is presented in details in Pélissier et al. (2011).
}
\source{
The basic dataset is available in \href{https://figshare.com/collections/Tree_demography_in_an_undisturbed_Dipterocarp_permanent_sample_plot_at_Uppangala_Western_Ghats_of_India/3304026}{Pélissier et al. (2011)}.
We only included the additional column \code{binome}.
}
\references{
Pélissier R., Pascal J.-P., Ayyappan N., Ramesh B. R., Aravajy S. & Ramalingam (2011). Twenty years tree demography in an undisturbed Dipterocarp permanent sample plot at Uppangala, Western Ghats of India. Ecology. 92: 1376.
}
\examples{
data(UPSP_Species)
str(UPSP_Species)
}
\keyword{datasets}
