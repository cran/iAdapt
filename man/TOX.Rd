\name{TOX}
\alias{TOX}
\docType{data}
\title{
Sample array of toxicity probabilities for 6 doses. Taken from Du et al.
}
\description{
This is a sample array of toxicity probabilities to be used for trial simulations with nTTP as the toxicity endpoint. In this example, we have 3 toxicity types, 6 test doses, and 5 AE grades (0-4). See the nTTP_simulation_example vignette for more details.
}
\usage{data("TOX")}
\format{
  The format is:
 num [1:6, 1:5, 1:3] 0.791 0.738 0.685 0.662 0.605 0.39 0.172 0.195 0.19 0.2 ...
}
\source{
https://pubmed.ncbi.nlm.nih.gov/30403559/
}
\examples{
data(TOX)
TOX
}
\keyword{datasets}
