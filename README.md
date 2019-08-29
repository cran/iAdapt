# iAdapt
R package for early phase adaptive two-stage clinical trial design for toxicity and immunologic outcomes in oncology.

## Installation

To install the package, run the following code:


```r
install.packages("devtools")
devtools::install_github("alyssamv/iAdapt")
```


## Overview

This package provides software based on the early phase trial design by [Chiuzan et al. (2018)](https://www.tandfonline.com/doi/abs/10.1080/19466315.2018.1462727). Stage 1 is safety-driven dose-escalation, and Stage 2 employs efficacy-driven randomization while continuing to monitor dose safety.

The design uses a likelihood paradigm, rather than rules. e.g. In Stage 1, when the likelihood ratio for a dose is greater than a prespecified threshold, the dose is considered acceptably safe and subsequent patients are enrolled on the next dose level. Conversely, if the likelihood ratio is less than or equal to the threshold, escalation is stopped. 

One hallmark of this design is its ability to identify the most effective dose in the presence of a non-monotone dose-response curve - a phenomenon common in immunotherapies. Additionally, it follows a frequentist framework, but allows for adaptive design components.

The function of this package is two-fold:

* Produce trial outcomes through simulation for an inputted scenario, and
* Implement the design in a real trial.

Vignettes are provided to walk the user through each function of the package. 

## Background

This design is relevant in the face of a non-monotonous dose-response relationship - a phenomenon that most often seen in immunologic therapies. Additionally, the probabilistic nature (as opposed to rule-based) provides an advantage in identifying the optimal dose to carry forward in development, by allowing more than one dose to be examined for efficacy. 

### Dose-response relationship 

Often, dose-finding designs rely on a monotone dose-response curve, meaning that as dose increases, we expect the drug's effectiveness to increase too. This is a convenient assumption, though not always accurate. Instead, a relationship may exist in which case some dose-escalation designs may falsely move a higher dose forward.