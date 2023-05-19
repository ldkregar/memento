
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Memento

<!-- badges: start -->
<!-- badges: end -->

Memento is a suite of functions to find change points in phenotype
states on a phylogeny.

Suppose you know the number of methylated and unmethylated read counts
at a CpG site across all tips of a phylogeny. The question Memento
addresses is whether the tip states are inherited from the root and
whether any of the internal nodes changed its state, which was then
inherited by its descendants.

Memento will take two vectors of counts as well as the phylogeny of tips
and return the most likely root state, which branches changed and what
type of change occurred on the branches. Memento relies on naive Bayes
as it generates all possible tip state vectors for up to m number of
independent changes on a phylogeny and then compares the observed data
to each possible expected vector of tip states. The vector with highest
similarity is selected and since we know how that vector was generated,
we can infer the root state and which branches changed.

Something about number of combinations

If the phylogeny is very big (say, contains more than 20 tips), Memento
allows the user to split the tree, infer most likely state changes on
each subtree and then combine the partial solutions to infer the most
likely solution for the full tree.

Memento also works on phased data and can infer whether a change
happened on maternal and/or paternal chromosome. This is referred to as
the 4 state model.

Since Memento computes probabilities for a large set of possible tip
states, it can also return second, third … best solution.

Memento also contains several visualisation functions.

Priors assume x.

## Installation

You can install the development version of Memento from
[GitHub](https://github.com/) with:

NB: This library is under active development so some unexpected
behaviour is to be expected. Drop me an email if you require help, I’m
at LD18 at sanger dot ac dot uk.

``` r
# install.packages("devtools")
devtools::install_github("ldkregar/memento")
```
