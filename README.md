# TreeGhats
Taxonomic and ecological database of trees of Western Ghats
[![Build Status](https://travis-ci.org/frmunoz/TreeGhats.svg?branch=master)](https://travis-ci.org/frmunoz/TreeGhats)

In this project, we aim to develop a comprehensive database on the taxon names and a tool to recognize and resolve taxon names for
the Western Ghats region in R software environment.
The database includes a comprehensive taxon list of shrubs and trees (those species with the girth of at least 10 cm or more at
1.3 m from the ground) of the Western Ghats with an up-to-date taxonomic reference. The tool contains all the taxon names that
have been commonly used in forest ecology research over the last decades in South India. For each taxon name, a taxonomic
status is issued with reference to The Flowering plants of the Western Ghats (TBGRI, Nayar, Beegam, and Sibi. 2014) and to
The Plant List (http://www.theplantlist.org/). Priority is given to the recent local flora in case of inconsistent taxonomic
status, but the database offers flexibility in the choice of the reference. When a taxon has a synonymous status, the correspondence
to the accepted name is given according to the selected reference. APGIII family names and standard form of
author’s names based on The International Plant Names Index are included.

A tool is also provided to analyze the taxonomic status of many species quickly and efficiently. It is available for use in the free
R software environment. The user provides list of species to the function, and the output contains the detailed taxonomic diagnostic
of species of the list, along with basic ecological information (habit, phenology…). The user can then easily analyze many lists
with heterogeneous taxonomic references, in order to perform comparative analyses of their patterns in terms of taxonomic and 
ecological diversity.

The database and the tool also come with a recent dated phylogeny resolved at genus level. It is provided in standard Newick format
that can be analyzed with many standard software. The user can then analyze patterns of phylogenetic diversity and test
evolutionary hypotheses on the origin of the patterns. 

Therefore, the overall package is intended to meet the need of modern scientists who wish to decipher the current and past drivers of
biodiversity from multiple inventory data over the whole Western Ghats. The current version of the database and tool will evolve
constantly in future with the acquired knowledge of new species and as the state of art changes. We hope it will thus foster extensive
and large-scale analyses for better understanding and preserving the biological heritage of WG forests.
