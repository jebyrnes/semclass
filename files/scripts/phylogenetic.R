### STRUCTURAL EQUATION MODELING: Phylogenetic Correction ###

### Author: Jon Lefcheck
### Last updated: 25 October 2017

#install.packages("caper")
library(caper)
library(ape)
library(nlme)
library(piecewiseSEM)

# Read in data
duffy <- read.csv("../data/duffy.csv")

# Read in tree
tree <- read.tree("../data/synalpheus_tree.txt")

plot(tree)

# Merge Synalpheus data and tree
duffy.merge <- comparative.data(tree, duffy, names.col = "Species")

# Run model without acknowledging phylogenetic non-independence
duffy.sem <- psem(
  
  lm(Host.Range ~ Eusociality.index + Total.Mass.Female, duffy.merge$data),
  
  lm(Rubble.Abundance ~ Eusociality.index + Host.Range, duffy.merge$data)
  
)

summary(duffy.sem)

# Repeat but include phylogenetic information
duffy.sem2 <- psem(
  
  pgls(Host.Range ~ Eusociality.index + Total.Mass.Female, duffy.merge),
  
  pgls(Rubble.Abundance ~ Eusociality.index + Host.Range, duffy.merge),
  
  duffy.merge
  
)

# Check fit
summary(duffy.sem2)