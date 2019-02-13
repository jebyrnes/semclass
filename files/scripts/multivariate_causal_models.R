#Libraries
library(dagitty)

#Build a graph
g <- dagitty("dag{
  x1 -> y1
  x1 -> y2
  x1 -> y3
  y3 -> y4
  y2 -> y4
  y1 -> y4
  y2 -> y1
             }")

#Show us the graph
plot(graphLayout(g))

#How do we shut the back door?
adjustmentSets(g, exposure = "y1", outcome = "y4")


# My Model
forest_mod <-  dagitty("dag{
  	waves -> kelp -> algae  
	  algae -> inverts
  	waves -> algae
}")


plot(graphLayout(forest_mod))

adjustmentSets(forest_mod, 
               exposure = "kelp", 
               outcome = "inverts")

#conditional independence relationships
impliedConditionalIndependencies(forest_mod)
