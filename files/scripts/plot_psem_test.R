
library(piecewiseSEM)
source("./plot_psem.R")


data(keeley)

mod <- psem(
  lm(rich ~ cover, data=keeley),
  lm(cover ~ firesev, data=keeley),
  lm(firesev ~ age, data=keeley),
  data = keeley
  
)

coefs(mod)

#
keeley_graph <- plot.psem(mod, return=TRUE)
#
plot.psem(mod)
#
plot.psem(mod, node_attrs = list(
  shape = "rectangle", color = "black",
  fillcolor = "orange", x = 3, y=1:4))
