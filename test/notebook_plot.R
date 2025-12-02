############################################################
## In this notebook, we provide template code to represent graphically
## the PCs obtained with the msPCA package.
## The code requires the dplyr, ggplot2, and ggforce packages.
############################################################
library(datasets)
df <- datasets::mtcars
S <- cor(df)

library(msPCA)
mspca_results <- msPCA::mspca(S, 3, c(6,6,6), verbose=TRUE)

library(dplyr)
library(ggplot2)
library(ggforce)

v <- mspca_results$x_best
row.names(v) <- row.names(S)

v

#Pick the two (original) dimensions on which to represent the PCs
dim1 <- row.names(v)[4]
dim2 <- row.names(v)[5]

data.frame(t(v)) %>%
  mutate(pc_id = paste("PC ", row_number(), sep=""), origin = 0) %>%
  filter( abs(!!sym(dim1)) + abs(!!sym(dim2)) > 0.01) %>%
  ggplot() +
  aes(x=!!sym(dim1), y=!!sym(dim2), label = pc_id) +
  geom_point(shape=4) + geom_text(hjust=-0.5,vjust=-0.5,size=6) +
  geom_segment(aes(x = origin, xend = !!sym(dim1), y = origin, yend = !!sym(dim2)),
               arrow = arrow(angle = 20, length = unit(0.15, "inches"))) +
  geom_point(aes(x=0, y=0), inherit.aes=FALSE) +
  geom_circle(aes(x0=0, y0=0, r=1), inherit.aes=FALSE) + coord_fixed()
