rm(list=ls())
require(potts); require(tidyverse)


set.seed(1012)
#==============================================================================-
# simulate data ----
#==============================================================================-
ncolor = as.integer(4)
tr.beta = log(1 + sqrt(ncolor))
theta = c(rep(0, ncolor), tr.beta)
nrow = 30
ncol = 30
x = matrix(sample(ncolor, nrow*ncol, replace=TRUE), nrow=nrow, ncol=ncol)
foo = packPotts(x, ncolor)
out = potts(foo, theta, nbatch = 100000)
x = unpackPotts(out$final)
stat = calc_t_full(x, ncolor)[ncolor+1]


data.plot = data.frame(x = rep(1:30, times = 30), y = rep(1:30, each = 30), Color = as.character(as.vector(x)))
potts.plot = data.plot %>% ggplot(aes(x = x, y = y, fill = Color)) +
  geom_tile(color = 'white') +
  scale_fill_manual(values = c(
    '1' = 'gray1',
    '2' = 'gray35',
    '3' = 'gray65',
    '4' = 'gray90'
  )) +
  theme_bw()
potts.plot


save(tr.beta, foo, x, stat, ncolor, file = "data/sim.RData")
