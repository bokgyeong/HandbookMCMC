rm(list=ls())
library(dplyr); library(ggplot2); library(egg); library(gridExtra); library(grid); library(scales)
library(batchmeans)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
squish_trans <- function(from, to, factor) {
  trans <- function(x) {
    if (any(is.na(x))) return(x)
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    return(x)
  }
  inv <- function(x) {
    if (any(is.na(x))) return(x)
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    return(x)
  }
  # return the transformation
  return(trans_new("squished", trans, inv))
}


al = 0.01
n = 200000
p = 2



# =============================================================================-
# DMH ----
# =============================================================================-
diagdmh = data.frame()
m = c(1:5)

for(i in m){
  
  ### ACD
  acddmhRep = timeacddmhRep = numeric(30)
  for(runID in 1:length(acddmhRep)){
    load(paste0('acdRep/dmh/m', i, '/acd', runID, '.RData'))
    acddmhRep[runID] = acddmh
    timeacddmhRep[runID] = timeacddmh
  }
  df = data.frame(Method = 'ACD', Algorithm = 'DMH', m = as.factor(i),
                  mean = mean(acddmhRep), lb = quantile(acddmhRep, 0.025), ub = quantile(acddmhRep, 0.975),
                  meantime = mean(timeacddmhRep), lbtime = quantile(timeacddmhRep, 0.025), ubtime = quantile(timeacddmhRep, 0.975),
                  meanBelow = ifelse(mean(acddmhRep) < qchisq(1-al, p*(p+1)/2), TRUE, FALSE),
                  ubBelow = ifelse(quantile(acddmhRep, 0.975) < qchisq(1-al, p*(p+1)/2), TRUE, FALSE))
  diagdmh = rbind(diagdmh, df)

  ### Computation time
  load(paste0('postSamp/dmh', i, '.RData'))
  df = data.frame(Method = 'Computation time', Algorithm = 'DMH', m = as.factor(i),
                  mean = rtime/nrow(postSamp) * n / 60,
                  lb = rtime/nrow(postSamp) * n / 60,
                  ub = rtime/nrow(postSamp) * n / 60,
                  meantime = NA, lbtime = NA, ubtime = NA, meanBelow = NA, ubBelow = NA)
  diagdmh = rbind(diagdmh, df)
}

diagdmh


# =============================================================================-
# ALR ----
# =============================================================================-
diagalr = data.frame()
# m = c(50, 200, 100, 300, 400, 500)
m = c(50, 100, 200, 300, 400)
cycle = 10

for(i in m){
  
  ### ACD
  acdalrRep = timeacdalrRep = numeric(30)
  for(runID in 1:length(acdalrRep)){
    load(paste0('acdRep/alr/d', i, 'm', cycle, '/acd', runID, '.RData'))
    acdalrRep[runID] = acdalr
    timeacdalrRep[runID] = timeacdalr
  }
  df = data.frame(Method = 'ACD', Algorithm = 'ALR', m = as.factor(i),
                  mean = mean(acdalrRep), lb = quantile(acdalrRep, 0.025), ub = quantile(acdalrRep, 0.975),
                  meantime = mean(timeacdalrRep), lbtime = quantile(timeacdalrRep, 0.025), ubtime = quantile(timeacdalrRep, 0.975),
                  meanBelow = ifelse(mean(acdalrRep) < qchisq(1-al, p*(p+1)/2), TRUE, FALSE),
                  ubBelow = ifelse(quantile(acdalrRep, 0.975) < qchisq(1-al, p*(p+1)/2), TRUE, FALSE))
  diagalr = rbind(diagalr, df)

  ### Computation time
  load(paste0('postSamp/alr', i, 'm', cycle, '.RData'))
  df = data.frame(Method = 'Computation time', Algorithm = 'ALR', m = as.factor(i),
                  mean = rtime/nrow(postSamp) * n / 60,
                  lb = rtime/nrow(postSamp) * n / 60,
                  ub = rtime/nrow(postSamp) * n / 60,
                  meantime = NA, lbtime = NA, ubtime = NA, meanBelow = NA, ubBelow = NA)
  diagalr = rbind(diagalr, df)
}

diagalr$m = factor(diagalr$m, levels = paste(m))


# =============================================================================-
# ABCMCMC ----
# =============================================================================-
diagabcmcmc = data.frame()
# eps_threshold = c(1, 1.4, 1.6, 1.7, 1.8, 2, 3, 4)
eps_threshold = c(1, 1.4, 1.6, 1.7, 1.8, 1.9, 2, 3)
# eps_threshold = c(1, 1.4, 1.6, seq(1.7, 1.8, by = 0.01), 2, 3)
cycle = 20

for(i in eps_threshold){
  
  ### ACD
  acdabcmcmcRep = timeacdabcmcmcRep = numeric(30)
  for(runID in 1:length(acdabcmcmcRep)){
    load(paste0('acdRep/abcmcmc/E', i, 'm', cycle, '/acd', runID, '.RData'))
    acdabcmcmcRep[runID] = acdabcmcmc
    timeacdabcmcmcRep[runID] = timeacdabcmcmc
  }
  df = data.frame(Method = 'ACD', Algorithm = 'ABC-MCMC', m = as.factor(i),
                  mean = mean(acdabcmcmcRep), lb = quantile(acdabcmcmcRep, 0.025), ub = quantile(acdabcmcmcRep, 0.975),
                  meantime = mean(timeacdabcmcmcRep), lbtime = quantile(timeacdabcmcmcRep, 0.025), ubtime = quantile(timeacdabcmcmcRep, 0.975),
                  meanBelow = ifelse(mean(acdabcmcmcRep) < qchisq(1-al, p*(p+1)/2), TRUE, FALSE),
                  ubBelow = ifelse(quantile(acdabcmcmcRep, 0.975) < qchisq(1-al, p*(p+1)/2), TRUE, FALSE))
  diagabcmcmc = rbind(diagabcmcmc, df)

  ### Computation time
  load(paste0('postSamp/abcmcmcE', i, 'm', cycle, '.RData'))
  df = data.frame(Method = 'Computation time', Algorithm = 'ABC-MCMC', m = as.factor(i),
                  mean = rtime/nrow(postSamp$theta) * n / 60,
                  lb = rtime/nrow(postSamp$theta) * n / 60,
                  ub = rtime/nrow(postSamp$theta) * n / 60,
                  meantime = NA, lbtime = NA, ubtime = NA, meanBelow = NA, ubBelow = NA)
  diagabcmcmc = rbind(diagabcmcmc, df)
}

# diagabcmcmc$m = factor(diagabcmcmc$m, levels = paste0(rev(eps_threshold)))
diagabcmcmc$m = factor(diagabcmcmc$m, levels = paste0(eps_threshold))




# =============================================================================-
# BSL ----
# =============================================================================-
diagbsl = data.frame()
tunpar = c(10, 30, 50, 100, 300)

for(i in tunpar){
  
  ### ACD
  acdbslRep = timeacdbslRep = numeric(30)
  for(runID in 1:length(acdbslRep)){
    load(paste0('acdRep/bsl/N', i, '/acd', runID, '.RData'))
    acdbslRep[runID] = acdbsl
    timeacdbslRep[runID] = timeacdbsl
  }
  df = data.frame(Method = 'ACD', Algorithm = 'BSL', m = as.factor(i),
                  mean = mean(acdbslRep), lb = quantile(acdbslRep, 0.025), ub = quantile(acdbslRep, 0.975),
                  meantime = mean(timeacdbslRep), lbtime = quantile(timeacdbslRep, 0.025), ubtime = quantile(timeacdbslRep, 0.975),
                  meanBelow = ifelse(mean(acdbslRep) < qchisq(1-al, p*(p+1)/2), TRUE, FALSE),
                  ubBelow = ifelse(quantile(acdbslRep, 0.975) < qchisq(1-al, p*(p+1)/2), TRUE, FALSE))
  diagbsl = rbind(diagbsl, df)
  
  ### Computation time
  load(paste0('postSamp/bsl', i, '.RData'))
  df = data.frame(Method = 'Computation time', Algorithm = 'BSL', m = as.factor(i),
                  mean = rtime/nrow(bslSamp) * n / 60,
                  lb = rtime/nrow(bslSamp) * n / 60,
                  ub = rtime/nrow(bslSamp) * n / 60,
                  meantime = NA, lbtime = NA, ubtime = NA, meanBelow = NA, ubBelow = NA)
  diagbsl = rbind(diagbsl, df)
}

diagbsl$m = factor(diagbsl$m, levels = paste0(tunpar))



# =============================================================================-
# Figures for ACD and time ----
# =============================================================================-

range.acd = c(
  0, 
  max(
    diagalr %>% filter(Method == 'ACD') %>% select(ub) %>% max,
    diagdmh %>% filter(Method == 'ACD') %>% select(ub) %>% max,
    diagabcmcmc %>% filter(Method == 'ACD') %>% select(ub) %>% max,
    diagbsl %>% filter(Method == 'ACD') %>% select(ub) %>% max
  ))

range.cost = range(
  c(
    # diagalr %>% filter(Method == 'Computation time') %>% select(mean) %>% range,
    diagdmh %>% filter(Method == 'Computation time') %>% select(mean) %>% range,
    diagabcmcmc %>% filter(Method == 'Computation time') %>% select(mean) %>% range,
    diagbsl %>% filter(Method == 'Computation time') %>% select(mean) %>% range
  )
)


plotacdalr = diagalr %>% 
  filter(Method == 'ACD') %>% 
  ggplot(aes(x = m, y = mean, group = 1, color = meanBelow)) +
  geom_point(aes(shape = meanBelow), size = 2) +
  geom_errorbar(aes(ymin = lb, ymax = ub), color = 'gray30', width = 0) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 15)) +
  scale_color_manual(values = c('FALSE' = '#F8766D', 'TRUE' = '#00BFC4')) +
  geom_hline(yintercept = qchisq(1-al, p*(p+1)/2), linetype = 2) +
  # scale_y_log10() +
  # coord_cartesian(ylim = range.acd) +
  labs(x = 'd', y = 'ACD', title = '(a) ALR') +
  theme_bw() +
  theme(legend.position = 'none')


plotacddmh = diagdmh %>% 
  filter(Method == 'ACD') %>% 
  ggplot(aes(x = m, y = mean, group = 1, color = meanBelow)) +
  geom_point(aes(shape = meanBelow), size = 2) +
  geom_errorbar(aes(ymin = lb, ymax = ub), color = 'gray30', width = 0) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 15)) +
  scale_color_manual(values = c('FALSE' = '#F8766D', 'TRUE' = '#00BFC4')) +
  geom_hline(yintercept = qchisq(1-al, p*(p+1)/2), linetype = 2) +
  # scale_y_log10() +
  # coord_cartesian(ylim = range.acd) +
  labs(x = 'm', y = 'ACD', title = '(b) DMH') +
  theme_bw() +
  theme(legend.position = 'none')


plotacdabcmcmc = diagabcmcmc %>% 
  filter(Method == 'ACD') %>% 
  ggplot(aes(x = m, y = mean, group = 1, color = meanBelow)) +
  geom_point(aes(shape = meanBelow), size = 2) +
  geom_errorbar(aes(ymin = lb, ymax = ub), color = 'gray30', width = 0) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 15)) +
  scale_color_manual(values = c('FALSE' = '#F8766D', 'TRUE' = '#00BFC4')) +
  geom_hline(yintercept = qchisq(1-al, p*(p+1)/2), linetype = 2) +
  # scale_y_log10() +
  # coord_cartesian(ylim = range.acd) +
  labs(x = expression(epsilon), y = 'ACD', title = '(c) ABC-MCMC') +
  theme_bw() +
  theme(legend.position = 'none')

plotacdbsl = diagbsl %>% 
  filter(Method == 'ACD') %>% 
  ggplot(aes(x = m, y = mean, group = 1, color = meanBelow)) +
  geom_point(aes(shape = meanBelow), size = 2) +
  geom_errorbar(aes(ymin = lb, ymax = ub), color = 'gray30', width = 0) +
  scale_shape_manual(values = c('FALSE' = 17, 'TRUE' = 15)) +
  scale_color_manual(values = c('FALSE' = '#F8766D', 'TRUE' = '#00BFC4')) +
  geom_hline(yintercept = qchisq(1-al, p*(p+1)/2), linetype = 2) +
  # scale_y_log10() +
  # coord_cartesian(ylim = range.acd) +
  labs(x = 'N', y = 'ACD', title = '(d) BSL') +
  theme_bw() +
  theme(legend.position = 'none')


plotcostalr = diagalr %>% 
  filter(Method == 'Computation time') %>% 
  ggplot(aes(m, mean)) +
  geom_point(size = 2, shape = 16) +
  theme_bw() +
  labs(x = 'd', y = 'Computation (min)')

plotcostdmh = diagdmh %>% 
  filter(Method == 'Computation time') %>% 
  ggplot(aes(m, mean)) +
  geom_point(size = 2, shape = 16) +
  coord_cartesian(ylim = range.cost) +
  theme_bw() +
  labs(x = 'm', y = 'Computation (min)')

plotcostabcmcmc = diagabcmcmc %>% 
  filter(Method == 'Computation time') %>% 
  ggplot(aes(m, mean)) +
  geom_point(size = 2, shape = 16) +
  coord_cartesian(ylim = range.cost) +
  theme_bw() +
  labs(x = expression(epsilon), y = 'Computation (min)')

plotcostbsl = diagbsl %>% 
  filter(Method == 'Computation time') %>% 
  ggplot(aes(m, mean)) +
  geom_point(size = 2, shape = 16) +
  coord_cartesian(ylim = range.cost) +
  theme_bw() +
  labs(x = 'N', y = 'Computation (min)')




# -----------------------------------------------------------------------------=
# Combine plots ----
# -----------------------------------------------------------------------------=

plotDiag = ggarrange(plotacdalr, plotacddmh, plotacdabcmcmc, # plotacdbsl,
                     # plotcostalr, plotcostdmh, plotcostabcmcmc,
                     byrow = T, nrow = 1)

ggsave(plot = plotDiag, filename = 'fig/ergmACD.pdf',
       width = 2.4*3, height = 2.5)



# =============================================================================-
# Summary statistics ----
# =============================================================================-
sampleAll = data.frame()


### ALR---gold standard
cycle = 10
for(i in c(200)){
  load(paste0('postSamp/alr', i, 'm', cycle, '.RData'))
  df = data.frame(Sample = 'ALR', m = as.factor(i), parameter = postSamp[(burn-1)+1:n,], time = rtime / nrow(postSamp) * n)
  sampleAll = rbind(sampleAll, df)
}


### DMH
for(i in c(3)){
  load(paste0('postSamp/dmh', i, '.RData'))
  df = data.frame(Sample = 'DMH', m = as.factor(i), parameter = postSamp[burn+1:n,], time = rtime / nrow(postSamp) * n)
  sampleAll = rbind(sampleAll, df)
}

### ABCMCMC
cycle = 20
for(i in c(1.9)){
  load(paste0('postSamp/abcmcmcE', i, 'm', cycle, '.RData'))
  df = data.frame(Sample = 'ABC-MCMC', m = as.factor(i), parameter = postSamp$theta[burn+1:n,], time = rtime / nrow(postSamp$theta) * n)
  sampleAll = rbind(sampleAll, df)
}


### VB
load("postSamp/vi.RData")
df = data.frame(Sample = 'VB', m = NA, parameter = postSamp, time = NA)
sampleAll = rbind(sampleAll, df)


# =============================================================================-
# Tail probabilities ----
# =============================================================================-

## theta1 ----
cut.1 = sampleAll %>% filter(Sample == 'ALR') %>% 
  summarise(Cut = quantile(parameter.1, c(0.05, 0.45, 0.5, 0.55, 0.95))) %>% unlist()

sumstatsAll.1 = sampleAll %>%
  group_by(Sample, m) %>%
  summarise(Mean = round(mean(parameter.1), 2),
            Median = round(median(parameter.1), 2),
            SD = round(sd(parameter.1), 2),
            Prob1 = round(mean(parameter.1 <= cut.1[1]), 2),
            Prob2 = round(mean(parameter.1 >= cut.1[2] & parameter.1 <= cut.1[3]), 2),
            Prob3 = round(mean(parameter.1 >= cut.1[3] & parameter.1 <= cut.1[4]), 2),
            Prob4 = round(mean(parameter.1 >= cut.1[5]), 2),
            Time = format(round(mean(time)/60/60, 2), nsmall = 2))

sumstatsAll.1 %>% select(Sample, m , Median, Prob1, Prob4, Time)


## theta2 ----
cut.2 = sampleAll %>% filter(Sample == 'ALR') %>% 
  summarise(Cut = quantile(parameter.2, c(0.05, 0.45, 0.5, 0.55, 0.95))) %>% unlist()

sumstatsAll.2 = sampleAll %>%
  group_by(Sample, m) %>%
  summarise(Mean = round(mean(parameter.2), 2),
            Median = round(median(parameter.2), 2),
            SD = round(sd(parameter.2), 2),
            Prob1 = round(mean(parameter.2 <= cut.2[1]), 2),
            Prob2 = round(mean(parameter.2 >= cut.2[2] & parameter.2 <= cut.2[3]), 2),
            Prob3 = round(mean(parameter.2 >= cut.2[3] & parameter.2 <= cut.2[4]), 2),
            Prob4 = round(mean(parameter.2 >= cut.2[5]), 2),
            Time = format(round(mean(time)/60/60, 2), nsmall = 2))

sumstatsAll.2 %>% select(Sample, m , Median, Prob1, Prob4, Time)



# =============================================================================-
# Density plots ----
# =============================================================================-

sampleAll$Sample = factor(
  sampleAll$Sample, 
  levels = c('ALR', 'DMH', 'ABC-MCMC', 'VB'),
  labels = c("'ALR with d = 200'", "'DMH with m = 3'", "'ABC-MCMC with '*epsilon*' = 1.9'", "'VB'")
)



## 2d density ----
p.den2d = sampleAll %>% 
  ggplot(aes(x = parameter.1, y = parameter.2)) +
  geom_density_2d(color = 'gray20') +
  # geom_point() +
  facet_wrap(~Sample, nrow = 1, labeller = label_parsed) +
  labs(x = expression(theta[1]), y = expression(theta[2])) +
  theme_bw()



ggsave(plot = p.den2d, filename = 'fig/ergmDen.pdf',
       width = 2.2*4, height = 2.5)



## theta1 ----
sampleAll %>% 
  ggplot(aes(x = parameter.1, color = Sample)) +
  geom_density()


## theta1 ----
sampleAll %>% 
  ggplot(aes(x = parameter.2, color = Sample)) +
  geom_density()




