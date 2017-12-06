###
### Calculate J and y
###

rm(list = ls())

capture <- read.delim("data/NCEI_NatalityCaptureData.csv", sep = ",", header = TRUE)
resight <- read.delim("data/NCEI_NatalityResightData.csv", sep = ",", header = TRUE)

library(lubridate)
resight$sitedate = mdy(resight$sitedate)
resight$sightyear = year(resight$sitedate)
years = sort(unique(resight$sightyear))
years = years[years > 1990]
## restrict to 1991

brands = unique(capture$brand[which(capture$cohort == 1991)])

y = t(sapply(brands, function(brand) {
  sapply(years, function(k) {
    as.numeric(sum(resight[resight$sightyear == k & resight$brand == brand, ]$withpup
        %in% c("Y", "y")) > 0)
  })
}))

z = t(sapply(brands, function(brand) {
  sapply(years, function(k) {
    as.numeric(sum(resight$sightyear == k & resight$brand == brand) > 0)
  })
}))

## view output

source("./hw7_mcmc_lessMH.R")

n.mcmc <- 5e4

out = hw7.mcmc2(y, z, n.mcmc = n.mcmc, r.tune = .5, d.tune = 3,
                alpha.p = 1, beta.p = 3, alpha.psi = 2, beta.psi = 2)

#out = readRDS("./data/results.2.rds")
n.burn <- .5*n.mcmc

out$d.mh.prop[1] / out$d.mh.prop[2]
out$r.mh.prop[1] / out$r.mh.prop[2]

plot(out$p.save[(n.burn+1):n.mcmc], type = "l")
plot(out$psi.save[(n.burn+1):n.mcmc], type = "l")
plot(out$r.save[100,(n.burn+1):n.mcmc], type = "l")
plot(out$d.save[4,(n.burn+1):n.mcmc], type = "l", main = "last seen year 0")
plot(out$d.save[2,(n.burn+1):n.mcmc], type = "l", main = 'last seen year 17')
plot(out$d.save[7,(n.burn+1):n.mcmc], type = "l", main = 'last seen year 6')
plot(out$d.save[38,(n.burn+1):n.mcmc], type = "l", main = 'last seen year 3')
plot(out$d.save[114,(n.burn+1):n.mcmc], type = "l", main = 'last seen year 2')
plot(out$d.save[100,(n.burn+1):n.mcmc], type = "l", main = "last seen year 1")

library(ggplot2)

plot.dat = data.frame(individual = sort(rep(1:nrow(y), n.mcmc - n.burn)),
                      r = as.vector(t(out$r.save)))
ggplot(plot.dat, aes(x=factor(individual), y=r)) + 
  geom_violin()


# still takes awhile to burn in - may be able to change that if we mess around with 
# the tuning parameters.

