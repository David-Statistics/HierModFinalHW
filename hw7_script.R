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
set.seed(1215)
n.mcmc <- 5e4

out = hw7.mcmc(y, z, n.mcmc = n.mcmc, r.tune = .5, d.tune = 3,
                alpha.p = 1.5, beta.p = 4.5, alpha.psi = 2, beta.psi = 2)

n.burn <- .1*n.mcmc

out$d.mh.prop[1] / out$d.mh.prop[2]
out$r.mh.prop[1] / out$r.mh.prop[2]

plot(out$p.save[(n.burn+1):n.mcmc], type = "l")
plot(out$psi.save[(n.burn+1):n.mcmc], type = "l")
plot(out$r.save[2,(n.burn+1):n.mcmc], type = "l")
plot(out$d.save[4,(n.burn+1):n.mcmc], type = "l", main = "last seen year 0")
plot(out$d.save[2,(n.burn+1):n.mcmc], type = "l", main = 'last seen year 17')
plot(out$d.save[7,(n.burn+1):n.mcmc], type = "l", main = 'last seen year 6')
plot(out$d.save[38,(n.burn+1):n.mcmc], type = "l", main = 'last seen year 3')
plot(out$d.save[114,(n.burn+1):n.mcmc], type = "l", main = 'last seen year 2')
plot(out$d.save[100,(n.burn+1):n.mcmc], type = "l", main = "last seen year 1")

cut.point = 6
supp = seq(0,40,.01)
dens = .2*dnorm(supp, 25, 3) + .8*dexp(supp, 1/8)
const = 1-.2*pnorm(cut.point, 25, 3) -.8*pexp(cut.point, 1/8)
hist(out$d.save[7,(n.burn+1):n.mcmc], freq = FALSE, breaks = 100)
lines(supp, dens, type = 'l', lwd = 2)

res = c(mean(out$p.save[(n.burn+1):n.mcmc]), 
        quantile(out$p.save[(n.burn+1):n.mcmc], c(.025, .975)))
res = rbind(res, 
            c(mean(out$psi.save[(n.burn+1):n.mcmc]), 
              quantile(out$psi.save[(n.burn+1):n.mcmc], c(.025, .975))))
rownames(res) = c("p", "psi")
colnames(res) = c("Posterior Mean", "95% CI Lower", "95% CI Upper")
res
