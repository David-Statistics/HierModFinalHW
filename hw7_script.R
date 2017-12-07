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

source("./hw7_mcmc.R")
set.seed(1215)
n.mcmc <- 5e4

out = hw7.mcmc(y, z, n.mcmc = n.mcmc, r.tune = .5, d.tune = 3,
                alpha.p = 1.5, beta.p = 4.5, alpha.psi = 2, beta.psi = 2)

n.burn <- .1*n.mcmc

out$d.mh.prop[1] / out$d.mh.prop[2]
out$r.mh.prop[1] / out$r.mh.prop[2]

plot(out$p.save[(n.burn+1):n.mcmc], type = "l")
plot(out$psi.save[(n.burn+1):n.mcmc], type = "l")
plot(out$r.save[4,(n.burn+1):n.mcmc], type = "l")
plot(out$d.save[4,(n.burn+1):n.mcmc], type = "l", main = "last seen year 0")
plot(out$d.save[2,(n.burn+1):n.mcmc], type = "l", main = 'last seen year 17')
plot(out$d.save[7,(n.burn+1):n.mcmc], type = "l", main = 'last seen year 6')
plot(out$d.save[38,(n.burn+1):n.mcmc], type = "l", main = 'last seen year 3')
plot(out$d.save[114,(n.burn+1):n.mcmc], type = "l", main = 'last seen year 2')
plot(out$d.save[100,(n.burn+1):n.mcmc], type = "l", main = "last seen year 1")

d = apply(z, 1, function(x) {
  if(sum(x) > 0) {
    return(max(which(x > 0)) - 1 + abs(rnorm(1,0,.1)))
  }
  return(0)
})
min.d = floor(d)

cut.points = c(0,1,2,3,4,5,10,15,17)
supp = seq(0,40,.01)
dens = .2*dnorm(supp, 25, 3) + .8*dexp(supp, 1/8)
par(mfrow = c(3,3), oma = c(2,0,0,0))
for(k in seq_along(cut.points)) {
  cut.point = cut.points[k]
  inds = which(min.d == cut.point)
  const = 1-.2*pnorm(cut.point, 25, 3) -.8*pexp(cut.point, 1/8)
  plot(density(as.vector(out$d.save[inds,(n.burn+1):n.mcmc])), breaks = 100, lwd = 2,
       main = paste("Individuals Last Seen at", cut.point, "years"), xlab = "d")
  lines(supp[supp > cut.point], dens[supp > cut.point]/const, type = 'l', lwd = 2,
        col = 'grey50')
}
par(xpd=NA)
legend(x = -117, y = -.5, legend=c("Posterior", "Prior"), lwd=c(2,2), col= c("black", "grey50"))

par(mfrow = c(3,3))
for(k in seq_along(cut.points)) {
  #cut.point = cut.points[k]
  inds = which(min.d == cut.point)
  #const = 1-.2*pnorm(cut.point, 25, 3) -.8*pexp(cut.point, 1/8)
  plot(density(as.vector(out$r.save[inds,(n.burn+1):n.mcmc])), breaks = 100, lwd = 2,
       main = paste("Individuals Last Seen at", cut.point, "years"), xlab = "r")
}


library(xtable)
res = c(mean(out$p.save[(n.burn+1):n.mcmc]), 
        quantile(out$p.save[(n.burn+1):n.mcmc], c(.025, .975)))
res = rbind(res, 
            c(mean(out$psi.save[(n.burn+1):n.mcmc]), 
              quantile(out$psi.save[(n.burn+1):n.mcmc], c(.025, .975))))
rownames(res) = c("p", "psi")
colnames(res) = c("Posterior Mean", "95% CI Lower", "95% CI Upper")
print(xtable(res))
