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

n.mcmc <- 1e5

out2 = hw7.mcmc2(y, z, n.mcmc = n.mcmc, r.tune = 1.5, d.tune = 1.5)

nburn <- .5*n.mcmc

plot(out2$p.save[-(1:nburn)], type = "l")
plot(out2$psi.save[-(1:nburn)], type = "l")
plot(out2$r.save[100,], type = "l")
plot(out2$lambda.save[-(1:nburn)], type = "l")



# still takes awhile to burn in - may be able to change that if we mess around with 
# the tuning parameters.

