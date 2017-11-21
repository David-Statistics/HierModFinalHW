rm(list = ls())

cap.dat = read.csv("./data/NCEI_SurvivalCaptureData.csv", stringsAsFactors = FALSE)
resight1.dat = read.csv("./data/NCEI_SurvivalResightData_1.csv", stringsAsFactors = FALSE)
resight2.dat = read.csv("./data/NCEI_SurvivalResightData_2.csv", stringsAsFactors = FALSE)

resight.dat = rbind(resight1.dat, resight2.dat)
rm(resight1.dat, resight2.dat)

## Who thought to themselves, 
## 'yeah sure, why not have our two data sets have different date formats'
## (Didn't install lubridate on this computer before take off so had to do this by hand)

cap.dat$year = unname(sapply(cap.dat$branddate, 
                                          function(s) {
                                            year.end = as.numeric(strsplit(s, "/")[[1]][3])
                                            if (year.end < 50) {
                                              return(year.end + 2000)
                                            } 
                                            return(year.end + 1900)
                                            }))
resight.dat$year = as.numeric(unname(sapply(resight.dat$sitedate, 
                                 function(s) strsplit(s, "/")[[1]][3])))

min.year = min(min(cap.dat$year), min(resight.dat$year))
max.year = max(max(cap.dat$year), max(resight.dat$year))

## Create a matrix of sighting occurances regardless of location
## Get an idea of time frame we're working with

sighting.mat = matrix(0, nrow = nrow(cap.dat), ncol = length(min.year:max.year))
rownames(sighting.mat) = cap.dat$brand
colnames(sighting.mat) = min.year:max.year
for(k in seq_along(cap.dat$brand)) {
  brand = cap.dat$brand[k]
  j = which(colnames(sighting.mat) == cap.dat$year[k])
  sighting.mat[k,j] = 1
  resight.years = resight.dat$year[which(resight.dat$brand == brand)]
  js = sapply(resight.years, function(r.year) which(colnames(sighting.mat) == r.year))
  for(j in js) {
    sighting.mat[k,j] = sighting.mat[k,j] + 1
  }
}
trimmed.sm = sighting.mat[rowSums(sighting.mat) > 1, ]
bern.tsm = trimmed.sm
bern.tsm[bern.tsm > 1] = 1
image(t(log(trimmed.sm+1)))
image(t(bern.tsm))

## Look at how frequent switching sites within a year is

inds = which(rowSums(sighting.mat) > 1)
res = sapply(inds, function(ind) {
  rd.inds = which(resight.dat$brand == cap.dat$brand[ind])
  mean(apply(table(resight.dat[rd.inds, c(3,9)]), 2, function(x) sum(x > 0) == 1))
}) 


