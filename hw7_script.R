####################
# survival data
capt <- read.delim("data/NCEI_SurvivalCaptureData.csv", sep = ",", header = TRUE)
head(capt)
resight1 <- read.delim("data/NCEI_SurvivalResightData_1.csv", sep = ",", header = TRUE)
resight2 <- read.delim("data/NCEI_SurvivalResightData_2.csv", sep = ",", header = TRUE)
head(resight1)
head(resight2)
resight <- rbind(resight1, resight2)
head(resight)

# we need to use the natality data set if we want to model abundance of parous females
  

# natality data
natcap <- read.delim("data/NCEI_NatalityCaptureData.csv", sep = ",", header = TRUE)
head(natcap)
natres <- read.delim("data/NCEI_NatalityResightData.csv", sep = ",", header = TRUE)
head(natres)

# change sex variable to character
which(natcap$sex!="FALSE") # all females
sex <- replace(natcap$sex, seq(1, length(natcap$sex)), "F")
natcap <- cbind(natcap[,1:2], sex, natcap[,4:6])
which(natres$adjsex!="FALSE") # all females
adjsex <- replace(natres$adjsex, seq(1, length(natres$adjsex)), "F")
natres <- cbind(natres[,1:7], adjsex, natres[,9:10])

# rename for clarity
capture <- natcap
resight <- natres

# get rid of unnecessary variables
capture.clean <- capture[,c(1,3,6)]
length(unique(capture.clean$cohort))
resight.clean <- resight[,c(1,7,8,10)]

# change withpup to binary
withpup <- as.matrix(resight.clean$withpup)
withpup <- apply(withpup, 1, function(x){ if(x == "Y" | x == "y") x <- 1 else x <- 0})

resight.clean <- cbind(resight.clean[,1:3], withpup)
colnames(resight.clean) <- c("resight_year", "brand", "adjsex", "withpup")

table(capture.clean$cohort)
# 1991 
# 2000

# just want the brands for year 1991 in the capture data
capture.91 <- capture.clean[which(capture.clean$cohort == 1991),]
# brands captured in 1991
brands.91 <- capture.91$brand
# just want the resight data for the brands captured in 1991
resight.91 <- resight.clean[which(resight.clean$brand %in% brands.91),]

head(capture.91)
head(resight.91)

# create data frame: col = brand, row = year 
# 0 = female sighted but not with pup
# 1 = female sighted with pup
# NA = female note sighted

# get rid of duplicates: same sighting for same individual in same year
dups <- duplicated(cbind(resight.91$resight_year, resight.91$brand, 
                         resight.91$withpup))
rs.91 <- resight.91[!dups,]
# order: withpup = 1 comes first 
rs.91.order <- rs.91[order(rs.91$withpup, decreasing = TRUE),, drop = FALSE]
# order: increasing with year
rs.91.o2 <- rs.91.order[order(rs.91.order$resight_year, decreasing = FALSE),, drop = FALSE]

y <- matrix(0, nrow = length(unique(rs.91.o2$brand)), ncol = length(1991:(max(rs.91.o2$resight_year))))
rownames(y) <- unique(rs.91.o2$brand)
colnames(y) <- 1991:(max(rs.91.o2$resight_year))
inds <- which(rs.91.o2$withpup == 1)
for(ind in inds) {
  i = which(rownames(y) == rs.91.o2$brand[ind])
  t = which(colnames(y) == rs.91.o2$resight_year[ind])
  y[i,t] = 1
}

resight.91.clean <- reshape(rs.91.o2, idvar = "brand", drop = "adjsex", v.names = "withpup",
        timevar = "resight_year", direction = "wide")


### Clean Data ###

capture.91
# data frame of the individuals (4-5 month old female pups) captured and 
# branded in 1991
resight.91.clean
# data.frame of resights for the individuals captured in 1991 and whether 
# or not they had a pup with them 


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

y = sapply(brands, function(brand) {
  sapply(years, function(k) {
    sum(resight[resight$sightyear == k & resight$brand == brand, ]$withpup
        %in% c("Y", "y"))
  })
})

J = sapply(brands, function(brand) {
  sapply(years, function(k) {
    sum(resight$sightyear == k & resight$brand == brand)
  })
})


