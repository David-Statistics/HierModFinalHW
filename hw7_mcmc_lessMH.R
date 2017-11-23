hw7.mcmc2 <- function(y, z, 
                     p.tune = .01, alpha.p = 1, beta.p = 1, 
                     psi.tune = .01, alpha.psi = 1, beta.psi = 1, 
                     d.tune = .8, r.tune = .8,
                     mu0 = 1.5, s0 = .2, a = 1, b = 10,
                     n.mcmc = 1e4) {
  
  ###
  ### Set up data
  ###
  
  # y = seen with pup 
  # z = seen at all 
  
  
  # starting values for number recruited
  r = apply(y, 1, function(x) {
    if(sum(x) > 0) {
      return(which(x > 0)[1] - abs(rnorm(1,0,r.tune)))
    } 
    return(18)
  })
  
  # starting values for number dead/alive?? 
  d = apply(z, 1, function(x) {
    if(sum(x) > 0) {
      return(max(which(x > 0)) + abs(rnorm(1,0,.1)))
    }
    return(1)
  })
  
  
  # how many were seen with a pup each year
  rs.y = rowSums(y)
  
  # how many were seen each year
  rs.z = rowSums(z)
  
  # first year we saw a pup for each individual
  max.r = ceiling(r)
  r[rowSums(y) == 0] = 1
  
  # last year we saw each individual
  min.d = floor(d)
  
  # total pups
  n.pups = sum(y)
  
  # total sightings
  n.sightings = sum(z)
  
  p = .01
  psi = .01
  
  p.save = numeric(n.mcmc)
  psi.save = numeric(n.mcmc)
  r.save = matrix(0, nrow = length(r), ncol = n.mcmc)
  d.save = matrix(0, nrow = length(d), ncol = n.mcmc)
  lambda.save = numeric(n.mcmc)
  
  
  # for mh acceptance ratio
  kp = 1
  kpsi = 1
  
  
  ###
  ### Start MCMC
  ###
  
  for(k in seq_len(n.mcmc)) {
    #print(k)
    ###
    ### update p
    ###
    
    n.avail.years = sapply(seq_along(r), function(i) {
      length(ceiling(r[i]):floor(d[i]))
    })
    
    p = rbeta(1, alpha.p + n.pups, beta.p + sum(n.avail.years) - n.pups)
    
    
    ###
    ###  update psi
    ###
    
    psi = rbeta(1, alpha.psi + n.sightings, beta.psi + sum(floor(d)) - n.sightings)
    
    ###
    ### update r
    ###
    
    r.prop = r + rnorm(length(r), 0, r.tune)
    valid.r = which(r.prop < max.r & r.prop > 0 & r.prop < d)
    n.avail.years.prop = n.avail.years
    n.avail.years.prop[valid.r] = sapply(valid.r, function(i) {
      length(ceiling(r.prop[i]):floor(d[i]))
    })
    mh.cuts = sapply(valid.r, function(i) {
      dbinom(rs.y[i], n.avail.years.prop[i], p, log = TRUE) +
        dlnorm(r.prop[i], mu0, s0, log = TRUE) -
        dbinom(rs.y[i], n.avail.years[i], p, log = TRUE) -
        dlnorm(r[i], mu0, s0, log = TRUE)
    })
    updates = valid.r[which(log(runif(length(valid.r))) < mh.cuts)]
    r[updates] = r.prop[updates]
    
    ###
    ### update lambda
    ###
    
    lambda = rgamma(1, a + length(d), b + sum(d))
    
    ###
    ### update d
    ###
    
    d.prop = d + rnorm(length(d), 0, d.tune)
    valid.d = which(d.prop > min.d & d.prop > 0 & d.prop > r)
    mh.cuts = sapply(valid.d, function(i) {
      dbinom(rs.z[i], floor(d.prop[i]), psi, log = TRUE) +
        dexp(d.prop[i], lambda, log = TRUE) -
        dbinom(rs.z[i], floor(d[i]), psi, log = TRUE) -
        dexp(d[i], lambda, log = TRUE)
    })
    updates = valid.d[which(log(runif(length(valid.d))) < mh.cuts)]
    d[updates] = d.prop[updates]
    
    p.save[k] = p
    psi.save[k] = psi
    r.save[,k] = r
    d.save[,k] = d
    lambda.save[k] = lambda
  }
  
  return(list(p.save = p.save,
              psi.save = psi.save,
              r.save = r.save,
              d.save = d.save,
              lambda.save = lambda.save,
              n.mcmc = n.mcmc))
  
}
