hw7.mcmc2 <- function(y, z, 
                     alpha.p = .5, beta.p = 1.5, 
                     alpha.psi = .01, beta.psi = .005, 
                     mu0 = 1.6, s20 = .0625,
                     a.sigma = 3, b.sigma = 1,
                     d.tune = .8, r.tune = .8,
                     a.lambda = 10, b.lambda = 250,
                     n.mcmc = 1e4) {
  
  ###
  ### Set up data
  ###
  
  # y = seen with pup 
  # z = seen at all 
  
  
  # starting values for time until recruitment
  r = apply(y, 1, function(x) {
    if(sum(x) > 0) {
      return(which(x > 0)[1] - abs(rnorm(1,0,r.tune)))
    } 
    return(18)
  })
  
  # starting values for time until death
  d = apply(z, 1, function(x) {
    if(sum(x) > 0) {
      return(max(which(x > 0)) + abs(rnorm(1,0,.1)))
    }
    return(1)
  })
  
  
  # how many years each individual was seen with a pup
  rs.y = rowSums(y)
  
  # how many years each individual was seen
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
  
  p = .1
  psi = .1
  n = length(d)
  mu_r = mu0
  s2_r = a.sigma/b.sigma
  s_r = sqrt(s2_r)
  
  p.save = numeric(n.mcmc)
  psi.save = numeric(n.mcmc)
  r.save = matrix(0, nrow = length(r), ncol = n.mcmc)
  d.save = matrix(0, nrow = length(d), ncol = n.mcmc)
  mu_r.save = numeric(n.mcmc)
  s2_r.save = numeric(n.mcmc)
  lambda.save = numeric(n.mcmc)
  r.mh.prop = c(0,0)
  d.mh.prop = c(0,0)
  
  ###
  ### Start MCMC
  ###
  t1 = Sys.time()
  for(k in seq_len(n.mcmc)) {
    if(k %% 1e4 == 0) {
      print(paste(k, Sys.time() - t1))
    }
    ###
    ### update p
    ###
    
    # count the number of years between r_i and d_i
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
      dbinom(rs.y[i], min(18, n.avail.years.prop[i]), p, log = TRUE) +
        dlnorm(r.prop[i], mu_r, s_r, log = TRUE) -
        dbinom(rs.y[i], min(18, n.avail.years[i]), p, log = TRUE) -
        dlnorm(r[i], mu_r, s_r, log = TRUE)
    })
    updates = valid.r[which(log(runif(length(valid.r))) < mh.cuts)]
    r[updates] = r.prop[updates]
    r.mh.prop[1] = r.mh.prop[1] + length(updates)
    r.mh.prop[2] = r.mh.prop[2] + sum(valid.r)
    
    ###
    ### update mu_r
    ###
    
    mu_r = rnorm(1, (s2_r * mu0 + s20*sum(log(r)))/(s2_r+n*s20), 
                 sqrt(s2_r*s20/(s2_r+n*s20)))
    
    ###
    ### update s2_r
    ###
    
    s2_r = 1/rgamma(1,n/2 + a.sigma, rate = sum((log(r) - mu_r)^2)/2 + 1/b.sigma)
    s_r = sqrt(s2_r)
    
    ###
    ### update lambda
    ###
    
    lambda = rgamma(1, a.lambda + length(d), b.lambda + sum(d))
    
    ###
    ### update d
    ###
    
    d.prop = d + rnorm(length(d), 0, d.tune)
    valid.d = which(d.prop > min.d & d.prop > 0 & d.prop > r)
    mh.cuts = sapply(valid.d, function(i) {
      dbinom(rs.z[i], min(18, floor(d.prop[i])), psi, log = TRUE) +
        dexp(d.prop[i], lambda, log = TRUE) -
        dbinom(rs.z[i], min(18, floor(d[i])), psi, log = TRUE) -
        dexp(d[i], lambda, log = TRUE)
    })
    updates = valid.d[which(log(runif(length(valid.d))) < mh.cuts)]
    d[updates] = d.prop[updates]
    d.mh.prop[1] = d.mh.prop[1] + length(updates)
    d.mh.prop[2] = d.mh.prop[2] + sum(valid.d)
    
    p.save[k] = p
    psi.save[k] = psi
    r.save[,k] = r
    d.save[,k] = d
    lambda.save[k] = lambda
    mu_r.save[k] = mu_r
    s2_r.save[k] = s2_r
  }
  
  return(list(p.save = p.save,
              psi.save = psi.save,
              r.save = r.save,
              d.save = d.save,
              lambda.save = lambda.save,
              mu_r.save = mu_r.save,
              s2_r.save = s2_r.save,
              r.mh.prop = r.mh.prop,
              d.mh.prop = d.mh.prop,
              n.mcmc = n.mcmc))
  
}
