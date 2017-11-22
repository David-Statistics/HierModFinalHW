hw7.mcmc <- function(y, z, 
                     p.tune = .05, alpha.p = 1, beta.p = 1, 
                     psi.tune = .05, alpha.psi = 1, beta.psi = 1, 
                     d.tune = .8, r.tune = .8,
                     mu0 = 1.5, s0 = .2, a = 1, b = 10,
                     n.mcmc = 1e4) {
  
  ###
  ### Set up data
  ###
  
  r = apply(y, 1, function(x) {
    if(sum(x) > 0) {
      return(which(x > 0)[1] - abs(rnorm(1,0,r.tune)))
    } 
    return(18)
  })
  
  d = apply(z, 1, function(x) {
    if(sum(x) > 0) {
      return(max(which(x > 0)) + abs(rnorm(1,0,.1)))
    }
    return(18)
  })
  rs.y = rowSums(y)
  rs.z = rowSums(z)
  max.r = ceiling(r)
  min.d = floor(d)
  
  p = .01
  psi = .01
  
  p.save = numeric(n.mcmc)
  psi.save = numeric(n.mcmc)
  r.save = matrix(0, nrow = length(r), ncol = n.mcmc)
  d.save = matrix(0, nrow = length(d), ncol = n.mcmc)
  lambda.save = numeric(n.mcmc)
  
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
    
    p.prop <- rnorm(1, p, p.tune)
    if(p.prop > 0 & p.prop < 1) {
      mh.num <- sum(dbinom(rs.y, n.avail.years, p.prop, log = TRUE)) +
        dbeta(p.prop, alpha.p, beta.p, log = TRUE)
      mh.den <- sum(dbinom(rs.y, n.avail.years, p, log = TRUE)) +
        dbeta(p, alpha.p, beta.p, log = TRUE)
      if(log(runif(1)) < mh.num - mh.den) {
        p = p.prop
      }
    }
    
    ###
    ###  update psi
    ###
    
    psi.prop <- rnorm(1, psi, psi.tune)
    if(psi.prop > 0 & psi.prop < 1) {
      mh.num <- sum(dbinom(rs.z, floor(d) - 1, psi.prop, log = TRUE)) +
        dbeta(psi.prop, alpha.psi, beta.psi, log = TRUE)
      mh.den <- sum(dbinom(rs.z, floor(d) - 1, psi, log = TRUE)) +
        dbeta(psi, alpha.psi, beta.psi, log = TRUE)
      if(log(runif(1)) < mh.num - mh.den) {
        psi = psi.prop
      }
    }
    
    
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
