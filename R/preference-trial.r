
setGeneric("preference_trial",
  function(sample_size, power, sigma2, phi, delta_nu, delta_pi, delta_tau, 
           alpha=0.05, theta=0.5, xi=1, nstrata=1L, k=1L) {
    standardGeneric("preference_trial")
  })

setMethod("preference_trial",
  signature(sample_size="integer", power="missing", sigma2="numeric", 
            phi="numeric", delta_nu="numeric", delta_pi="numeric", 
            delta_tau="numeric", alpha="numeric", theta="numeric",
            xi="numeric", nstrata="integer", k="integer"),
  function(sample_size, sigma2, phi, delta_nu, delta_pi, delta_tau, 
           alpha, theta, xi, nstrata, k) {
    ret <- NULL
    max_len <- max(c(length(sample_size), length(sigma2), length(phi), 
                     length(delta_nu), length(delta_pi), length(delta_tau), 
                     length(alpha), length(theta), length(xi), length(nstrata), 
                     length(k)))
    for (i in 1:length(max_len)) {
      ret <- rbind(ret, 
        data.frame(sample_size=sample_size[i %% length(sample_size)], 
                   sigma2=sigma2[i %% length(sigma2)], 
                   phi=phi[i %% length(phi)], 
                   delta_nu=delta_nu[i %% length(delta_nu)], 
                   delta_pi=delta_pi[i %% length(delta_pi)], 
                   delta_tau=delta_tau[i %% length(delta_tau)], 
                   alpha=alpha[i %% length(alpha)], 
                   theta=theta[i %% length(theta)], 
                   xi=xi[i %% length(xi)], 
                   nstrata=nstrata[i %% length(nstrata)],
                   k=k[i %% length(k)]))
    }
    class(ret) <- c("preference_trial", class(ret))
    ret
  })

setMethod("preference_trial",
  signature(sample_size="missing", power="numeric", sigma2="numeric", 
            phi="numeric", delta_nu="numeric", delta_pi="numeric", 
            delta_tau="numeric", alpha="numeric", theta="numeric", 
            xi="numeric", nstrata="integer", k="integer"),
  function(sample_size, sigma2, phi, delta_nu, delta_pi, delta_tau, 
           alpha, theta, xi, nstrata, k) {
    # We'll set sample size to power to create the ret variable. Then
    # we'll fill it in with the actual sample size.
    sample_size <- power
    ret <- NULL
    max_len <- max(length(sample_size), length(sigma2), length(phi), 
                   length(delta_nu), length(delta_pi), length(delta_tau), 
                   length(alpha), length(theta), length(xi), 
                   length(nstrata), length(k))
    for (i in 1:length(max_len)) {
      ret <- rbind(ret, 
        data.frame(sample_size=sample_size[i %% length(sample_size)], 
                   sigma2=sigma2[i %% length(sigma2)], 
                   phi=phi[i %% length(phi)], 
                   delta_nu=delta_nu[i %% length(delta_nu)], 
                   delta_pi=delta_pi[i %% length(delta_pi)], 
                   delta_tau=delta_tau[i %% length(delta_tau)], 
                   alpha=alpha[i %% length(alpha)], 
                   theta=theta[i %% length(theta)], 
                   xi=xi[i %% length(xi)], 
                   nstrata=nstrata[i %% length(nstrata)],
                   k=k[i %% length(k)]))
    }
    for (i in 1:nrow(ret)) {
      ret$sample_size[i] <- max(unlist(
        overall_sample_size(power[i %% max_len], ret$phi[i], ret$sigma2[i],
                            ret$delta_pi[i], ret$delta_nu[i], ret$delta_tau[i],
                            ret$alpha[i], ret$theta[i], ret$xi[i], 
                            ret$nstrata[i])))
    }
    class(ret) <- c("preference_trial", class(ret))
    ret
  })
