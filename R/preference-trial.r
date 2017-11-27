
# Circular indexing to recycle values.
cind <- function(i, vec) {
  (i-1) %% length(vec) + 1
}

setGeneric("preference.trial",
  function(sample_size, power, sigma2, phi, delta_nu, delta_pi, delta_tau, 
           alpha=0.05, theta=0.5, xi=1, nstrata=1L, k=1L) {
    standardGeneric("preference.trial")
  })

setMethod("preference.trial",
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
        data.frame(sample_size=sample_size[cind(i, sample_size)],
                   sigma2=sigma2[cind(i, sigma2)], 
                   phi=phi[cind(i, phi)],
                   delta_nu=delta_nu[cind(i, delta_nu)],
                   delta_pi=delta_pi[cind(i, delta_pi)],
                   delta_tau=delta_tau[cind(i, delta_tau)],
                   alpha=alpha[cind(i, alpha)], 
                   theta=theta[cind(i, theta)],
                   xi=xi[cind(i, xi)],
                   nstrata=nstrata[cind(i, nstrata)],
                   k=k[cind(i, k)]))
    }
    class(ret) <- c("preference.trial", class(ret))
    ret
  })

setMethod("preference.trial",
  signature(sample_size="missing", power="numeric", sigma2="numeric", 
            phi="numeric", delta_nu="numeric", delta_pi="numeric", 
            delta_tau="numeric", alpha="numeric", theta="numeric", 
            xi="numeric", nstrata="integer", k="integer"),
  function(power, sigma2, phi, delta_nu, delta_pi, delta_tau, 
           alpha, theta, xi, nstrata, k) {
    # We'll set sample size to power to create the ret variable. Then
    # we'll fill it in with the actual sample size.
    sample_size <- power
    ret <- NULL
    max_len <- max(length(power), length(sigma2), length(phi), 
                   length(delta_nu), length(delta_pi), length(delta_tau), 
                   length(alpha), length(theta), length(xi), 
                   length(nstrata), length(k))
    for (i in 1:length(max_len)) {
      ret <- rbind(ret, 
        data.frame(sample_size=power[cind(i, power)], 
                   sigma2=sigma2[cind(i, sigma2)], 
                   phi=phi[cind(i, phi)], 
                   delta_nu=delta_nu[cind(i, delta_nu)],
                   delta_pi=delta_pi[cind(i, delta_pi)],
                   delta_tau=delta_tau[cind(i, delta_tau)], 
                   alpha=alpha[cind(i, alpha)], 
                   theta=theta[cind(i, theta)], 
                   xi=xi[cind(i, xi)],
                   nstrata=nstrata[cind(i, nstrata)],
                   k=k[cind(i, k)]))
    }
    for (i in 1:nrow(ret)) {
      ret$sample_size[i] <- max(unlist(
        overall_sample_size(power[cind(i, seq_len(max_len))], ret$phi[i], 
                            ret$sigma2[i], ret$delta_pi[i], ret$delta_nu[i], 
                            ret$delta_tau[i], ret$alpha[i], ret$theta[i], 
                            ret$xi[i], ret$nstrata[i])))
    }
    class(ret) <- c("preference.trial", class(ret))
    ret
  })
