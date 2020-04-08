abg_eqns <- function(t, vars, parms) {
    E <- vars[1]
    I <- vars[2]
    a <- parms[1]
    b <- parms[2]
    g <- parms[3]
    dE <- -a*E + b*I
    dI <- a*E - g*I
    return(list(c(dE, dI)))
}

v0 <- c(E=0, I=1)
parms <- c(alpha=1, beta=2, gamma=1, epsilon=0)
out <- deSolve::ode(v0, c(0,25,50), abg_eqns, parms)

remp <- function(mat, col=3) {
  t2 <- mat[3,]
  t1 <- mat[2,]
  
  (log(t2[3]) - log(t1[3])) / (t2[1] - t1[1])
}

lplus <- function(parms) {
  a <- parms[1]
  b <- parms[2]
  g <- parms[3]
  apg <- a+g
  
  0.5*(-apg + sqrt(apg^2 + 4*a*(b-g)))
}

print(remp(out))
print(lplus(parms))

u0 <- c(S=1e6, E=0, I=1, Is=0, R=0)
out2 <- deSolve::ode(u0, c(0, 10, 20), seir_equations, parms)

print(remp(out2,4))
print(lplus(parms))

beta <- calcbeta_approx(c(T0=5, D0=14, A0=3))
parms2 <- c(alpha=1/3, beta=beta, gamma=1/14, epsilon=0)
out2a <- deSolve::ode(u0, c(0,10,20), seir_equations, parms2)
print(remp(out2a, 4))
print(lplus(parms2))
print(calclplus(parms2[1], parms2[2], parms2[3]))

parms3 <- c(alpha=1/3, beta=beta, gamma=1/14, epsilon=1/5)
out3 <- deSolve::ode(u0, c(0,10,20), seir_equations, parms3)
itot <- out3[,'I'] + out3[,'Is']
print((log(itot[3]) - log(itot[2])) / (out3[3,1] - out3[2,1]))
print(lplus(parms3))
print(calclplus(parms3[1], parms3[2], parms3[3]))

modparms <- c(A0=3, T0=6.35339, D0=14, Ts=5)
calcbeta_test <- function(parms)
{
  beta_guess <- calcbeta_approx(parms)
  alpha <- calcalpha(parms)
  gamma <- calcgamma(parms)
  targfun <- function(beta) {
    parms[['T0']] - calct0(alpha, beta, gamma)
  }
  solv <- nleqslv::nleqslv(beta_guess, targfun)
  solv$x
}
print(calcbeta_test(modparms))

mp <- CovMitigation:::param_defaults
for (n in names(modparms)){
  mp[[n]] <- modparms[[n]]
}
