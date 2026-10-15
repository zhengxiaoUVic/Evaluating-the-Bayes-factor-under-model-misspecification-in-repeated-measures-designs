Q <- function(a) {
  #' Input -
  #' a:  number of levels
  #' Details - 
  #' This approach is equivalent to constructing Helmert contrasts.
  #' `apply(contr.helmert(a), 2, function(col) col / sqrt(sum(col^2)))`
  #' Output - the projection matrix of size a × (a-1)

  if (a==2) return(c(1, -1)) # For a 2 × 2 repeated-measures design, we set a=b=1 for easier interpretation.
  S <- diag(a) - matrix(1, nrow=a, ncol=a) / a
  e <- qr(S)
  Q <- qr.Q(e) %*% diag(sign(diag(qr.R(e))))
  Q[, which(abs(e$qraux) > 1e-3), drop=F]
}

mVCOV <- function(a=2, b=2, sde=1, sd0=0.5, sd1=0.25, sd2=0.25, sd12=0.25) {
  #' Input -
  #' a:         number of levels of factor A
  #' b:         number of levels of factor B
  #' sde:       standard deviation for residuals
  #' sd0:       standard deviation for random intercepts
  #' sd1:       standard deviation for random slopes for projected factor A
  #' sd2:       standard deviation for random slopes for projected factor B
  #' sd12:      standard deviation for random slopes for the projected interaction AB
  #' Output - the marginal variance-covariance matrix for a two-way repeated-measures design

  sde^2 * diag(a*b) +
    sd0^2 * matrix(1, nrow=a*b, ncol=a*b) +
    sd1^2 * kronecker(matrix(1, nrow=b, ncol=b), Q(a) %*% t(Q(a))) +
    sd2^2 * kronecker(Q(b) %*% t(Q(b)), matrix(1, nrow=a, ncol=a)) +
    sd12^2 * kronecker(Q(b) %*% t(Q(b)), Q(a) %*% t(Q(a)))
}


# Verify unique components
length(unique(c(round(mVCOV(a=4, b=4, sd1=0.2, sd2=0.3, sd12=0.1), 3))))

mVCOV(a=2, b=2, sde=4, sd0=sqrt(5), sd1=1, sd2=1, sd12=1)
mVCOV(a=2, b=2, sde=2, sd0=sqrt(8), sd1=2, sd2=2, sd12=2) # same marginal variance-covariance matrix (partial identifiability)

mVCOV(a=2, b=2, sde=2, sd0=sqrt(3), sd1=0, sd2=1, sd12=sqrt(2))
mVCOV(a=2, b=2, sde=sqrt(2), sd0=sqrt(3.5), sd1=sqrt(0.5), sd2=sqrt(1.5), sd12=sqrt(2.5))

options(digits=3)
# Examples
mVCOV()
mVCOV(a=2, b=3)
mVCOV(a=2, b=3, sd1=0, sd2=0, sd12=0)
mVCOV(a=3, b=3)
mVCOV(a=3, b=3, sd1=0, sd2=0, sd12=0)
mVCOV(a=3, b=3, sd12=0)

# Compound symmetry holds in the presence of random slopes (sd1^2 - sd12^2 / b = sd2^2 - sd12^2 / a = 0)
mVCOV(a=3, b=3, sd1=1, sd2=1, sd12=sqrt(3))
mVCOV(a=3, b=3, sd1=0.5, sd2=0.5, sd12=sqrt(3*0.5^2))
mVCOV(a=2, b=3, sd1=1, sd2=sqrt(3), sd12=sqrt(3))


# Verify the coefficients
mVCOV(sde=0, sd0=0, sd1=1, sd2=0, sd12=0)
mVCOV(a=2, b=3, sde=0, sd0=0, sd1=1, sd2=0, sd12=0)
mVCOV(a=2, b=3, sde=0, sd0=0, sd1=0, sd2=1, sd12=0)
mVCOV(a=2, b=3, sde=0, sd0=0, sd1=0, sd2=0, sd12=1)
mVCOV(a=3, b=3, sde=0, sd0=0, sd1=1, sd2=0, sd12=0)
mVCOV(a=3, b=3, sde=0, sd0=0, sd1=1, sd2=1, sd12=0)
mVCOV(a=3, b=3, sde=0, sd0=0, sd1=0, sd2=0, sd12=1)





C <- function(a, b, orthonormal=T) {
  #' Input -
  #' a:  number of levels of factor A
  #' b:  number of levels of factor B
  #' Details - 
  #' https://doi.org/10.5964/meth.8415
  #' Output - the orthonormal / orthogonal contrast matrix of size ab × (ab-1)
  #'          for y₁₁, ..., yₐ₁, y₁₂, ..., yₐb
  
  B_matrix <- model.matrix(~A*B,
                           expand.grid(A=paste0("A",1:a), B=paste0("B",1:b)),
                           contrasts.arg=list(A="contr.poly", B="contr.poly")) # try also "contr.helmert"
  C_matrix <- solve(B_matrix)
  if(orthonormal) {
    apply(C_matrix, 1, function(row) row / sqrt(sum(row^2)))[,-1]
  } else {
    t(C_matrix)[-1,]
  }
}

# Sphericity
a <- 3; b <- 3
Sigma <- mVCOV(a=a, b=b, sd1=1, sd2=2, sd12=3)
Cmat <- C(a,b)
round(t(Cmat) %*% Sigma %*% Cmat,2)

Qmat <- cbind(kronecker(rep(1,b), Q(a)),
              kronecker(Q(b), rep(1,a)),
              kronecker(Q(b), Q(a)))
round(t(Qmat) %*% Sigma %*% Qmat,2)

