find_cf <- function(x, outrange, N, nu, sg, g) {
  # N: single integer value
  # nu, sg, g: single numeric values
  # this function finds the characteristic function of given data x
  # theoutput will be in range of [output(1) output(2)] and N number of cell
  # in each dimension, nu,sg,g are the learning parameters
  # Teta is te real part of characteristic function
  # TetaI is the imaginary part of characteristic function
  # Alfa is the domain of the output generated

  n <- dim(x)[1]  ;  d <- dim(x)[2]
  alfax <- Rfast::rep_row(seq(from = outrange[1], to = outrange[2], length.out = N), d)
  iX <- matrix(0, n, d)

  for ( i in 1:d ) {
    out <- outer(x[, i], alfax[i, ], "-")^2
    iX[, i] <- Rfast::rowMins(out, value = FALSE)
  }

  N.d <- N^d
  # teta and W for real part network
  Teta <- Rfast::Rnorm(N.d)
  W <- Rfast::Rnorm(N.d)

  # tetea and W for imaginary part network
  TetaI <- Rfast::Rnorm(N.d)
  WI <- Rfast::Rnorm(N.d)

  Alfa <- createalfa(alfax)
  K <- matrix(0, nrow(Alfa), N.d)
  # calculate kernel once
  for ( j in 1:N.d ) {
    K[, j] <- exp( -Rfast::rowsums(Rfast::eachrow(Alfa, Alfa[j, ], oper = "-")^2 ) / ( 2 * sg^2) )
    K[, j] <- K[, j] / sum(K[, j])
  }
  K[K < 1e-14] <- 0

  for ( i in 1:n)  {
    #print(i)
    ind <- iX[i, ]
    if ( d > 1 )  ind <- (ind[-d] - 1) * N + ind[d]
    # real part learning
    Teta <- Teta + g * ( cos( Alfa %*% t( x[i, , drop = FALSE] ) ) - Teta)
    tmp <- sum( Teta - Rfast::rowsums( (K * W) * as.vector(Teta) ) )
    W <- as.vector(W + nu * tmp * Teta * K[, ind])
    # imaginary part learning
    TetaI <- TetaI + g * ( sin( Alfa %*% t( x[i, , drop=FALSE] ) ) - TetaI )
    tmp <- sum( TetaI - Rfast::rowsums( K * WI * as.vector(TetaI) ) )
    WI <- as.vector(WI + nu * tmp * TetaI * K[, ind])
  } ## end  for ( i in 1:n)

  list(Teta = Teta, TetaI = TetaI, Alfa = Alfa)
}
