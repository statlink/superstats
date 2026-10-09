main <- function(N = 20, d = 2, n = 500, nu = 1.8, sg = 0.25, g = 0.00036, dist = 1) {
  # N: number of resolution for characteristic function
  # d: number of dimension of generated data
  # n: number of data to be generated
  # nu, sg, g: learning parameter
  # dist: data distribution 0:gaussian, 1:uniform

  if ( dist == 0 ) { # generate d dimension of gaussian
    Mu <- 1 * matrix(1, 1, d)
    Sg2 <- 2 * diag(d)
    x <- Rfast::rmvnorm(n, Mu, Sg2)
    # define output range
    outrange <- c(-5, 5)
  } else {
    # or generate d dimension of uniform data or uniform in range of [a b]
    a <-  -2
    b <- 3
    x <- a + (b - a) * matrix( Rfast2::Runif(n * d), n, d )
    # define output range
    outrange <- c(a - 1, b + 1)
  }

  dev.new()
  # plot generated data
  if ( d == 2 ) {
    plot(x[, 1], x[, 2], col = "blue", pch = 19, xlab = "1st dimension", ylab = "2nd dimension", cex.axis = 1.3, cex.lab = 1.3)
    title("Generated data", cex.main = 1.8)
    par(cex.lab = 1.8)
  } else {
    plot(x, col = "red", pch = 19, xlab = "Data number", ylab = "1st dimension", cex.axis = 1.3, cex.lab = 1.3)
    title("Generated data")
  }

  # calculate characteristic function by rbf network described by paper
  res <- find_cf(x, outrange, N, nu, sg, g)

  dev.new()
  # plot real and imaginary part
  if ( d == 2 ) {
    par( mfrow = c(1, 2) )
    persp(x = 1:N, y = 1:N, z = matrix(res$Teta, nrow = N), main = "Peal part of C.F. from network", col = "lightblue",
          xlab = "x", ylab = "y", zlab = "z", cex.axis = 1.3, cex.lab = 1.3)
    persp(x = 1:N, y = 1:N, z = matrix(res$TetaI, nrow = N), main = "Imaginary part of C.F. from network", col = "lightblue",
          xlab = "x", ylab = "y", zlab = "z", cex.axis = 1.3, cex.lab = 1.3)

  } else if ( d == 1 ) {
    par( mfrow = c(1, 2) )
    plot(res$Teta, main = "Peal part of C.F. from network", xlab = "Index", ylab = "Value",
         cex.main = 1.8, cex.axis = 1.3, cex.lab = 1.3)
    plot(res$TetaI, main = "Imaginary part of C.F. from network", xlab = "Index", ylab = "Value",
         cex.main = 1.8, cex.axis = 1.3, cex.lab = 1.3)
  }

  # theoretic characteristic function  plot from wikipedia
  # https://en.wikipedia.org/wiki/Characteristic_function_(probability_theory)
  # look table of theoretic characteristic functions of certain distributions
  # and also we plotted the error between theoretic characteristic and
  # calculated by network optimization

  dev.new()
  if ( dist == 1 ) { # uniform theoric result
    #fr <- ( exp(i * b * res$Alfa) - exp(i * a * res$Alfa) ) / ( i * (b - a) * res$Alfa )
    fr <- res$Alfa
    fr[, 1] <- ( exp( complex(imaginary = b * res$Alfa[, 1]) ) - exp( complex( imaginary = a * res$Alfa[, 1]) ) ) /
                 complex( imaginary = (b - a) * res$Alfa[, 1] )
    fr[, 2] <- ( exp( complex(imaginary = b * res$Alfa[, 2]) ) - exp( complex( imaginary = a * res$Alfa[, 2]) ) ) /
                 complex( imaginary = (b - a) * res$Alfa[, 2] )

    if ( d == 2 ) {
      # Calculate fr
      fr <- fr[, 2] * fr[, 1]
      par( mfrow = c(2, 2) )
      persp(x = 1:N, y = 1:N, z = matrix(Re(fr), nrow = N), col = "lightblue", main = "Real part of C.F. from theory",
            xlab = "x", ylab = "y", zlab = "z", cex.axis = 1.3, cex.lab = 1.3)
      persp(x = 1:N, y = 1:N, z = matrix(Im(fr), nrow = N), col = "lightblue", main = "Imaginary part of C.F. from theory",
            xlab = "x", ylab = "y", zlab = "z", cex.axis = 1.3, cex.lab = 1.3)
      persp(x = 1:N, y = 1:N, z = matrix(Re(fr) - res$Teta, nrow = N), col = "lightblue", main = "Error of real part of C.F.",
            xlab = "x", ylab = "y", zlab = "z", cex.axis = 1.3, cex.lab = 1.3)
      persp(x = 1:N, y = 1:N, z = matrix(Im(fr) - res$TetaI, nrow = N), col = "lightblue", main = "Error of imaginary part of C.F.",
            xlab = "x", ylab = "y", zlab = "z", cex.axis = 1.3, cex.lab = 1.3)

    } else if ( d == 1 ) {
      par( mfrow = c(2, 2) )
      plot(Re(fr), main = "Real part of C.F. from theory", cex.main = 1.8, cex.axis = 1.3, cex.lab = 1.3)
      plot(Im(fr), main = "Imaginary part of C.F. from theory", cex.main = 1.8, cex.axis = 1.3, cex.lab = 1.3)
      plot(Re(fr) - res$Teta, main = "Error of real part of C.F.", cex.main = 1.8, cex.axis = 1.3, cex.lab = 1.3)
      plot(Im(fr) - res$TetaI, main = "Error of imaginary part of C.F.", cex.main = 1.8, cex.axis = 1.3, cex.lab = 1.3)
    }

  } else { # normal theoretic result from wikipedia
    fr <- exp( complex( imaginary = res$Alfa %*% t(Mu), real = - 0.5 * diag( res$Alfa %*% Sg2 %*% t(res$Alfa) ) ) )

    if ( d == 2 ) {
      par( mfrow = c(2, 2) )
      persp(z = matrix(Re(fr), nrow = N), theta = 30, phi = 30, main = "Real part of C.F. from theory", col = "lightblue",
            xlab = "x", ylab = "y", zlab = "z", cex.axis = 1.3, cex.lab = 1.3)
      persp(z = matrix(Im(fr), nrow = N), theta = 30, phi = 30, main = "Imaginary part of C.F. from theory", col = "lightblue",
            xlab = "x", ylab = "y", zlab = "z", cex.axis = 1.3, cex.lab = 1.3)
      persp(z = matrix(Re(fr) - res$Teta, nrow = N), theta = 30, phi = 30, main = "Error of real part of C.F.", col = "lightblue",
            xlab = "x", ylab = "y", zlab = "z", cex.axis = 1.3, cex.lab = 1.3)
      persp(z = matrix(Im(fr) - res$TetaI, nrow = N), theta = 30, phi = 30, main = "Error of imaginary part of C.F.", col = "lightblue",
            xlab = "x", ylab = "y", zlab = "z", cex.axis = 1.3, cex.lab = 1.3)

    } else if ( d == 1 ) {
      par( mfrow = c(2, 2) )
      plot(Re(fr), main = "Real part of C.F. from theory", cex.main = 1.8, cex.axis = 1.3, cex.lab = 1.3)
      title(main = "Real part of C.F. from theory", cex.main = 1.8)
      plot(Im(fr), main = "Imaginary part of C.F. from theory", cex.main = 1.8, cex.axis = 1.3, cex.lab = 1.3)
      title(main = "Imaginary part of C.F. from theory", cex.main = 1.8)
      plot(Re(fr) - res$Teta, main = "Error of real part of C.F.", cex.main = 1.8, cex.axis = 1.3, cex.lab = 1.3)
      title(main = "Error of real part of C.F.", cex.main = 1.8)
      plot(Im(fr) - res$TetaI, main = "Error of imaginary part of C.F.", cex.main = 1.8, cex.axis = 1.3, cex.lab = 1.3)
      title(main = "Error of imaginary part of C.F.", cex.main = 1.8)
    }
  }

  # error of the real and imaginary
  real_error <- sum( ( res$Teta - Re(fr) )^2 )
  imag_error <- sum( ( res$TetaI - Im(fr) )^2 )
  res <- c(real_error, imag_error)
  names(res) <- c("Re error", "Im error")
  res
}
