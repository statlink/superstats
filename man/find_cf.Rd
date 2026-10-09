\name{find_cf}
\alias{find_cf}
\title{
Characteristic function estimation with a learning algorithm
}

\description{
Estimates the characteristic function of multivariate data over a grid of points
using a learning (neural network type) algorithm.
}

\usage{
find_cf(x, outrange, N, nu, sg, g)
}

\arguments{
\item{x}{
A numeric matrix with the data. The rows are the observations (n) and the columns
are the variables (d).
}
\item{outrange}{
A vector with two values, the lower and upper limits of the grid at which the
characteristic function is evaluated. The same range is used in every dimension.
}
\item{N}{
A single integer value, the number of grid points (cells) in each dimension. The
characteristic function is evaluated at \eqn{N^d} points in total.
}
\item{nu}{
A single numeric value, the learning rate for the weights of the networks.
}
\item{sg}{
A single positive numeric value, the standard deviation of the Gaussian kernel.
}
\item{g}{
A single numeric value, the learning rate used when updating the estimate of the
characteristic function with each observation.
}
}

\details{
This function finds the characteristic function of the given data \code{x}. The
characteristic function is evaluated on a regular grid of \eqn{N^d} points, with \code{N}
equally spaced values between \code{outrange[1]} and \code{outrange[2]} in each of the \eqn{d}
dimensions. The grid is created with \code{\link{createalfa}}. The parameters \code{nu},
\code{sg} and \code{g} are the learning parameters of the algorithm.

Two networks are used, one for the real part and one for the imaginary part of the
characteristic function. Their values are initialised randomly and then updated
once for each observation (row of \code{x}), using the cosine and sine of the inner
product between the grid points and the observation. A Gaussian kernel with standard deviation
\code{sg}, evaluated between all pairs of grid points, is used in the updating of the weights.

Since the starting values are random, the results depend on the random number
generator, so use \code{set.seed} for reproducibility. The kernel matrix has
\eqn{N^d \times N^d} elements, hence the computational cost and memory requirements grow
very fast with \code{N} and with the number of variables.
}

\value{
A list including:
\item{Teta}{
The real part of the characteristic function, evaluated at the grid points
(one row for each row of \code{Alfa}).
}
\item{TetaI}{
The imaginary part of the characteristic function, evaluated at the grid points.
}
\item{Alfa}{
The domain of the output generated, i.e. a matrix with \eqn{N^d} rows and \eqn{d} columns
containing the grid points at which the characteristic function was evaluated.
}
}

\references{
Alzeley O. (2020). On a Statistical Learning Algorithm of Superstatistics.
Computational and Mathematical Methods, 2(3), e1085.
}

\author{
Original Matlab code by Omar Alzeley.

R implementation and documentation: Manos Papadakis and Michail Tsagris.
}

\seealso{
\code{\link{main}}, \code{\link{createalfa}}
}

\examples{
set.seed(1)
## univariate data
x <- matrix( rnorm(100), ncol = 1 )
a <- find_cf(x, outrange = c(-2, 2), N = 11, nu = 0.1, sg = 0.5, g = 0.1)
cbind(a$Alfa, a$Teta, a$TetaI)

## bivariate data
y <- matrix( rnorm(100 * 2), ncol = 2 )
b <- find_cf(y, outrange = c(-1, 1), N = 5, nu = 0.1, sg = 0.5, g = 0.1)
dim(b$Alfa)  ## 25 rows, 2 columns
}
