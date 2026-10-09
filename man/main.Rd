\name{main}
\alias{main}
\title{
Characteristic function estimation with a learning algorithm: simulation example
}

\description{
Generates data from a Gaussian or a uniform distribution, estimates their
characteristic function with a learning algorithm and compares it with the
theoretical characteristic function.
}

\usage{
main(N = 20, d = 2, n = 500, nu = 1.8, sg = 0.25, g = 0.00036, dist = 1)
}

\arguments{
\item{N}{
Number of grid points (resolution) in each dimension at which the characteristic
function is evaluated.
}
\item{d}{
Number of dimensions of the generated data. Plots are produced only for \code{d = 1}
and \code{d = 2}.
}
\item{n}{
Number of data to be generated.
}
\item{nu}{
The value of the learning parameter (learning rate of the weights of the networks),
see \code{\link{find_cf}}.
}
\item{sg}{
The value of the learning parameter (standard deviation of the Gaussian kernel),
see \code{\link{find_cf}}.
}
\item{g}{
The value of the learning parameter (learning rate for updating the estimate of the
characteristic function), see \code{\link{find_cf}}.
}
\item{dist}{
The data distribution, 0 is for Gaussian and 1 is for uniform.
}
}

\details{
If \code{dist = 0}, the data are generated from a \code{d}-dimensional Gaussian distribution
with mean vector \eqn{(1,\ldots,1)} and covariance matrix \eqn{2 I_d}, and the characteristic
function is evaluated over the range \eqn{[-5, 5]} in each dimension. Otherwise, the data are
generated from the uniform distribution on \eqn{[-2, 3]} in each dimension, and the characteristic
function is evaluated over the range \eqn{[-3, 4]}. The characteristic function is then estimated
with \code{\link{find_cf}} and compared with the theoretical characteristic function of the
distribution that generated the data. The theoretical characteristic function is used to
compute the error of the estimated one.

Three graphics devices are opened (via \code{dev.new}). They contain the generated data,
the real and the imaginary part of the estimated characteristic function, and the real and
imaginary part of the theoretical characteristic function together with the errors (differences
between the theoretical and the estimated values).

The following is the theory behind the use of uniform distributions in superstatistics.
Let \eqn{U[0,N]} be a uniform variable on the interval \eqn{[0,N]}.
Observe that
\deqn{\mathrm{E}\, U[0,N]=\frac{1}{N}\int_0^N x \, dx = \frac{N}{2}.}
Write
\deqn{
\begin{array}{rcl}
\mathrm{E}\, e^{-hU[0,N]} & = & \displaystyle \frac{1}{N}\int_0^N e^{-hx} dx =
e^{-hN/2}\frac{1}{N}\int_0^N e^{-h(x-N/2)} dx \\
& = & \displaystyle e^{-hN/2}\int_0^N e^{-hN(x/N-1/2)} \frac{dx}{N}.
\end{array}
}
Now, by considering a new variable of integration \eqn{t=\frac{x-N/2}{N}},
\deqn{\mathrm{E}\, e^{-hU[0,N]}=e^{-hN/2}\int_{-1/2}^{1/2}e^{-hN t} dt.}
Assume that \eqn{c=hN} is small. Then, by Taylor expansion,
\eqn{e^{-c t}=1-c t+ c^2t^2/2-\cdots}
and so
\deqn{
\mathrm{E}\, e^{-hU[0,N]}=e^{-hN/2}\int_{-1/2}^{1/2}
\left(1-c t+ \frac{c^2t^2}{2}-\cdots\right)dt = e^{-hN/2}\left(1 +\frac{c^2}{24} +\cdots\right).
}
Consider now the uniform distribution case, \eqn{U[a,a+b]}.

Let the random variable \eqn{\beta} be distributed uniformly in \eqn{[a, a + b]} for some
\eqn{a \geq 0}. Then, its probability density function reads as
\deqn{
f(\beta)=
\left\lbrace
\begin{array}{ll}
\displaystyle \frac{1}{b}, & \mbox{if } a < \beta < a+b, \\
0, & \mbox{otherwise.}
\end{array}
\right.
}
The mean and the variance of this distribution are therefore
\deqn{
\beta_0 = a + \frac{b}{2}, \qquad \sigma^2 = \frac{b^2}{12}. \qquad (1)
}

The superstatistics for this more general case can be found by the same argument.
First, one should note that the generalised Boltzmann factor is defined by
Beck and Cohen (2003) as
\deqn{B = \int_0^\infty f(\beta) \exp(-\beta E) d\beta,}
direct integration of which gives
\deqn{
\int_0^\infty f(\beta) \exp(-\beta E) d\beta = \int_a^{a + b} \frac{1}{b} \exp(-\beta E) d\beta =
\frac{\exp(-a E)}{b E} \left[1 - \exp(-b E)\right].
}

Substituting \eqn{a} from (1) into the final expression of the integral, we have
\deqn{
B = \exp(-E \beta_0) \frac{\exp\left(E \frac{b}{2}\right) - \exp\left(-E \frac{b}{2}\right)}{b E}.
}

Taking into account the following series expansions for small \eqn{bE}:
\deqn{
\frac{\exp\left(\pm E \frac{b}{2}\right)}{b E} \approx \frac{1}{bE} \pm \frac{1}{2} +
\frac{bE}{8} \pm \frac{(bE)^2}{48} + \cdots
}
we have
\deqn{
\frac{\exp\left(E \frac{b}{2}\right) - \exp\left(- E \frac{b}{2}\right)}{b E} \approx
1 + \frac{(bE)^2}{24} + \frac{(bE)^4}{1920} + \cdots
}
Substituting this into the expression for \eqn{B}, we arrive at the expansion
\deqn{
B \approx \exp(-\beta_0 E) \left[1 + \frac{(b E)^2}{24} + \frac{(b E)^4}{1920} + \cdots\right].
}

Let \eqn{x} be a discrete random variable with a discrete uniform distribution with parameter
\eqn{N}, i.e. \eqn{x \in \lbrace a + 1, a + 2, \ldots, a + N \rbrace}. Then, the probability
mass function is given by \eqn{f(x) = \frac{1}{N}} and the mean is given by
\deqn{
\begin{array}{rcl}
x_0 & = & \displaystyle \sum x f(x) = \frac{1}{N}\sum_{x = a + 1}^{a + N} x =
\frac{1}{N}\left( (a + 1) N + \sum_{j = 0}^{N - 1} j\right) \\
& = & \displaystyle \frac{1}{N}\left( (a + 1)N + \frac{N (N - 1)}{2}\right) =
a + 1 + \frac{N - 1}{2} = a + \frac{N + 1}{2}.
\end{array}
}

Let \eqn{I_{z, N} = \frac{1}{N}\sum_{j = a + 1}^{a + N} z^{j}}. Then, taking into account that
\eqn{a = x_0 - \frac{N + 1}{2}}, we have
\deqn{
\begin{array}{rcl}
I_{z, N} & = & \displaystyle \frac{z^{a + 1} + z^{a + 2} + \cdots + z^{a + N}}{N} \\
& = & \displaystyle \frac{z^{x_0 - \frac{N + 1}{2} + 1} + z^{x_0 - \frac{N + 1}{2} + 2} +
\cdots + z^{x_0 - \frac{N + 1}{2} + N}}{N} \\
& = & \displaystyle z^{x_0} \frac{z^{- \frac{N + 1}{2}} \left(z + z^2 + \cdots + z^N\right)}{N} \\
& = & \displaystyle z^{x_0} \frac{z^{- \frac{N + 1}{2}} z \left(z^N - 1\right)}{N \left(z - 1\right)} \\
& = & \displaystyle \frac{z^{x_0}}{z - 1} \frac{z^{\frac{1 + N}{2}} - z^{\frac{1 - N}{2}}}{N}.
\end{array}
}

Note that we can also use the exponential form \eqn{z^{\frac{N}{2}} = \exp\left[\frac{N}{2} \ln z\right]}
to write
\deqn{
I_{z, N} = \frac{z^{x_0 + \frac{1}{2}}}{z - 1}
\frac{\exp\left[\frac{N}{2} \ln z\right] - \exp\left[-\frac{N}{2} \ln z\right]}{N}.
}

Taking into account the expansions above, we obtain
\deqn{
\begin{array}{rcl}
\exp\left[\frac{N}{2} \ln z\right] - \exp\left[-\frac{N}{2} \ln z\right] & = &
\displaystyle \sum_{n = 0}^\infty \frac{1}{n!} \left[\left(\frac{N}{2}\right)^n \ln^n z -
\left(-\frac{N}{2}\right)^n \ln^n z \right] \\
& = & \displaystyle \sum_{n = 0}^\infty \frac{1}{n!} \left[\left(\frac{N}{2}\right)^n -
\left(-1\right)^n \left(\frac{N}{2}\right)^n \right] \ln^n z.
\end{array}
}

Evidently,
\deqn{
\left(\frac{N}{2}\right)^n - \left(-1\right)^n \left(\frac{N}{2}\right)^n =
\left\lbrace
\begin{array}{ll}
\displaystyle 2\left(\frac{N}{2}\right)^{2k + 1}, & n = 2k + 1, \\
0, & n = 2k.
\end{array}
\right.
}

Therefore, finally, we arrive at
\deqn{
I_{z, N} = \frac{2z^{x_0 + \frac{1}{2}}}{N\left(z - 1\right)} \sum_{k = 0}^\infty
\frac{1}{\left(2k + 1\right)!} \left(\frac{N}{2}\right)^{2k + 1} \ln^{2k + 1} z.
}

The Boltzmann factor for discrete uniform distributions reads as
\deqn{B = \frac{1}{N} \sum_{j = a + 1}^{a + N} \exp(-j E).}
Then, taking into account that \eqn{a = x_0 - \frac{N + 1}{2}}, we have
\deqn{
\begin{array}{rcl}
B & = & \displaystyle \frac{\exp\left[\left(\frac{1}{2} - x_0\right) E\right]}{\exp(E) - 1}
\frac{\exp\left(\frac{1}{2} N E\right) - \exp\left(-\frac{1}{2} N E\right)}{N} \\
& = & \displaystyle \frac{2}{N} \frac{\exp\left[\left(\frac{1}{2} - x_0\right) E\right]}{\exp(E) - 1}
\sinh\left(\frac{1}{2} N E\right).
\end{array}
}

Taking into account that
\eqn{\sinh(x) = x + \frac{x^3}{3!} + \cdots = \sum_{n = 0}^\infty \frac{x^{2 n + 1}}{(2n + 1)!}},
we obtain
\deqn{
\begin{array}{rcl}
B & = & \displaystyle \frac{\exp\left[\left(\frac{1}{2} - x_0\right) E\right]}{\exp(E) - 1} \cdot
\frac{2}{N} \left[ \frac{1}{2} N E + \frac{(N E)^3}{2^3 \cdot 3!} +
\frac{(N E)^5}{2^5 \cdot 5!} + \cdots \right] \\
& = & \displaystyle \frac{E \exp\left[\left(\frac{1}{2} - x_0\right) E\right]}{\exp(E) - 1}
\left[ 1 + \frac{(N E)^2}{24} + \frac{(N E)^4}{1920} + \cdots \right].
\end{array}
}
}

\value{
A named numeric vector with two values: the error of the real part (\code{"Re error"}) and the
error of the imaginary part (\code{"Im error"}) of the estimated characteristic function. Each
error is the sum of the squared differences between the estimated and the theoretical values
over the grid points.

As a side effect, three plots (graphics devices) are produced when \code{d} is 1 or 2, see the
details.
}

\references{
Alzeley O. (2020). On a Statistical Learning Algorithm of Superstatistics.
Computational and Mathematical Methods, 2(3), e1085.

Beck C. and Cohen E.G.D. (2003). Superstatistics. Physica A, 322: 267--275.
}

\author{
Original Matlab code by Omar Alzeley.

R implementation and documentation: Manos Papadakis and Michail Tsagris.
}

\seealso{
\code{\link{find_cf}}, \code{\link{createalfa}}
}

\examples{
## bivariate uniform data
main(N = 10, d = 2, n = 100)

## univariate Gaussian data
main(N = 10, d = 1, n = 100, dist = 0)
}


