\name{createalfa}
\alias{createalfa}
\title{
Grid of all combinations of parameter values
}

\description{
Creates a matrix with all possible combinations of the candidate values of several parameters.
}

\usage{
createalfa(alfax)
}

\arguments{
\item{alfax}{
A numeric matrix with one row per parameter. The i-th row contains the candidate
values of the i-th parameter. All rows must have the same number of columns.
}
}

\details{
The function creates all possible combinations (the Cartesian product, or full grid)
of the values stored in the rows of \code{alfax}, by picking exactly one value from
each row. It works recursively: if \code{alfax} has a single row, the transposed
row is returned; otherwise, each value of the first row is repeated and bound
(\code{cbind}) with the combinations obtained from the remaining rows.

If \code{alfax} has \eqn{p} rows and \eqn{m} columns, the output has \eqn{m^p} rows and
\eqn{p} columns. The first column varies slowest and the last column varies fastest.
Note that the number of rows grows exponentially with the number of rows of
\code{alfax}, so large inputs can require a lot of memory.
}

\value{
A matrix with \code{ncol(alfax)^nrow(alfax)} rows and \code{nrow(alfax)} columns. Each
row of the matrix is one combination of values, with the j-th element taken from
the j-th row of \code{alfax}.
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
\code{\link{main}}, \code{\link{find_cf}}
}

\examples{
## two parameters with three candidate values each
x <- rbind( c(1, 2, 3), c(10, 20, 30) )
createalfa(x)

## three parameters with two candidate values each
y <- rbind( c(0.1, 0.5), c(1, 2), c(5, 6) )
dim( createalfa(y) )  ## 8 rows, 3 columns

## a single row is simply transposed
createalfa( rbind(1:3) )
}
