\name{powell}
\alias{powell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{powell}
\description{Optimizes a function using Powell's UObyQA algorithm.}
\usage{
powell(par, fn, control = powell.control(), check.hessian = TRUE, ...)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{Starting values for objective function}
  \item{fn}{A function to be optimized. The function takes the
    parameters (\code{par}) as its first argument.}
  \item{control}{A list of control parameters}
  \item{check.hessian}{logical; if \code{TRUE} an eigenvalue
    decomposition is used to check the hessian for positive
    definiteness.}
  \item{\dots}{Additional arguments to be passed to \code{fn}}
}
\details{
  This function seeks the least value of a function of many variables,
  by a trust region method that forms quadratic models by
  interpolation. The algorithm is described in "UOBYQA: unconstrained
  optimization by quadratic approximation" by M.J.D. Powell, Report
  DAMTP 2000/NA14, University of Cambridge.
}
\value{
  A list with components
  \item{par}{The final values of the parameters.}
  \item{value}{The final value of the function being optimized.}
  \item{counts}{The number of times the function is called.}
  \item{hessian}{A symmetric matrix of the estimated Hessian.}
  \item{eigen.hessian}{If \code{check.hessian} is \code{TRUE} the
    eigenvalues and eigenvectors; otherwise \code{NULL}.}
  \item{convergence}{0 if converged, 1 otherwise.}
  \item{control}{The input control parameters.}
  \item{message}{Information about the model fit. This will be non-null
    only if \code{check.hessian} is \code{TRUE} and the resulting
    hessian is not positive definite.}
  \item{call}{The original call to the optimizer.}
}
\references{http://plato.asu.edu/topics/problems/nlounres.html}
\author{Sundar Dorai-Raj}

\seealso{\code{\link{optim}}}
\examples{
set.seed(1)
fn <- function(beta, y, x, w) \{
  # binomial deviance using double log link
  mu <- exp(x \%*\% beta)
  logLik <-  - y * mu + (w - y) * log(1 - exp(-mu))
  -2 * sum(logLik)
\}
n <- 1000
beta <- c(-1, 0.5)
w <- rpois(n, 100)
x <- rep(c("A", "B"), length = n)
X <- model.matrix(~ x, data.frame(x))
y <- rbinom(n, w, exp(-exp(X \%*\% beta)))

x1 <- powell(beta, fn, y = y, x = X, w = w)
x2 <- optim(beta, fn, y = y, x = X, w = w, hessian = TRUE)
x3 <- glm(1 - y/w ~ x, data = data.frame(x, y, w),
          family = binomial("cloglog"), weights = w)

# compare coefficients from 3 fits
rbind(powell = x1$par, optim = x2$par, glm = coef(x3))
# compare standard errors from 3 fits
rbind(powell = sqrt(diag(2 * solve(x1$hessian))), 
      optim = sqrt(diag(2 * solve(x2$hessian))), 
      glm = sqrt(diag(vcov(x3))))
}
\keyword{optimize}% at least one, from doc/KEYWORDS
