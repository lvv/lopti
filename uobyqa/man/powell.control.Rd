\name{powell.control}
\alias{powell.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{powell.control}
\description{Controls convergence in \code{powell}.}
\usage{powell.control(trace = 0, rhobeg = log(100)/10, rhoend = 1e-04,
                      maxit = 1e+05, fnscale = 1, parscale = 1, ...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trace}{an integer that determines the level of information
    produced at each step of the optimizer. The value of \code{trace}
    should be set to 0, 1, 2 or 3, which controls the amount of
    printing. Specifically, there is no output if \code{trace=0} and
    there is output only at the return if \code{trace=1}. Otherwise,
    each new value of RHO is printed, with the best vector of variables
    so far and the corresponding value of the objective
    function. Further, each new value of F with its variables are output
    if \code{trace=3}.}  
  \item{rhobeg}{About one tenth of the greatest expected change to a variable.}
  \item{rhoend}{Indicates the accuracy that is required in the final
    values of the variables.}
  \item{maxit}{The upper bound on the number of times \code{fn} is to
    be called.}
  \item{fnscale}{A scalar to divide the function value
    by. \code{fnscale} less than zero implies a maximization rather
    than a minimization.}
  \item{parscale}{A vector to divide the parameter values by. If only
    one value is given it is replicated to the length of \code{par}.}
  \item{\dots}{Additional arguments that will be ignored!}
}
\details{\code{rhobeg} and \code{rhoend} must be set to the initial and
  final values of a trust region radius, so both must be positive with
  \code{rhoend<=rhobeg}.}
\author{Sundar Dorai-Raj}

\seealso{\code{\link{powell}}}
\keyword{optimize}
