\name{SpliceNet}
\alias{SpliceNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SpliceNet: Differential Gene Co-expression Netwok from Level 3 RNASeq data at Splicing Isoform Resolution 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SpliceNet(IP, OP, A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IP}{
%%     ~~Describe \code{IP} here~~
}
  \item{OP}{
%%     ~~Describe \code{OP} here~~
}
  \item{A}{
%%     ~~Describe \code{A} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (IP, OP, A) 
{
    input <- IP
    output <- OP
    outpval <- paste(output, "Pval.txt", sep = "-")
    alpha <- as.numeric(A)
    library(pracma)
    library(MSBVAR)
    library(tseries)
    library(base)
    IT <- read.table(input, header = FALSE, sep = "\t")
    IM <- as.matrix(IT)
    IM <- IM[, -2]
    ID <- IM[, 1]
    IDU <- unique(ID)
    IM <- IM[, -1]
    D <- length(ID)
    RM <- matrix(data = NA, nrow = (length(IDU)), ncol = (length(IDU)))
    PM <- matrix(data = NA, nrow = length(IDU), ncol = length(IDU))
    for (i in 1:length(IDU)) {
        RM[i, i] = 1
        PM[i, i] = 0
    }
    sp1 = 1
    ep1 = 1
    ii = 0
    for (i in sp1:D) {
        if ((ID[i] != ID[i + 1]) || (i == D)) {
            ii = ii + 1
            ep1 = i
            X1 <- IM[sp1:ep1, ]
            X1 <- t(X1)
            if (sp1 == ep1) {
                X1 <- t(X1)
            }
            sp2 = 1
            ep2 = 1
            {
                jj = 0
                for (j in sp2:D) {
                  if ((ID[j] != ID[j + 1]) || (j == D)) {
                    jj = jj + 1
                    ep2 = j
                    {
                      X2 <- IM[sp2:ep2, ]
                      X2 <- t(X2)
                      if (sp2 == ep2) {
                        X2 <- t(X2)
                      }
                      n1 = nrow(X1)
                      n2 = nrow(X2)
                      p1 = ncol(X1)
                      p2 = ncol(X2)
                      N1 = n1 - 1
                      N2 = n2 - 1
                      X3 <- cbind(X1, X2)
                      X3 <- matrix(X3, nrow = nrow(X3), ncol = ncol(X3))
                      mode(X3) <- "numeric"
                      Sigma = cov(X3)
                      A11 = N1 * Sigma[1:p1, 1:p1]
                      A22 = N1 * Sigma[(p1 + 1):(p1 + p2), (p1 + 
                        1):(p1 + p2)]
                      A12 = N1 * Sigma[1:p1, (p1 + 1):(p1 + p2)]
                      A21 = N1 * Sigma[(p1 + 1):(p1 + p2), 1:p1]
                      A11 = as.matrix(A11)
                      A22 = as.matrix(A22)
                      A12 = matrix(A12, p1, p2)
                      A21 = matrix(A21, p2, p1)
                      M = A21 \%*\% inv(A11) \%*\% A12 \%*\% inv(A22)
                      r1 = p2/p1
                      r2 = p2/(n1 - 1 - p1)
                      hn = sqrt(r1 + r2 - r1 * r2)
                      eg = r2/(r1 + r2)
                      vg = 2 * (hn^2) * (r1^2) * (r2^2)/((r1 + 
                        r2)^4)
                      TS = 0
                      TS = (vg^(-0.5)) * (sum(diag(M)) - p2 * 
                        eg)
                      if ((is.infinite(TS)) | (is.nan(TS)) | 
                        (is.null(TS))) {
                        TS = -10000
                      }
                      result = NULL
                      result[1] = TS
                      pvalue = 1 - pnorm(TS)
                      cvi = 1 - alpha
                      CV = qnorm(cvi)
                      result[2] = CV
                      if (TS >= CV) {
                        result[3] = 1
                      }
                      else {
                        result[3] = 0
                      }
                      RM[ii, jj] <- result[3]
                      PM[ii, jj] <- pvalue
                      sp2 = j + 1
                    }
                  }
                }
            }
            sp1 = i + 1
        }
    }
    RM <- rbind(IDU, RM)
    PM <- rbind(IDU, PM)
    IDU <- c("I-MATRIX", IDU)
    RM <- cbind(IDU, RM)
    PM <- cbind(IDU, PM)
    write(RM, file = output, ncolumns = ncol(RM), append = FALSE, 
        sep = "\t")
    write(PM, file = outpval, ncolumns = ncol(PM), append = FALSE, 
        sep = "\t")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
