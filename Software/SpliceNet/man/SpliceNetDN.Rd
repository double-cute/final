\name{SpliceNetDN}
\alias{SpliceNetDN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SpliceNet: Differential Gene Co-expression Netwok from Level 3 RNASeq data at Splicing Isoform Resolution
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SpliceNetDN(NN, CN)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NN}{
%%     ~~Describe \code{NN} here~~
}
  \item{CN}{
%%     ~~Describe \code{CN} here~~
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
function (NN, CN) 
{
    input1 <- NN
    output2 <- CN
    IT1 <- read.table(input1, header = FALSE, sep = "\t")
    IM1 <- as.matrix(IT1)
    ID1 <- IM1[, 1]
    IM1 <- IM1[, -1]
    IM1 <- IM1[-1, ]
    ID1 <- ID1[-1]
    IT2 <- read.table(input2, header = FALSE, sep = "\t")
    IM2 <- as.matrix(IT2)
    ID2 <- IM2[, 1]
    IM2 <- IM2[, -1]
    IM2 <- IM2[-1, ]
    ID2 <- ID2[-1]
    D <- (length(ID1))
    key <- c("", "")
    temp <- c("")
    for (i in 1:D) {
        for (j in 1:D) {
            if ((is.na(IM1[i, j])) && (is.na(IM2[i, j]))) {
            }
            else {
                if (IM1[i, j] != IM2[i, j]) {
                  key <- sort(c(ID1[i], ID1[j]))
                  val <- paste(key[1], "\t", key[2], "\t", IM1[i, 
                    j], "\t", IM2[i, j], "\n", sep = "")
                  temp <- c(temp, val)
                }
            }
        }
    }
    temp <- unique(temp)
    for (i in 1:length(temp)) {
        cat(temp[i])
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
