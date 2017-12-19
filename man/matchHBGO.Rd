\name{matchHBGO}
\alias{matchHBGO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
matchHBGO(sampleGenes, GOref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleGenes}{
%%     ~~Describe \code{sampleGenes} here~~
}
  \item{GOref}{
%%     ~~Describe \code{GOref} here~~
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
function (sampleGenes, GOref) 
{
    library(reshape2)
    library(dplyr)
    col_num = ncol(sampleGenes)
    data_long <- melt(sampleGenes, id.vars = 0, variable.name = "geneID", 
        value.name = "GO_ID")
    colnames(GOref)[1] <- "GO_ID"
    merged_data <- merge(data_long, GOref, by = "GO_ID")
    merged_data <- arrange(merged_data, geneID)
    colnames(merged_data)[3] <- "function"
    colnames(merged_data)[4] <- "differential_regulation"
    return(merged_data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
