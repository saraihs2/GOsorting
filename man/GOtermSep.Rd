\name{GOtermSep}
\alias{GOtermSep}

\title{
Separate GO terms for each Gene 
}
\description{
This function separates the single string of mutiple GO terms into a vectorized list of separate GO terms per gene. The list is then stored as a matrix so that all GO terms are characterized for each gene. 
}
\usage{
GOtermSep(TrimmedData)
}

\arguments{
  \item{TrimmedData}{
Trimmed data output from readTrinOutput function
}
}

\value{
Returns a data frame of GO terms listed per gene id. 

Output is argument 1 for matchHBGO function
}

\author{
Sarai Stuart
}

\seealso{
\code{\link{readTrinOutput}}
\code{\link{matchHBGO}}
}

\examples

#Path to test file
path <- system.file("extdata", 
                       "New_GO_test_file.csv",
                       package = "GOsorting")
                       
test.csv <- read.csv(path)

#Insert "path" argument into function
test <- readTrinOutput(test.csv)


#Insert "myTrimData" argument into function
GOtermSep(test)


