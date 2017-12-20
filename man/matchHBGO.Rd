\name{matchHBGO}
\alias{matchHBGO}

\title{
Match my GO terms to honey bee GO terms
}
\description{
Match the GO terms found on my data list to those found to differentially regulated in the honey bee. Matches indicated conserved GO functions between data sets. Matched GO terms are returned in addition to the input gene id, function of GO term, and the differential regualtion that was found in the honey bee. Function requires the installed packages reshape2 and diplyr.   
}
\usage{
matchHBGO(sampleGenes, GOref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleGenes}{
Argument is the data frame of separated GO terms per gene id, can be created from \code{\link{GOtermSep}}
}
  \item{GOref}{
Argument is the reference data frame of GO terms from the honey bee, can be input from \code{\link{RefCSV}}
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
https://github.com/hadley/reshape
https://github.com/tidyverse/dplyr
}
\author{
Sarai
}


\seealso{
\code{\link{GOtermSep}}
\code{\link{RefCSV}}

}
\examples

# Path my GO term test file
myGOterm <- system.file("extdata",
                          "./myGOterm.csv",
                          package = "GOsorting")
                          
# Path reference honey bee file
HbGOtable <- system.file("extdata",
                          "./HbGOtable.csv",
                          package = "GOsorting")
                          
# Create a data frame of matched GO terms between my data and the honey bee                           
matchHBGO(myGOterm,HbGOtable)



