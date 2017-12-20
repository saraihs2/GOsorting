\name{RefCSV}
\alias{RefCSV}

\title{
Read in honey bee GO reference file 
}
\description{
Read in the .CSV file containing honey bee GO terms, GO functions, and whether or not the gene was up or down regulated. Reference data file is from the honey bee alarm pheromone experiment. Differential genes were measure in the brain and regulation response is measure as alarm pheromone : control.
}
\usage{
RefCSV(HBrefPath)
}

\arguments{
  \item{HBrefPath}{
Argument file from the honey bee alarm pheromone experiment
}
}

}
\value{
Returns a data frame containing GO terms, functions, and regulation (up or down)
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Alaux et al. 2009
}
\author{
Sarai Stuart
}

\seealso{
\code{\link{matchHBGO}}
\code{\link{GOtermSep}}
}

\examples

# Path to Honey Bee reference file
HBrefPath <- system.file("extdata",
                          "./HB_alarm_GOterms.csv",
                          package = "GOsorting")
                          
# Create the reference table to compare to input data                           
HbGOtable <- RefCSV(HBrefPath)


