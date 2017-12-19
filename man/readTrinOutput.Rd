\name{readTrinOutput}
\alias{readTrinOutput}

\title{
Read, process, and trim Trinotate output 
}
\description{
The function processes differentially expressed RNA data that has been analyzed for Gene Ontology enrichment. Function reads in data, processes out data with unidentified genes, and trims data columns needed for GO matching. Function is compatible with the the .CSV-converted output file format from Trinotate:Transcriptome Functional Annotation and Analysis suite.
}
\usage{
readTrinOutput(path)
}

\arguments{
  \item{path}{
User provided path to .CSV Trinotate-formated file 
}
}

\value{
Returns a data frame with the columns of the following headers : gene_id, transcript_id, gene_ontology_blast, gene_ontology_pfam

Rows without a gene_id are filtered out

Output is the argument for GOtermSep Function
}
\references{
https://trinotate.github.io/
}
\author{
Sarai Stuart
}

\seealso{
\code{\link{GOtermSep}}
}
\examples{
# Path to test file
 path <- system.file("extdata",
                       "New_GO_test_file.csv",
                       package = "GOsorting")

# Insert "path" argument into function
readTrinOutput(path)
