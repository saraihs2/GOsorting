###!!!!must convert Trinotate output Excel doc into a .csv file input to GOsorting.R!!!####


####FUNCTION 1: read in and trim data from Trinotate output ####

#read and extract data from CSV file
#define read function to import the data file
readTrinOutput <- function(path){
  
  #all data from Trinotate output
  fullTable <- read.csv(path,
                        
                        #retain header lables from csv file
                       header = TRUE,
                       
                       #separate strings by commas
                       sep = ",",
                       
                       #lable missing data: convert "." to "NA" in R environment
                       na.strings = ".",
                       
                       #column headers
                       col.names = c("gene_id",
                                     "transcript_id",
                                     "sprot_Top_BLASTX_hit",
                                     "RNAMMER",
                                     "prot_id",
                                     "prot_coords",
                                     "sprot_Top_BLASTP_hit",
                                     "custom_pombe_pep_BLASTX",
                                     "custom_pombe_pep_BLASTP",
                                     "Pfam",
                                     "SignalP",
                                     "TmHMM",
                                     "eggnog",
                                     "Kegg",
                                     "gene_ontology_blast",
                                     "gene_ontology_pfam",
                                     "transcript",
                                     "peptide"),
                       
                       #import all data as strings
                       stringsAsFactors = FALSE)
 
  #select columns of interest from fullTable and put output to a data frame called dataTrim 
  #while retaining headers
  dataTrim <- data.frame(fullTable$gene_id,
                         fullTable$transcript_id,
                         fullTable$gene_ontology_blast,
                         fullTable$gene_ontology_pfam, 
                         
                         #keep data as strings during data trim
                         stringsAsFactors = FALSE)
  
  #covert trimmed dataframe headers to original headers from Trinotate output file 
  #i.e. "fullTable$gene_id" will simply be "gene_id" 
  colnames(dataTrim) <- c("gene_id",
                        "transcript_id",
                        "gene_ontology_blast",
                        "gene_ontology_pfam")
  
  #filter out rows of trimmed data in which no gene id is available
  myTrimData <- subset(dataTrim,(!is.na(dataTrim$gene_id)))
  
  return(myTrimData)  
  
}



####FUNCTION 2: Separate individual GO terms listed per gene ####

#separate individual GO terms listed for each gene in the column gene_ontology_blast
GOtermSep <- function(TrimmedData){
  
  #size for preallocation of memory is the length of gene_id column because we are matching 
  #gene id and it's associated GO terms
  nrows <- length(TrimmedData$gene_id)
  
  #choose the GO:####### pattern retain only this pattern (remove biological process text, etc.)
  pattern <- 'GO:[:digit:]{7}'
  
  #create a variable called "myList" with which to populate a list of GO terms
  myList <- list()
  
  #make the length of "myList" that is equal to the length of "nrows"; gene id list will be equivalent in 
  #length to list of GO terms for each gene id
  length(myList) <- length(nrows)
  
  #create list of the gene ids
  #gene ids are not included if not identified by gene name
  geneNames <- character(nrows) 
  
  #there is now an output list for all of the identified gene ids
  
  #loop through the gene_ontology_blast to identify and retain only the GO terms 
  for(i in 1:nrows){
    
    #move through gene id column and output names into geneNames list
    geneNames[i] <- TrimmedData$gene_id[i]
    
    #move through column that contains the GO terms we want to parse output items in list
    GOblastCol <- myTrimData$gene_ontology_blast[i] 
    
    #output a vector of strings that matches the GO:####### pattern from GOblastCol
    GOtermVec <- as.vector(c(stringr::str_match_all(GOblastCol, pattern)[[1]]))
    
    #populated myList with the vector of GO terms listed for each gene
    myList[[i]] <- GOtermVec
  }
  names(myList) <- geneNames
  return(as.data.frame(myList))
}



####FUNCTION 3: read in the reference honey bee GO file from Alaux et al. 2009 ####

#read in modified csv file  
RefCSV <- function(HBrefPath){
  
  HbGOtable <- read.csv(HBrefPath,
                      
                      #keep column headers
                      header = TRUE,
                      
                      #separate strings at the commas
                      sep = ",",
                      
                      #convert missing data (.) to NAs
                      na.strings = ".",
                      
                      #keep all strings as strings
                      stringsAsFactors = FALSE,
                      
                      #column headers
                      col.names = c("GO #",
                                    "GO name",
                                    "change in expression: alarm pher : control"))
  
  return(HbGOtable)
}

HBrefPath <- "./HB_alarm_GOterms.csv"
HbGOtable <- RefCSV(HBrefPath)

####FUNCTION 4: match my GO terms with the differentially enriched GO terms in the brain of the
#honey be defensive response, and subset the list of terms ####

#define function with two arguments, i.e. both my data matrix and Hb matrix must be specified in 
#order to match terms present in both sets of data
matchHBGO<- function(sampleGenes, GOref){
  
  library(reshape2)
  library(dplyr)

  #create variable that is the length of my sample data
  col_num=ncol(sampleGenes)
  
  #transform my data format from wide to long for comparison
  data_long <- melt(sampleGenes, id.vars =0 , variable.name = "geneID", value.name="GO_ID")
  
  #label column 1 of honey bee data set GO_ID
  colnames(GOref)[1] <- "GO_ID"
  
  #merge my data set and honey bee data set based on the GO IDs that match
  merged_data<-merge(data_long, GOref, by="GO_ID")
  
  #sort merged data frame using the geneID header
  merged_data<-arrange(merged_data, geneID)
  
  #add the header "GO function" to column 3
  colnames(merged_data)[3]<-"function"
  
  #add the header "differential regulation" to column 4 
  colnames(merged_data)[4]<-"differential_regulation"
  
  return(merged_data)
}













