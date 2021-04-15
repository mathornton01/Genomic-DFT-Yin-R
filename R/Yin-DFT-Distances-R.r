# Author: Micah A. Thornton 
# Date: April 2021 
# Summary: Translation of much of the DFT Distance Framework into an R-Package



#' Encode a Single Genome
#' 
#' Encodes a single genomic signal accepted as a character string into either 
#' a four by signal length, or two by signal length matrix depending on whether
#' the user specifies 2D or 4D, that contains the encoded genomic signal. 
#' @param stringGenome The genomic signal the user would like to work with 
#' @param dimension Either '2D' or '4D' character vector, '4D' by default, this 
#' specifies whether to use the binary four-signal encoding originally described
#' in 2014 by Yin and Yau, or three-vauled two-signal encoding. 
#' @param strategy This is only used for the 2D encoding, it determines which of 
#' the possible 2D encodings is used, this is determined for ensembles by 
#' considering the concentration of nucleotides which are Adonine and Cytosine or
#' Adonine and Guanine. 
#' @return Returns a matrix of dimension 2xSignalLength, or 4xSignalLength where 
#' each row is determined by the encoding strategy, for 4D encoding, each row 
#' represents the presence or absence of a specific nucleotide, for the 2D 
#' encoding, the two signals of three values -1, 0, 1 are determined according 
#' to the strategy choosen. 
#' @examples 
#' EncodedSignal2D <- encodeGenome('ACCACTTGAAGAGAGCCCGGGAT', '4D'); 
#' EncodedSignal4DAC <- encodeGenome('ACCACTTGAAGAGACCCGGGAT', '2D', 'AC'); 
#' EncodedSignal4DAG <- encodeGenome('ACCACTTGAAGAGACCCGGGAT','2D','AG'); 
#' @export
encodeGenome <- function(stringGenome,dimension='4D',strategy="AC"){
  if (dimension == "4D"){
    strUpper <- toupper(stringGenome);
    encA <- gsub('A','1', gsub('[^A]', '0',strUpper));
    encC <- gsub('C','1', gsub('[^C]', '0',strUpper)); 
    encG <- gsub('G','1', gsub('[^G]', '0',strUpper)); 
    encT <- gsub('T','1', gsub('[^T]', '0',strUpper)); 
    encodedSignal <- rbind(as.numeric(unlist(strsplit(encA,''))),
                           as.numeric(unlist(strsplit(encC,''))),
                           as.numeric(unlist(strsplit(encG,''))),
                           as.numeric(unlist(strsplit(encT,''))));
    return(encodedSignal); 
  }
  if (dimension == "2D"){
    encodedSignal <- matrix(0,nrow=2,ncol=nchar(stringGenome))
    Alocs <- as.vector(gregexpr('A',toupper(stringGenome))[[1]]);
    Clocs <- as.vector(gregexpr('C',toupper(stringGenome))[[1]]); 
    Glocs <- as.vector(gregexpr('G',toupper(stringGenome))[[1]]); 
    Tlocs <- as.vector(gregexpr('T',toupper(stringGenome))[[1]]); 
    if (strategy == "AC"){
      encodedSignal[1,Alocs] <- 0; 
      encodedSignal[2,Alocs] <- -1;
      encodedSignal[1,Clocs] <- 0; 
      encodedSignal[2,Clocs] <- 1; 
      encodedSignal[1,Glocs] <- 1; 
      encodedSignal[2,Glocs] <- 0; 
      encodedSignal[1,Tlocs] <- -1; 
      encodedSignal[2,Tlocs] <- 0; 
    }
    if (strategy == "AG"){
      encodedSignal[1,Alocs] <- 0; 
      encodedSignal[2,Alocs] <- -1;
      encodedSignal[1,Clocs] <- 1; 
      encodedSignal[2,Clocs] <- 0; 
      encodedSignal[1,Glocs] <- 0; 
      encodedSignal[2,Glocs] <- 1; 
      encodedSignal[1,Tlocs] <- -1; 
      encodedSignal[2,Tlocs] <- 0; 
    }
    return(encodedSignal)
  }
}

#' Determine Optimal 2D encoding
#' 
#' Uses the ensemble of strings provided by the user as a list of genomic 
#' signals and determines the optimal encoding strategy by the algorithim 
#' proposed in Yin and Yau 2015, where the distance from A,C concentration 
#' and A,G concentration to one-half are used to determine which encoding will 
#' provide the best balance among all the possible encodings. 
#' @param stringGenomes A list of strings representing genomic signals 
#' @return A character vector containing the name of the optimal strategy, either
#' 'AC' or 'AG'. 
#' @examples
#' getOptimal2DStrategy(list('ACCCTTAACCAAGGAGGGAGAGTTTCCCCGGGGGAGG',
#'                            'CCCCCCACACAAACCCTTGGGGAAAACCCGGAAGGCCCCCC'));
#' @export
getOptimal2DStrategy <- function(stringGenomes){
  total <- 0; 
  AC <- 0; 
  AG <- 0; 
  for (stringGenome in stringGenomes){
    total <- total + nchar(stringGenome); 
    AC <- AC + sum(table(strsplit(toupper(stringGenome),''))[c('A','C')]);
    AG <- AG + sum(table(strsplit(toupper(stringGenome),''))[c('A','G')]);
  }
  Rac <- abs(AC/total - 0.5); 
  Rag <- abs(AG/total - 0.5); 
  if (Rac <= Rag) {return('AC')}
  else {return('AG')}
} 

#' Encode an ensemble of genomes into a numerical form
#' 
#' Will produce a list of either the 4D or the 2D representation of a set of 
#' signals and return them in a list. 
#' @param stringGenomes is a list containing the genomes in a string format 
#' @param dimension is a character string either '2D' or '4D'.
#' @return A list containing matrixes of various column length, but either 2 or 
#' 4 rows. 
#' @examples
#' encodeGenomes(list('ACCCCAATTAGAGGGACTTTGGGAACCGGAGAT','GGCCCAGAGGGAAACCGGT'),
#'               '2D'); 
#' encodeGenomes(list('ACCCCAATTAGAGGGACTTTGGGAACCGGAGAT','GGCCCAGAGGGAAACCGGT'), 
#'               '4D'); 
#' @export
encodeGenomes <- function(stringGenomes,dimension='2D') {
  if (dimension=='2D'){
  strat <- getOptimal2DStrategy(stringGenomes);
  }
  encodedSignals <- list();
  return(lapply(stringGenomes, function(x){
    return(encodeGenome(x,dimension,strat));
  }))
}

tstForMultiGenomes <- function(numStrings = 100, avLength, deviation){
  genomeStrings <- list(numStrings); 
  lengths <- round(rnorm(numStrings,avLength,deviation)); 
  for (i in 1:numStrings){
    genomeStrings[i] <- paste(sample(c('A','C','G','T'),numStrings,replace=T),sep='',collapse='')
  }
  return(encodeGenomes(genomeStrings));
}