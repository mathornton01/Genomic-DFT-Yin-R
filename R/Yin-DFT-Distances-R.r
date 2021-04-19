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
      #encodedSignal[1,Alocs] <- 0; 
      encodedSignal[2,Alocs] <- -1;
      #encodedSignal[1,Clocs] <- 0; 
      encodedSignal[2,Clocs] <- 1; 
      encodedSignal[1,Glocs] <- 1; 
      #encodedSignal[2,Glocs] <- 0; 
      encodedSignal[1,Tlocs] <- -1; 
      #encodedSignal[2,Tlocs] <- 0; 
    }
    if (strategy == "AG"){
      #encodedSignal[1,Alocs] <- 0; 
      encodedSignal[2,Alocs] <- -1;
      encodedSignal[1,Clocs] <- 1; 
      #encodedSignal[2,Clocs] <- 0; 
      #encodedSignal[1,Glocs] <- 0; 
      encodedSignal[2,Glocs] <- 1; 
      encodedSignal[1,Tlocs] <- -1; 
      #encodedSignal[2,Tlocs] <- 0; 
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

#' Get the Fourier Power Spectrum for a single encoded Genomic Signal 
#' 
#' This function produces the power spectrum of the Fourier transform for a 
#' single genomic signal that has been encoded using either the 2D or 4D 
#' representation, the function will produce an error if it is not supplied with
#' a matrix of values that has a number of rows equal to 2 or 4. 
#' @param encodedSignal is a genomic signal of interest, for which the average 
#' power spectrum will be computed and returned to the user, you may encode 
#' your genomic character strings by using the encodeGenomes or encodeGenome 
#' function in this package. 
#' @return A vector of values indicating the average power spectral density 
#' (according to the Fourier Transform) for the encoded genomes four, or two 
#' constituent signals. 
#' @examples 
#' MTHFR100 <- "TGGCCAGGTATAGTGGCTCATACCTGTAATCCCAGCACTCTGGGAGACCGAAGCAGTATCACCTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATG"; 
#' encMTHFR100 <- encodeGenome(MTHFR100); 
#' psMTHFR100 <- getPowerSpectraSingle(encMTHFR100); 
#' plot(psMTHFR100, type='l',xlab='Frequency/Sequency',
#'      ylab='Power Spectral Density', main="Power Spectrum of First 100 nucleotides of MTHFR"); 
#' @export
getPowerSpectraSingle <- function(encodedSignal){
  PS <- function(x) {return(abs(x)^2)}
  if (!nrow(encodedSignal) %in% c(2,4)){
    print("Please encode genomic string prior to attempting to get the power spectrum"); 
    return(); 
  } 
  fc <- t(apply(encodedSignal,1,fft))
  ps <- t(apply(fc,1,PS)); 
  return(colMeans(ps))
}

#' Get the Fourier Power Spectra for an ensemble of encoded genomic signals 
#' 
#' This is a wrapper for the getPowerSpectraSingle function, that allows for 
#' the direct return of power spectra for each of the signals in an ensemble. 
#' @param encodedEnsemble a list of encoded genomes that are produced using 
#' one of the encoding methods programmed in this package. 
#' @return A list of Power spectra for the corresponding genomic sequences. 
#' @examples 
#' genStrings1 <- list('ACCAAGGATATTAGGACCC','CCCCAGGGAGATTTAGG','CCCGGGAGAGATTTAG'); 
#' encStrings1 <- encodeGenomes(genStrings1); 
#' psStrings1 <- getPowerSpectraEnsemble(encStrings1); 
#' @export
getPowerSpectraEnsemble <- function(encodedEnsemble){
  return(lapply(encodedEnsemble,getPowerSpectraSingle));
}

#' Evenly Scale Signals from their initial size of 'n' to size 'm'. 
#' 
#' This function will scale a power spectrum from an initial size of n to a 
#' size m.  Note that the new signal size m must be larger than the original 
#' signal size n, but cannot be too large (larger than twich the original size).
#' 
#' @param genomicPS The Power Spectrum of the genetic sequence that is being 
#' scaled, in general this could be any arbitrary real sequence, however in 
#' keeping with the spirit of the packages overall usability, this function is 
#' defined in terms that relate to the genomic nature of the data used. 
#' @param scaleTo The length to which it is desired to scale the original 
#' spectrum.  This is the same as the value m in the original paper by Yin et 
#' al. (2015). 
#' 
#' @return The scaled power spectrum, which is composed of the original sequence
#' of length n, and is scaled to the new size of m or 'scaleTo'. 
#' @examples 
#' tg <- "ACCAGGAGATTAGAGCCCCAGAGTAGAGCCCCAGAGATTAGAGCCAGAGTGAGAGCCGANNNAGAGC"; 
#' pstg <- getPowerSpectraSingle(encodeGenome(tg,'2D')); 
#' scaled <- evenlyScaleSingle(pstg,80); 
#' @export
evenlyScaleSingle <- function(genomicPS, scaleTo){
  n <- length(genomicPS); 
  m <- scaleTo; 
  if (m <= n){
    print(paste("Cannot scale sequence (or does not make sense to) 
                of size ", n, " to size ", m)); 
  }
  if (m >= 2*n){
    print("Scaling sequence to more than two times it's original size will 
          highly dilute the initial characteristics of the sequence.");
  }
  Tm <- numeric(m); 
  Tn <- genomicPS; 
  Tm[1] <- Tn[1]; 
  for (k in 2:m){
    Q <- k*n/m; 
    R <- floor(Q); 
    if (R == 0){
      R <- 1;
    }
    if (Q-R == 0){
      Tm[k] <- Tn[Q]; 
    }
    else {
      Tm[k] <- Tn[R] + (Q-R)*(Tn[R+1]-Tn[R]);
    }
  }
  return(Tm);
}

#'  Evenly scale an ensemble of spectra such that their
#'
#'  This function will take an ensemble of genomic power spectra, and scale them
#'  all such that they are of a length equivalent to the maximum length of a 
#'  sequence in the ensemble. 
#'  @param spectraList A list containing the power spectra of genomic signals 
#'  for sequences of interest. 
#'  @return A list of scaled spectra all of which should be the same length. 
#'  @examples 
#'  tg <- c("ACCCAAGAGAGAGCCCCCGAGAGAGAGAGAGAGAGCCCCGAGAGAGCGAGACGAGAC","TAGAGCCGAGATAGAGCCGAGAGTTAGAC","CGGAGAGNNGGAGAGCCCGAGAGTTTGAGNN")
#'  eg <- encodeGenomes(tg); 
#'  ps <- getPowerSpectraEnsemble(eg); 
#'  sps <- evenlyScaleEnsemble(ps);
#'  @export
evenlyScaleEnsemble <- function(spectraList){
  spectraLengths <- lapply(spectraList,length); 
  maxLength <- max(unlist(spectraLengths));
  maxIndex <- which(unlist(lapply(spectraList,length)) == max(unlist(lapply(spectraList,length))));
  scaledSpectra <- list(length(spectraList)); 
  scaledSpectra[[maxIndex]] <- spectraList[[maxIndex]];
  for (i in 1:length(spectraList)){
    if (i == maxIndex){
      next;
    }
    scaledSpectra[[i]] <- evenlyScaleSingle(spectraList[[i]], maxLength); 
  }
  return(scaledSpectra); 
}

#' Get the Power Spectra Distance
#' 
#' Computes the Euclidean distances among all of the sequences for all of the 
#' power spectra, applying a standard Euclidean distance measure to the entire 
#' computed spectrum.  The result is returned as a standard pairwise distances 
#' matrix. 
#' @param genomeList Genetic strings expected in a list. 
#' @return pair-wise distances matrix computed as the euclidean distance among 
#' the various power spectra for the sequences provided. 
#' @examples 
#' tg <- c("ACCCAAGAGAGAGCCCCCGAGAGAGAGAGAGAGAGCCCCGAGAGAGCGAGACGAGAC","TAGAGCCGAGATAGAGCCGAGAGTTAGAC","CGGAGAGNNGGAGAGCCCGAGAGTTTGAGNN")
#' dm <- getPowerSpectraDistances(tg); 
#' @export
getPowerSpectraDistances <- function(genomeList){
  eg <- encodeGenomes(genomeList);
  ps <- getPowerSpectraEnsemble(eg); 
  sps <- evenlyScaleEnsemble(ps); 
  dmat <- matrix(0,length(sps),length(sps)); 
  for (i in 1:length(sps)){
    for (j in i:length(sps)){
      dmat[i,j] <- sqrt(sum((sps[[i]]-sps[[j]])^2))
    }
  }
  return(dmat + t(dmat)); 
}

tstForMultiGenomes <- function(numStrings = 100, avLength, deviation){
  genomeStrings <- list(numStrings); 
  lengths <- round(rnorm(numStrings,avLength,deviation)); 
  for (i in 1:numStrings){
    genomeStrings[i] <- paste(sample(c('A','C','G','T'),lengths[i],replace=T),sep='',collapse='')
  }
  dmat <- getPowerSpectraDistances(genomeStrings);
  print(dmat); 
}

tstForMultiGenomes(avLength=300,deviation=5);