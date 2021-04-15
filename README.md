# Genomic-DFT-Yin-R Package
This is a simple R-Package that allows for the user to implement the DFT distance calculations and phylogeny building 
procedures that are initially discussed in Yin and Yau 2014, and extended in Yin and Yau 2015.  To install and load 
the R Library see the installation section. 

## Installation 

In order to install R-libraries directly from Github, the user must first install the devtools library with: 
```r
install.packages('devtools'); 
```
Once the devtools package has been installed, the user can load the library using the standard 'library' or 'require' functions
then from the devtools package the user can install the Genomic-DFT-Yin-R-Package by using the 'install_github' function. 
```r
library(devtools);
install_github('mathornton01/Genomic-DFT-Yin-R'); 
```
Now R will display several messages that give various updates of the status of the download and installation.  Once the installation 
has been completed, the user can load the library with the typical 'library' and 'require' functions. 
```r
library(YinGenomicDFTDistances); 
```

## Checking the Documentation 

Users proficient in R already know that the man command or '?' prior to a function will allow the user to pull up a reference on 
a particular function in a package if that documentation exists.  To check the documentation for the functions in the Yin DFT 
package, the same ? notation may be used.  For instance, to check the documentation on how to encode an ensemble of genome 
strings (represented as a list of character strings) one may use the following command: 

```r
?encodeGenomes
``` 

## Doing a Simple Analysis 

