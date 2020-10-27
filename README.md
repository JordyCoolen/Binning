# Binning
For bayesian integration studies genes of each dataset needs to be
grouped in bins based on their value compared to the positive predictive geneset and/or the negative predictive geneset.

To have a optimized bayesian integration.
This script automates the binning by optimizing for a high likelihood score in a automated fashion.
Possible to optimize for optimal likelihood or centered around zero.
Currently only a global optimalization is possible, which is in most cases the preferred method of optimalization.

The bin borders can be globally optimized for a 2,3,4 or 5 bins.
Or the more preferred method let the binning tool calculate the most optimal number of bins and borders (--bins 0)

The binning process uses a negative and postive geneset/identifier which all should be present in the dataset. If missing values are present or there is overlap in identifiers between the negative and positive set the program will give a error indicating what the problem is.

## Algorithm
                reads data file
                        |
                extract genes and values
                        |
                checks double genes
                        |
                get genes from positive and negative
                        |
                check if genes in positive and negative
                are present in input data file
                        |
                check if there is no overlap between
                positive and negative list
                        |
                sort values while maintaining gene order
                        |
                extract sorted genes to list
                        |
                extract sorted values to list
                        |
                get MIN and MAX values
                        |
               set MIN and MAX borders
                        |
       sweeps all combinations of each bin and remembers the score
                        |
         highest scoring option will be given and plotted

## Install
download or pull this repository then:

```bash
conda create -n <environmentname> python=2.7
conda activate <environmentname>
conda install -c conda-forge matplotlib
```

## Usage
```python
python Binning.py -d <datafile> -n <negativelist> -p <positivelist> -o <outputfilename> [--bins 0/2/3/4] [--type local/global] [--score 0/1]
```
## Parameters
--bins    0: calculates the optimum number of bins with a max of 5 bins, with a high likelihood score (preferred method)
--bins    2/3/4/5 force the code to optimize for this many bins
--score   0: |log2((P/totalP)/(N/totalN))| for optimal likelikhood
--score.  1: log2((P/totalP)/(N/totalN)) for closest to 0 optimalization

P = number of positives in bin
totalP = number of total postives in positivelist
N = number of negatives in bin
totalN = number of total negatives in negativelist

## Known problems
```python
--type local
```
Is currently not working.

## Version
27 October 2020

## Authors
J.P.M. Coolen
D.R. Garza

## License
[MIT](https://choosealicense.com/licenses/mit/)
