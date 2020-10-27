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

## Example results

--score 0 (optimal likelihood method)


Number of Gene/identifiers for the positive and negative set

totalP 53
totalN 294

Result:
overall score = 8.13212179365

0: (0, 11, 2.734786296106959, 1.0, 101.0)

1: (11, 18, 5.056714390994322, 101.0, 162.0)

2: (18, 347, 0.3406211065510634, 162.0, 5572.0)

Number of positive and negatives in optimized bins

0 pos  6
0 neg  5
1 pos  6
1 neg  1
2 pos  41
2 neg  288

plot with the results:

![alt tag](https://github.com/JordyCoolen/Binning/blob/main/example/example_opt_likelihood.png)

--score 1 (closest to zero --> 0 <--)


Number of Gene/identifiers for the positive and negative set

totalP 53
totalN 294


Result:
overall score = 3.44931450115

0: (0, 104, 1.1009141949048564, 1.0, 1243.0)

1: (104, 173, 0.6237549837182151, 1243.0, 2441.0)

2: (173, 347, -1.7246453225303384, 2441.0, 5572.0)

Number of positive and negatives in optimized bins

0 pos  29
0 neg  75
1 pos  15
1 neg  54
2 pos  9
2 neg  165

plot with the results:

![alt tag](https://github.com/JordyCoolen/Binning/blob/main/example/example_closest_to_zero2.png)

## Known problems
```python
--type local
```
Is currently not working.

## Version
27 October 2020

## Authors
J.P.M. Coolen and D.R. Garza

## License
[MIT](https://choosealicense.com/licenses/mit/)
