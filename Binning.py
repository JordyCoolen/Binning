#!/usr/bin/env python

import sys, os, argparse
from time import sleep
from math import log
import csv
import collections
import cProfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INFO = 'Tool to bin values using positive and negative sets'
USAGE = 'python Binning.py -d <datafile> -n <negativelist> -p <positivelist> -o <outputfilename> [--bins 0/2/3/4] [--type local/global] [--score 0/1]'

__version__="3.2"

# TODO:
# code has alot of redundancy, but because of time constrains we keep it this way
# make variable bins (in while), is hard to implement
# loops could be made smaller for global
# for i in xrange(1, len(rankedPN) - args.bins + (1 * args.minbinsize) - (args.minbinsize))
# for i in xrange(1, len(rankedPN) - (args.bins*args.minbinsize) + 1)

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO+'/n'+USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--data", type=str, required=True,
                        help="input data file (as full path)")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="output name (as full path)")
    parser.add_argument("-p", "--positive", type=str, required=True,
                        help="positive set (as full path)")
    parser.add_argument("-n", "--negative", type=str, required=True,
                        help="negative set (as full path)")
    parser.add_argument("--bins", type=int, default=4, choices=[0, 2, 3, 4, 5],
                        help="number of bins, should be either 0,2, 3, 4 or 5, "
                        "0 means loop over 2 - 5 bins to find optimal bin default = 4")
    parser.add_argument("--type", type=str, default='global', choices=['local','global'],
                        help="use either local or global optimization, default global, currently local is not working")
    parser.add_argument("--minbinsize", type=int, default=10,
                        help="minimum size per bin "
                        "default = 10")
    parser.add_argument("--score", type=int, choices=[0, 1], default=0,
                        help="scoring type "
                        "0 = |log2((P/totalP)/(N/totalN))|"
                        " method optimal likelihood"
                        "############################"
                        " 1 = log2((P/totalP)/(N/totalN))"
                        " method --> 0 <--, "
                        " default = 0")

    # parse all arguments
    args = parser.parse_args()
    
    # check if file exists
    if not os.path.exists(args.data):
        sys.exit('File is not present {}'.format(args.data))
    if not os.path.exists(args.positive):
        sys.exit('File is not present {}'.format(args.positive))
    if not os.path.exists(args.negative):
        sys.exit('File is not present {}'.format(args.negative))
    
    return args

def run(args):
    '''
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
                
    '''
    print 'Start finding best bins'
     
    read = open(args.data, 'r')
     
    genes = []
    allvalues = []
    removed = []
     
    data = csv.DictReader(read)
    for row in read:
        splitrow = row.strip().split('\t')
        if splitrow[0] == 'ID':
            continue
        if splitrow[1] == 'NA':
            print 'NA value found: {}'.format(row)
            removed.append(splitrow[0])
        else:
            genes.append(splitrow[0])
            allvalues.append(float(splitrow[1]))
        
# TEST
#     genes = genes[:1000]
#     allvalues = allvalues[:1000]
    
    # calculate min distance and get smaller distance for setting border
    dis = [i - j for i in allvalues for j in allvalues]
    mindis = 0.9 * min([i for i in dis if i > 0])
    
#     # add randomization to avoid grouping of P or N's
#     allvalues = [i+np.random.uniform(0,mindis) for i in allvalues] 
#     
#     # calculate min distance and get smaller distance for setting border
#     dis = [i - j for i in allvalues for j in allvalues]
#     mindis = 0.9 * min([i for i in dis if i > 0])
    
    # check for double ID's
    check_for_double_genes(genes)
    
    # get values of positive set
    print 'Positives'
    poslist = to_list(args.positive)
    poslist = remove(poslist, removed)
    check_for_double_genes(poslist)
    totalP = len(poslist)
    check_presence(genes, poslist)
            
    # get values of negative set
    print 'Negatives'
    neglist = to_list(args.negative)
    neglist = remove(neglist, removed)
    check_for_double_genes(neglist)
    totalN = len(neglist)
    check_presence(genes, neglist)
    
    print 'Number of P:{}'.format(totalP)
    print 'Number of N:{}'.format(totalN)
    ratioPN = totalP / float(totalN)
    print 'Ratio P/N:{}'.format(ratioPN)
    
    # check ID's of positive should not be in the negative list
    print 'Test if pos and neg list do not share genes'
    incorrect = False
    for i in poslist:
        if i in neglist:
            incorrect = True
            print 'Gene present in both list: {}'.format(i)
    if incorrect:
        sys.exit('Positive and Negative share genes')
        
    # sort values with corresponding genes to maintain order
    sortedRes = sorted(zip(allvalues, genes), key=lambda x: x[0])

    # extract genes
    sortedgenes = [x[1] for x in sortedRes]
    
    # extract values
    sortedvalues = [x[0] for x in sortedRes]
    
    MIN = sortedvalues[0]  # first element is minimal
    MAX = sortedvalues[-1]  # last element is maximum
 
    print 'Minimum value: {}'.format(MIN)
    print 'Maximum value: {}'.format(MAX)
    
    # change create list of genes by changing gene ID's with P or N
    # according to if its a positive or Negative
    
    rankedvaluesPN = []
    rankedPN = []
    indexPN = []
    
    for i, g in enumerate(sortedgenes):
        if g in poslist:
            indexPN.append(i)
            rankedvaluesPN.append(sortedvalues[i])
            rankedPN.append(1)
        if g in neglist:
            indexPN.append(i)
            rankedvaluesPN.append(sortedvalues[i])
            rankedPN.append(0)
    
    print rankedvaluesPN
    print rankedPN
    print indexPN
    
    #TODO this was a test to avoid border and PN problems:
    # make sortedvalues unique
    # match indexes to last value
    unique_indexPN = unique_indexing(rankedvaluesPN)
    print unique_indexPN
    
    # check if total P and total N matches rankedPN
    if rankedPN.count(1) != totalP:
        print 'Total P does not match, debug'
    if rankedPN.count(0) != totalN:
        print 'Total N does not match, debug'
    
    # select type of bin optimization
    if args.type == 'local':
        res_dict = main_loop(np.array(rankedPN), totalP, totalN)
        print res_dict
        make_plot(rankedPN, res_dict)
    if args.type == 'global':
        print 'binnummer:(indexPN,indexPN,score,leftborder,rightborder)'
        if args.bins == 2:
            MAXscore, res_dict = two_binning(args, rankedPN, indexPN, rankedvaluesPN, totalP, totalN, MIN, MAX, mindis)
        if args.bins == 3:
            MAXscore, res_dict = three_binning(args, rankedPN, indexPN, rankedvaluesPN, totalP, totalN, MIN, MAX, mindis)
        if args.bins == 4:
            MAXscore, res_dict = four_binning(args, rankedPN, indexPN, rankedvaluesPN, totalP, totalN, MIN, MAX, mindis)
        if args.bins == 5:
            MAXscore, res_dict = five_binning(args, rankedPN, indexPN, rankedvaluesPN, totalP, totalN, MIN, MAX, mindis)
        
        # optimizing 2 - 5 bins
        if args.bins == 0:
            MAXscore = 0
            args.bins = 2
            score, res = two_binning(args, rankedPN, indexPN, rankedvaluesPN, totalP, totalN, MIN, MAX, mindis)
            if 0.6*score > MAXscore:
                MAXscore = score
                res_dict = res
                args.bins = 3
                print score, res
                score, res = three_binning(args, rankedPN, indexPN, rankedvaluesPN, totalP, totalN, MIN, MAX, mindis)
                if 0.6*score > MAXscore:
                    MAXscore = score
                    res_dict = res
                    args.bins = 4
                    print score, res
                    score, res = four_binning(args, rankedPN, indexPN, rankedvaluesPN, totalP, totalN, MIN, MAX, mindis)
                    if 0.6*score > MAXscore:
                        MAXscore = score
                        res_dict = res
                        args.bins = 5
#                         print score, res
#                         score, res = five_binning(args, rankedPN, indexPN, rankedvaluesPN, totalP, totalN, MIN, MAX, mindis)
#                         if 0.9*score > MAXscore:
#                             MAXscore = score
#                             res_dict = res
#                             print score, res
  
        if results_check(res_dict, rankedPN, totalP, totalN):
            print '#'*20
            print 'Result:'
            print MAXscore
            print res_dict
            # write output to file
            with open(args.output+'.txt', 'w') as out:
                out.write('Parameters\n')
                out.write(str(args)+'\n')
                out.write('Results:\n')
                out.write('Positives:{}\n'.format(totalP))
                out.write('Negatives:{}\n'.format(totalN))
                out.write(str(MAXscore)+'\n')
                out.write('binnumber:(indexPN,indexPN,score,leftborder,rightborder)\n')
                out.write(str(res_dict))
            make_plot(rankedPN, res_dict, outputname=args.output)
              
    print 'Finished binning'
    
def two_binning(args, rankedPN, indexPN, rankedvaluesPN, totalP, totalN, MIN, MAX, mindis):
    
    # init
    if args.score == 0:
        MAXSUM = 0
    else:
        OPTIMAL = 200
    
    # loop over rankedPN
    for i in xrange(1, len(rankedPN) - args.bins + 1):
        bin1 = get_bin(0, i+1, rankedPN)
        
        # count number of P and N
        P = bin1.count(1)
        N = bin1.count(0)
        print("bin1: P:%s" %(P))
        print("bin1: N:%s" %(N))
        # calculate objective function to get score
        bin1score = objective_function(P, N, totalP, totalN, args.score)
        print'bin1score: %s' %(bin1score)
        # size of the bin
        if (indexPN[i] < args.minbinsize) and (rankedvaluesPN[i]-rankedvaluesPN[i+1] != 0):
            continue
        
        # to end
        bin2 = get_bin(i+1, len(rankedPN), rankedPN)
        
        # count number of P and N
        P = bin2.count(1)
        N = bin2.count(0)
        print("bin2: P:%s" %(P))
        print("bin2: N:%s" %(N))
        # calculate objective function to get score
        bin2score = objective_function(P, N, totalP, totalN, args.score)
        print'bin2score: %s' %(bin2score)
        # size of the bin
        if (indexPN[-1] - indexPN[i]) < args.minbinsize:
            continue
            
        SUM = abs(sum([bin1score, bin2score]))
        print'SUMscore: %s' %(SUM)
        
        # check if borders are different
        if check_borders(MIN, rankedvaluesPN[i+1]):
            if check_borders(rankedvaluesPN[i+1], MAX):
                
                if args.score == 0:
                    if SUM > MAXSUM:
                        MAXSUM = SUM
                        # binnummer:(indexPN,indexPN,score,leftborder,rightborder)
                        res_dict = {0:(0,i+1,bin1score,MIN, rankedvaluesPN[i+1]),
                                    1:(i+1,len(rankedPN),bin2score,rankedvaluesPN[i+1],MAX)}
                else:
                    if abs(0-SUM) < OPTIMAL:
                        OPTIMAL = abs(0-SUM)
                        MAXSUM = sum([abs(bin1score),abs(bin2score)])
                        # binnummer:(indexPN,indexPN,score,border)
                        res_dict = {0:(0,i+1,bin1score,MIN, rankedvaluesPN[i+1]),
                                    1:(i+1,len(rankedPN),bin2score,rankedvaluesPN[i+1],MAX)}
    
    return MAXSUM, res_dict

def three_binning(args, rankedPN, indexPN, rankedvaluesPN, totalP, totalN, MIN, MAX, mindis):
    
    # init
    if args.score == 0:
        MAXSUM = 0
    else:
        OPTIMAL = 200
    
    # loop over rankedPN
    for i in xrange(1, len(rankedPN) - args.bins + 1):
        bin1 = get_bin(0, i+1, rankedPN)
        
        # count number of P and N
        P = bin1.count(1)
        N = bin1.count(0)
        
        # calculate objective function to get score
        bin1score = objective_function(P, N, totalP, totalN, args.score)
        
        # size of the bin
        if (indexPN[i] < args.minbinsize) and (rankedvaluesPN[i]-rankedvaluesPN[i+1] != 0):
            continue
            
        for j in xrange(i+1, len(rankedPN) - args.bins + 2):
            bin2 = get_bin(i+1, j+1, rankedPN)
            
            # count number of P and N
            P = bin2.count(1)
            N = bin2.count(0)
            
            # calculate objective function to get score
            bin2score = objective_function(P, N, totalP, totalN, args.score)
            
            # size of the bin
            if (indexPN[j] < args.minbinsize) and (rankedvaluesPN[j]-rankedvaluesPN[j+1] != 0):
                continue
        
            # to end
            bin3 = get_bin(j+1, len(rankedPN), rankedPN)
            
            # count number of P and N
            P = bin3.count(1)
            N = bin3.count(0)
            
            # calculate objective function to get score
            bin3score = objective_function(P, N, totalP, totalN, args.score)
            
            # size of the bin
            if (indexPN[-1] - indexPN[j+1]) < args.minbinsize:
                continue
            
            SUM = abs(sum([bin1score, bin2score, bin3score]))
            
            # check if borders are different
            if check_borders(MIN, rankedvaluesPN[i+1]):
                if check_borders(rankedvaluesPN[i+1], rankedvaluesPN[j+1]):
                    if check_borders(rankedvaluesPN[j+1], MAX):
            
                        if args.score == 0:
                            if SUM > MAXSUM:
                                MAXSUM = SUM
                                # binnummer:(indexPN,indexPN,score,border)
                                res_dict = {0:(0,i+1,bin1score,MIN, rankedvaluesPN[i+1]),
                                            1:(i+1,j+1,bin2score,rankedvaluesPN[i+1],rankedvaluesPN[j+1]),
                                            2:(j+1,len(rankedPN),bin3score,rankedvaluesPN[j+1],MAX)}
                        else:
                            if abs(0-SUM) < OPTIMAL:
                                OPTIMAL = abs(0-SUM)
                                MAXSUM = sum([abs(bin1score),abs(bin2score),abs(bin3score)])
                                # binnummer:(indexPN,indexPN,score,border)
                                res_dict = {0:(0,i+1,bin1score,MIN, rankedvaluesPN[i+1]),
                                            1:(i+1,j+1,bin2score,rankedvaluesPN[i+1],rankedvaluesPN[j+1]),
                                            2:(j+1,len(rankedPN),bin3score,rankedvaluesPN[j+1],MAX)}
    
    return MAXSUM, res_dict
   
def four_binning(args, rankedPN, indexPN, rankedvaluesPN, totalP, totalN, MIN, MAX, mindis):
    
    # init
    if args.score == 0:
        MAXSUM = 0
    else:
        OPTIMAL = 200
    
    # loop over rankedPN
    for i in xrange(1, len(rankedPN) - args.bins + 1):
        bin1 = get_bin(0, i+1, rankedPN)
        
        # count number of P and N
        P = bin1.count(1)
        N = bin1.count(0)
        
        # calculate objective function to get score
        bin1score = objective_function(P, N, totalP, totalN, args.score)
        
        # size of the bin
        if (indexPN[i] < args.minbinsize) and (rankedvaluesPN[i]-rankedvaluesPN[i+1] != 0):
            continue
        
        # next loop
        for j in xrange(i+1, len(rankedPN) - args.bins + 2):
            bin2 = get_bin(i+1, j+1, rankedPN)
            
            # count number of P and N
            P = bin2.count(1)
            N = bin2.count(0)
            
            # calculate objective function to get score
            bin2score = objective_function(P, N, totalP, totalN, args.score)
            
            # size of the bin
            if ((indexPN[j] - indexPN[i]) < args.minbinsize) and (rankedvaluesPN[j]-rankedvaluesPN[j+1] != 0):
                continue
            
            for k in xrange(j+1, len(rankedPN) - args.bins + 3):
                bin3 = get_bin(j+1, k+1, rankedPN)
                
                # count number of P and N
                P = bin3.count(1)
                N = bin3.count(0)
                
                # calculate objective function to get score
                bin3score = objective_function(P, N, totalP, totalN, args.score)
                
                # size of the bin
                
                if ((indexPN[k] - indexPN[j]) < args.minbinsize) and (rankedvaluesPN[k]-rankedvaluesPN[k+1] != 0):
                    continue
            
                # to end
                bin4 = get_bin(k+1, len(rankedPN), rankedPN)
                
                # count number of P and N
                P = bin4.count(1)
                N = bin4.count(0)
                
                # calculate objective function to get score
                bin4score = objective_function(P, N, totalP, totalN, args.score)
                
                # size of the bin
                if (indexPN[-1] - indexPN[k]) < args.minbinsize:
                    continue
                
                SUM = abs(sum([bin1score, bin2score, bin3score, bin4score]))
                
                # check if borders are different
                if check_borders(MIN, rankedvaluesPN[i+1]):
                    if check_borders(rankedvaluesPN[i+1], rankedvaluesPN[j+1]):
                        if check_borders(rankedvaluesPN[j+1], rankedvaluesPN[k+1]):
                            if check_borders(rankedvaluesPN[k+1], MAX):
                
                                if args.score == 0:
                                    if SUM > MAXSUM:
                                        MAXSUM = SUM
                                        # binnummer:(indexPN,indexPN,score,border)
                                        res_dict = {0:(0,i+1,bin1score,MIN, rankedvaluesPN[i+1]),
                                                    1:(i+1,j+1,bin2score,rankedvaluesPN[i+1],rankedvaluesPN[j+1]),
                                                    2:(j+1,k+1,bin3score,rankedvaluesPN[j+1],rankedvaluesPN[k+1]),
                                                    3:(k+1,len(rankedPN),bin4score,rankedvaluesPN[k+1],MAX)}
                                else:
                                    if abs(0-SUM) < OPTIMAL:
                                        OPTIMAL = abs(0-SUM)
                                        MAXSUM = sum([abs(bin1score),abs(bin2score),abs(bin3score),abs(bin4score)])
                                        # binnummer:(indexPN,indexPN,score,border)
                                        res_dict = {0:(0,i+1,bin1score,MIN, rankedvaluesPN[i+1]),
                                                    1:(i+1,j+1,bin2score,rankedvaluesPN[i+1],rankedvaluesPN[j+1]),
                                                    2:(j+1,k+1,bin3score,rankedvaluesPN[j+1],rankedvaluesPN[k+1]),
                                                    3:(k+1,len(rankedPN),bin4score,rankedvaluesPN[k+1],MAX)}
    
    return MAXSUM, res_dict

def five_binning(args, rankedPN, indexPN, rankedvaluesPN, totalP, totalN, MIN, MAX, mindis):
    
    # init
    if args.score == 0:
        MAXSUM = 0
    else:
        OPTIMAL = 200
    
    # loop over rankedPN
    for i in xrange(1, len(rankedPN) - args.bins + 1):
        bin1 = get_bin(0, i+1, rankedPN)
        
        # count number of P and N
        P = bin1.count(1)
        N = bin1.count(0)
        
        # calculate objective function to get score
        bin1score = objective_function(P, N, totalP, totalN, args.score)
        
        # size of the bin
        if (indexPN[i] < args.minbinsize) and (rankedvaluesPN[i]-rankedvaluesPN[i+1] != 0):
            continue
        
        # next loop
        for j in xrange(i+1, len(rankedPN) - args.bins + 2):
            bin2 = get_bin(i+1, j+1, rankedPN)
            
            # count number of P and N
            P = bin2.count(1)
            N = bin2.count(0)
            
            # calculate objective function to get score
            bin2score = objective_function(P, N, totalP, totalN, args.score)
            
            # size of the bin
            if ((indexPN[j] - indexPN[i]) < args.minbinsize) and (rankedvaluesPN[j]-rankedvaluesPN[j+1] != 0):
                continue
            
            for k in xrange(j+1, len(rankedPN) - args.bins + 3):
                bin3 = get_bin(j+1, k+1, rankedPN)
                
                # count number of P and N
                P = bin3.count(1)
                N = bin3.count(0)
                
                # calculate objective function to get score
                bin3score = objective_function(P, N, totalP, totalN, args.score)
                
                # size of the bin
                if ((indexPN[k] - indexPN[j]) < args.minbinsize) and (rankedvaluesPN[k]-rankedvaluesPN[k+1] != 0):
                    continue
            
                for l in xrange(k+1, len(rankedPN) - args.bins + 4):
                    bin4 = get_bin(k+1, l+1, rankedPN)
                    
                    # count number of P and N
                    P = bin4.count(1)
                    N = bin4.count(0)
                    
                    # calculate objective function to get score
                    bin4score = objective_function(P, N, totalP, totalN, args.score)
                    
                    # size of the bin
                    if ((indexPN[l] - indexPN[k]) < args.minbinsize) and (rankedvaluesPN[l]-rankedvaluesPN[l+1] != 0):
                        continue
                    
                    # to end
                    bin5 = get_bin(l+1, len(rankedPN), rankedPN)
                    
                    # count number of P and N
                    P = bin5.count(1)
                    N = bin5.count(0)
                    
                    # calculate objective function to get score
                    bin5score = objective_function(P, N, totalP, totalN, args.score)
                    
                    # size of the bin
                    if (indexPN[-1] - indexPN[l]) < args.minbinsize:
                        continue
                
                    SUM = abs(sum([bin1score, bin2score, bin3score, bin4score, bin5score]))
                    
                    # check if borders are different
                    if check_borders(MIN, rankedvaluesPN[i+1]):
                        if check_borders(rankedvaluesPN[i+1], rankedvaluesPN[j+1]):
                            if check_borders(rankedvaluesPN[j+1], rankedvaluesPN[k+1]):
                                if check_borders(rankedvaluesPN[k+1], rankedvaluesPN[l+1]):
                                    if check_borders(rankedvaluesPN[l+1], MAX):
                    
                                        if args.score == 0:
                                            if SUM > MAXSUM:
                                                MAXSUM = SUM
                                                # binnummer:(indexPN,indexPN,score,border)
                                                res_dict = {0:(0,i+1,bin1score,MIN, rankedvaluesPN[i+1]),
                                                            1:(i+1,j+1,bin2score,rankedvaluesPN[i+1],rankedvaluesPN[j+1]),
                                                            2:(j+1,k+1,bin3score,rankedvaluesPN[j+1],rankedvaluesPN[k+1]),
                                                            3:(k+1,l+1,bin4score,rankedvaluesPN[k+1],rankedvaluesPN[l+1]),
                                                            4:(l+1,len(rankedPN),bin5score,rankedvaluesPN[l+1],MAX)}
                                        else:
                                            if abs(0-SUM) < OPTIMAL:
                                                OPTIMAL = abs(0-SUM)
                                                MAXSUM = sum([abs(bin1score),abs(bin2score),abs(bin3score),abs(bin4score), abs(bin5score)])
                                                # binnummer:(indexPN,indexPN,score,border)
                                                res_dict = {0:(0,i+1,bin1score,MIN, rankedvaluesPN[i+1]),
                                                            1:(i+1,j+1,bin2score,rankedvaluesPN[i+1],rankedvaluesPN[j+1]),
                                                            2:(j+1,k+1,bin3score,rankedvaluesPN[j+1],rankedvaluesPN[k+1]),
                                                            3:(k+1,l+1,bin4score,rankedvaluesPN[k+1],rankedvaluesPN[l+1]),
                                                            4:(l+1,len(rankedPN),bin5score,rankedvaluesPN[l+1],MAX)}
    
    return MAXSUM, res_dict

def check_borders(leftvalue, rightvalue):
    if leftvalue - rightvalue == 0:
        res = False
    else:
        res= True
    return res

def remove(all, removed):
    '''
        Remove removed from all
    '''
    return [i for i in all if i not in removed]

def check_for_double_genes(genes):   
    '''
        Checks if there are double genes
        present in list. And reports
        the double values.
    '''
    if len(set(genes)) != len(genes):
        print 'Double genes present please correct:'
        print 'Unique list: {}'.format(len(set(genes)))
        print 'All in list: {}'.format(len(genes))
        print '{} doubles'.format(len(genes) - len(set(genes)))
        checked = []
        for v in genes:
            if v in checked:
                print 'Please remove: {}'.format(v)
            checked.append(v)
        sys.exit('EXIT')

def to_list(input):
    '''
        Single column to list
    '''
    list = []
    with open(input, 'r') as input:
        for line in input:
            list.append(line.strip())
    return list

def check_presence(all, selection):
    '''
        Check if data in selection is also in all
    '''
    print 'Check if data for present for selection'
    incorrect = False
    for i in selection:
        if i not in all:
            incorrect = True
            print 'Not present: {}'.format(i)
    if incorrect:
        sys.exit('Genes in Positive/Negative not in data file')

def get_bin(lowerindex, upperindex, rankedPN):

    bin = rankedPN[lowerindex:upperindex]
    return bin

def objective_function(P, N, totalP, totalN, score):
    ''' score = 0 is likelihood
        score = 1 is closest to 0 ( --> 0 <-- )
    '''
    if (P==0) or (N==0):
        return 0
    if score == 0:
        return abs(log((float(P)/totalP)/(float(N)/totalN),2))
    if score == 1:
        return log((float(P)/totalP)/(float(N)/totalN),2)
    
def loop(d1array, totalP, totalN):
    
    stop = 2
    objective = 0
    for i in xrange(4, len(d1array)+1):
        obj = objective_function(sum(d1array[np.arange(i)]),len(d1array[np.arange(i)])-sum(d1array[np.arange(i)]),totalP, totalN)
        if obj>objective:
            objective=obj
            stop = i
    return (objective, stop)


def right_slider(partition1, partition2, totalP, totalN):
    lik = objective_function(sum(partition1),len(partition1)-sum(partition1),totalP, totalN)+objective_function(sum(partition2),len(partition2)-sum(partition2),totalP, totalN)
    best = len(partition1)
    for i in xrange(10, len(partition1)):
        p1 = partition1[0:i]
        p2 = np.concatenate([partition2, partition1[i+1::]]) 
        lik2 = objective_function(sum(p1),len(p1)-sum(p1),totalP, totalN)+objective_function(sum(p2),len(p2)-sum(p2),totalP, totalN)
        print lik, lik2
        if lik2>lik:
            best=i-1
            lik=lik2
    return best

def main_loop(d1array, totalP, totalN):
    bins={}
    c = len(d1array)
    x=-1
    counter=0
    while x<c:
        x+=1
        array = d1array[x::]
        best_partition =loop(array, totalP, totalN)
        next_partition = loop(array[(x+best_partition[1]-1)::], totalP, totalN)
        array2 = array[(x+best_partition[1]-1):next_partition[1]]
        rs = right_slider(array, array2, totalP, totalN)

        best_partition = (best_partition[0],rs)
        check_last = len(d1array[(best_partition[1]+x)::])
        if check_last<4:
            bins[counter] = (x, best_partition[1]+x+check_last, best_partition[0])
            x+=best_partition[1]+x+check_last
        else:
            bins[counter] = (x, best_partition[1]+x, best_partition[0])
            x+=best_partition[1]-1
        print x
        counter+=1
    return bins

def make_plot(data, result_dict,outputname, v=None):
    counter=0
    tp=sum(data)
    tn =len(data)-tp
    xv=[]
    xl=[]
    plt.style.use('ggplot')
    for i in result_dict:
        data_range=data[result_dict[i][0]:result_dict[i][1]]
        pos = sum(data_range)
        neg = len(data_range)-sum(data_range)
        plt.bar(counter+1, pos/float(tp), color='b')
        plt.bar(counter, neg/float(tn), color='r')
        xv.append(counter+0.5)
        xl.append('[%.3f'%result_dict[i][3]+': %.3f]:'%result_dict[i][4]+'%.4f'%result_dict[i][2])
        counter+=4
    plt.xticks(xv,xl,fontsize=8)
    if not outputname.endswith('.png'):
        plt.savefig(outputname+'.png')
    else:
        plt.savefig(outputname)

def results_check(res_dict, rankedPN, totalP, totalN):
    P = 0
    N = 0
    for i in res_dict:
        bin = rankedPN[res_dict[i][0]:res_dict[i][1]]
        print i, 'pos ', bin.count(1)
        print i, 'neg ', bin.count(0)
        P += bin.count(1)
        N += bin.count(0)
    
    print 'results_check'
    print 'P',P
    print 'N',N
    print 'totalP',totalP
    print 'totalN',totalN
    
    if (P+N) == (totalP + totalN):
        return True

#TODO: test not implemented
def unique_indexing(values):
    ''' takes a list and outputs unique list of values
        and a matched dictionary with indexes matching the original
        index of the last value of reoccuring values.
        
        Example:
        [0, 1, 2, 2, 3, 4, 9, 9, 9, 9, 10, 11, 12]
        Gives:
        dic {0: 0, 1: 1, 2: 3, 3: 4, 4: 5, 9: 9, 10: 10, 11: 11, 12: 12}
        
        Input:
        - values, list of values
    '''
    count=0
    dic={}
    for i in values:
        if i in dic.keys():
            dic[i]=count
        else:
                dic[i]=count
        count+=1
    return dic

if __name__ == "__main__":
#     # load arguments set global workdir
    args = parse_args()
    run(args)

