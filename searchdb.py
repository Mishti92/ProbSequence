# Constants.
wmer_size = 7

#import modules
import numpy as np
import math
import pickle
import sys


# Measure the time spent in the algorithm.
def tic():
    #Homemade version of matlab tic and toc functions
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print "Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds."
    else:
        print "Toc: start time not set"


def findhits( searchsequence , numhits):
    for index in xrange( len(searchsequence) - (wmer_size-1)  ):
        currentseq = searchsequence[ index : index + wmer_size ]
        yield currentseq, pointers[ currentseq ][0:numhits], pointers_idx[ currentseq ][0:numhits]

if __name__ == '__main__':

#    print "Loading database.."
    tic()
    pointers = pickle.load( open( "pointers.p", "rb" ) )
    pointers_idx = pickle.load( open( "pointers_idx.p", "rb" ) )
    toc()


    if len(sys.argv) == 2:
        searchsequence = sys.argv[1]
    else:
        #default
        searchsequence = 'CTGAGGGGCCGGGGAAGGGGGCCAATCCTACCTTGGGGTGACTGAGAGGCTCAGGGAGGTGGCTCACTCCTGGGCTGCCATCAGGCTTTCCCAGGGAGCGTTCTGGCACCTGATGACTGAGGCCTGGGTGTTAACCCTGGCACTCCACTTCCTGCCTCCCTGATTTCTCCCCCTTTTTGGGGTGTCGTATATGAAGTGAGGTAACAGACGTGATGGTGCTGGCTGAGCCGCTGGCCCTGGTGTGGGGTGGGTGGAGTTGGGCCTGAAATGTGGGCATCCAGGCTGAGTGAGGCCCTCCCCGCCACCCCCAGGGCCTGGGCTGACGCCTCTCAGCAGCCCATTCTAGCCCCTCCCCCAGGAGCCAGGACCTACGTGTAGACTTTGTCACCAAT'
        
    
    count = 0
    hits_list=[]
    prob_hits_list=[]
    index_hits_list=[]
    for k1,k2,k3 in findhits(searchsequence,10):
        hits_list.append([count, k1])
        prob_hits_list.append(k2)
        index_hits_list.append(k3)
        count =count +1
    
    pickle.dump( searchsequence, open( 'picklefiles/input_sequence.p', "wb" ) )
    pickle.dump( hits_list, open( 'picklefiles/hits_list.p', "wb" ) )
    pickle.dump( prob_hits_list, open( 'picklefiles/prob_hits_list.p', "wb" ) )
    pickle.dump( index_hits_list, open( 'picklefiles/index_hits_list.p', "wb" ) )
    
