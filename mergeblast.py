# Constants.

numparts = 60 #how many subparts of db? (at /blast folder, if 0...60, choose 60, not 61)

#import modules
import numpy as np
import math
import pickle


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

if __name__ == '__main__':


    pointers = pickle.load( open( "blast/blast_0.p", "rb" ) )
    pointers_idx = pickle.load( open( "blast/blastidx_0.p", "rb" ) )

    print "Loading database entries (takes 1 minute in my computer).."
    tic()
    pointers = pickle.load( open( "blast/blast_0.p", "rb" ) )
    pointers_idx = pickle.load( open( "blast/blastidx_0.p", "rb" ) )

    for id in range(1, numparts + 1):
        x = pickle.load( open( "blast/blast_%d.p" % id, "rb" ) )
        y = pickle.load( open( "blast/blastidx_%d.p" % id, "rb" ) )

        for key in pointers:
            if key in x:
                pointers[key] = np.concatenate([pointers[key],x[key]])
                pointers_idx[key] = np.concatenate([pointers_idx[key],y[key]])

        for key in x:
            if key not in pointers:
                pointers[key] = x[key]
                pointers_idx[key] = y[key]

    print "Database loading complete.. "
    toc()



    print "Sorting started.. "
    tic()
    for key in pointers:
        sortedindices = sorted(range(len(pointers[key])),key=lambda x : pointers[key][x],reverse=True)
        pointers[key] = pointers[key][sortedindices]
        pointers_idx[key] = pointers_idx[key][sortedindices]
    print "Sorting complete.."
    toc()

    pickle.dump( pointers, open( 'pointers.p' , "wb" ) )
    pickle.dump( pointers_idx, open( 'pointers_idx.p', "wb" ) )