# Constants.
wmer_size = 7 # Blast uses w = 7, but 11, 15 are also possible choices.

# Only consider sequences above certain threshold.
countrelevant = 0 # Number of relevant sequences (relevancy determined by probthreshold).
probthreshold = 4 # Only consider w-mer sequences if they have larger probability than this threshold.
                    # Purpose is to not bloat the database with very unlikely "hits".

#import modules
import numpy as np
import math


from sys import getsizeof

# Convert base 10 number to base 4.
def base4(n):
    if n == 0: return [str(0)]
    digits = []
    while n:
        digits.append(str(n % 4))
        n /= 4
    return digits[::-1]

# Convert base 4 number to base 10.
def base10(s):
    num = 0
    s = str(s); s_len = len(s);
    for x in xrange(s_len):
        num += (4**x)*int(s[s_len-x-1])
    return num


# Convert nucleotides to numbers.
# I'm using number representation instead of nucleotides. [a,c,g,t] => [0,1,2,3]
# So, for example, ACTCGAA => [0131200].  nuc2num('ACTCGAA')
# CAUTION: OUTPUT STILL STRING. USE int(output) to convert to integers.
mapnucleotides = {'A': '0', 'C': '1', 'G': '2', 'T': '3'}
def nuc2num(s):
    num = [mapnucleotides[x] if x in mapnucleotides else -1 for x in s]
    return "".join(num)

# Convert numbers to nucleotides.
mapnums = {'0': 'A', '1': 'C', '2': 'G', '3': 'T'}
def num2nuc(no):
    s = [mapnums[x] if x in mapnums else '-1' for x in no]
    return s

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

    # READ IN THE DATABASE OF THE NUCLEOTIDES.
    seq = open('chr22.maf.ancestors.42000000.complete.boreo.fa', 'r').read()

    # READ IN THE PROBABILITY FILE.
    probn = np.zeros(len(seq)) #prob of each nucleotide in the database
    with open('chr22.maf.ancestors.42000000.complete.boreo.conf') as f:
        index = 0
        for l in f:
            for prob in l.split():
                probn[index] = float( prob )
                index += 1


    # CREATE AND FILL THE DATABASE.

    # Generate list of lists, which will hold the locations at which the
    # w-mer occurs in the database.


    num_wmers = 4**wmer_size
    pointers = [[] for x in xrange(num_wmers)]


    # Scan whole database.
    tic()
    wmer_str = ['*']*wmer_size


    #for index in xrange( len(seq) - (wmer_size-1)  ):
    for index in xrange( 100  ):
        print index
        probw = probn[index : index + wmer_size] # Current w-mer.

        currentseq = seq[index : index + wmer_size] # Current nucleotide sequence in the database.
        currentseq_num = nuc2num(seq[index : index + wmer_size]) # Numerical version of current seq.

        wildcards = [i for i, prob in enumerate(probw) if prob < 1.00 ] # If sequence is ACGTAAA but the 3rd and 7th positions are not 100% G and A: AC*TAA*
        deterministicpos = list( set(xrange(wmer_size)) - set(wildcards) ) # All the positions with probability = 1.00

        numpermutations = 4**len(wildcards) # Number of possibilities that can fit to this subsequence. Using above example, there are 2 positions that can be filled with any of the 4, hence 4^2
        permutationlength = len( base4(numpermutations - 1) ) # Assuming we have 64 permutations, numbering will be 0, .., 63. 0 is "0" in base4, but 63 is "333", we want to make "0"=>"000" for easier processing.


        # Place deterministic nucleotides into the current w-mer string. If we have
        # CAACTAA with only A's are fully deterministic,
        # now we'll have *AA**AA. (Actually *00**00, number representation of nucleotides)

        if len(deterministicpos)>0:
            map(wmer_str.__setitem__, deterministicpos, [ nuc2num(currentseq[i]) for i in deterministicpos])

        # Probability of current w-mer.
        wmer_prob_org = np.zeros(wmer_size)
        if len(deterministicpos)>0:
            wmer_prob_org[deterministicpos] = probn[deterministicpos]


        # Try all permutations and add to database.
        for currentpermutation in xrange(numpermutations):
            wmer_prob = wmer_prob_org

            permutestr = base4(currentpermutation)
            permutestr = ['0'] * ( permutationlength - len(permutestr) ) + permutestr

            if len(wildcards)>0:
                map(wmer_str.__setitem__, wildcards, permutestr) # This is the number encoding of the w-mer. hash-index=0123211 => encoded w-mer=ACGTGCC.
            hashindex = base10( int( ('').join(wmer_str) ) ) # wmer_str is on base4, to use it as hash-index, must convert to base 10.

            # W-mer probability. It is same for all values because we assign (for this project)
            # the same probability for each of the remaining nucleotides..
            for i,k in enumerate(wildcards):
                if currentseq_num[k] == permutestr[i]: # We try each combination that can fit to current w-mer. If the combination's kth character/nucleotide is actually the character given in the original sequence, then use probability directly.
                    wmer_prob[k] = probw[k]
                else: # If not, just add 1/3th of 1-prob
                    wmer_prob[k] = (1-probw[k])/3

            # Log sum trick to avoid vanishing probabilities.
            wmer_prob = [math.log(k, 2) for k in wmer_prob] # convert all probs to log. scale
            b = max( wmer_prob )
            probsum = [math.exp(k) for k in \
                                            [p-b for p in wmer_prob]]
            probsum = b + np.sum(probsum)

            if probsum > probthreshold: # If the probability of current sequence is above certain threshold, we put it to the database.
                countrelevant += 1
            #pointers[hashindex].append([index,wmer_prob]) # We add an entry to the database, in the following form: (Index that the w-mer occurs/begins, Probability of this w-mer occuring in this position w.r.t. the probabilities provided in database)

    toc()
    print countrelevant

