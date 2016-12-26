
# coding: utf-8

# In[285]:

import numpy as np
import matplotlib.pyplot as plt
import pickle
from operator import itemgetter
import math

# In[400]:

seq = open("chr22.maf.ancestors.42000000.complete.boreo.fa").read()
prob =(open("chr22.maf.ancestors.42000000.complete.boreo.conf").read()).split()

m= len(seq)

#If using log and sum
# match_score = 4
# match_low_score = 4
# mismatch_score = 1

match_score = 2
match_low_score = 1
mismatch_score = 0
wmer_size=7


#Loads previous results from searchdb.py file 
hits_list = pickle.load( open( "picklefiles/hits_list.p", "rb" ) )
prob_hits_list = pickle.load( open( "picklefiles/prob_hits_list.p", "rb" ) )
index_hits_list = pickle.load( open( "picklefiles/index_hits_list.p", "rb" ) )
inp_query = pickle.load( open( "picklefiles/input_sequence.p", "rb" ) )



# In[525]:

def match_nuc(nuc1,nuc2):
    #Matches nucleotide and returns true or false
    if nuc1.lower() == nuc2.lower():
        return True
    else:
        return False

# In[526]:

#If using transition, transversion scores
#def check_nuc(nuc1,nuc2):
##     assert len(a)==len(b)==1
#    if (nuc1 in ['a','A'] and nuc2 in ['g','G']) or (nuc2 in ['a','A'] and nuc1 in ['g','G']):
#        return "transition"
#    elif (nuc1 in ['c','C'] and nuc2 in ['t','T']) or (nuc2 in ['c','C'] and nuc1 in ['t','T']):
#        return "transition"
#    else:
#        return "transversion"


# In[406]:


def check_prob(seqind):
    #Returns probability of the nucleotide in the database using index
    return float(prob[seqind])  
    

# In[408]:

def high_prob(nuc_ind):
    #Threshold set to 0.6 for probability. Checked with 0.5. Results better with 0.6
    if check_prob(nuc_ind) > 0.6:
        return True
    else:
#         print check_prob(nuc_ind)
        return False

def find_evalue(m,n,score):
    k=0.1
    lambdaS=0.2
    return k*m*n*math.exp(-lambdaS*score)
    

# In[410]:

def return_score(inp_nuc,hit_nuc,hit_nuc_ind):
    #Returns the score for the two nucleotides being compared
    if match_nuc(inp_nuc,hit_nuc):
        if high_prob(hit_nuc_ind):
#            print 'Match with high prob'
            return match_score*check_prob(hit_nuc_ind)
#             return np.log(match_score*check_prob(hit_nuc_ind))
        else:
#            print 'Match with low prob'
            return match_low_score*check_prob(hit_nuc_ind)
#             return np.log(match_low_score*check_prob(hit_nuc_ind))
    else:
        if high_prob(hit_nuc_ind):
#            print 'Mis-match with high prob'
#             return np.log(mismatch_score)
            return mismatch_score
        else:
            match_prob = (1 - check_prob(hit_nuc_ind))/3
#             print match_prob
#            print 'Mis-match with low prob'
            return match_low_score*match_prob
#             return np.log(match_low_score*match_prob)


# In[413]:

def score_hits(inp_ind,hit_ind):
    #Function to score hits from left and right of the hit with respect to the input query
    inp_ind_left=inp_ind
    inp_ind_right=inp_ind + wmer_size
    hit_ind_left= hit_ind
    hit_ind_right= hit_ind + wmer_size
    
    
    #Initiliazing lists and variables
    score_left = 1
    score_list_left=[]
    score_right = 1
    score_list_right=[]
    prob_list_left=[]
    prob_list_right=[]
    
    #Left of hit wmer
    for i in range(len(inp_query)-wmer_size):
        if hit_ind_left == 0 or inp_ind_left == 0:
            break
        inp_ind_left = inp_ind_left -1
        hit_ind_left = hit_ind_left -1        
        r_score = return_score(inp_query[inp_ind_left], seq[hit_ind_left], hit_ind_left)
#         print 'r_score', r_score
        score_left = score_left + r_score
        prob_list_left.append(r_score)
        score_list_left.append(score_left)
        
#To visualize the score towards left of the hit
#    plt.plot(score_list_left)
#    plt.ylabel('Score from Left of Hit')
#    plt.xlabel('Towards Left of Hit Index')
#    plt.show()

    #Right of hit wmer
    for i in range(len(inp_query)-wmer_size):
        if hit_ind_right == len(seq) or inp_ind_right == len(inp_query):
            break
        r_score = return_score(inp_query[inp_ind_right], seq[hit_ind_right], hit_ind_right)
        score_right = score_right + r_score
#         print 'rscore',r_score
        score_list_right.append(score_right)
        prob_list_right.append(r_score)
        inp_ind_right = inp_ind_right + 1
        hit_ind_right = hit_ind_right + 1
        
            
#To visualize the score towards right of the hit
#    plt.plot(score_list_right)
#    plt.ylabel('Score from Right of Hit')
#    plt.xlabel('Towards Right of Hit Index')
#    plt.show()
#    
    return score_list_left, score_list_right, prob_list_left, prob_list_right


if __name__ == '__main__':
    
    matched_list=[]    
    print 'Input Query: '
    print inp_query
    n = len(inp_query)
    print '-------------------------------------'
    print '-------------------------------------'
    for i in range(len(hits_list)):
        inp_ind = hits_list[i][0]
        for h in range(len(index_hits_list[i])):
            hit_ind = index_hits_list[i][h]
            prob_hit = prob_hits_list[i][h]
        #function to calculate score and find alignments
            score_list_left, score_list_right, prob_list_left, prob_list_right= score_hits(inp_ind,hit_ind)

            #Towards right of the hit
            threshold = 6
            if len(score_list_right) < threshold:
                threshold = 1
            if not score_list_right:
                i_right = 0
                prob_right = 0
            else:
                for i1 in range(len(score_list_right)):
                    if i1+threshold >=len(score_list_right):
                        threshold = threshold -1
                    if score_list_right[i1] == score_list_right[i1+threshold]:
                        i_right = i1-1
                        prob_right = score_list_right[i1]
                        break


            # In[535]:

            #Towards left of the hit
            threshold = 6
            if len(score_list_left) < threshold:
                threshold = 1
            if not score_list_left:
                i_left = 0
                prob_left = 0
            else:
                for i2 in range(len(score_list_left)):
                    if i2+threshold >=len(score_list_left):
                        threshold = threshold -1              
                    if score_list_left[i2] == score_list_left[i2+threshold]:
                        i_left = i2+1
                        prob_left = score_list_left[i2]
                        break
        


            prob_final = prob_hit
            for item in prob_list_right[:i_right]:
                if item != 0:
                    prob_final = prob_final+item
            for item in prob_list_left[:i_left]:
                if item != 0:
                    prob_final = prob_final+item
            inp_query_string = inp_query[inp_ind-i_left:inp_ind+wmer_size+i_right]
            
#            print 'Index in database: '
            hit_ind_db=hit_ind-i_left
            
#            print 'Matched Alignment in Database: '
            matched_alignment = seq[hit_ind-i_left:hit_ind+wmer_size+i_right]
    
            prob_alignment = prob[hit_ind-i_left:hit_ind+wmer_size+i_right]
            
            temp_list = [inp_query_string,hit_ind_db,matched_alignment,prob_final,len(matched_alignment)]
    
            if temp_list not in matched_list:
                matched_list.append(temp_list)
            
            
    matched_list = sorted(matched_list, key=itemgetter(4),reverse=True)
    
    #Prints output
    for item in matched_list:

        print 'Input Query SubString:'
        print item[0]
        print 'Index of Matched Alignment in DB: ', item[1]
        print 'Matched Alignment: '
        print item[2]
        print 'Score: '
        print item[3]
        print 'E Value: '
        print find_evalue(m,n,item[3])
        print '--------------------------------'
