
# coding: utf-8

# In[285]:

import numpy as np
import matplotlib.pyplot as plt
import pickle


# In[400]:

seq = open("seq.txt").read()
prob =(open("prob.txt").read()).split()
# match_score = 4
# match_low_score = 4
# mismatch_score = 1
match_score = 1
match_low_score = 1
mismatch_score = 0
wmer_size=7


# In[523]:

# print len(seq)
# print len(prob)


# In[524]:

# print seq[2]
# print prob[2]


# In[525]:

def match_nuc(nuc1,nuc2):
    if nuc1.lower() == nuc2.lower():
        #Match - prob - high or low
        return True
    else:
        return False
#         if check_nuc(nuc1,nuc2) == 'transition':
#             #do computation for transition
#             print 'transition'
#         elif check_nuc(nuc1,nuc2) == 'transversion':
#             #do computation for transversion
#             print 'transversion'


# In[526]:

def check_nuc(nuc1,nuc2):
#     assert len(a)==len(b)==1
    if (nuc1 in ['a','A'] and nuc2 in ['g','G']) or (nuc2 in ['a','A'] and nuc1 in ['g','G']):
        return "transition"
    elif (nuc1 in ['c','C'] and nuc2 in ['t','T']) or (nuc2 in ['c','C'] and nuc1 in ['t','T']):
        return "transition"
    else:
        return "transversion"


# In[527]:

nuc1='a'
nuc2='A'
match_nuc(nuc1,nuc2)


# In[406]:

def check_prob(seqind):
    #high and low
#     print prob[seqind]
    return float(prob[seqind])  
    


# In[407]:

check_prob(2)


# In[408]:

def high_prob(nuc_ind):
    if check_prob(nuc_ind) > 0.6:
        return True
    else:
#         print check_prob(nuc_ind)
        return False


# In[409]:

high_prob(2)


# In[410]:

def return_score(inp_nuc,hit_nuc,hit_nuc_ind):
#     print inp_nuc, hit_nuc
    #compares 
    if match_nuc(inp_nuc,hit_nuc):
        if high_prob(hit_nuc_ind):
            print 'Match with high prob'
            return match_score*check_prob(hit_nuc_ind)
#             return np.log(match_score*check_prob(hit_nuc_ind))
        else:
            print 'Match with low prob'
            return match_low_score*check_prob(hit_nuc_ind)
#             return np.log(match_low_score*check_prob(hit_nuc_ind))
    else:
        if high_prob(hit_nuc_ind):
            print 'Mis-match with high prob'
#             return np.log(mismatch_score)
            return mismatch_score
        else:
            match_prob = (1 - check_prob(hit_nuc_ind))/3
#             print match_prob
            print 'Mis-match with low prob'
            return match_low_score*match_prob
#             return np.log(match_low_score*match_prob)


# In[536]:

hits_list = pickle.load( open( "hits_list.p", "rb" ) )
prob_hits_list = pickle.load( open( "prob_hits_list.p", "rb" ) )
index_hits_list = pickle.load( open( "index_hits_list.p", "rb" ) )


# In[537]:

prob_hits_list[0]
index_hits_list[0]


# In[413]:

def score_hits(inp_ind,hit_ind):
    #Function to score hits from left and right of the wmer with respect to the input query
    inp_ind_left=inp_ind
    inp_ind_right=inp_ind + wmer_size
    hit_ind_left= hit_ind
    hit_ind_right= hit_ind + wmer_size
    score_left = 1
    score_list_left=[]
    score_right = 1
    score_list_right=[]
    prob_list_left=[]
    prob_list_right=[]

    # test_list1=[]
    # test_list2=[]

    # print len(inp_query)

    #Left of hit wmer
    for i in range(len(inp_query)-wmer_size):
        if hit_ind_left == 0 or inp_ind_left == 0:
            break
        inp_ind_left = inp_ind_left -1
        hit_ind_left = hit_ind_left -1
#         print inp_query[inp_ind_left], seq[hit_ind_left], hit_ind_left
#         print check_prob(hit_ind_left)
        r_score = return_score(inp_query[inp_ind_left], seq[hit_ind_left], hit_ind_left)
#         print 'r_score', r_score
        score_left = score_left + r_score
        prob_list_left.append(r_score)
        score_list_left.append(score_left)
#         print "score", score_left
        

#     print score_list_left
#     print prob_list_left
    plt.plot(score_list_left)
    plt.ylabel('Score from Left of Hit')
    plt.xlabel('Towards Left of Hit Index')
    plt.show()

    #Right of hit wmer

    for i in range(len(inp_query)-wmer_size):
#         print inp_query[inp_ind_right], seq[hit_ind_right], hit_ind_right
#         print check_prob(hit_ind_right)
        r_score = return_score(inp_query[inp_ind_right], seq[hit_ind_right], hit_ind_right)
        score_right = score_right + r_score
#         print 'rscore',r_score
        score_list_right.append(score_right)
        prob_list_right.append(r_score)
#         print "score", score_right
        inp_ind_right = inp_ind_right + 1
        hit_ind_right = hit_ind_right + 1
        if hit_ind_right == len(seq) or inp_ind_right == len(inp_query):
            break
            
#     print score_list_right
#     print prob_list_right
    plt.plot(score_list_right)
    plt.ylabel('Score from Right of Hit')
    plt.xlabel('Towards Right of Hit Index')
    plt.show()
    
    return score_list_left, score_list_right, prob_list_left, prob_list_right


    #     test_list1.append(inp_query[inp_ind])
    #     test_list2.append(seq[hit_ind])

    


# In[414]:

print hits_list[2]


# In[528]:

# ACCTCAGACCCCTCCCT-TCTCCAC-TCGGGCTTTTCTCCCATCTTTG input query
# AACTAACCACCACCCCTG-TCTCCAC-TCACCGGAACAGAGACTCCCC sequence
# inp_query = 'ACCTCAGACCCCTCCCTTCTCCACTCGGGCTTTTCTCCCATCTTTG'
# inp_query = 'AACTAACCACCACCCCTGTCTCCACTCACCGGAACAGAGACTCCCC'
inp_query="CTGAGGGGCCGGGGAAGGGGGCCAATCCTACCTTGGGGTGACTGAGAGGCTCAGGGAGGTGGCTCACTCCTGGGCTGCCATCAGGCTTTCCCAGGGAGCGTTCTGGCACCTGATGACTGAGGCCTGGGTGTTAACCCTGGCACTCCACTTCCTGCCTCCCTGATTTCTCCCCCTTTTTGGGGTGTCGTATATGAAGTGAGGTAACAGACGTGATGGTGCTGGCTGAGCCGCTGGCCCTGGTGTGGGGTGGGTGGAGTTGGGCCTGAAATGTGGGCATCCAGGCTGAGTGAGGCCCTCCCCGCCACCCCCAGGGCCTGGGCTGACGCCTCTCAGCAGCCCATTCTAGCCCCTCCCCCAGGAGCCAGGACCTACGTGTAGACTTTGTCACCAAT"
inp_ind = hits_list[6][0]
# inp_ind = 17
print index_hits_list[6][2]
hit_ind = index_hits_list[6][2]
# hit_ind = 19
prob_hit = prob_hits_list[6][2]
print prob_hit


# In[529]:

score_list_left, score_list_right, prob_list_left, prob_list_right= score_hits(inp_ind,hit_ind)


# In[530]:

_len=len(inp_query)
print _len
# _len = 120
# threshold = _len/10
threshold = 8
print threshold


# In[531]:

#Towards right of the hit
threshold = 8
if len(score_list_right) < threshold:
    threshold = 1
    print threshold
if not score_list_right:
    i_right = 0
    prob_right = 0
else:
    for i in range(len(score_list_right)):
        print i, i+threshold
#         print score_list_right[i]
        if score_list_right[i] == score_list_right[i+threshold]:
            print 'Will stop'
            i_right = i-1
            prob_right = score_list_right[i]
            break
        else:
            print 'We are good'


# In[535]:

#Towards left of the hit
threshold = 8
print len(score_list_left)
if len(score_list_left) < threshold:
    threshold = 1
#     print threshold
# print score_list_left
if not score_list_left:
    i_left = 0
    prob_left = 0
else:
    for i in range(len(score_list_left)):
        print i, i_left+threshold
#         print score_list_left[i]
        if score_list_left[i] == score_list_left[i+threshold]:
            print 'Stop now.'
            i_left = i+1
            prob_left = score_list_left[i]
            break
        else:
            print 'We are good'


# In[533]:

print inp_query[inp_ind-i_left:inp_ind+wmer_size+i_right]
print seq[hit_ind-i_left:hit_ind+wmer_size+i_right]
print prob[hit_ind-i_left:hit_ind+wmer_size+i_right]
print prob_list_right[:i_right]


# In[534]:

print prob_hit
print prob_list_left
print prob_list_right


# In[522]:

prob_final = 1
for item in prob_list_right[:i_right]:
    if item != 0:
        prob_final = prob_final*item
print prob_final


# In[ ]:

if __name__ == '__main__':


# In[ ]:




# In[ ]:



