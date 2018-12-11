from __future__ import division
from Bio.Data import IUPACData
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_protein
import numpy as np
import random
from random import randint
import math
import copy

import charge_calculator

def random_seq(l):
    seq = ""
    for p in l:
        r_num = randint(0, len(p) - 1)
        seq = seq + p[r_num]
    return MutableSeq((seq), IUPACData.protein_letters)


def point_mut(seqin, l):
    seq = copy.deepcopy(seqin)
    pos_num = randint(0, len(seq) - 1)
    p = l[pos_num]
    q = randint(0, len(p) - 1)
    seq[pos_num] = p[q]
    return seq


def score_func(seq):
    chargeobj = charge_calculator.Charge_Calculator(seq)
    pH_x = [5.8, 9]
    charges = chargeobj.charge_values(pH_x)

    c1 = charges[0]
    c2 = charges[1]

    s1 = c1**2
    s2 = c2**2

    w1 = 1
    w2 = 1

    total_score = -(1/(w1*s1)) - w2*s2

    return total_score

def accept_probability(old_score, new_score, T):
       if new_score < old_score:
              return 1.0
       else:
              return math.exp((old_score - new_score)/T)

def anneal(seq, score_func_, mut_func_, perm_res):
    original = copy.deepcopy(seq)
    old_score = score_func_(seq)
    T = 1.0
    T_min = 0.00001
    alpha = 0.99
    i = 1
    scores = []
    iterations = []
    while T > T_min:
        new_sol = mut_func_(seq, perm_res)
        #print(sol, new_sol)
        new_score = score_func_(new_sol)
        ap = accept_probability(old_score, new_score, T)
        #print(ap, old_score, new_score, T)
        randnum = random.uniform(0,1)
        #print(ap, randnum)
        if ap > randnum:
            print("accepted")
            seq = new_sol
            old_score = new_score
        else:
            print("rejected")
        print(i, T, seq, score_func_(seq)) # old_score)
        scores.append(score_func_(seq))
        iterations.append(i)
        T = T*alpha
        i +=1
    np.savetxt('start_' + str(original) + '_end_' + str(seq) + '.dat', np.column_stack((iterations, scores)), fmt='%f')
    return seq, score_func_(seq)










