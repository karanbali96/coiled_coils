#!/usr/bin/python

import math, getopt, sys

import simulated_annealing
import mut_functions
import score_functions

perm_res = [['I'],
       ['K','H', 'S', 'T', 'Y', 'R', 'D', 'E'],
       ['E', 'K','H', 'S', 'T', 'Y','R', 'D'],
       ['L'],
       ['E'],
       ['H','S', 'T', 'Y', 'R', 'D', 'E', 'K'],
       ['K']]

perm_res_concat = perm_res + perm_res + perm_res + perm_res


seq = mut_functions.random_seq(perm_res_concat)
simulated_annealing.anneal(seq, score_functions.score_func, mut_functions.point_mut, perm_res_concat)






