#!/usr/bin/python

import math, getopt, sys

import simulated_annealing
import mut_functions

perm_res = [['I'],
       ['K','H', 'S', 'T', 'Y', 'C', 'R', 'D', 'E'],
       ['E', 'K','H', 'S', 'T', 'Y', 'C', 'R', 'D'],
       ['L'],
       ['E'],
       ['H','S', 'T', 'Y', 'C', 'R', 'D', 'E', 'K'],
       ['K']]


seq = mut_functions.random_seq(perm_res)
simulated_annealing.anneal(seq, simulated_annealing.score_func, mut_functions.point_mut, perm_res)






