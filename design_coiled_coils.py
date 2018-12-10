#!/usr/bin/python

import math, getopt, sys

import simulated_annealing

perm_res = [['I'],
       ['K','H', 'S', 'T', 'Y', 'C', 'R', 'D', 'E'],
       ['E', 'K','H', 'S', 'T', 'Y', 'C', 'R', 'D'],
       ['L'],
       ['E'],
       ['H','S', 'T', 'Y', 'C', 'R', 'D', 'E', 'K'],
       ['K']]


seq = simulated_annealing.random_seq(perm_res)
simulated_annealing.anneal(seq, simulated_annealing.score_func, simulated_annealing.point_mut, perm_res)






