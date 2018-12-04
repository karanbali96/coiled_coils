from __future__ import division
from Bio.Data import IUPACData
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_protein
import numpy as np
import random
from random import randint
import math
import copy


positive_pKs = {'Nterm': 7.5, 'K': 10.0, 'R': 12.0, 'H': 5.98}
negative_pKs = {'Cterm': 3.55, 'D': 4.05, 'E': 4.45, 'C': 9.0, 'Y': 10.0}
pKcterminal = {'D': 4.55, 'E': 4.75}
pKnterminal = {'A': 7.59, 'M': 7.0, 'S': 6.93, 'P': 8.36, 'T': 6.82, 'V': 7.44,
               'E': 7.7}
charged_aas = ('K', 'R', 'H', 'D', 'E', 'C', 'Y')

l = [['I'],
       ['K','H', 'S', 'T', 'Y', 'C', 'R', 'D', 'E'],
       ['E', 'K','H', 'S', 'T', 'Y', 'C', 'R', 'D'],
       ['L'],
       ['E'],
       ['H','S', 'T', 'Y', 'C', 'R', 'D', 'E', 'K'],
       ['K']]
class Charge(object):

    def __init__(self, ProteinSequence):
        """Initialize the class."""
        self.sequence = ProteinSequence
        self.amino_acids_content = self.count_amino_acids(self.sequence)
        self.charged_aas_content = self._select_charged(self.amino_acids_content)

    def count_amino_acids(self, seq):
        """Count standard amino acids, return a dict.

        Counts the number times each amino acid is in the protein
        sequence. Returns a dictionary {AminoAcid:Number}.

        The return value is cached in self.amino_acids_content.
        It is not recalculated upon subsequent calls.
        """
        #if self.amino_acids_content is None:
        prot_dic = dict((k, 0) for k in IUPACData.protein_letters)
        for aa in prot_dic:
            prot_dic[aa] = seq.count(aa)

            #self.amino_acids_content = prot_dic

        return prot_dic


    # This function creates a dictionary with the contents of each charged aa,
    # plus Cterm and Nterm.
    def _select_charged(self, AminoAcidsContent):
        charged = {}
        for aa in charged_aas:
            charged[aa] = float(AminoAcidsContent[aa])
        charged['Nterm'] = 1.0
        charged['Cterm'] = 1.0
        return charged

    # This function calculates the total charge of the protein at a given pH.
    def _chargeR(self, pH, pos_pKs, neg_pKs):
        PositiveCharge = 0.0
        for aa, pK in pos_pKs.items():
            CR = 10 ** (pK - pH)
            partial_charge = CR / (CR + 1.0)
            PositiveCharge += self.charged_aas_content[aa] * partial_charge

        NegativeCharge = 0.0
        for aa, pK in neg_pKs.items():
            CR = 10 ** (pH - pK)
            partial_charge = CR / (CR + 1.0)
            NegativeCharge += self.charged_aas_content[aa] * partial_charge

        return PositiveCharge - NegativeCharge

    def charge_values(self, pH_x):

        pos_pKs = dict(positive_pKs)
        neg_pKs = dict(negative_pKs)
        nterm = self.sequence[0]
        cterm = self.sequence[-1]
        if nterm in pKnterminal:
            pos_pKs['Nterm'] = pKnterminal[nterm]
        if cterm in pKcterminal:
            neg_pKs['Cterm'] = pKcterminal[cterm]

        #pH_x = [5.8, 9]
        charge_y = []


        for pH in pH_x:
            charge = self._chargeR(pH, pos_pKs, neg_pKs)
            charge_y.append(charge)


        #print('charge at 5.8 =', charge_y[0], 'charge at 9 =', charge_y[1])



        #print(s1)
        #print(s2)

        return charge_y

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
    chargeobj = Charge(seq)
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

def anneal(sol):
    original = copy.deepcopy(sol)
    old_score = score_func(sol)
    T = 1.0
    T_min = 0.00001
    alpha = 0.99
    i = 1
    scores = []
    iterations = []
    while T > T_min:
        new_sol = point_mut(sol, l)
        #print(sol, new_sol)
        new_score = score_func(new_sol)
        ap = accept_probability(old_score, new_score, T)
        #print(ap, old_score, new_score, T)
        randnum = random.uniform(0,1)
        #print(ap, randnum)
        if ap > randnum:
            print("accepted")
            sol = new_sol
            old_score = new_score
        else:
            print("rejected")
        print(i, T, sol, score_func(sol)) # old_score)
        scores.append(score_func(sol))
        iterations.append(i)
        T = T*alpha
        i +=1
    np.savetxt('start_' + str(original) + '_end_' + str(sol) + '.dat', np.column_stack((iterations, scores)), fmt='%f')
    return sol, score_func(sol)










