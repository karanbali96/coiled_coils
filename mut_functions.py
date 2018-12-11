from __future__ import division
from Bio.Data import IUPACData
from Bio.Seq import MutableSeq
from random import randint
import copy



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

