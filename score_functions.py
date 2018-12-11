from __future__ import division

import charge_calculator

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

