import numpy as np
import scipy.special as sp


"""
Fucntion that retrns the number of steps required for a chain of ON/OFF states
to all go ON starting from OFF
"""

def ptree(nodes, p):
    if nodes == 0:
        return []

    weights = [(1-p)**i * p**(nodes-i)*sp.binom(nodes, i) for i in range(nodes + 1)]
    conf = np.arange(nodes, -1, -1)

    nodes_on = np.random.choice(conf, p=weights)

    return [nodes_on] + ptree(nodes - nodes_on, p)
