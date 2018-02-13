import numpy as np
import matplotlib.pyplot as plt

def heaviside(k,x):
    return 2.0*(1.0/(1.0+pow(np.e,-k*x)) - 0.5 )

def mirrorHeaviside(k,x,s):
    return 1.0/(1.0+pow(np.e,k*(x-s)))
