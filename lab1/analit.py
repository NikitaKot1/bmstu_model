import numpy as np

def analit_1(u):
    return 3*np.exp(u) - u**2 - 2*u - 2

def analit_2(u):
    return np.exp(u**2) - (u**2 + 1) / 2