
import pandas as pd 
import numpy as np 

from hdna import *
from tqdm import tqdm

from scipy.optimize import dual_annealing

EXPNAME = 'sumofinchandpk55'

notes = """
Trying some last times
"""

# Import experimental data from Hertel 
expdata = pd.read_csv('./data/herteldata.csv', names=['sequences', 'experimental'])
# Clean the dataframe 
expdata = expdata.drop(0)
expdata['experimental'] = ['{:e}'.format(float(e)) for e in expdata['experimental']]

limit = len(expdata)
data = expdata.copy().iloc[:limit]
data['index'] = data.index 
data.set_index(data['sequences'], inplace=True)

MOD = Model(
    stacking='nostacking',
    Na=0.15,
    min_nucleation=1)
MOD.setgeometry(theta=120, phi = 120)
 
MOD.alpha = 1
MOD.gamma = 0
MOD.kappa = 1

OPT = Options(Nsim=5000)

H = HDNA(data, EXPNAME, model=MOD, options=OPT)

with open(f'results/{EXPNAME}/notes.txt', 'w') as savenote:
    savenote.write(notes)
    savenote.close()


# zipping = 7.7e7
# sliding = 1.35e5 #zipping*np.exp(-(2)/(CONST.R*MOD.kelvin))

zipping = 8.8e7
sliding = 4e5

H.run([zipping, sliding])