
import pandas as pd 
import numpy as np 

from hdna import *
from tqdm import tqdm

from scipy.optimize import dual_annealing

EXPNAME = 'sumofinchandpk44'

notes = """

back to sum of rates for slidings. now try to put a cap on fwd to not let it go higher than the initially proposed value
removed /100 on length 2 slidings
Changed zipping and sliding slightly
sliding up to 2 
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
    min_nucleation=1)
MOD.setgeometry(theta=90, phi = 120)
 
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

zipping = 7.7e7
sliding = 2e5

H.run([zipping, sliding])