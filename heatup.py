from glob import glob
import planet_prop as pp
import os
import SH_mesa_build as mb
import numpy as np
import sys

Nmass = int(sys.argv[1])
NXenv = int(sys.argv[2])
model = glob(os.getcwd()+"/models/*Mass_%02d_Xenv_%02d_L_00.mod"%(NXenv,Nmass))[0]
models = glob(os.getcwd()+"/models/*Mass_%02d_Xenv_%02d_L_*.mod"%(NXenv,Nmass))
import re 
maxdone = np.max(np.array(map( lambda x: re.search(".*_L_(\d+)",x).groups() , models ) ,dtype=int) )
print maxdone

NL = 25

M1,M2 = 10.,20.
X1,X2 = 0.01,0.5
NM,NX = 3,3
L_M_whittle=6e-8
L_M_min=6e-8
L_M_max=0.6e2


sep=0.1
Tstar= 5780 
Rstar= 1.0 * pp.rsun 
Teq=Tstar*(Rstar/(2.0*sep*pp.au_to_cm))**0.5


mb.heat_up_model(model,L_M_min,L_M_max,NL,Teq,Xbondi=0.11) 
