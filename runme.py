import os
import planet_prop as pp
import SH_mesa_build as mb
import numpy as np

#M1,M2 = 5.,20.
M1,M2 = 1.84210526,20.
X1,X2 = 0.001,0.5
NM,NX = 24,20
L_M_whittle=6e-8
L_M_min=6e-8
L_M_max=1e2

sep=0.1
Tstar= 5780 
Rstar= 1.0 * pp.rsun 
Teq=Tstar*(Rstar/(2.0*sep*pp.au_to_cm))**0.5

#mb.build_starting_Mcore_Xenv_grid(M1,M2,NM,X1,X2,NX,Teq,L_M_min,os.getcwd()+"/models",restart=False) 
mb.build_starting_Mcore_Xenv_grid(M1,M2,NM,X1,X2,NX,Teq,L_M_min,os.getcwd()+"/models",restart=True) 
