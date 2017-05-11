import planet_prop as pp
import mesa_build as mb
import numpy as np
import mesa_namelist as mn
import os
import shutil
import tempfile

from argparse import ArgumentParser
parser = ArgumentParser(description='Evolve MESA planet models')

# Model specification
parser.add_argument('-M','--mass', metavar='N', type=int, default=0, help='core mass model number')
parser.add_argument('-X','--Xenv', metavar='N', type=int, default=-1, help='envelope fraction model number')
parser.add_argument('-L','--luminosity', metavar='N', type=int, default=-1, help='luminosity model number')

# Job type
parser.add_argument('--job_type', default=None, help='Specify job type')

# Input parameters
parser.add_argument('-I','--input_directory',metavar="DIR",default="./",help="directory of input models")
parser.add_argument('-O','--output_directory',metavar="DIR",default="./",help="directory for output models")
parser.add_argument('--separation', metavar='FLOAT', type=float, default=0.1, help='planet-star separation')
parser.add_argument('--max_age', metavar='FLOAT', type=float, default=.7E10, help='maximum age to evolve model')
parser.add_argument('--cadence', metavar='FLOAT', type=float, default=2E8, help='output cadence')
parser.add_argument('--Sigma', metavar='FLOAT', type=float, default=250, help='Column depth for "F-Sigma" irradiation')
parser.add_argument('--no_core', default=False, action='store_true', help='Turn off contribution of core\'s heat capacity')
parser.add_argument('--redo', default=False, action='store_true', help='Re-compute models that already have existing output files')
parser.add_argument('--hottest', default=False, action='store_true', help='evolve only the hottest model')


args = parser.parse_args()

print "Job type: %s"%args.job_type

INTIALIZE=(args.job_type=='initialize')
EVOLVE=(args.job_type=='evolve')
EVAPORATE=(args.job_type=='evaporate')
HEATUP=(args.job_type=='heatup')


NOJOB = not np.any([INTIALIZE,EVAPORATE,EVOLVE,HEATUP])

if NOJOB:
	print "No job type specified."

NMASS = args.mass
NXENV = args.Xenv
NLUM = args.luminosity
All_Xenv=(NXENV<0)
All_Lum=(NLUM<0) or HEATUP
	

indir=os.path.abspath(args.input_directory)+"/"
outdir=os.path.abspath(args.output_directory)+"/"
print "Input directory set to: %s"%indir
print "Output directory set to: %s"%outdir

# Variable controls
Tmax=args.max_age  # Maximum Age
cadence=args.cadence # Output cadence
Sigma=args.Sigma # heating depth in F-Sigma routine
sep=args.separation # Planet-star separation

# Constants
sb=5.6704e-5 # Steffan-Boltzmann Constant 
Tstar= 5780 # Stellar Temperature 
Rstar= 1.0 * pp.rsun # Stellar Radius
Teq=Tstar*(Rstar/(2.0*sep*pp.au_to_cm))**0.5 # Equilibrium Temprature
Flux=4.0*sb*Teq**4.0 # Incident flux
Thold=1e14 # holding time

# create temporary run directory
curr_dir=os.getcwd()

Mcore=np.linspace(1.84210526,20.,24)[-(NMASS+1)] # Planet Core Mass
#Mcore=np.linspace(5.,20.,20)[-(NMASS+1)] # Planet Core Mass
min_mass=(1+1e-5)*Mcore*pp.mearth/pp.msun # Minimum mass envelope stopping condition
L_M_min=6e-8
L_M_max=0.6e2
Lgrid = np.logspace(np.log10(L_M_min),np.log10(L_M_max),25)

# return to current directory and delete temporary run directory
import SH_mesa_build as shmb
import SH_mesa_namelist as shnm
import glob
import re


if All_Xenv:
	N_Xenv_N_L_list = np.array( [re.search("Xenv_(\d+)_L_(\d+)",f).groups() for f in glob.glob(indir+"/model_Mass_%02d_*"%NMASS) ]   ,dtype=int)
elif All_Lum:
	N_Xenv_N_L_list = np.array( [re.search("Xenv_(\d+)_L_(\d+)",f).groups() for f in glob.glob(indir+"/model_Mass_%02d_Xenv_%02d_*"%(NMASS,NXENV)) ] ,dtype=int)
else:
	N_Xenv_N_L_list = np.array( [re.search("Xenv_(\d+)_L_(\d+)",f).groups() for f in glob.glob(indir+"/model_Mass_%02d_Xenv_%02d_L_%02d*"%(NMASS,NXENV,NLUM)) ] ,dtype=int)


if HEATUP: 

	dirpath=shmb.get_temporary_work_directory()
	print "temporary run directory:"
	print dirpath
	os.chdir(dirpath)


	history_temp_out = "./temp_hist.dat"
	cs=(1.38e-16*Teq/(2.35*1.67e-24))**0.5


	for NXENV in set(N_Xenv_N_L_list[:,0]):
		msk = N_Xenv_N_L_list[:,0]==NXENV
		max_L_done = np.max(N_Xenv_N_L_list[msk][:,1])
		mod_file = indir+"model_Mass_%02d_Xenv_%02d_L_%02d.mod"%(NMASS,NXENV,max_L_done)
		for NLUM in range(max_L_done+1,25):
			print "Heating model X=%02d L=%02d"%(NXENV,NLUM)
			mass,radius=shmb.get_mass_radius(mod_file)
			Rbondi=6.67e-8*mass/(2.0*cs**2.0)
			XB = radius / Rbondi
			print "XB of last model: %g / %g = %g"%(radius,Rbondi,XB)
			if XB > 0.11:
				#os.remove(mod_file)
				break
			mod_out = outdir+"model_Mass_%02d_Xenv_%02d_L_%02d.mod"%(NMASS,NXENV,NLUM)
			# heat
			L_M_want = Lgrid[NLUM]
			nml=shnm.build_namelist_mp_hold(Thold,L_M_want,Flux,Sigma)
			shmb.run_mesa_namelist(nml,mod_file,mod_out,history_temp_out)

			# finalize
			nml=shnm.build_namelist_ev_Fsigma(1e4,Flux,Sigma,False)
			shmb.run_mesa_namelist(nml,mod_out,mod_out,history_temp_out)

			mod_file = mod_out

	os.chdir(curr_dir)
	shutil.rmtree(dirpath)

elif not NOJOB:
	print "%d input files found"%len(N_Xenv_N_L_list)
	dirpath=shmb.get_temporary_work_directory()
	print "temporary run directory:"
	print dirpath
	os.chdir(dirpath)

	for NXENV,NLUM in N_Xenv_N_L_list:
		mod_file = indir+"model_Mass_%02d_Xenv_%02d_L_%02d.mod"%(NMASS,NXENV,NLUM)
		L_M_ratio = Lgrid[NLUM]
		
		# Relax initial model to specified current irradiation
		if INTIALIZE:
			nml=shnm.build_namelist_mp_hold(Thold,L_M_ratio,Flux,Sigma,log_directory=dirpath)
			hist_out = outdir+"../delete_me.dat" #outdir+"hist_Mass_%02d_Xenv_%02d_L_%02d.mod"%(NMASS,NXENV,NLUM)
			mod_out = outdir+"model_Mass_%02d_Xenv_%02d_L_%02d.mod"%(NMASS,NXENV,NLUM)
		# Simulate evolution
		if EVOLVE:
			if args.hottest:
				msk = N_Xenv_N_L_list[:,0]==NXENV
				max_L_done = np.max(N_Xenv_N_L_list[msk][:,1])
				if NLUM != max_L_done:
					continue
			nml=shnm.build_namelist_ev_Fsigma(Tmax,Flux,Sigma,False,log_directory=dirpath,HIST_CADENCE=cadence,include_core=(not args.no_core))
			hist_out=outdir+"ne_hist_Mass_%02d_Xenv_%02d_L_%02d.dat"%(NMASS,NXENV,NLUM)
			mod_out=outdir+"ne_model_Mass_%02d_Xenv_%02d_L_%02d.mod"%(NMASS,NXENV,NLUM)
		# Simulate photoevaporation
		if EVAPORATE:
			nml=shnm.build_namelist_ev_Fsigma_we(Tmax,Flux,Sigma,Teq,sep,min_mass,False,log_directory=dirpath,HIST_CADENCE=cadence,include_core=(not args.no_core))
			hist_out=outdir+"evap_hist_Mass_%02d_Xenv_%02d_L_%02d.dat"%(NMASS,NXENV,NLUM)
			mod_out=outdir+"evap_model_Mass_%02d_Xenv_%02d_L_%02d.mod"%(NMASS,NXENV,NLUM)
			
		if (not os.path.isfile(mod_out)) or args.redo:
			print "running model %s"%mod_out
			shmb.run_mesa_namelist(nml,mod_file,mod_out,hist_out)
			print "saved model file %s"%(mod_out)
		else:
			print "%s exists, skipping to next model..."%(mod_out)
	
	os.chdir(curr_dir)
	shutil.rmtree(dirpath)
