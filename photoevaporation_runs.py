import planet_prop as pp
import mesa_build as mb
import numpy as np
import mesa_namelist as mn
import os
import shutil
import tempfile
import sys

All_Xenv=True
NMASS=int(sys.argv[1])
if len(sys.argv)>2:
	NXENV=int(sys.argv[2])
	All_Xenv=False
if len(sys.argv)>3:
	NLUM=int(sys.argv[3])
	All_Lum=False

topdir="/projects/p20783/sjh890/01_MESA_Projects/planet_mesa/05/models/"
outdir="/projects/p20783/sjh890/01_MESA_Projects/planet_mesa/05/evap_models/"

Tmax=8.e9 # Maximum Age
sb=5.6704e-5 # Steffan-Boltzmann Constant 
Sigma=250 # heating depth in F-Sigma routine
sep=0.1 # Planet-star separation
Tstar= 5780 # Stellar Temperature 
Rstar= 1.0 * pp.rsun # Stellar Radius
Teq=Tstar*(Rstar/(2.0*sep*pp.au_to_cm))**0.5 # Equilibrium Temprature
Flux=4.0*sb*Teq**4.0 # Incident flux
HIST_CADENCE = 5.e7 # history output cadence

def my_build_namelist_ev_Fsigma_we(Tmax,Flux,Sigma,Temp,sep,stop_mass,ES,log_directory="LOGS"):

	import f90nml #fortran namelist tools
	from f90nml.namelist import Namelist as NmlDict #namelist class

	#build namelist
	nml=NmlDict()
	sj=NmlDict() #star_job
	cn=NmlDict() #controls
	pg=NmlDict() #pgplot

	#declare some xctrl arrays
	x_ctrl=[0.0]*30
	x_integer_ctrl=[0]*5
	x_logical_ctrl=[False]*10

	x_ctrl[8]=1 #safety factor
	x_ctrl[9]=1 #safety factor
	x_ctrl[11]=sep # sep AU
	x_ctrl[12]=1 # stellar mass Msun
	x_ctrl[18]=Temp
	x_ctrl[23]=6.9183e7 #saturation for X-rays
	x_ctrl[24]=10**-3.6
	x_ctrl[25]=1.19


	x_logical_ctrl[0]=True # radio-active heating
	x_logical_ctrl[1]=True # core heat capacity
	x_logical_ctrl[5]=False # core/envelope heating
	x_logical_ctrl[6]=True #evaporation
	x_logical_ctrl[7]=ES # include evolution of star or not

	x_logical_ctrl[8]=True # uniform history data
	x_ctrl[0]= HIST_CADENCE # history cadence
	cn['max_years_for_timestep'] = HIST_CADENCE # history cadence

	pg['read_extra_pgstar_inlist1']= False 
	sj['read_extra_star_job_inlist1']= False 

	sj['show_log_description_at_start'] = False
	sj['load_saved_model'] = True
	sj['saved_model_name'] = 'model_in.mod'
	sj['save_model_when_terminate'] = True
	sj['save_model_filename'] = 'model_out.mod'
	sj['set_initial_age'] = True
	sj['initial_age'] = 0
	sj['set_initial_dt']=True
	sj['years_for_initial_dt']=1e5
	sj['pgstar_flag'] = False
	sj['change_v_flag'] = True
	sj['new_v_flag'] = True
	sj['eos_file_prefix'] = 'mesa'
	sj['kappa_file_prefix'] = 'gs98'
	sj['kappa_lowT_prefix'] = 'lowT_Freedman11'
	sj['change_lnPgas_flag'] = True
	sj['new_lnPgas_flag'] = True
	sj['set_initial_model_number'] = True
	sj['initial_model_number'] = 0
	sj['set_irradiation'] = True
	sj['set_to_this_irrad_flux']= Flux
	sj['irrad_col_depth'] = Sigma

	cn['log_directory']=log_directory
	cn['photo_directory']=log_directory
	cn['photostep']=1000

	cn['T_mix_limit'] = -1
	cn['logQ_limit'] = 10
	cn['mixing_length_alpha'] = 1.89
	cn['MLT_option'] = 'Henyey'
	cn['write_header_frequency'] = 10
	cn['terminal_cnt'] = 100
	cn['profile_interval'] = 10000
	cn['varcontrol_target'] = 3e-4
	cn['mesh_delta_coeff'] = 1
	cn['max_age'] = Tmax
	cn['max_years_for_timestep']=5e7
	cn['use_other_energy']=True
	cn['use_other_wind']=True
	cn['star_mass_min_limit']=stop_mass
	cn['x_ctrl']=x_ctrl
	cn['x_logical_ctrl']=x_logical_ctrl
	cn['x_integer_ctrl']=x_integer_ctrl
        cn['min_q_for_k_below_const_q'] = 0.9
        cn['min_q_for_k_const_mass'] = 0.5

	nml['star_job']=sj
	nml['controls']=cn
	nml['pgstar']=pg

	return nml;	




# create temporary run directory
curr_dir=os.getcwd()

Mcore=np.linspace(5.,20.,20)[-(NMASS+1)] # Planet Core Mass
min_mass=(1+1e-5)*Mcore*pp.mearth/pp.msun # Minimum mass envelope stopping condition

# return to current directory and delete temporary run directory
import SH_mesa_build as shmb
dirpath=shmb.get_temporary_work_directory()
print "temporary run directory:"
print dirpath
os.chdir(dirpath)
import glob
import re


if All_Xenv:
	N_Xenv_N_L_list = np.array( [re.search("Xenv_(\d+)_L_(\d+)",f).groups() for f in glob.glob(topdir+"model_Mass_%02d_Xenv_*.mod"%NMASS) ]   ,dtype=int)
elif All_Lum:
	N_Xenv_N_L_list = np.array( [re.search("Xenv_(\d+)_L_(\d+)",f).groups() for f in glob.glob(topdir+"model_Mass_%02d_Xenv_%02d_*.mod"%(NMASS,NXENV)) ] ,dtype=int)
else:
	N_Xenv_N_L_list = np.array( [re.search("Xenv_(\d+)_L_(\d+)",f).groups() for f in glob.glob(topdir+"model_Mass_%02d_Xenv_%02d_L_%02d*"%(NMASS,NXENV,NLUM)) ] ,dtype=int)

for NXENV,NLUM in N_Xenv_N_L_list:
	mod_file = topdir+"model_Mass_%02d_Xenv_%02d_L_%02d.mod"%(NMASS,NXENV,NLUM)
	print mod_file
	# make the namelist
	nml=my_build_namelist_ev_Fsigma_we(Tmax,Flux,Sigma,Teq,sep,min_mass,False,log_directory=dirpath)
	hist_out = outdir+"evap_hist_Mass_%02d_Xenv_%02d_L_%02d.mod"%(NMASS,NXENV,NLUM)
	mod_out = outdir+"evap_model_Mass_%02d_Xenv_%02d_L_%02d.mod"%(NMASS,NXENV,NLUM)
	shmb.run_mesa_namelist(nml,mod_file,mod_out,hist_out)
	## return to current directory and delete temporary run directory
os.chdir(curr_dir)
shutil.rmtree(dirpath)

