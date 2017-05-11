import shutil
import numpy as np
import planet_prop as pp   # tools for planetary properties
import os

# collection of mesa functions for building and running the code
SRCDIR="/projects/p20783/sjh890/01_MESA_Projects/planet_mesa/00_MyScripts"

def fortney_rRockIron(rmf, logm):
	return 0.7932 + logm**2*(0.0975 + 0.0592*rmf) + logm*(0.4938 + 0.2337*rmf) + 0.3102*rmf

def get_core_density(Mcore):
	Rcore= fortney_rRockIron( 2./3. , np.log10(Mcore) )   
	rho_core=Mcore*pp.mearth/(4.0/3.0*np.pi*(Rcore*pp.rearth)**3.0)
	rho_core=(np.float64(rho_core).item())
	return rho_core

# create temporary run directory
def get_temporary_work_directory():
	import tempfile
	
	dirpath = tempfile.mkdtemp()
	print "temporary run directory:"
	print dirpath
	#os.chdir(dirpath)
	os.system("ln -s %s/Owen_Rates %s"%(SRCDIR,dirpath)) 
	os.system("ln -s %s/make_planets/create_planet %s"%(SRCDIR,dirpath)) 
	return dirpath

def run_mesa_namelist(nml,model_in,model_out,hist_out,SETUP_ONLY=False):
	import shutil
	shutil.copyfile(model_in,"./model_in.mod")
	nml.write('inlist',force=True)
	if not SETUP_ONLY:
		os.system('./create_planet')
		shutil.copyfile("model_out.mod",model_out)
		shutil.copyfile("history.data",hist_out)

def build_starting_Mcore_Xenv_grid(McoreLow,McoreHigh,Ncore,XLow,XHigh,Nenv,Teq,L_M_start,output_dir,start_file=SRCDIR+"/3.00_Mj.mod",restart=True):
	import f90nml #fortran namelist tools
	import SH_mesa_namelist as mn # for generating mesa namelist
	#usful constants
	sb=5.6704e-5
	kb=1.38e-16
	curr_mass=5.695799690399999e+30 # exact starting mass of 3Mjup model
	Flux=4.0*sb*Teq**4.0
	ES=False

	Sigma=250 # heating depth in F-Sigma routine
	Thold=1e14 # holding time
	Twhittle=1e14 # mass removal time-scale

	Mcore_Grid=np.linspace(McoreHigh,McoreLow,Ncore)
	Xenv_Grid=np.logspace(np.log10(XHigh),np.log10(XLow),Nenv)

	if restart:
		import glob
		import re
		donefiles = glob.glob("%s/model*.mod"%output_dir) 	
		doneindicies = np.array(map(lambda x: re.search("model_Mass_(\d+)_Xenv_(\d+)",x).groups(),donefiles),dtype=int)
		istart=np.max(doneindicies[:,0])
		print "starting core mass: %02d"%istart
		jstart=np.max(doneindicies[doneindicies[:,0]==istart][:,1])+1	
		if jstart==Nenv:
			jstart=-1
			istart=istart+1
		print "Restarting from model (%d,%d)"%(istart,jstart)
	else:
		istart=0
		jstart=-1
	
	curr_dir=os.getcwd()
	work_dir = get_temporary_work_directory()
	os.chdir(work_dir)
	for i,Mcore in enumerate(Mcore_Grid):
		if i<istart:
			print "skipping core mass %d"%i
			continue
		for j,Xenv in enumerate(Xenv_Grid):
			if j<=jstart:
				print "skipping Xenv %d"%j
				continue

			model_out = "%s/model_Mass_%02d_Xenv_%02d_L_%02d.mod"%(output_dir,i,j,0)
			history_out = "%s/history_Mass_%02d_Xenv_%02d_L_%02d.dat"%(output_dir,i,j,0)
			model_temp_out = "%s/temp_output.mod"%work_dir
			history_temp_out = "%s/temp_hist.dat"%work_dir

			if i==0 and j==0:
				init_file=start_file
			elif i==0:
				# use smallest existing envelope fraction
				init_file="%s/model_Mass_%02d_Xenv_%02d_L_%02d.mod"%(output_dir,i,j-1,0)
			else:
				init_file="%s/model_Mass_%02d_Xenv_%02d_L_%02d.mod"%(output_dir,i-1,j,0)

			# make a fresh model from the start file
			rho_core = get_core_density(Mcore)
			
			# Adjust core
			print "Adjusting core mass to %g"%Mcore
			nml = mn.build_namelist_pc( rho_core,Mcore*pp.mearth / pp.msun ,log_directory=work_dir)
			run_mesa_namelist(nml,init_file,model_temp_out,history_temp_out)

			# set the envelope fraction
			print "Setting envelope mass to %g"%Xenv
			target_mass = (1.0+Xenv) * Mcore * pp.mearth
			mdot=(target_mass-curr_mass)/Twhittle/pp.msun
			target_mass=(np.float64(target_mass).item())
			mdot=(np.float64(mdot).item())
			nml = mn.build_namelist_mp_whittle(Twhittle,L_M_start,target_mass/pp.msun,mdot,Flux,Sigma)
			run_mesa_namelist(nml,model_temp_out,model_temp_out,history_temp_out)
			
			
			# evolve to steady state holding L_M fixed
			print "Relaxing thermal state"
			nml=mn.build_namelist_mp_hold(Thold,L_M_start,Flux,Sigma)
			run_mesa_namelist(nml,model_temp_out,model_temp_out,history_temp_out)

			# evolve short time under irradiation
			print "finalizing initial model"
			nml=mn.build_namelist_ev_Fsigma(1e4,Flux,Sigma,False)
			run_mesa_namelist(nml,model_temp_out,model_out,history_out)

				#
	shutil.rmtree(work_dir)

def get_mass_radius(modelfile):
	from plot_models import ParseMESAModelFile
	modeldict=ParseMESAModelFile(modelfile,HLINE_N=15)
	radius = np.max(np.exp(modeldict['lnR']) )
	mass= modeldict['M/Msun'] * pp.msun
	return mass,radius
	

def heat_up_model(model_file,L_M_start,L_M_stop,N_L_M,Teq,Xbondi=0.1):
	# evolve to steady state holding L_M fixed
	import f90nml #fortran namelist tools
	import SH_mesa_namelist as mn # for generating mesa namelist
	import mesa as ms
	import re
	

	sb=5.6704e-5
	Flux=4.0*sb*Teq**4.0
	Sigma=250 # heating depth in F-Sigma routine	

	cs=(1.38e-16*Teq/(2.35*1.67e-24))**0.5
	Thold=1e14 # holding time
	
	curr_dir=os.getcwd()
	work_dir = get_temporary_work_directory()
	os.chdir(work_dir)

	model_temp_out = "%s/temp_output.mod"%work_dir
	history_temp_out = "%s/temp_hist.dat"%work_dir


	L_M_grid=np.logspace(np.log10(L_M_start),np.log10(L_M_stop),N_L_M)
	for i,L_M_want in enumerate(L_M_grid):
		
		model_out = re.sub("L_\d+","L_%02d"%i,model_file)
		# skip existing files
		if os.path.isfile(model_out):
			shutil.copyfile(model_out,model_temp_out)
			continue

		# heat up to target L/M
		print "relaxing to specified thermal state"
		nml=mn.build_namelist_mp_hold(Thold,L_M_want,Flux,Sigma)
		if i==0:
			run_mesa_namelist(nml,model_file,model_temp_out,history_temp_out)
		else:
			run_mesa_namelist(nml,model_temp_out,model_temp_out,history_temp_out)

		# evolve short time under irradiation
		print "finalizing  model"
		nml=mn.build_namelist_ev_Fsigma(1e4,Flux,Sigma,False)
		run_mesa_namelist(nml,model_temp_out,model_out,history_temp_out)

		mass,radius = get_mass_radius(model_out)

		Rbondi=6.67e-8*mass/(2.0*cs**2.0)
		XB=radius/Rbondi
		print "model %02d Xbondi %g"%(i,XB)
		if XB > Xbondi:
			os.remove(model_out)
			break
		else:
			print "saving model: %s"%model_out
	shutil.rmtree(work_dir)
	os.chdir(curr_dir)
