import re
import planet_prop as pp
import glob
import numpy as np
import matplotlib.pyplot as plt
import sys
toNumString= lambda x: re.sub("D","E",x)
def StringToNum(string,numtype=float):
	try:
		return numtype(re.sub("D","E",string))
	except:
		return string

def read_file(fi):
	with open(fi,"r") as f:
		return [ l.split() for l in f.readlines() ]


def make_history_dictionary(histfi):
	histlines = read_file(histfi)
	hist_heads=histlines[5]
	hist_data=np.array(histlines[6:],dtype=float).T
	histdict = {}
	for i,k in enumerate(hist_heads):
		histdict.update({k:hist_data[i]})
	return histdict
	

class profile():
	def __init__(self,fi):
		profilelines = read_file(fi)

		self.headers = profilelines[5]
		self._data = np.array(profilelines[6:],dtype=float).T
		
		self.metadata_headers = profilelines[1]
		self._metadata = np.array(profilelines[2],dtype=float)
		
		self.age = self._metadata[self.metadata_headers.index("star_age")]
	def get(self,string):
		return self._data[self.headers.index(string)]

def ParseMESAModelFile(finame,HLINE_N = 14):
	with open(finame,'r') as fi:
		lines = [l.split() for l in fi.readlines()  ] # remove comment lines: if not re.match("^\!",l)
	

	headers = lines[HLINE_N]
	
	scalar_properties_lines = lines[4:HLINE_N-2]
	scalar_properties = [ l[:2] for l in lines[4:HLINE_N-2] ]

	properties = dict()

	for prop in scalar_properties:
		entry = StringToNum(prop[1])
		properties.update( {prop[0]:entry})
	
	datalines = lines[HLINE_N+1:-8]
	datalines = [map(toNumString,y) for y in datalines]
	
	ndata = np.array(datalines,dtype=float)
	ndata = ndata[:,1:]
	

	for i,h in enumerate(headers):
		properties.update({h:ndata[:,i]})

	return properties
if __name__=="__main__":
	histfile = sys.argv[1]
	histdict = make_history_dictionary(histfile)

	fig,axarr=plt.subplots(3,sharex=True)
	
	t_cen =0.5 * (histdict['star_age'][1:] + histdict['star_age'][:-1] )
	E_cen = 0.5 * (histdict['internal_energy'][1:] + histdict['internal_energy'][:-1])
	dt = histdict['star_age'][1:] - histdict['star_age'][:-1]
	dE = histdict['internal_energy'][1:] - histdict['internal_energy'][:-1]
	Edot = dE/dt
	tkh = -1*E_cen/Edot 
	
	# KH timescale
	axarr[0].loglog(histdict['star_age'],np.abs(histdict['kh_timescale_fixed']),label='kh timescale')
	axarr[0].loglog(t_cen,tkh,label='E/Edot')
	axarr[0].loglog(histdict['star_age'],3*histdict['star_age'],label='3x age')
	axarr[0].legend(loc=2)
	axarr[0].set_xlim(xmin=1e6)
	# Internal Energy
	axarr[1].loglog(histdict['star_age'],np.abs(histdict['internal_energy']),label='internal energy')
	axarr[1].set_ylabel("internal energy")
	axarr[1].set_xlim(xmin=1e6)
	# Luminosity
	dt = histdict['star_age'][1:] - histdict['star_age'][:-1]
	t_cen =0.5 * (histdict['star_age'][1:] + histdict['star_age'][:-1] )
	dE = histdict['internal_energy'][1:] - histdict['internal_energy'][:-1]
	
	axarr[2].plot(histdict['star_age'],np.log10(histdict['internal_luminosity']),label='L_int')
	axarr[2].plot(t_cen ,np.log10(np.abs( Edot / pp.lsun / pp.secyer ) ),label='dE/dt')
	axarr[2].set_ylabel("log10(internal luminosity)")
	axarr[2].legend()
	
	axarr[2].set_xlim(xmin=1e6)
	
	plt.show()
