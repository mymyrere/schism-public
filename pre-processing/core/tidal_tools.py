import os
import numpy as np
from scipy import interpolate
from scipy.spatial import Delaunay
from netCDF4 import Dataset
from ttide.t_getconsts import t_getconsts
from ttide.t_vuf import t_vuf
#from matplotlib import delaunay
from ttide.t_predic import t_predic

def extrapolate_amp(lon,lat,lonr, latr, maskr, arr):
	x, y = np.ma.array(lon), np.ma.array(lat)
	x.mask, y.mask = arr.mask, arr.mask

	tri = Delaunay(x.compressed(), y.compressed())

	interp = tri.nn_extrapolator(arr.compressed())
	arr = interp(lonr, latr)
	arr = np.ma.masked_where(maskr==0, arr)

	return arr
def extrapolate_pha(lon,lat,lonr, latr, maskr, arr):
	x, y = np.ma.array(lon), np.ma.array(lat)
	x.mask, y.mask = arr.mask, arr.mask
	tri = Delaunay(x.compressed(), y.compressed())

	cx, cy = np.cos(arr), np.sin(arr)
	interp = tri.nn_extrapolator(cx.compressed())
	cx = interp(lonr, latr)
	interp = tri.nn_extrapolator(cy.compressed())
	cy = interp(lonr, latr)
	arr = np.arctan2(cy, cx)
	arr = np.ma.masked_where(maskr==0, arr)
	return arr

def extract_HC(modfile,vars, lon, lat, conlist=None, logger=None):
	"""
	Extract harmonic constituents and interpolate onto points in lon,lat
	set "z" to specifiy depth for transport to velocity conversion
	set "constituents" in conlist
	Returns:
	    u_re, u_im, v_re, v_im, h_re, h_im, omega, conlist
	"""

	###
	lon = np.asarray(lon)
	lat = np.asarray(lat)

	# Make sure the longitude is between 0 and 360
	lon = np.mod(lon,360.0)

	###
	# Read the filenames from the model file
	pathfile = os.path.split(modfile)
	path = pathfile[0]


	f = Dataset(modfile,'r')

	###
	# Read the grid file
	X=f.variables['lon'][:]
	Y=f.variables['lat'][:]
	if len(X.shape)==1:
		X, Y = np.meshgrid(X, Y)

	if 'dep' in f.variables:
		depth=f.variables['dep'][:]
	else:
		depth=np.zeros((X.shape))

	
	###
	# Check that the constituents are in the file
	conList = []
	conIDX=[]
	if conlist is not None:
		for ncon in range(0,len(f.variables['cons'])):
			if f.variables['cons'][ncon][:] in conlist:
				x=''
				conList.append(''.join([x+n for n in f.variables['cons'][ncon]]))
				conIDX.append(ncon)
	else:
		for ncon in range(0,len(f.variables['cons'])):
			x=''
			conList.append(''.join([x+n for n in f.variables['cons'][ncon]]))
			conIDX.append(ncon)

	const = t_getconsts(np.array([]))
	consindex = [np.where(const[0]['name']==con.ljust(4))[0][0] for con in conList]
	tfreq = (2*np.pi)*const[0]['freq'][consindex]/3600.

	###
	# Now go through and read the data for each

	
	maskr = []
	var={}
	# interpolating to ROMS grid requires special care with phases!!
	#    this is subjected to parallelization via muitiprocessing.Process - do it!
	#    there is a sandbox in roms repo with a starting exercise
	for var0 in vars:
		var[var0]={}
		N=-1
		for con,ncon in zip(const[0]['name'][consindex],conIDX):
			N=N+1
			logger.info("Interpolating %s for %s" %(var0, con))
			
			if "amp" in var0:
				var[var0][con] = extrapolate_amp(X, Y, lon, lat, maskr, f.variables[var0][ncon,:,:]) 
			else:
				var[var0][con] = extrapolate_pha(X, Y, lon, lat, maskr, f.variables[var0][ncon,:,:])



	return var,tfreq,const[0]['name'][consindex]

def get_tide(cons,freq,amp,pha,t,lat0):
	nodes=len(amp[amp.keys()[0]])

	tide=np.ones((nodes,1))

	for node in range(0,nodes):
		tidecon=np.ones((len(cons),4))
		for i,con in enumerate(cons):
			tidecon[i,0]=amp[con][node]
			tidecon[i,1]=1e8
			tidecon[i,2]=pha[con][node]
			tidecon[i,3]=1e8

		tide[node,:]=t_predic(t,np.array(cons),np.array(freq),tidecon,lat=lat0)

	return tide
