""" 
Script to make a psf model for a set of HST frames in a given filter.  The
script will make raw psf models using TinyTim for each frame and simulate one
star in the middle (i.e. on the pixel position of the galaxy) of each frame.
PSF defocus will be read from the header for each frame.

The frames are then drizzled and the background subtracted. Finally the psf
model of the drizzled image is cut out from the subtracted image. 

Parameters:
imlist      -- ASCII file with list of flt images, one per line. 
loclist     -- ASCII file with a list of extension ID, and x/y pixel
		coordinates for the raw psf model. One row per frame.
tinybase    -- Filename of the tiny1 base parameter file. NB, filter
		and instrument needs to be correctly set.
simpos      -- ASCII file with simulation position (RA,DEC). World 
		coordinates in degrees. 
drizconfig  -- Filename of drizzlepac configuration file.
--outpsf    -- Root of final output PSF model (optional, default is 'PSF'
--psfscale X    -- Scale the psf model with the factor X before co-adding. 
		The model is also scaled with the exposure time of the frame.
		(optional, default is 3).
--order O       -- The order of the interpolation used when shifting psf model
		to sub-pixel accurate positions (optional, default is 1).
--defkey    -- FITS header keyword for defocus values (optional, default is
		DFOC_MOD).

All or part of the parameters can be supplied as a text file with
one parameter per row and with the file name preceded by @. Note
that the correct order of the positional arguments must be kept.

"""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from astropy.io import fits
import argparse
import pyfits as pyf
import numpy as np
from subprocess import call
from scipy import rand
from pyraf import iraf
from scipy.ndimage.interpolation import geometric_transform as geotr
import configparser
#from drizzlepac import skytopix
from astropy import wcs 

from zap import Zap as zap
from numpy import ceil, loadtxt, c_, zeros, modf, array
from numpy import round as rr
import sys
import fileinput
from drizzlepac import astrodrizzle
adriz=astrodrizzle.AstroDrizzle



filt = ['775', '782']

from utils import basic_params
import os
#-----------------------------------------------------------------------------
# IRAF/tinytim setup
#-----------------------------------------------------------------------------
hsel=iraf.hselect
#rd2xy=skytopix.rd2xy
imcent=iraf.imcentroid
# Dirty fix, specify directory where tinytim can be found

#tinydir='/home/jens/tools/tinytim-7.5/'

#tinydir='/home/sourabh/miniconda3/lib/python3.5/site-packages/tinytim-7.5/'  
#-----------------------------------------------------------------------------
# Functions/classes
#-----------------------------------------------------------------------------
def makeconf(conffile,defoc,loc,root,offset):
	if ((defoc/0.011-offset)<0.): 
		sign='m'
	else:
		sign='p'
	defocstr=sign+str(abs(int(100*np.round((defoc/0.011-offset),2))))
	locstr= str(int(loc[0]))+'_'+str(int(loc[1]))
	
	with open(root+'_foc'+defocstr+'_loc'+locstr+'.conf','w') as fout:
		print (conffile)
		for modline in fileinput.input(conffile):
			outline=modline
			if modline.rstrip('\n').endswith('chip')|modline.rstrip('\n').endswith('= Focus'):
				outline = str(defoc)+' # Z4 = Focus\n'
			if modline.rstrip('\n').endswith('Position 1'):
				outline = str(int(loc[0]))+' '+\
					str(int(loc[1]))+'  # Position 1\n'
			if modline.rstrip('\n').endswith('rootname'):
				outline = root+'_foc'+defocstr+'_loc'+locstr+\
					' # Output PSF file rootname\n'
			fout.write(outline)

	return root+'_foc'+defocstr+'_loc'+locstr+'.conf',root+'_foc'+defocstr+'_loc'+locstr

def runtiny(conf,params,root,offset):
	#rpsfnames=runtiny(tinybase,dllist,psfbase,0.)

	# params should contain [defocus, xloc, yloc]
	psfoutnames=[]
	if params.ndim==1:
		pars=params
		defoc=(pars[0]+offset)*0.011
		loc=pars[1:]
		# Create config file
		conffile,psfoutname=makeconf(conf,defoc,loc,root,offset)
		psfoutnames=[psfoutname]
		# Run tiny2 + tiny3
		call(['./tiny2',conffile])
		call(['./tiny3',conffile])
	else:
		for pars in params:
			defoc=(pars[0]+offset)*0.011
			loc=pars[1:]
			# Create config file
			conffile,psfoutname=makeconf(conf,defoc,loc,root,offset)
			print(conffile,psfoutname)
			psfoutnames.append(psfoutname)
			# Run tiny2 + tiny3
			call(['./tiny2',conffile])
			call(['./tiny3',conffile])
	return psfoutnames

def shift_func(out_coo,dx,dy):
	"""Function used for geometric transformation of images. This
	implements a simple x,y shift."""
	return(out_coo[0]-dx,out_coo[1]-dy)

def hstsim(ind, images,psfs,simwcs,exts,psfscale,order):
	outlist=[fname.replace('flc.fits','sim.fits') for fname in images]
	#outlist=[fname.replace('single.fits','sim.fits') for fname in outlist]

	# Go through the list of images and add sources
	for ii,fname in enumerate(images):
		print('*******************************************')
		print('Adding simulated source to '+fname)
		print('in extension '+str(exts[ii]))
		print('PSF file(s) used: '+psfs[ii])
		print('=> output file: '+outlist[ii])
		print('')
		hdu_psf = fits.open(psfs[ii])
		psf = hdu_psf[0].data
		extyes=exts[ii]

		hdu_gal= fits.open(fname)
		exptime=hdu_gal[0].header['EXPTIME']
		hdu_gal.writeto(outlist[ii], overwrite = True)
		#hdr=hdul[extyes].header
		# Read coordinate positions from wcs coord file
		# and use WCS conversion in fits file.
		hdu1 = fits.open(fname)
		if ind ==4:
			hdr1 = hdu1[5].header  #### for ULIRG 5 it is extension 5 others its 2
		else:
			hdr1 = hdu1[2].header
		w = wcs.WCS(hdr1)

		world =np.array([simwcs[ii][0], simwcs[ii][1]]  )
		print (world[0],world[1])
		simcoo = w.wcs_world2pix(world[0],world[1], 1)#165.5588688980612, 38.04309910403013, 1)

		simcoo = modf(simcoo)

		print(simcoo[0][0], simcoo[0][1], simcoo[1][0], simcoo[1][1])
		#simcoo=rd2xy(fname+'['+str(extyes)+']',ra=simwcs[0],
		#		dec=simwcs[1],verbose=True)
		#simcoo=list(int (simcoo[x1]) for x1 in range(2))#np.round(simcoo) 
		#print(simcoo)

		# Get cutout region from frame, subtract one
		# from all values to adjust to python indices.
		
		xlo=simcoo[1][0]-psf.shape[0]/2.-1
		xhi=simcoo[1][0]+psf.shape[0]/2.-1
		ylo=simcoo[1][1]-psf.shape[1]/2.-1
		yhi=simcoo[1][1]+psf.shape[1]/2.-1    

		boxcoo=ceil(c_[xlo,xhi,ylo,yhi])[0].astype(int)
		print (xlo, yhi, ylo, yhi)
		print (boxcoo)
		sim=zeros(hdu_gal[extyes].data.shape)
		# work on a copy of the frame
		print ("extensions are these ==========================>", extyes)
		sim=hdu_gal[extyes].data
		psf_shifted=zeros(psf.shape)
		print('Co-adding the PSF model to the input frame')
		shifts=(simcoo[0][1],simcoo[0][0])
		print (exptime, psfscale, shifts, order, extyes)
		psf_shifted=geotr(psf,shift_func, extra_arguments=shifts,order=order)*exptime*psfscale#*10**5
		#psf_shifted = psf
		sim1 = np.copy(sim)
		fits.writeto( "psfshift.fits", data = psf_shifted, overwrite = True)
		fits.writeto( "psf.fits", data = psf,  overwrite = True)
		fits.writeto( "before_adding.fits", data = sim,  overwrite = True)
		
		#print (psf_shifted.shape)
		#if boxcoo[2]<boxcoo[0]:
		a1 = boxcoo[2]
		a2 = boxcoo[3]
		a3 = boxcoo[0]
		a4 = boxcoo[1]
		sim1[boxcoo[2]:boxcoo[3], boxcoo[0]:boxcoo[1]] = sim[boxcoo[2]:boxcoo[3], boxcoo[0]:boxcoo[1]] + psf_shifted
		#sim1[a1:a2, a3:a4] = sim[a1:a2, a3:a4] + psf_shifted
		#for i1 in range (psf_shifted.shape[0]):
		#	for j1 in range(psf_shifted.shape[0]):
		#		sim1[a1+i1, a3+j1] = sim[a1+i1, a3+j1] +psf_shifted[i1,j1] 
		#else:
			#sim1[boxcoo[0]:boxcoo[1], boxcoo[2]:boxcoo[3]] = sim[boxcoo[0]:boxcoo[1], boxcoo[2]:boxcoo[3]] + psf_shifted

		print ("psfscale -->", psfscale)
		print ("exptime -->" , exptime)

		#fits.writeto("magic.fits", data = sim1-sim, overwrite = True)
		fits.writeto("magic.fits", data = sim1, overwrite = True)

		# Write to output fits file
		hdu_gal[extyes].data = sim1
		hdu_gal.writeto(outlist[ii], overwrite = True)
		#pyf.update(outlist[ii],sim,hdr,extyes)
		#elif (inst=='SBC')|(inst=='HRC'):
		#    # For ACS:SBC/HRC there is only one chip
		#    # which corresponds to extension 1.
		#    extyes=1
		#    hdul=pyf.open(fname)
		#    exptime=hdul[0].header['EXPTIME']
		#    hdul.writeto(outlist[ii])
		#    hdr=hdul[extyes].header
		#    # Read coordinate positions from wcs coord file
		#    # and use WCS conversion in fits file.
		#    simcoo=rd2xy(fname+'['+str(extyes)+']',ra=simwcs[0],
		#            dec=simwcs[1],verbose=False)
		#    simcoo=modf(simcoo) 
		#    # Get cutout region from frame, subtract one
		#    # from all values to adjust to python indices.
		#    xlo=simcoo[1][0]-psf.shape[0]/2.-1
		#    xhi=simcoo[1][0]+psf.shape[0]/2.-1
		#    ylo=simcoo[1][1]-psf.shape[1]/2.-1
		#    yhi=simcoo[1][1]+psf.shape[1]/2.-1          
		#    boxcoo=ceil(c_[xlo,xhi,ylo,yhi])[0].astype(int)
		#    print(boxcoo)
		#    sim=zeros(hdul[extyes].data.shape)
		#    # work on a copy of the frame
		#    sim[:,:]=hdul[extyes].data
		#    psf_shifted=zeros(psf.shape)
		#    print('Co-adding the PSF model to the input frame')
		#    shifts=(simcoo[0][1],simcoo[0][0])
		#    psf_shifted=geotr(psf,shift_func,
		#        extra_arguments=shifts,order=order)*\
		#        exptime*psfscale
		#    sim[boxcoo[2]:boxcoo[3],boxcoo[0]:boxcoo[1]]+=psf_shifted
		#    # Write to output fits file
		#    pyf.update(outlist[ii],sim,hdr,extyes)

		#elif (inst=='UVIS'):
		#    # For WFC3/UVIS the galaxy is always in chip2 (UVIS2)
		#    # which corresponds to extension 1.
		#    # For Tololo 1247, the galaxy is in chip1...
		#    extyes=1
		#    hdul=pyf.open(fname)
		#    exptime=hdul[0].header['EXPTIME']
		#    hdul.writeto(outlist[ii])
		#    hdr=hdul[extyes].header
		#    # Read coordinate positions from wcs coord file
		#    # and use WCS conversion in fits file.
		#    simcoo=rd2xy(fname+'['+str(extyes)+']',ra=simwcs[0],
		#            dec=simwcs[1],verbose=False)
		#    simcoo=modf(simcoo) 
		#    # Get cutout region from frame, subtract one
		#    # from all values to adjust to python indices.
		#    xlo=simcoo[1][0]-psf.shape[0]/2.-1
		#    xhi=simcoo[1][0]+psf.shape[0]/2.-1
		#    ylo=simcoo[1][1]-psf.shape[1]/2.-1
		#    yhi=simcoo[1][1]+psf.shape[1]/2.-1          
		#    boxcoo=ceil(c_[xlo,xhi,ylo,yhi])[0].astype(int)
		#    sim=zeros(hdul[extyes].data.shape)
		#    # work on a copy of the frame
		#    sim[:,:]=hdul[extyes].data
		#    psf_shifted=zeros(psf.shape)
		#    print('Co-adding the PSF model to the input frame')
		#    shifts=(simcoo[0][1],simcoo[0][0])
		#    psf_shifted=geotr(psf,shift_func,
		#        extra_arguments=shifts,order=order)*\
		#        exptime*psfscale
		#    sim[boxcoo[2]:boxcoo[3],boxcoo[0]:boxcoo[1]]+=psf_shifted
		#    # Write to output fits file
		#    pyf.update(outlist[ii],sim,hdr,extyes)

		hdu_gal.close()
	return outlist

def wcstopixfile(wcscoords,refima,imcoords):
	"""Function that converts coordinates from wcs pixel. Input and output
	are ascii files with one pair of coordinates per row. The image refima
	contains the reference system to use."""
	wcscoo=np.loadtxt(wcscoords)
	imcoo=np.zeros(wcscoo.shape)
	if len(wcscoo.shape)==1: 
		wcscoo=np.array([wcscoo])
		imcoo=np.array([imcoo])
	for ii,radec in enumerate(wcscoo):
		

		#xytmp=rd2xy(refima+'[0]',ra=radec[0],dec=radec[1],verbose=False)
		fits_file = fits.open(refima)
		w = wcs.WCS(fits_file[2].header)

		#world = np.array([simwcs[ii][0], simwcs[ii][1]]  )
		#print (world[0],world[1])
		print(radec)
		xytmp = w.wcs_world2pix(radec[0],radec[1], 1)#165.5588688980612, 38.04309910403013, 1)
		print(xytmp)

		imcoo[ii,:]=np.array([xytmp[0],xytmp[1]])
	print (imcoords)
	np.savetxt(imcoords,imcoo,fmt='%g')
	return imcoo

#-----------------------------------------------------------------------------
# Main
#-----------------------------------------------------------------------------
def ha_psf(im_drz, ind, imlist,loclist,tinybase,simpos,drzcfg,outpsf,psfscale,order,defkey):
	"""Main script""" 
	# Get inputs and perform some sanity checks
	print('*******************************************')
	print('Reading input files.')
	imin=open(imlist)
	images=imin.read().split()
	imin.close()
	simwcs=loadtxt(simpos)
	locin=loadtxt(loclist)
	
	if len(locin.shape)>1:
		n_l=locin.shape[0]
	else:
		n_l=1  

	# Abort if files mismatch
	if (n_l)!=len(images):
		print('Number of frames is not equal to number of defocus  \
			values given. Aborting!')
		sys.exit()
	# Check that all frames come from the same detector and filter
	print('Detector and filters for all frames:')
	notfirst=False
	dets=[]
	filts=[]
	for ima in images:
		t=hsel(images=ima+'[0]', fields="DETECTOR", expr="yes", Stdout=1)
		dets.append(t[0])
		t=hsel(images=ima+'[0]', fields="FILTER,FILTER1,FILTER2", expr="yes", Stdout=1)
		filts.append(t[0].replace("N/A", "").replace("\t", "").replace("CLEAR1L", "").replace("CLEAR2L", "").replace("CLEAR1S",""))
	if (len(set(dets))>1)|(len(set(filts))>1):
		print('Image list contains data from multiple configurations.\
			Aborting!')
		sys.exit()

	detector=dets[0]
	print (detector)
	filt1=filts[0]
	print (filt1)

	# STEP 1 First create the raw psf models
	print('********************************************************************************') 
	print('Creating raw psf models, one per frame ...')
	print('********************************************************************************') 
	
	# Start by reading defocus values from headers
	# Create tinytim input list, dllist
	if (n_l==1):
		xylocs = locin[1:]
		exts = [locin.astype(int)[0]]
		hdr = pyf.getheader(images[0])
		dllist=np.r_[hdr[defkey],xylocs]
		meandef = dllist[0]
	else:
		xylocs = locin[:,1:]
		exts = locin[:,0].astype(int) 
		defocii = np.zeros(n_l)
		for ii in range(n_l):
			hdr = pyf.getheader(images[ii])
			defocii[ii] = hdr[defkey]
		meandef = defocii.mean()
		dllist=np.c_[defocii,xylocs]
	if meandef>0:
		finalout=outpsf+'_defp'+str(abs(int(rr(meandef,2)*100)))+'.fits'
	else:
		finalout=outpsf+'_defm'+str(abs(int(rr(meandef,2)*100)))+'.fits'

	# Temporary file names for the raw psfs
	psfbase='psf_'+str(int(rand()*1e8))
	print ('PSF base is give by ', psfbase)
	# Offset is set to zero, make sure to correct input
	rpsfnames=runtiny(tinybase,dllist,psfbase,0.)
	rpsfnames=[name+'.fits' for name in rpsfnames]
	# Rename the raw models, remove the rest.
	call(['rename s/00.fits/.fits/ *00.fits'],shell=True)
	call(['rm -f '+psfbase+'*tt3 '+psfbase+'*00_psf.fits '+psfbase+'*conf'],shell=True)
	print('The following raw psf models have been created:')
	for name in rpsfnames: print(name)
	
	# STEP 2 Add raw psfs to the flt frames
	print('********************************************************************************') 
	print('Adding simulated stars to flt frames.')
	print('********************************************************************************')   
	simimas=hstsim(ind, images,rpsfnames,simwcs,exts,psfscale,order)
	print('The following simulated images have been created: ')
	if len(simimas)>1:
		for name in simimas: print(name)
	else:
		print(simimas)

	# STEP 3 Drizzle flt and sim frames
	# First drizzle the non-simulated data
	print('********************************************************************************') 
	print('Drizzling unsimulated data.')
	print('********************************************************************************') 
	# Setup drizzlepac
	tmpcfg = "tmp"+str(int(rand()*1e8))+".cfg"
	sbccfg  = "sbc_driz.cfg"
	uviscfg = "uvis_driz.cfg"
	config_dir = '/home/sourabh/ULIRG_package/config/'
	wfccfg  = config_dir+ "astrodrizzle_WFC_conf.cfg"
	hrccfg  = "hrc_driz.cfg"
	with open(drzcfg) as ff: oblines=ff.read().splitlines()
	if detector == "SBC" : 
		with open(tmpcfg, 'w') as fout:
			for modline in fileinput.input(sbccfg):
				outline=modline
				for obline in oblines:
					if modline.startswith(obline.split()[0]):
						outline=obline+'\n'
				fout.write(outline)
	elif detector == "UVIS":
		with open(tmpcfg, 'w') as fout:
			for modline in fileinput.input(uviscfg):
				outline=modline
				for obline in oblines:
					if modline.startswith(obline.split()[0]):
						outline=obline+'\n'
				fout.write(outline)
	elif detector == "WFC": 
		with open(tmpcfg, 'w') as fout:
			for modline in fileinput.input(wfccfg):
				outline=modline
				
				#print (outline , "chhhhhhhhhhhh")
				'''
				for obline in oblines:
					if modline.startswith(obline.split()[0]):
						outline=obline+'\n'
				'''
				fout.write(outline)
	elif detector == "HRC": 
		with open(tmpcfg, 'w') as fout:
			for modline in fileinput.input(hrccfg):
				outline=modline
				for obline in oblines:
					if modline.startswith(obline.split()[0]):
						outline=obline+'\n'
				fout.write(outline)
	else: 
		print ("!!!!error --"), detector, ("not one of [ SBC | UVIS | WFC | HRC ]")
		sys.exit()  
	print ("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW")
	 
	nosima= os.path.dirname(images[0])+'/'+filt1+'_nosim.fits'
	tmpout='tmp'+str(int(rand()*1e8))
	adriz(input=images,output=tmpout,runfile='astrodrizzle.log',final_wht_type = 'ERR',
		configobj=tmpcfg)
	print ("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW")
	if os.path.exists(tmpout+'_drc.fits'):
		call(['cp',tmpout+'_drc.fits',nosima])
	else:
		call(['cp',tmpout+'_drz.fits',nosima])


	#zap(tmpout+'_drc_sci.fits')
	#zap(tmpout+'_drc_wht.fits')
	#zap(tmpout+'_drc_ctx.fits')
	
	# Then drizzle the simulated images
	print('********************************************************************************') 
	print('Drizzling simulated data.')
	print('********************************************************************************') 
	sima=os.path.dirname(images[0])+'/'+ filt1+'_sim.fits'
	tmpout='tmp'+str(int(rand()*1e8))
	adriz(input=simimas,output=tmpout,runfile='astrodrizzle.log',final_wht_type = 'ERR',
		configobj=tmpcfg)
	if os.path.exists(tmpout+'_drc.fits'):
		call(['cp',tmpout+'_drc.fits',sima])
	else:
		call(['cp',tmpout+'_drz.fits',sima])


	#zap(tmpout+'_drc_sci.fits')
	#zap(tmpout+'_drc_wht.fits')
	#zap(tmpout+'_drc_ctx.fits')
	
	# STEP 4 Subtract model frame and cutout the psf model
	# determine size of the frame
	print('********************************************************************************') 
	print('Removing background and cutting out PSF model.')
	print('********************************************************************************') 
	xsize=int(hsel(images=sima+ '[1]',fields='NAXIS1',expr='yes',Stdout=1)[0])
	ysize=int(hsel(images=sima+'[1]',fields='NAXIS2',expr='yes',Stdout=1)[0])
	impos_approx=np.array([xsize/2,ysize/2])
	print ("position of image approximately ", impos_approx)
	imreg='['+str(int((impos_approx[0]-65)))+':'+\
		str(int((impos_approx[0]+65)))+','+\
		str(int((impos_approx[1]-65)))+':'+\
		str(int((impos_approx[1]+65)))+']'

	################################

	##### two ways to find positin f simulated psf
	##### and then find shift between them ###
	# First find the exact position of the simulated psf
	sima_sub=sima.replace('sim','sim_sub')
	hdu_sim = fits.open(sima)
	hdu_nosim = fits.open(nosima)
	dat1 = hdu_sim[1].data - hdu_nosim[1].data
	hdu_sim[1].data = dat1
	#dat1[np.isnan(dat1)] = 0.0
	#iraf.imarith(sima,'-',nosima, sima.replace('sim','sim_sub'))
	hdu_sim.writeto(sima.replace('sim', 'sim_sub'), overwrite = True)
	hdu_sim.close()
	
	simimfile='tmp'+str(int(rand()*1e8))+'_im.coords'
	

	imcoords = wcstopixfile(simpos,sima,simimfile)
	
	iraf.imcent.unlearn()

	imcent(sima_sub+'[1]', coords=simimfile, niterate = 20,Stdout=tmpout+'.dat')
	
	
	with open(tmpout+'.dat', 'r') as f: cooline=f.readlines()[1].split()[1:4]
	cooline=[float(cooline[ii]) for ii in [0,2]]
	
	shift=(imcoords[0][1]-cooline[1],imcoords[0][0]-cooline[0])
	
	print (shift, cooline, imcoords )
	size = 300
	
	x1 = int(imcoords[0][1]-size/2)
	x2 = int(imcoords[0][1]+size/2)
	y1 = int(imcoords[0][0]-size/2)
	y2 = int(imcoords[0][0]+size/2)
	'''
	x1 = int(im_drz[ind][1]-size/2)
	x2 = int(im_drz[ind][1]+size/2)
	y1 = int(im_drz[ind][0]-size/2)
	y2 = int(im_drz[ind][0]+size/2)
	'''
	psfstamp = dat1[x1:x2, y1:y2]
	#psf_s1=geotr(psfstamp,shift_func,extra_arguments=shift,order=order)
	#psf_s2 = psf_s1/psf_s1.sum()
	psf_s2=psfstamp/psfstamp.sum()

	#iraf.imcopy(sima_sub+imreg,finalout)

	
	hdu1 = fits.PrimaryHDU(data = psf_s2)
	header1 = hdu1.header
	header1["CREATOR"] = ("SSC", "Sept 17 2018")
	hdu1.writeto("%s.fits"%(outpsf), overwrite = True)
	'''
	psfstamp=pyf.getdata(finalout)
	hdr=pyf.getheader(finalout)
	# Shift the psf to the center

	pyf.update(finalout,psf_s,hdr, output_verify = "ignore")
	'''
#-----------------------------------------------------------------------------
# Script I/O
#-----------------------------------------------------------------------------
def main(config_file):
	#im_drz =[(1320.2164578776137, 4232.254568212932),(3133.44270371664,	5998.802594324208),(1077.9077757159448, 1782.7589540070876),(4144.009171783876, 5989.1062047339765),(5859.227208853262, 4226.265089816052)]
	im_drz =[(1320.2164578776137, 4232.254568212932),(3133.44270371664,	5998.802594324208),(1078.92343744, 1784.99832295),(4144.009171783876, 5989.1062047339765),(5859.227208853262, 4226.265089816052)]

	for i in range(5):
		if i ==2:
			section_gal = 'NULIRG%s' % (int(i + 1))
			params, params_gal = basic_params(config_file, 'basic', section_gal)
			drizconfig = params['wfc_config']
			psfscale = float(params['psfscale'])
			order = int(params['order'])
			defkey = params['defkey']
			gal_name = params_gal['name']
			primary_dir = params['data_dir'] + gal_name + '/'
			txt_dir = primary_dir + 'HA_PSF/TXT/' 
			call('./tiny1')
			for j in range(2):
				outpsf = params['outpsf']+'_%s_gal%s'%(filt[j], i+1)

				imlist = txt_dir + 'imlist%s_gal%s.txt' % (filt[j], i + 1)
				loclist = txt_dir + 'loclist%s_gal%s.txt' % (filt[j], i + 1)
				tinybase = params['config_dir'] + 'tinybase%s_gal%s.txt' % (filt[j], i + 1)
				simpos = txt_dir + 'simpos%s_gal%s.txt' % (filt[j], i + 1)
				ha_psf(im_drz, i, imlist,loclist,tinybase,simpos,drizconfig,
				outpsf, psfscale,order,defkey)



if __name__ == '__main__':
	main(config_file)
