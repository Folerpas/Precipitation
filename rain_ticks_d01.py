#!/usr/bin/env python
# -⁻- coding: UTF-8 -*-

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import os
import time
from mpl_toolkits.basemap import Basemap, addcyclic
from datetime import datetime,timedelta
from scipy.ndimage.filters import minimum_filter, maximum_filter

#--------------------------------------------------------
####################### FUNCTIONS  #######################
def georeferencio():
	"""compute latitudes, longitudes, x, y"""
	x_dim = len(ncfile.dimensions['west_east'])
	y_dim = len(ncfile.dimensions['south_north'])
	dx = float(ncfile.DX)
	dy = float(ncfile.DY)
	width_meters = dx * (x_dim - 1)
	height_meters = dy * (y_dim - 1)
	cen_lat = float(ncfile.CEN_LAT)
	cen_lon = float(ncfile.CEN_LON)
	truelat1 = float(ncfile.TRUELAT1)
	truelat2 = float(ncfile.TRUELAT2)
	standlon = float(ncfile.STAND_LON)
	truelat1 = float(ncfile.TRUELAT1)
	truelat2 = float(ncfile.TRUELAT2)
	mp = Basemap(resolution='i',projection='lcc',\
	    width=width_meters,height=height_meters,\
	    lat_0=cen_lat,lon_0=cen_lon,lat_1=truelat1,\
	    lat_2=truelat2)
	lons1 = ncfile.variables['XLONG'][0]
	lats = ncfile.variables['XLAT'][0]
	nlats = len(lats)
	nlons = len(lons1)
	x,y = mp(lons1,lats)
	return x,y, mp,nlats,nlons,lons1,lats;

def dailyrain(t):
	"""compute daily rain """
	rain = (ncfile.variables['RAINNC'][t,:,:]+ncfile.variables['RAINC'][t,:,:])-(ncfile.variables['RAINNC'][t-3,:,:]+ncfile.variables['RAINC'][t-3,:,:])
	return rain;

def pressure(t):
	''' Convert Surface Pressure to Sea Level Pressure'''	
	t2 =  ncfile.variables['T2'][t,:,:]
	u10 = ncfile.variables['U10'][t,:,:]
	v10 = ncfile.variables['V10'][t,:,:]
	hgt = ncfile.variables['HGT'][t,:,:]
	psfc = ncfile.variables['PSFC'][t,:,:]
	stemps = t2+6.5*hgt/1000.
	mslp = psfc*np.exp(9.81/(287.0*stemps)*hgt)*0.01 + (6.7 * hgt / 1000)
	return mslp;


def ticks(matriz,mode='wrap',window=10):
    '''find min and max values from an inputa matrix'''
    baixa = minimum_filter(matriz, mode=mode, size=window)
    alta = maximum_filter(matriz, mode=mode, size=window)
    return np.nonzero(matriz == baixa), np.nonzero(matriz == alta)

#--------------------------------------------------------
####################### CONSTANTS  #######################

yyyy = str(time.localtime()[0])
mm = (str(time.localtime()[1])).zfill(2)
dd = (str(time.localtime()[2])).zfill(2)
today=str(yyyy)+str(mm)+str(dd)
date=datetime(int(yyyy),int(mm),int(dd),00) #comezo no dia no que estamos
delta=timedelta(hours=3)


#----------------------------------------------------------------------------------------
####################### Download from Meteogalicia & open netcdf  #######################

thredds_meteogalicia='https://mandeo.meteogalicia.es/thredds/fileServer/wrf_2d_36km/fmrc/files/'+str(today)+'/'
netcdf='wrf_arw_det_history_d01_'+str(today)+'_0000.nc4'
PATH_DESCARGA='/Users/sabelasanfiz'

if os.path.isfile(PATH_DESCARGA+netcdf):
	ncfile=Dataset(PATH_DESCARGA+netcdf)
	print('archivo abierto '+str(netcdf))
else:
	os.system('wget '+thredds_meteogalicia+netcdf+' '+PATH_DESCARGA)
	print('descarga completada ')
	ncfile=Dataset(PATH_DESCARGA+netcdf)
	print('archivo abierto '+str(netcdf))


#--------------------------------------------
################## PLOT  ####################

for t in range(ncfile.variables['Times'].shape[0]): 
	x=georeferencio()[0]
	y=georeferencio()[1]
	mp=georeferencio()[2]
	nlats=georeferencio()[3]
	nlons=georeferencio()[4]
	lons1=georeferencio()[5]
	lats=georeferencio()[6]
	#define plot
	fig = plt.figure(figsize=(10,8),frameon=False)
	ax1 = fig.add_subplot(111)
	ax1.set_title('Presion (hPa), Rain 3hr (l/m^2) '+str(date),fontsize=14,bbox=dict(facecolor='white', alpha=0.65),x=0.5,y=1.05,weight = 'demibold')

	#plot pressure
	clevs = np.array(np.arange(940,1050,5))
	mslp = pressure(t)
	cs1 = mp.contour(x,y,mslp,clevs,linewidths=1.5,animated=True)

	# plot rain
	#clevs=[ 0.2, 0.5, 1, 2, 5,10, 15, 20, 25, 30, 35, 40, 45,50]
	#colors=['#A0FBE8','#00FEFE','#00C8FE','#0096FE','#0064FE','#0032FE','#3200FE','#6400FE','#9600FE','#C800FE','#FA00FE','#C800C8','#960096']
	#rain = dailyrain(t)
	#cs2 = mp.contourf(x,y,rain,clevs,colors=colors,extend='both')
	#cs2.cmap.set_over('#640064')
	#cs2.cmap.set_under((1,1,1,0))
	#cbar=plt.colorbar(cs2,orientation='vertical',extend='both',shrink=0.5)

	#basemap properties
	mp.fillcontinents(color='#F6E4D8',lake_color='#C9D3E8',zorder=0)
	mp.drawcoastlines(linewidth=0.2)
	mp.drawcoastlines(linewidth=0.2)
	mp.drawparallels(np.arange(10.,70.,5.),labels=[True,False,False,False],color = '0.25', linewidth = 0.1)
	mp.drawmeridians(np.arange(-40.,20.,5.),labels=[False,False,False,True],color = '0.25', linewidth = 0.1)

	#plot high's & low's. 
	mslp=pressure(t)
	local_min, local_max = ticks(mslp, mode='wrap', window=50)
	mslp, lons = addcyclic(mslp, lons1[0,:])

	xlows = x[local_min]; ylows = y[local_min]; lowvals = mslp[local_min];
	xhighs = x[local_max];yhighs = y[local_max]; highvals = mslp[local_max];
	 
	# plot lows as blue L's, with min pressure value underneath.
	xyplotted = []
	# don't plot if there is already a L or H within dmin meters.
	yoffset = 0.022*(mp.ymax-mp.ymin)
	dmin = yoffset
	for x,y,p in zip(xlows, ylows, lowvals):
	    if x < mp.xmax and x > mp.xmin and y < mp.ymax and y > mp.ymin:
		 dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
		 if not dist or min(dist) > dmin:
		     plt.text(x,y,'B',fontsize=16,fontweight='bold',ha='center',va='center',color='b')
		     plt.text(x,y-yoffset,repr(int(p)),fontsize=9,ha='center',va='top',color='b',bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))
		     xyplotted.append((x,y))
	# plot highs as red H's, with max pressure value underneath.
	xyplotted = []
	for x,y,p in zip(xhighs, yhighs, highvals):
	    if x < mp.xmax and x > mp.xmin and y < mp.ymax and y > mp.ymin:
		 dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
		 if not dist or min(dist) > dmin:
		     plt.text(x,y,'A',fontsize=16,fontweight='bold',ha='center',va='center',color='r')
		     plt.text(x,y-yoffset,repr(int(p)),fontsize=9,ha='center',va='top',color='r',bbox = dict(boxstyle="square",ec='None',fc=(1,1,1,0.5)))
		     xyplotted.append((x,y))

	plt.savefig('./'+str(t)+'_slp.png',dpi=300,bbox_inches='tight')
	print ('dia: '+str(date)+' feito')
	plt.close()
	date=date+delta

del t;

print('imaxes feitas :)')
#creo un gif
os.system('convert -delay 100 -loop 0 *.png slp'+str(today)+'.gif')
#eliminamos todo, agas o gif
os.system('rm '+str(archivo))
#os.system('rm *png')

