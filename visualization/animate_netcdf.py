#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : animate_netcdf.py
# Author            : Tom Durrant <t.durrant@metocean.co.nz>
# Date              : 28.04.2016
# Last Modified Date: 28.08.2018
# Last Modified By  : Tom Durrant <t.durrant@metocean.co.nz>
# encoding: utf-8


from pymo.core.lib import UDSQuery
import urllib,time
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import numpy.ma as ma
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import logging
import datetime
import os
import pandas as pd
import argparse
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from glob import glob

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class Plot(object):

    def __init__(self, dset, var,nframes=None, proj=None, startdate=None,
                 enddate=None, domain=None):
        self.bmap=None
        self.var=var
        self._set_cblabel()
        self.nframes = nframes
        self.proj = proj or 'PlateCarree'
        if isinstance(dset,str):
            self.dset = xr.open_dataset(dset)
        elif isinstance(dset,list):
            self.dset = xr.open_mfdataset(dset)
        else:
            self.dset = dset
        if startdate:
            sd = datetime.datetime.strptime(startdate, "%Y%m%d_%Hz")
        else:
            sd = None
        if enddate:
            ed = datetime.datetime.strptime(enddate, "%Y%m%d_%Hz")
        else:
            ed = None
        self.dset = self.dset.sel(time=slice(sd, ed))
        if domain:
            lat1, lon1, lat2, lon2 = map(float, domain.split(','))
            self.dset = self.dset.sel(latitude=slice(lat1, lat2),
                                      longitude=slice(lon1, lon2))
        if self.var=='wndsp':
            try:
                self.dset[self.var]=np.sqrt(self.dset['ugrd10m']**2+self.dset['vgrd10m']**2)
            except Exception:
                self.dset[self.var]=np.sqrt(self.dset['uwnd']**2+self.dset['vwnd']**2)
        if self.var=='current':
            self.dset[self.var]=np.sqrt(self.dset['ucur']**2+self.dset['vcur']**2)
        if self.var=='tp':
            if not 'tp' in self.dset.variables:
                self.dset[self.var]=1/self.dset['fp']

    def _set_cblabel(self):
        if self.var=='hs':
            self.cblabel="$H_{s} [m]$"
        if self.var=='ugrdsfc':
            self.cblabel="$U_{10} [ms^{-1}]$"
        if self.var=='vgrdsfc':
            self.cblabel="$V_{10} [ms^{-1}]$"
        if self.var=='tp':
            self.cblabel="$T_{p} [sec]$"
        if self.var=='wndsp':
            self.cblabel="$U_{10} [ms^{-1}]$"


    def animate(self,**kwargs):
        fig=plt.figure(figsize=(12,6))
        ims=[]
        ims.append(plotMap(self.dset[self.var][0], **kwargs))
        kwargs.update({'add_colorbar': False})
        kwargs.update({'add_labels': False})
        kwargs.update({'decorate': False})
        vmin,vmax = ims[0][0].get_clim()
        kwargs.update({'vmin': vmin})
        kwargs.update({'vmax': vmax})
        nframes = self.nframes or self.dset[self.var].shape[0]
        for ii in np.arange(0,nframes):
            logger.info('Plotting %s'%str(ii))
            ims.append(plotMap(self.dset[self.var][ii], **kwargs))
        ani = animation.ArtistAnimation(fig, ims, interval=100, repeat_delay=3000, blit=True)
        return ani

    def plot_index_full(self,index,**kwargs):
        self.plot_index(index,**kwargs)
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "3%", pad="3%")
        plt.colorbar(self.im, cax=cax)

def plotMap(da, ax=None, ptype='pcolormesh',cblabel=None,decorate=True,
            title=False, proj='PlateCarree', **kwargs):
    if not ax:
        clon = (da.longitude.values.min() + da.longitude.values.max())/2.
        clon = (( clon + 180. ) % 360. ) - 180.
        clat = (da.latitude.values.min() + da.latitude.values.max())/2.
        if proj == 'PlateCarree':
            projection = ccrs.PlateCarree(central_longitude=clon)
        elif proj == 'Orthographic':
            projection = ccrs.Orthographic(central_longitude=clon,central_latitude=clat)
        else:
            raise Exception("Projection %s not supported" % proj)
        ax = plt.axes(projection=projection)
    kwargs.update({'ax': ax})
    if ptype=='pcolormesh':
        mpl = da.plot.pcolormesh(transform=ccrs.PlateCarree(), **kwargs)
    elif ptype=='contour':
        mpl = da.plot.contour(transform=ccrs.PlateCarree(),zorder=2, **kwargs)
        mpl = da.plot.contourf(transform=ccrs.PlateCarree(),zorder=1, **kwargs)
    else:
        raise Exception("Plot type %s not recognised" % ptype)
    plt.title('')
    if decorate:
        coast = NaturalEarthFeature(category='physical', scale='50m',
                                    facecolor='gray',
                                    edgecolor='black',name='coastline')
        ax.add_feature(coast)
        plt.axes(ax)
        if proj != 'Orthographic':
            gl = ax.gridlines(draw_labels=True)
            gl.xlabels_top = gl.ylabels_right = False
            gl.yformatter = LATITUDE_FORMATTER
            gl.xformatter = LONGITUDE_FORMATTER
    ax.set_extent([da.longitude.values.min(), da.longitude.values.max(), da.latitude.values.min(), da.latitude.values.max()])

    ts = pd.to_datetime(da.time.values)
    ttl = plt.text(0.5, 1.01, "Valid Time: %s" % ts.strftime('%Y%m%d %Hz'),
            horizontalalignment='center', verticalalignment='bottom',
            transform=ax.transAxes)
    return mpl, ttl


def test_anim():

    tmp=Plot(glob('/data/ww3/tinyapp/latest/bob*nc')[0],'hs',nframes=4)
    ani=tmp.animate()
    #ani.save('test.gif',writer='imagemagick',fps=5,bitrate=-1)
    plt.show()

def main():
    parser = argparse.ArgumentParser(
                        epilog="Example: \npython wwtimesteps test.nc")
    parser.add_argument('ncfile', type=str,
                    help='netcdf file')
    parser.add_argument('-v', '--variable', type=str, default='hs',
                    help='variable')
    parser.add_argument('-o', '--output', type=str, default='out.gif',
                    help='output')
    parser.add_argument('-p','--proj', type=str, default='PlateCarree',
                    help='Map projection')
    parser.add_argument('-f','--fps', type=int, default=5,
                    help='Frames per second')
    parser.add_argument('-r','--range', type=float, nargs='+', default=[0,10],
                    help='Range')
    parser.add_argument('-c','--cmap', type=str, default='jet',
                    help='cmap')
    parser.add_argument('-nf','--nframes', type=int, default=None,
                    help='cmap')
    parser.add_argument('-s','--startdate', type=str, default=None,
                    help='start date ,%Y%m%d_%Hz')
    parser.add_argument('-e','--enddate', type=str, default=None,
                    help='end date ,%Y%m%d_%Hz')
    parser.add_argument('-d','--domain', type=str, default=None,
                    help='Domain lat1,lon1,lat2,lon2')
    args = parser.parse_args()
    tmp=Plot(args.ncfile, args.variable, nframes=args.nframes,
            startdate=args.startdate, enddate=args.enddate, domain=args.domain)
    ani=tmp.animate(vmin=args.range[0],vmax=args.range[1],cmap=plt.get_cmap(args.cmap),
                    proj=args.proj)
    ani.save(args.output,writer='imagemagick',fps=args.fps,bitrate=-1)

if __name__ == '__main__':
    # test_anim()
    main()


