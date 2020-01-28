from bctide import develop_bnd
from trimesh import OPEN_BOUNDARY
import numpy as np
from matplotlib.dates import num2date,date2num
from tidal_tools import extract_HC,get_tide
from res_tools import get_file_interpolator,vertical_interpolation
from filetype import create_ncTH






class OpenBoundaries(object):
    def __init__(self,obc,hgrid,vgrid,t0,t1,z0=0.001,logger=None):
        '''Docstring'''  

        if logger:
            self.logger = logger

        self.obc=obc
        self.hgrid=hgrid
        self.vgrid=vgrid
        self.t0=t0
        self.t1=t1
        self.bnd_nodes=[bnd.nodes for bnd in self.hgrid.mesh.boundaries if bnd.btype == OPEN_BOUNDARY]

        self.tidal=False
        self.residual=False
        self.i23d=2
        self.ivs=1
        self.z0=z0
        self.lat0=np.mean(self.hgrid.latitude)

        bnd=obc.get('bnd',None)
        if bnd: 
            self.llat,self.llon,self.zz=self.get_all_open_bnd(bnd)
            self.zz=np.array(self.zz)
        else:
            self.llat=self.hgrid.latitude
            self.llon=self.hgrid.longitude
            self.zz=self.hgrid.depth



    def get_all_open_bnd(self,bnd):
        if  type(bnd[0])==str:
            bnd=[int(x) for x in bnd[0].split()]

        Lat_ocean=[]
        Lon_ocean=[]
        Zlayer_ocean=[]
        for nn in bnd:
            ocean_boundary = self.bnd_nodes[nn-1] 
            Lat_oce=[]
            Lon_oce=[]
            Zlayer_oce=[]
            for node_i in ocean_boundary:
                Lat_ocean.append(self.hgrid.latitude[node_i])
                Lon_ocean.append(self.hgrid.longitude[node_i])
                Zlayer_ocean.append(self.vgrid.sigma_to_zlayer(node_i,self.hgrid.h[node_i],0.,0.1))

        return Lat_ocean,Lon_ocean,Zlayer_ocean        


    def add_res(self,res):
        self.f_out,self.LONR,self.LATR=get_file_interpolator(res['filename'],res['vars'],self.llon,self.llat)
        self.residual=True

    def add_tide(self,tidal):
    	self.HC,self.tfreq,self.const=extract_HC(tidal['filename'],tidal['vars'],self.llon,self.llat, logger=self.logger)
    	self.tidal=True



    def create_Dthnc(self,fileout,TimeSeries):      
        # create file
        if self.i23d==3:
            Nlev=self.zz.shape[1]
        else:
            Nlev=1

        time_Series,nc=create_ncTH(fileout,len(self.llon),Nlev,self.ivs,np.round((TimeSeries-TimeSeries[0])*24*3600))



        for n in range(0,len(TimeSeries)):
            total=np.zeros(shape=(self.ivs,len(self.llon),Nlev))

            # get tide
            if self.tidal:
                var=self.HC.keys()
                var=list(set([x.replace('_amp','').replace('_pha','') for x in var]))

                for i,v in enumerate(var):
                    # horizontal interpolation
                    tmp=get_tide(self.const,self.tfreq,self.HC[v+'_amp'],self.HC[v+'_pha'],np.array(TimeSeries[n]),self.lat0)

                    # vertical interpolation
                    total[i,:,:]=vertical_interpolation(tmp,self.zz,z0=self.z0)





    def create_th(self):
        pass

    def make_boundary(self,filename,dt=3600):
        if self.logger:
            self.logger.info("  Writing %s" %filename)
        TimeSeries=np.arange(date2num(self.t0),date2num(self.t1)+1,dt/(24.*3600.))
       
        if filename.endswith('.th.nc'):
            self.create_Dthnc(filename,TimeSeries)
        elif filename.endswith('.th'):
            self.create_th()







