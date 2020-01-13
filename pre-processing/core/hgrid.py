import sys,os
sys.path.append(os.path.join('/home/remy/Software/schism-public/','BayDeltaSCHISM','pyschism')) 
from base_io import BaseIO
from schism_mesh import SchismMeshGr3Reader
from trimesh import OPEN_BOUNDARY  # OPEN_BOUNDARY
from pyproj import Proj, transform
from schism_setup import SchismSetup

import argparse
# -------------------------------------------- BASE MODEL CLASSES -----------------------------
WGS84 = 4326 # Spatial Reference System
#-------------------------------------------- GRID METHODS -----------------------------
def transform_proj(x, y, in_espg=2193, out_espg=4326):
    ''' Performs cartographic transformations (converts from
        longitude,latitude to native map projection x,y coordinates and
        vice versa) using proj (http://trac.osgeo.org/proj/).
    '''
    in_proj = Proj(init='epsg:%s'%in_espg)
    out_proj = Proj(init='epsg:%s'%out_espg)
    x2, y2 = transform(in_proj, out_proj, x, y)
    return x2, y2

class HorizontalGrid(BaseIO):
    """A class that manages I/O of HGRID files
    - Transformation of grid files
    - Make available all grid data & metadata
    """

    def __init__(self,path_gr3,format='gr3',epsg=2193, logger=None):
        """ Constructor"""
        super(HorizontalGrid, self).__init__()

        if logger:
            self.logger = logger

        self.path_gr3=path_gr3
            
        # ----------------------------------------------------------- grid properties -----------       
        self.hgrid = self.load()
        self.mesh = self.hgrid
        self.nopen = self.mesh.n_boundaries(OPEN_BOUNDARY) # number of open boundaries segments    
        self.nnode = [bnd.n_nodes() for bnd in self.mesh.boundaries if bnd.btype == OPEN_BOUNDARY]
        self.h = self.mesh.nodes[:, 2]
        self.epsg=epsg

        if epsg == WGS84:
            self.longitude = self.mesh.nodes[:, 0]
            self.latitude = self.mesh.nodes[:, 1]
        else:
            self.path_ll=self.path_gr3[:-3]+'ll'
            easting = self.mesh.nodes[:, 0]
            northing = self.mesh.nodes[:, 1]
            self.longitude, self.latitude = transform_proj(easting, northing, in_espg=epsg, out_espg=WGS84)
            if not os.path.isfile(self.path_ll):
                sc=SchismSetup()
                sc.input.mesh=self.hgrid
                self.logger.info('\tCreating: %s' % (self.path_gr3[:-3]+'ll'))
                sc.write_hgrid_ll(self.path_gr3[:-3]+'ll',input_epsg=epsg, output_epsg=WGS84)


 

 
    def __repr__(self):
        rep = "\n%s\n" % self.__class__
        rep += "\ninfo: %s\n\t" % self.info.__dict__
        rep += "\nGrid types:\n \t path_ll: '%s' \n\t path_gr3: '%s'\n" % (self.path_ll, self.path_gr3)
        return rep

    def get_properties(self):
        '''Docstring'''
        pass

    def load(self):
        '''Docstring'''
        return SchismMeshGr3Reader().read(fpath_mesh=self.path_gr3,read_boundary=True)

    def close(self):
        '''Docstring'''
        self.hgrid.close()

    def link_to_root(self,rootdir):
        if not os.path.isfile(os.path.join(rootdir,'hgrid.gr3') ):
            os.symlink(self.path_gr3, os.path.join(rootdir,'hgrid.gr3') )

        if self.epsg != WGS84 and not os.path.isfile(os.path.join(rootdir,'hgrid.ll') ):    
            os.symlink(self.path_ll, os.path.join(rootdir,'hgrid.ll') )

def run_as_script():
    parser = argparse.ArgumentParser(description="SCHISM grid reader")
    parser.add_argument('-g', '--hgrid', type=str,
                        help='grid file')
    parser.add_argument('-f', '--format', type=str,default='gr3',
                        help='grid format')
    parser.add_argument('-e', '--epsg', type=int, default=2193,
                        help="grid coordinate reference")

    args = parser.parse_args()
    
    
    gr=HorizontalGrid(args.hgrid,format=args.format,epsg=args.epsg,logger=None)

if __name__ == '__main__':
    run_as_script()