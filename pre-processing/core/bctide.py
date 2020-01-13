from matplotlib.dates import date2num

from ttide.t_getconsts import t_getconsts
from ttide.t_vuf import t_vuf

class BCinputs(object):
    """
    Class that manages input file generation for a model simulation.
    """

    def __init__(self,obc,lat0,t0, logger=None):
        '''Docstring'''  

        if logger:
            self.logger = logger

        self.lat0 = lat0
        self.t0= t0
        self.t0dec=date2num(t0)
        self.obc=obc

    def Calculate(lat, t0, cons):
        const = t_getconsts(np.array([]))
        consindex = [np.where(const[0]['name']==con.ljust(4))[0][0] for con in cons]
        # V: astronomical phase, U: nodal phase modulation, F: nodal amplitude correction
        v,u,f = t_vuf('nodal', np.array([t0]), consindex, lat)
        tear = 360.*(v+u)
        tfreq = (2*np.pi)*const[0]['freq'][consindex]/3600.
        talpha = const[0]['name'][consindex]
        return talpha, tfreq, tear, f

    def _write_bctides(self,filename, talpha=[], tfreq=[], tear=[], tnf=[]):
        '''
           ----------------------------------- Elevation b.c. section ------
            iettype(j): (DEFAULT: 5)

                0: elevations are not specified for this boundary (in this case the velocity must be specified).

                1: time history of elevation on this boundary
                  > no input in bctides.in; time history of elevation is read in from elev.th (ASCII);

                2: forced by a constant elevation
                  > ethconst (constant elevation value for this segment)

                3: forced by tides
                  > for k=1, nbfr
                  >   alpha(k) !tidal constituent name
                  >   for i=1, nond(j) !loop over all open boundary nodes on this segment
                  >       emo((j,i,k) efa (j,i,k) !amplitude and phase for each node on this open boundary
                  >   end for i
                  > end for k

                4: space- and time-varying input
                  > no input in this file; time history of elevation is read in from elev2D.th(binary)

                5: combination of (3) and (4):
                   > time history of elevation is read in from elev2D.th, and then added to tidal B.C. specified below
                   > for k = 1, nbfr
                   >     alpha(k) !tidal constituent name
                   >     for i = 1, nond(j) !loop over all open boundary nodes on this segment
                   >         emo((j, i, k) efa(j, i, k) !amplitude and phase for each node on this open boundary
                   >     end for i
                   > end for k

            ----------------------------------- Velocity b.c. section -----------------------
            ifltype(j): (DEFAULT: 5)

                0: velocity not specified (no input needed)

                1: time history of discharge on this boundary
                  > no input in this file; time history of discharge is read in from flux.th (ASCII)

                2: this boundary is forced by a constant discharge
                  > vthconst (constant discharge; note that a negative number means inflow)

                3: vel. (not discharge!) is forced in frequency domain
                  > for k=1, nbfr
                  >    vmo(j,k) vfa(j,k) (uniform amplitude and phase along each boundary segment)
                  >    end for;

                4 or -4: 3D input
                        > time history of velocity (not discharge!) is read in from uv3D.th (binary)

                5: combination of (4) and tides
                  > time history of velocity (not discharge!) is read in from uv3D.th (binary) and then added to tidal velocity specified below
                  >    for k=1, nbfr
                  >       alpha(k) (tidal constituent name)
                  >       for i=1, nond(j) (loop over all open boundary nodes on this segment)
                  >           umo(j,i,k) ufa(j,i,k) vmo(j,i,k) vfa(j,i,k) !amplitude and phase for (u,v) at each node on this open boundary
                  >       end for i
                  >    end for k

                -1: Flanther type radiation b.c. (iettype must be 0 in this case)
                   > eta_mean (mean elevation below)
                   > for i=1,nond(j) (loop over all nodes)
                   >     eta_m0(i) !mean elev at each node
                   > end for i
                   > vn_mean (mean normal velocity)
                   > for i=1,nond(j)
                   >    qthcon(1:Nz,i,j) (mean normal velocity at the node; at all levels)
                   > end for i

            ----------------------------- Temperature b.c. section ---------------------------
            itetype: (DEFAULT: ?)

                0: temperature not specified (no input needed)

                1: time history of temperature on this boundary
                  > tobc !nudging factor (between 0 and 1 with 1 being strongest nudging) for inflow; time history of temperature will be read in from TEM_1.th (ASCII)

                2: this boundary is forced by a constant temperature
                  > tthconst (constant temperature on this segment)
                  > tobc (nudging factor (between 0 and 1) for inflow)

                3: initial temperature profile for inflow
                  > tobc (nudging factor (between 0 and 1) for inflow)

                4: 3D input
                  > tobc (nudging factor (between 0 and 1); time history of temperature is read in from TEM_3D.th (binary)

            ---------------------------- Salintiy B.C. section ------------------------------
            isatype: (DEFAULT: ?)

                the structure is similar to temperature

            ---------------------------- Tracers B.C. section -------------------------------
            If any tracer module is invoked, you also need the corresponding B.C. part for each tracer module

            isatype: (DEFAULT: 0)

                the structure is similar to temperature and salinity

            end for !j: open boundary segment
            *** Notes: the tidal amplitudes and phases can be generated using utility scripts shown on the web.'''


        self.bctides = open(self.bctides_fname, 'w')
        # 48-character start time info string (only used for visualization with xmvis6)
        self.bctides.write(self.nest.timerun.t0.strftime('%d/%m/%Y %H:%M:%S\n'))
        # ntip tip_dp > # of constituents used in earth tidal potential; cut-off depth for applying tidal potential (i.e., it is not calculated when depth < tip_dp).
        self.bctides.write("0 1000\n")
        # nbfr > # of tidal boundary forcing frequencies
        self.bctides.write("%d\n"%len(talpha))
        for k in range(0, len(talpha)):
            # tidal constituent name
        	self.bctides.write("%s\n"%(talpha[k])) 
            # angular frequency (rad/s), nodal factor, earth equilibrium argument (deg)
        	self.bctides.write("%.6f %.6f %.6f\n" % (tfreq[k], tnf[k], np.mod(tear[k],360)))
        #number of open boundary segments
        self.bctides.write('%s\n'%self.nest.obc.nopen)


        # looping through all open boundaries
        for btype, bdict in self.nest.obc.btype.items():

            iettype = bdict['iettype']
            ifltype = bdict['ifltype']
            itetype = bdict['itetype']
            isatype = bdict['isatype']

            # open boundary types
            bctypes = '{nnode} {iettype} {ifltype} {itetype} {isatype}'.format(nnode=bdict['nnode'],
                                                                               iettype=iettype['type'],
                                                                               ifltype=ifltype['type'],
                                                                               itetype=itetype['type'],
                                                                               isatype=isatype['type'])

            # need to make the tracer bnd info. call consistent to other bnd types (i.e ocean, rivers)        
            if self.nest.tracers:
                bctypes = bctypes.__add__('{tracers}') 
                trctypes = []

                if btype in self.nest.tracers.obc.keys():
                    for tid, tdict in self.nest.tracers.obc[btype].items():
                        for tmodule in self.nest.tracers.obc['modules']: # to preserve module order in string
                            if tid == tmodule:
                                trctypes.append(' {itrtype}'.format(itrtype=tdict['itrtype']['type']))
                else:
                    for tmodule in self.nest.tracers.obc['modules']:
                        if tmodule in self.nest.tracers.obc['active_modules']:
                            trctypes.append(' {itrtype}'.format(itrtype=0))
                
                bctypes = bctypes.format(tracers=''.join(trctypes))
            
            self.bctides.write(bctypes)
            self.bctides.write("\n")

            #----------------------------------------------------- Elevation b.c section ---------

            #elevations are not specified (velocity must be specified).
            if iettype['type'] == 0:
                pass
                # Check if velocity is specified or raise exception.

            #time history of elevation (no input needed; elevation is read from elev.th)
            elif iettype['type'] == 1:
                pass
                # Check if file exists or raise exception.

            #forced by a constant elevation
            elif iettype['type'] == 2:
		    	self.bctides.write("{:.2f}\n".format(iettype['const']))

		    #space and time-varying input (no input needed; elevation is read from elev2D.th)
            elif iettype['type'] == 4:
		    	pass
                # Check if file exists or raise exception.

            #iettype=3: forced by tides
            #iettype=5: combination of 3 and 4 (elevation is read in from elev2D.th', and then added to tidal B.C.)
            elif iettype['type'] in (3, 5):
                et = self.get_tide(bdict['lons'], bdict['lats'], outfmt='constituents')['et']
                for k in range(len(talpha)):
                    # tidal constituent name
                    self.bctides.write("%s\n"%(talpha[k])) 
                    # loop over all open boundary nodes on this segment
                    for i in range(len(bdict['nodes'])):
                        # amplitude and phase for each node on this open boundary
                        self.bctides.write("{:.6f} {:.6f}\n".format(et.amp[k,i],np.mod(degrees(et.pha[k,i]),360)))

            #----------------------------------------------------- Velocity b.c section ----------

            #velocity not specified
            if ifltype['type'] == 0:
                pass

            #time history of discharge (read from flux.th)
            elif ifltype['type'] == 1:
                pass
                # Check if file exist or raise exception.

            #forced by a constant discharge (negative value means inflow)
            elif ifltype['type'] == 2:
		        self.bctides.write("{:.4f}\n".format(ifltype['const']))

            #velocity forced in frequency domain
            elif ifltype['type'] == 3:
		        print "not implemented yet"

            # 3D input (time history of velocity is read in from uv3D.th)
            elif ifltype['type'] in (4,-4):
                # Check if file exist or raise exception.
                if ifltype['type'] == -4:
		        	self.bctides.write("{:.2f} {:.2f}\n".format(ifltype['inflow'],ifltype['outflow']))

            # combination of 4 and tides
            elif ifltype['type'] in (5,-5):
                # Check if file exist or raise exception.
                
                if ifltype['type'] == -5:
		        	self.bctides.write("{:.2f} {:.2f}\n".format(ifltype['inflow'],ifltype['outflow']))

                ut = self.get_tide(bdict['lons'], bdict['lats'], outfmt='constituents')['ut']
                vt = self.get_tide(bdict['lons'], bdict['lats'], outfmt='constituents')['vt']
                for k in range(len(talpha)):
                    # tidal constituent name
                    self.bctides.write("%s\n"%(talpha[k])) 
                    # loop over all open boundary nodes on this segment
                    for i in range(len(bdict['nodes'])):
                        # amplitude and phase for (u,v) at each node on this open boundary
                        self.bctides.write("{:.6f} {:.6f} {:.6f} {:.6f}\n".format(ut.amp[k,i],np.mod(degrees(ut.pha[k,i]),360),
                                                                                  vt.amp[k,i],np.mod(degrees(vt.pha[k,i]),360)))

	        # radiation condition - Flather (iettype must be 0)
            elif ifltype['type'] == -1:
		    	self.bctides.write("{:.2f} {:.2f}\n".format(ifltype['eta_m0'],ifltype['qthcon']))

            #----------------------------------------------------- Temperature b.c section -------

            # temperature not specified
            if itetype['type'] == 0:
                pass

            # time history or initial profile for inflow or 3D input
            if itetype['type'] == 1 or itetype['type'] == 3 or itetype['type'] == 4:
		    	self.bctides.write("{:.2f}\n".format(itetype['tobc']))

            # forced by a constant temperature
            elif itetype['type'] == 2:
		    	# constant temperature on this segment
		    	self.bctides.write("{:.2f}\n".format(itetype['const']))
		    	# nudging factor (between 0 and 1) for inflow
		    	self.bctides.write("{:.2f}\n".format(itetype['tobc']))

            #----------------------------------------------------- Salinity b.c section ----------

            # salinity not specified
            if isatype['type'] == 0:
                pass 

            # time history or initial profile for inflow or 3D input
            elif isatype['type'] == 1 or isatype['type'] == 3 or isatype['type'] == 4:
		    	self.bctides.write("{:.2f}\n".format(isatype['tobc']))

            # forced by a constant salinity
            elif isatype['type'] == 2:
		    	# constant salinity on this segment
		    	self.bctides.write("{:.2f}\n".format(isatype['const']))
		    	# nudging factor (between 0 and 1) for inflow
		    	self.bctides.write("{:.2f}\n".format(isatype['tobc']))

            #----------------------------------------------------- Tracers b.c section -----------
            if self.nest.tracers:
  
                if btype in self.nest.tracers.obc.keys():
                    for tid, tdict in self.nest.tracers.obc[btype].items():
                        itrtype = tdict['itrtype']
                        for tmodule in self.nest.tracers.obc['modules']: # to preserve module order in string
                            if tid == tmodule:
                                tthconst = self.nest.tracers.obc['bctides'][tid][btype]['tthconst'] 
                                tobc = self.nest.tracers.obc['bctides'][tid][btype]['tobc'] 

                                # tracer not specified
                                if itrtype['type'] == 0:
                                    pass

                                # time history or initial profile for inflow or 3D input
                                elif itrtype['type'] == 1 or itrtype['type'] == 3 or itrtype['type'] == 4:
		                        	self.bctides.write("{tobc}\n".format(tobc=tobc))
                                
                                # forced by a constant tracer
                                elif itrtype['type'] == 2:
                                    # constant tracer on this segment
                                    self.bctides.write("{tthconst}\n".format(tthconst=tthconst))
                                    # nudging factor (between 0 and 1) for inflow
                                    self.bctides.write("{tobc}\n".format(tobc=tobc))

        self.bctides.flush()
        self.bctides.close()


    def make_bctides(self,filename):
        '''Docstring'''
        
        if self.logger:
            self.logger.info("Writing Bctides.in")


        talpha, tfreq, tear, tnf = Calculate(self.lat0, self.t0dec, self.cons)
        self._write_bctides(filename,talpha, tfreq, tear, tnf)

