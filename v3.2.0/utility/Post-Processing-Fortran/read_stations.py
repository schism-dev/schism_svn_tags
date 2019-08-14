import numpy
import pylab # optional, used only for plotting at the end of the script
import time,os
import extract_mod as rog

# the python module extract_mod includes a fortran module which contains
# the global variables of the extract_mod module
# this shallow copy is simply a coding convenience
a = rog.extract_mod
# specify a bad value to be used in the extracted data
a.fill_in = numpy.nan

def pass_xyz(a,x,y,z=None):
   '''simple function to assign x,y, and z locations for data extraction in the fortran global variables. If no z values are provided, set global variable specifying profile/transect output'''
   a.nxy = len(x)
   a.x00 = x
   a.y00 = y
   a.z00 = z
   if z is None:
     a.itransect = 1
   else: 
     a.itransect = 0 # 0 for depth profile at specified point, 1 for timeseries at single depth

class outputrecord_fixedstation():
    '''creates an object for storing model output outside of the fortran global variables'''
    def __init__(self,a):
       self.time = numpy.array([])
       if a.itransect == 1 and a.i23d == 3:
         # allocate empty arrays shaped appropriately to append profile data
         self.depths = numpy.ones([a.nvrt,a.nxy,0])
         self.values1 = numpy.ones([a.nvrt,a.nxy,0])
         if a.ivs == 2:
            self.values2 = numpy.ones([a.nvrt,a.nxy,0])
       else:
         # allocate empty arrays shaped appropriately to append timeseries at x,y,z points
         self.depths = numpy.ones([a.nxy,0])
         self.values1 = numpy.ones([a.nxy,0])
         if a.ivs == 2:
            self.values2 = numpy.ones([a.nxy,0])

# examples of manual assignment of locations, comment out either (A) or (B)
# (A)
# manually assign locations at which to extract single depth stations, the 3rd array is the depths
# pass_xyz(a,numpy.array([350000, 344000]),numpy.array([289600,292000]),numpy.array([-2,-10]))

# (B)
# manually assign locations at which to extract a profile, with no third argument, all depths will be extracted
pass_xyz(a,numpy.array([350000, 344000]),numpy.array([289600,292000]))

# specify coordinate system of binary model output
a.ics =1 # 1 for cartesian coordinates, 2 for latlon in binary output

# specify coordinate system of requested z coordinates 
a.ifs = 0 # 0 (input z are z coordinates) or 1 (input z are relative to free surface;

# there is also a function in the module to read station locations from a file

# rog.read_station('station.ts.sample')

# read_station handles location files in station, transect, xyt, and xyzt formats
# the station location file format specifies ics, ifs, itransect 
# as well as the locations and times to extract
 
# set base file name
ftype = 'salt.63'
basedir = 'outputs'
# specify range of output files
ran = [1,2]

first = True
# loop through time period of interest
for i in ran:#range(1,ndays):
  # construct file name of binary output
  fname = os.path.join(basedir,'%d_%s' % (i,ftype))
  # extract header information from first binary output file
  if first:
    rog.readheader(fname)
    # find parent element for each location at which data will be extracted
    rog.find_parents()
    # preallocate storage for extracted data in the fortran module
    a.outtime = numpy.ones(a.nrec)*numpy.nan
    if a.itransect == 1 and a.i23d == 3:
       a.varout = numpy.ones([3,a.nvrt,a.nxy,a.nrec],order = 'F')*numpy.nan
    else:
       # allocate a degenerate dimension for depth to be consistent with the profile reader
       a.varout = numpy.ones([3,1,a.nxy,a.nrec],order = 'F')*numpy.nan
    # convert start time specified in header of output file to unix time
    base_time = time.mktime(time.strptime(''.join(a.start_time).strip(),'%m/%d/%Y %H:%M:%S %Z'))
    # create a storage structure for the extracted data, this one is designed for fixed stations
    output = outputrecord_fixedstation(a)
    first = False

  # read data from file
  rog.readdata(fname)
#       At this point, for time series input (itransect=0), 
#       varout(1:ivs,1:1,1:nxy,1:nrec) is the output with
#       times given by outtime(1:nrec) (in sec), where nrec is # of time steps within 
#       each stack, nxy is the # of points in bp_file, ivs=1 indicates
#       scalar (1) output; ivs=2 indicates vector outputs (e.g. 1=u; 2=v). 
#
#       For transect input (itransect=1; for 3D variables only), 
#       varout(1:ivs,1:nvrt,1:nxy,1:nrec) 
#       is the final output with times given by outtime(1:nrec) (in sec). 
#       The vertical structure (i.e. z coordinates) is given by 
#       varout(3,1:nvrt,1:nxy,1:nrec), where nvrt is the total # of vertical levels.
#
# make deep copies of the output data to prevent data from being overwritten 
# during the next read in the loop
  output.time = numpy.append(output.time,a.outtime[:a.nrec]+ base_time)
# single level, so explicitly extract single level to reduce singleton dimension
  if a.itransect == 1 and a.i23d == 3:
     output.depths = numpy.append(output.depths,a.varout[2,:,:,:],2) 
     output.values1 = numpy.append(output.values1,a.varout[0,:,:,:],2)
     if a.ivs == 2:
       output.values2 = numpy.append(output.values2,a.varout[1,:,:,:],2)
  else:
     # output array for single depth stations has nan depth data, since depth was an argument
     output.depths = a.z00
     output.values1 = numpy.append(output.values1,a.varout[0,0,:,:],1)
     if a.ivs == 2:
       output.values2 = numpy.append(output.values2,a.varout[1,0,:,:],1)
  
for i in range(0,a.nxy):
  pylab.figure()
  if a.itransect == 1 and a.i23d == 3:
    pylab.contourf(numpy.tile(output.time,(a.nvrt,1)),output.depths[:,i,:],output.values1[:,i,:])
    pylab.savefig('test_%s.png' % i, format = 'png')
  else:
    pylab.plot(output.time,output.values1[i,:])
    pylab.savefig('test_ts_%s.png' % i, format = 'png')
