"""
Demo of using functions in ezTestNc.py.
To run, do 'python ezTestNcDemo.py'.

20111021 rsz Created.
20111027 rsz Add limits tests.
"""
import sys,os
from ezTestNc import *
import numpy as np
from datetime import datetime, timedelta
from netCDF4 import num2date, date2num
from ezTest import *

#######################################################################
class ezTestNcDemo(object):
  """ A test suite with single setup/teardown and many tests. """
  
  def __init__(self):
    self.name = 'ezTestNc Example Suite'
    
    self.tests = [
      [self.testFileExists, 'Mock data files exist.'],
      [self.testAllGlobalAttsEqual, 'Global attributes are equal for 2 files.'],
      [self.testAttEquals, 'Attribute values (global and variable).'],
      [self.testDimEquals, 'Dim sizes.'],
      [self.testDimsEqual, 'All dims equal for 2 files.'],
      [self.testFilesEqual, 'All data and metadata equal for 2 files.'],
      [self.testFormatEquals, 'File formats.'],
      [self.testFormatsEqual, 'File formats match for 2 files.'],
      [self.testHasAtt, 'Has specific attributes.'],
      [self.testHasDim, 'Has specific dims.'],
      [self.testHasSameVars, 'Two files has same list of vars.'],
      [self.testHasShape, 'Var has specific shape.'],
      [self.testHasVar, 'File has specific variable.'],
      [self.testIsRec, 'Test record dimensions.'],
      [self.testMinValueEquals, 'Minimum var values.'],
      [self.testMaxValueEquals, 'Maximum var values.'],
      [self.testShapesEqual, 'Verify var shapes.'],
      [self.testValueEquals, 'Check var values.'],
      [self.testValuesEqual, 'Compare var slices between 2 files.'],
      [self.testVarsEqualMeta, 'Compare var metadata for 2 files.'],
      [self.testVarsEqualData, 'Compare var data for 2 files.'],
      [self.testVarsEqual, 'Compare vars for 2 files.'],
      [self.testVarTypeEquals, 'Check specific var type.'],
      [self.testVarTypesEqual, 'Compare var types for 2 files.'],
      [self.testAllVarsEqualMeta, 'Compare metadata for all vars for 2 files.'],
      [self.testAllVarsEqualData, 'Compare data for all vars for 2 files.'],
      [self.testAllVarsEqual, 'Compare all vars for 2 files.'],
      [self.testValueLimits, 'Check if any or all values are min,max,nan,pinf,ninf,-0,denorm,finite.'],
    ]
  
  #--------------------------------------------------------------------  
  def setUp(self):
    """ Create mock netcdf data files. """
    self.filenames = ezTmpFilenames(suffix='.nc', n=5)
    self.CreateFile1()
    self.CreateFile2()
    self.CreateFile3()
    self.CreateFile4()
    self.CopyFile1()
  
  #--------------------------------------------------------------------  
  def tearDown(self):
    """ Delete mock netcdf data files. """
    #ezDeleteFiles(self.filenames)

  #--------------------------------------------------------------------  
  def CreateFile1(self): 
    nc = Dataset(self.filenames[0], 'w', format='NETCDF3_CLASSIC')
    nc.author = 'ezTestNcDemo'
    nc.place = '1 Palmer Square'
    
    nt,nz,ny,nx = 5,50,91,361
    dt = nc.createDimension('time', None)
    dz = nc.createDimension('z', nz)
    dy = nc.createDimension('y', ny)
    dx = nc.createDimension('x', nx)
    
    vt = nc.createVariable('time', 'f8', ('time',))
    vz = nc.createVariable('level', 'f4', ('z',))
    vy = nc.createVariable('latitude', 'i1', ('y',))
    vx = nc.createVariable('longitude', 'i2', ('x',))    
    temp = nc.createVariable('temp', 'f4', ('time','z','y','x',))
    
    temp.units = 'Kelvin'
    vt.units = 'hours since 0001-02-03 12:34:56.7'
    vt.calendar = 'gregorian'
    vy.units = 'degrees north'
    vx.units = 'degrees east'
    
    vx[:] = np.linspace(-180.,180,nx)
    vy[:] = np.linspace(0,ny-1,ny)
    vz[:] = np.linspace(1,nz,nz)**2
    dates = [datetime(2001,3,1)+n*timedelta(hours=12) for n in range(nt)]
    vt[:] = date2num(dates, units=vt.units, calendar=vt.calendar)
    temp[:,:,:,:] = np.reshape(np.linspace(1, nt*nz*ny*nx, nt*nz*ny*nx), [nt,nz,ny,nx])
    
    nc.close()
    
  #--------------------------------------------------------------------  
  def CreateFile2(self):
    nc = Dataset(self.filenames[1], 'w', format='NETCDF3_64BIT')
    nc.author = 'ezTestNcDemo'
    nc.place = '1 Palmer Square'
    
    nt,nz,ny,nx = 5,50,91,361
    dt = nc.createDimension('time', None)
    dz = nc.createDimension('z', nz)
    dy = nc.createDimension('y', ny)
    dx = nc.createDimension('x', nx)
    
    vt = nc.createVariable('time', 'f8', ('time',))
    vz = nc.createVariable('level', 'f4', ('z',))
    vy = nc.createVariable('latitude', 'i1', ('y',))
    vx = nc.createVariable('longitude', 'i2', ('x',))    
    temp = nc.createVariable('temp', 'f4', ('time','z','y','x',))

    vt.units = 'hours since 0001-02-03 12:34:56.7'
    vt.calendar = 'gregorian'

    vx[:] = np.linspace(-180.,180,nx)
    vy[:] = np.linspace(0,ny-1,ny)
    vz[:] = np.linspace(1,nz,nz)**2
    dates = [datetime(2001,3,1)+n*timedelta(hours=12) for n in range(nt)]
    vt[:] = date2num(dates, units=vt.units, calendar=vt.calendar)
    # Reverse of temp data in first file.
    temp[:,:,:,:] = np.reshape(np.linspace(nt*nz*ny*nx, 1, nt*nz*ny*nx), [nt,nz,ny,nx])

    nc.close()

  #--------------------------------------------------------------------  
  def CreateFile3(self): 
    nc = Dataset(self.filenames[2], 'w', format='NETCDF4')
    nc.method = 'CreateFile3'
    
    nt,n1,n2,n3 = 5,10,100,200
    dt = nc.createDimension('t', None)
    d1 = nc.createDimension('dim1', n1)
    d2 = nc.createDimension('dim2', n2)
    d3 = nc.createDimension('dim3', n3)

    var1 = nc.createVariable('var1', 'i1', ('t','dim1','dim2','dim3',))
    var1[...] = np.reshape(np.ones((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    
    nc.close()

  #--------------------------------------------------------------------  
  def CopyFile1(self): 
    cmd = 'cp -f %s %s' % (self.filenames[0], self.filenames[3])
    ezShell(cmd)

  #--------------------------------------------------------------------  
  def CreateFile4(self): 
    """ Create file with limits values per datatype. """
    nc = Dataset(self.filenames[4], 'w', format='NETCDF3_CLASSIC')
    nc.method = 'CreateFile4'
    
    nt,n1,n2,n3 = 5,10,20,30
    dt = nc.createDimension('t', None)
    d1 = nc.createDimension('d1', n1)
    d2 = nc.createDimension('d2', n2)
    d3 = nc.createDimension('d3', n3)
    
    v = nc.createVariable('int8anymin', 'int8', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    limits = np.iinfo(np.int8)
    v[1,2,3,4] = limits.min

    v = nc.createVariable('int8allmin', 'int8', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = limits.min

    v = nc.createVariable('int8anymax', 'int8', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = limits.max
    
    v = nc.createVariable('int8allmax', 'int8', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = limits.max

    # - - - - - - - - - - - - - - - - - - - - - - - -
    v = nc.createVariable('int16anymin', 'int16', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    limits = np.iinfo(np.int16)
    v[1,2,3,4] = limits.min

    v = nc.createVariable('int16allmin', 'int16', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = limits.min

    v = nc.createVariable('int16anymax', 'int16', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = limits.max
    
    v = nc.createVariable('int16allmax', 'int16', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = limits.max

    # - - - - - - - - - - - - - - - - - - - - - - - -
    v = nc.createVariable('int32anymin', 'int32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    limits = np.iinfo(np.int32)
    v[1,2,3,4] = limits.min

    v = nc.createVariable('int32allmin', 'int32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = limits.min

    v = nc.createVariable('int32anymax', 'int32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = limits.max
    
    v = nc.createVariable('int32allmax', 'int32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = limits.max

    # - - - - - - - - - - - - - - - - - - - - - - - -
    v = nc.createVariable('float32anymin', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    limits = np.finfo(np.float32)
    v[1,2,3,4] = limits.min

    v = nc.createVariable('float32allmin', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = limits.min

    v = nc.createVariable('float32anymax', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = limits.max
    
    v = nc.createVariable('float32allmax', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = limits.max

    # - - - - - - - - - - - - - - - - - - - - - - - -
    v = nc.createVariable('float32anynan', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = np.nan

    v = nc.createVariable('float32allnan', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = np.nan

    v = nc.createVariable('float32anypinf', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = np.inf

    v = nc.createVariable('float32allpinf', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = np.inf

    v = nc.createVariable('float32anyninf', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = -np.inf

    v = nc.createVariable('float32allninf', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = -np.inf

    nzero = np.arctan2(-0.0, 0.0)
    v = nc.createVariable('float32anynzero', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = nzero
    
    v = nc.createVariable('float32allnzero', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = nzero
    
    denorm = np.finfo(np.float32).tiny * np.finfo(np.float32).eps
    v = nc.createVariable('float32anydenorm', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = denorm

    v = nc.createVariable('float32alldenorm', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = denorm
 
    v = nc.createVariable('float32anyfinite', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = np.inf

    v = nc.createVariable('float32allfinite', 'float32', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.random.random((1,nt*n1*n2*n3)), [nt,n1,n2,n3])

    # - - - - - - - - - - - - - - - - - - - - - - - -
    v = nc.createVariable('float64anymin', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    limits = np.finfo(np.float64)
    v[1,2,3,4] = limits.min

    v = nc.createVariable('float64allmin', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = limits.min

    v = nc.createVariable('float64anymax', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = limits.max
    
    v = nc.createVariable('float64allmax', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = limits.max

    # - - - - - - - - - - - - - - - - - - - - - - - -
    v = nc.createVariable('float64anynan', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = np.nan

    v = nc.createVariable('float64allnan', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = np.nan

    v = nc.createVariable('float64anypinf', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = np.inf

    v = nc.createVariable('float64allpinf', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = np.inf

    v = nc.createVariable('float64anyninf', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = -np.inf

    v = nc.createVariable('float64allninf', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = -np.inf

    v = nc.createVariable('float64anynzero', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = nzero
    
    v = nc.createVariable('float64allnzero', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = nzero
    
    denorm = np.finfo(np.float64).tiny * np.finfo(np.float64).eps
    v = nc.createVariable('float64anydenorm', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = denorm

    v = nc.createVariable('float64alldenorm', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[...] = denorm
 
    v = nc.createVariable('float64anyfinite', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.zeros((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    v[1,2,3,4] = np.inf

    v = nc.createVariable('float64allfinite', 'float64', ('t','d1','d2','d3',))
    v[...] = np.reshape(np.random.random((1,nt*n1*n2*n3)), [nt,n1,n2,n3])
    
    nc.close()

  #--------------------------------------------------------------------
  def testFileExists(self):
    for f in self.filenames:
      if not os.path.exists(f):
        print 'ERROR: File "%s" doesn\'t exist.' % (f)
        return False
        
    return True
  
  #--------------------------------------------------------------------
  def testAllGlobalAttsEqual(self): 
    msg = ezMsg()
    res = ncAllGlobalAttsEqual(self.filenames[0], self.filenames[1], msg=msg)
    if not res: print msg.msg
    return res
    
  #--------------------------------------------------------------------
  def testAttEquals(self): 
    msg = ezMsg()
    res = ncAttEquals(self.filenames[0], 'author', var=None, value='ezTestNcDemo', msg=msg)
    if not res: print msg.msg
    return res

  #--------------------------------------------------------------------
  def testDimEquals(self):
    msg = ezMsg()
    res = ncDimEquals(self.filenames[0], 'x', 361, msg=msg)
    if not res: print msg.msg
    return res  

  #--------------------------------------------------------------------
  def testDimsEqual(self): 
    msg = ezMsg()
    res = ncDimsEqual(self.filenames[0], self.filenames[3], msg=msg)
    if not res: print msg.msg
    return res  
  
  #--------------------------------------------------------------------
  def testFilesEqual(self): 
    msg = ezMsg()
    res = ncFilesEqual(self.filenames[0], self.filenames[3], msg=msg)
    if not res: print msg.msg
    return res  

  #--------------------------------------------------------------------
  def testFormatEquals(self): 
    msg = ezMsg()
    res1 = ncFormatEquals(self.filenames[0], '32', msg=msg)
    if not res1: print msg.msg
    
    res2 = ncFormatEquals(self.filenames[1], '64', msg=msg)
    if not res2: print msg.msg

    res3 = ncFormatEquals(self.filenames[2], 'hdf', msg=msg)
    if not res3: print msg.msg
    
    return res1 and res2 and res3

  #--------------------------------------------------------------------
  def testFormatsEqual(self): 
    msg = ezMsg()
    res = ncFormatsEqual(self.filenames[0], self.filenames[3], msg=msg)
    if not res: print msg.msg
    return res
 
  #--------------------------------------------------------------------
  def testHasAtt(self):
    msg = ezMsg()
    res1 = ncHasAtt(self.filenames[0], 'author', msg=msg)
    if not res1: print msg.msg
    
    res2 = ncHasAtt(self.filenames[0], 'units', var='time', msg=msg)
    if not res2: print msg.msg

    res3 = ncHasAtt(self.filenames[0], 'offset', var='time', msg=msg)

    return res1 and res2 and not res3

  #--------------------------------------------------------------------
  def testHasDim(self):
    msg = ezMsg()
    res1 = ncHasDim(self.filenames[0], 'x', msg=msg)
    if not res1: print msg.msg
    
    res2 = ncHasDim(self.filenames[0], 'w', msg=msg)

    return res1 and not res2

  #--------------------------------------------------------------------
  def testHasSameVars(self):
    msg = ezMsg()
    res1 = ncHasSameVars(self.filenames[0], self.filenames[1], msg=msg)
    if not res1: print msg.msg

    res2 = ncHasSameVars(self.filenames[0], self.filenames[2], msg=msg)

    return res1 and not res2
    
  #--------------------------------------------------------------------
  def testHasShape(self):
    msg = ezMsg()
    res1 = ncHasShape(self.filenames[0], 'temp', shape=[5,50,91,361], msg=msg)
    if not res1: print msg.msg

    res2 = ncHasShape(self.filenames[0], 'temp', shape=[15,150,191,61], msg=msg)

    return res1 and not res2

  #--------------------------------------------------------------------
  def testHasVar(self):
    msg = ezMsg()
    res1 = ncHasVar(self.filenames[0], 'temp', msg=msg)
    if not res1: print msg.msg

    res2 = ncHasVar(self.filenames[0], 'tmp', msg=msg)

    return res1 and not res2

  #--------------------------------------------------------------------
  def testIsRec(self):
    msg = ezMsg()
    res1 = ncIsRec(self.filenames[0], 'time', msg=msg)
    if not res1: print msg.msg

    res2 = ncIsRec(self.filenames[0], 'temp', msg=msg)

    return res1 and not res2

  #--------------------------------------------------------------------
  def testMinValueEquals(self):
    msg = ezMsg()
    res1 = ncMinValueEquals(self.filenames[0], 'temp', '...', 1, msg=msg)
    if not res1: print msg.msg

    res2 = ncMinValueEquals(self.filenames[0], 'temp', ':,:,:,:', 1, msg=msg)
    if not res2: print msg.msg

    res3 = ncMinValueEquals(self.filenames[0], 'temp', '...', -1, msg=msg)
    
    return res1 and res2 and not res3

  #--------------------------------------------------------------------
  def testMaxValueEquals(self):
    msg = ezMsg()
    mx = 5*50*91*361
    res1 = ncMaxValueEquals(self.filenames[0], 'temp', '...', mx, msg=msg)
    if not res1: print msg.msg

    res2 = ncMaxValueEquals(self.filenames[0], 'temp', ':,:,:,:', mx, msg=msg)
    if not res2: print msg.msg

    res3 = ncMaxValueEquals(self.filenames[0], 'temp', '...', mx-1, msg=msg)
    
    return res1 and res2 and not res3

  #--------------------------------------------------------------------
  def testShapesEqual(self):
    msg = ezMsg()
    res1 = ncShapesEqual(self.filenames[0], 'temp', self.filenames[1], 'temp', msg=msg)
    if not res1: print msg.msg

    res2 = ncShapesEqual(self.filenames[0], 'temp', self.filenames[1], 'time', msg=msg)

    return res1 and not res2

  #--------------------------------------------------------------------
  def testValueEquals(self):
    msg = ezMsg()
    res1 = ncValueEquals(self.filenames[0], 'temp', [0,0,0,0], 1, msg=msg)
    if not res1: print msg.msg

    res2 = ncValueEquals(self.filenames[0], 'temp',  [-1,-1,-1,-1], 5*50*91*361, msg=msg)
    if not res2: print msg.msg
    
    res3 = ncValueEquals(self.filenames[0], 'temp',  [0,0,0,0], 0, msg=msg)
    
    return res1 and res2 and not res3

  #--------------------------------------------------------------------
  def testValuesEqual(self):
    msg = ezMsg()
    res1 = ncValuesEqual(self.filenames[0], 'temp', self.filenames[3], 'temp', '...', msg=msg)
    if not res1: print msg.msg

    res2 = ncValuesEqual(self.filenames[0], 'temp', self.filenames[1], 'temp', '...', msg=msg)

    return res1 and not res2  
  
  #--------------------------------------------------------------------
  def testVarsEqualMeta(self):
    msg = ezMsg()
    res1 = ncVarsEqual(self.filenames[0], 'temp', self.filenames[3], 'temp', atts=True, data=False, msg=msg)
    if not res1: print msg.msg

    res2 = ncVarsEqual(self.filenames[0], 'temp', self.filenames[1], 'temp', atts=True, data=False, msg=msg)

    return res1 and not res2  
  
  #--------------------------------------------------------------------
  def testVarsEqualData(self):
    msg = ezMsg()
    res1 = ncVarsEqual(self.filenames[0], 'temp', self.filenames[3], 'temp', atts=False, data=True, msg=msg)
    if not res1: print msg.msg

    res2 = ncVarsEqual(self.filenames[0], 'temp', self.filenames[1], 'temp', atts=False, data=True, msg=msg)

    return res1 and not res2  

  #--------------------------------------------------------------------
  def testVarsEqual(self):
    msg = ezMsg()
    res1 = ncVarsEqual(self.filenames[0], 'temp', self.filenames[3], 'temp', atts=True, data=True, msg=msg)
    if not res1: print msg.msg

    res2 = ncVarsEqual(self.filenames[0], 'temp', self.filenames[1], 'temp', atts=True, data=True, msg=msg)

    return res1 and not res2  
  
  #--------------------------------------------------------------------
  def testVarTypeEquals(self):
    msg = ezMsg()
    res1 = ncVarTypeEquals(self.filenames[0], 'temp', 'float32', msg=msg)
    if not res1: print msg.msg
    
    res2 = ncVarTypeEquals(self.filenames[0], 'temp', 'float64', msg=msg)

    return res1 and not res2  

  #--------------------------------------------------------------------
  def testVarTypesEqual(self):
    msg = ezMsg()
    res1 = ncVarTypesEqual(self.filenames[0], 'temp', self.filenames[3], 'temp', msg=msg)
    if not res1: print msg.msg

    res2 = ncVarTypesEqual(self.filenames[0], 'temp', self.filenames[1], 'time', msg=msg)

    return res1 and not res2  
 
  #--------------------------------------------------------------------
  def testAllVarsEqualMeta(self):
    msg = ezMsg()
    res1 = ncAllVarsEqual(self.filenames[0], self.filenames[3], atts=True, data=False, msg=msg)
    if not res1: print msg.msg

    res2 = ncAllVarsEqual(self.filenames[0], self.filenames[1], atts=True, data=False, msg=msg)

    return res1 and not res2  
  
  #--------------------------------------------------------------------
  def testAllVarsEqualData(self):
    msg = ezMsg()
    res1 = ncAllVarsEqual(self.filenames[0], self.filenames[3], atts=0, data=1, msg=msg)
    if not res1: print msg.msg

    res2 = ncAllVarsEqual(self.filenames[0], self.filenames[1], atts=0, data=1, msg=msg)

    return res1 and not res2  

  #--------------------------------------------------------------------
  def testAllVarsEqual(self):
    msg = ezMsg()
    res1 = ncAllVarsEqual(self.filenames[0], self.filenames[3], atts=1, data=1, msg=msg)
    if not res1: print msg.msg

    res2 = ncAllVarsEqual(self.filenames[0], self.filenames[1], atts=1, data=1, msg=msg)

    return res1 and not res2  
 
  #--------------------------------------------------------------------
  def testValueLimits(self):
    msg = ezMsg()
    res = True
    fn = self.filenames[4]
    
    for any in ['any', 'all']:
      for size in ['8','16','32']:
        for value in ['min', 'max']:
          var = 'int%s%s%s' % (size, any, value)
          res1 = ncValuesAre(fn, var, '...', any, value, msg=msg)
          res = res and res1
          if not res1 and msg.msg: print msg.msg
   
      for size in ['32','64']:
        for value in ['min', 'max', 'nan', 'pinf', 'ninf', 'denorm']:
          var = 'float%s%s%s' % (size, any, value)
          res1 = ncValuesAre(fn, var, '...', any, value, msg=msg)
          res = res and res1
          if not res1 and msg.msg: print msg.msg

        var = 'float%s%snzero' % (size, any)
        res1 = ncValuesAre(fn, var, '...', any, '-0', msg=msg)
        res = res and res1
        if not res1 and msg.msg: print msg.msg
          
        var = 'float%s%sfinite' % (size, any)
        res1 = ncValuesAre(fn, var, '...', any, 'finite', msg=msg)
        res = res and res1
        if not res1 and msg.msg: print msg.msg
        
    return res
  
#######################################################################
if __name__=='__main__':
  ezTestRunner([ezTestNcDemo()], opt=ezTestOpt()).run()
