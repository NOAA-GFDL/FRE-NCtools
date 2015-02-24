"""
Utility functions for testing Netcdf programs.
To use, just copy the file to your project path and do 'from ezTestNc import *'.
For demo, do 'python ezTestNcDemo.py'.

20111021 rsz Created.
20111026 rsz Added nan,ninf,pinf,-0,denorm,finite value compare function.
"""
from netCDF4 import Dataset
from netCDF4 import numpy as np

##################################################################
def ncAllGlobalAttsEqual(filename1, filename2, msg=None):
  res = True
  nc1 = Dataset(filename1)
  nc2 = Dataset(filename2)
  
  atts1 = nc1.ncattrs()
  atts2 = nc1.ncattrs()
  
  n1 = len(atts1)
  n2 = len(atts2)
  
  if n1 <> n2:
    if msg:
      msg.msg = 'Number of global attributes %d <> %d in files "%s" and "%s".' % (n1, n2, filename1, filename2)
    res = False

  if atts1 <> atts2:
    if msg:
      msg.msg = 'Global attribute names %s <> %s in files "%s" and "%s".' % (str(atts1), str(atts2), filename1, filename2)
    res = False

  if res:
    for i in xrange(n1):
      v1 = atts1[i]
      v2 = atts2[i]
      if v1 <> v2:
        if msg:
          msg.msg = 'Global attribute %s %s <> %s in files "%s" and "%s".' % (name, str(v1), str(v2), filename1, filename2)
        res = False
        break
    
  nc1.close()
  nc2.close()

  return res
  
##################################################################
def ncAllValuesEqual(filename1, var1, filename2, var2, msg=None):
  """ Test all values between 2 files and 2 variables, capable of handling large variables. """
  if not ncHasVar(filename1, var1, msg=msg): return False
  if not ncHasVar(filename2, var2, msg=msg): return False
  
  if not ncShapesEqual(filename1, var1, filename2, var2, msg) : return False
    
  res = True
  nc1 = Dataset(filename1)
  nc2 = Dataset(filename2)
  dims1 = nc1.dimensions[var1]
  dims2 = nc1.dimensions[var2]
  if ncIsRec(filename1, dims1[0]) or ncIsRec(filename2, dims2[0]):
    nt = len(nc1.dimensions[var1])
    for t in xrange(nt):
      res = (nc1.variables[var1][t,...] == nc2.variables[var2][t,...]).all()
      if not res: break
  else:
    res = (nc1.variables[var1][...] == nc2.variables[var2][...]).all()
  
  nc1.close()
  nc2.close()

  if (not res) and msg:
    msg.msg = '%s <> %s in files "%s" and "%s".' % (var1, var2, filename1, filename2)

  return res

##################################################################
def ncAttEquals(filename, att, var=None, value=None, msg=None):
  """ Test if var or global att equals a value. """
  if not ncHasAtt(filename, att, var, msg=msg): return False
  
  nc = Dataset(filename)
  if var == None:
    res = (getattr(nc,att) == value)
  else:
    res = (getattr(nc.variables[var],att) == value)
    
  # Handle list of values.
  if hasattr(res,'all'): 
    res = res.all() 
    
  nc.close()

  if (not res) and msg:
    if var:
      msg.msg = '%s.%s <> %s in file "%s".' % (var, att, str(value), filename)
    else:
      msg.msg = 'Global attribute %s <> %s in file "%s".' % (att, str(value), filename)

  return res
  
##################################################################
def ncDimEquals(filename, dim, size, msg=None):
  """ Test if size equals a value. """
  if not ncHasDim(filename, dim, msg=msg): return False
  
  nc = Dataset(filename)
  dsize = len(nc.dimensions[dim])
  res = (dsize == size)
  nc.close()

  if (not res) and msg:
    msg.msg = 'Dimension %s size %d <> %d in file "%s".' % (dim, dsize, size, filename)
    
  return res

##################################################################
def ncDimsEqual(filename1, filename2, msg=None):
  """ Compares all dims between two files. """
  nc1 = Dataset(filename1)
  nc2 = Dataset(filename2)
  
  n1 = len(nc1.dimensions)
  n2 = len(nc2.dimensions)
  
  res = True
  if n1 <> n2:
    if msg:
      msg.msg = 'Number of dimensions %d <> %d in files "%s" and "%s".' % (n1, n2, filename1, filename2)
    res = False
  
  if res:
    dims1 = nc1.dimensions
    dims2 = nc2.dimensions
    names1,names2 = [],[]
    for d in dims1: names1.append(d)
    for d in dims2: names2.append(d)
    
    if names1 <> names2:
      if msg:
        msg.msg = 'Dimensions %s <> %s in files "%s" and "%s".' % (str(dims1), str(dims2), filename1, filename2)
      res = False
    
    for d in dims1:
      len1 = len(dims1[d])
      len2 = len(dims2[d])
      if len1 <> len2:
        if msg:
          msg.msg = 'Dimension %s length %d <> %d in files "%s" and "%s".' % (d, len1, len2, filename1, filename2)
        res = False
        break
    
  nc1.close()
  nc2.close()
  
  return res
  
##################################################################
def ncFilesEqual(filename1, filename2, formats=True, meta=True, data=True, msg=None):
  """ Full comparison of files. """
  if formats:
    if not ncFormatsEqual(filename1, filename2, msg=msg): return False
    
  if meta:
    if not ncDimsEqual(filename1, filename2, msg=msg): return False
    if not ncAllGlobalAttsEqual(filename1, filename2, msg=msg): return False
  
  return ncAllVarsEqual(filename1, filename2, atts=meta, data=data, msg=msg)
  
##################################################################
def ncFormatEquals(filename, fmt, msg=None):
  """ fmt = 32|64|hdf """
  c = fmt[0]
  if c == '3': format = 'NETCDF3_CLASSIC'
  elif c == '6': format = 'NETCDF3_64BIT'
  else: format = 'NETCDF4'
  
  nc = Dataset(filename)
  theformat = nc.file_format
  res = theformat.count(format)
  nc.close()
  
  if (not res) and msg:
    msg.msg = 'Format %s <> %s in file "%s".' % (format, theformat, filename)
  
  return res

##################################################################
def ncFormatsEqual(filename1, filename2, msg=None):
  nc1 = Dataset(filename1)
  nc2 = Dataset(filename2)
  format1 = nc1.file_format
  format2 = nc2.file_format
  nc1.close()
  nc2.close()
  
  res = (format1 == format2)
  if (not res) and msg:
    msg.msg = 'Format %s <> %s for files "%s" and "%s".' % (format1, format2, filename1, filename2)

  return res

##################################################################
def ncGetAtt(filename, att, var=None):
  """ Returns var or global att value or None. """
  if not ncHasAtt(filename, att, var): return None
  
  nc = Dataset(filename)
  if var == None:
    res = getattr(nc,att)
  else:
    res = getattr(nc.variables[var],att)

  return res
  
##################################################################
def ncHasAtt(filename, name, var=None, msg=None):
  nc = Dataset(filename)
  if var == None:
    res = (nc.ncattrs().count(name) <> 0)
  else:
    res = nc.variables.has_key(var) and (nc.variables[var].ncattrs().count(name) <> 0)
    
  nc.close()
  
  if (not res) and msg:
    if var:
      msg.msg = 'Variable "%s" doesn\'t have attribute "%s" in file "%s".' % (var, name, filename)
    else:
      msg.msg = 'Global attribute "%s" doesn\'t exist in file "%s".' % (name, filename)
    
  return res

##################################################################
def ncHasDim(filename, name, msg=None):
  nc = Dataset(filename)
  res = nc.dimensions.has_key(name)
  nc.close()

  if (not res) and msg:
    msg.msg = 'Dimension %s doesn\'t exist in file "%s".' % (name, filename)
  
  return res

##################################################################
def ncHasSameVars(filename1, filename2, msg=None):
  """ Checks if both files have same list of variables. """
  res = True
  nc1 = Dataset(filename1)
  nc2 = Dataset(filename2)
  vars1 = nc1.variables
  vars2 = nc2.variables
  nv1 = len(vars1)
  nv2 = len(vars2)
  
  if nv1 <> nv2:
    if msg:
      msg.msg = 'Number of variables %d <> %d in files "%s" and "%s".' % (nv1, nv2, filename1, filename2)
      
    res = False
    
  if res:
    names1, names2 = [],[]
    for name in vars1: names1.append(name)
    for name in vars2: names2.append(name)
    if names1 <> names2:
      if msg:
        msg.msg = 'Variable names %s <> %s in files "%s" and "%s".' % (str(names1), str(names2), filename1, filename2)
      res = False
      
  nc1.close()
  nc2.close()

  return res
  
##################################################################
def ncHasShape(filename, var, shape=(), msg=None):
  if not ncHasVar(filename, var): return False
  
  if type(shape) <> type(()):
    shape = eval('('+str(shape)[1:-1]+')')
    
  nc = Dataset(filename)
  sh = nc.variables[var].shape
  res = (sh == shape)
  nc.close()

  if (not res) and msg:
    msg.msg = 'Variable %s shape %s <> %s in file "%s".' % (var, str(sh), str(shape),filename)

  return res

##################################################################
def ncHasVar(filename, name, msg=None):
  nc = Dataset(filename)
  res = nc.variables.has_key(name)
  nc.close()

  if (not res) and msg:
    msg.msg = 'Variable %s doesn\'t exist in file "%s".' % (name, filename)

  return res
  
##################################################################
def ncIsRec(filename, name, msg=None):
  if not ncHasDim(filename, name): return False
  
  nc = Dataset(filename)
  res = nc.dimensions[name].isunlimited()
  nc.close()

  if (not res) and msg:
    msg.msg = '%s is not a record dimension in file "%s".' % (name, filename)

  return res

##################################################################
def ncMinValueEquals(filename, var, slicestr, value, fill=False, msg=None):
  """
  Test min value for slice. Slice can be '...', ':,:', '0,2,:'. 
  If fill, then ignores 'value' and checks min value for flag.
  """
  if not ncHasVar(filename, var): return False
  
  nc = Dataset(filename)
  minv = eval('nc.variables[var][%s].min()' % (slicestr))
  
  if fill:
    res = (str(minv)=='--')
  else:
    res = (np.array(value).astype(minv.dtype) == minv)
    
  nc.close()

  if (not res) and msg:
    msg.msg = 'Variable %s[%s] min value %s <> %s in file "%s".' % (var, slicestr, str(minv), str(value), filename)

  return res

##################################################################
def ncMaxValueEquals(filename, var, slicestr, value, fill=False, msg=None):
  """ Test min value for slice, such as :,:. """
  if not ncHasVar(filename, var): return False
  
  nc = Dataset(filename)
  maxv = eval('nc.variables[var][%s].max()' % (slicestr))
  if fill:
    res = (str(maxv)=='--')
  else:
    res = (np.array(value).astype(maxv.dtype) == maxv)
    
  nc.close()

  if (not res) and msg:
    msg.msg = 'Variable %s[%s] max value %s <> %s in file "%s".' % (var, slicestr, str(maxv), str(value), filename)

  return res

##################################################################
def ncShapesEqual(filename1, var1, filename2, var2, msg=None):
  """ Compare between two variables in two files. """
  has1 = ncHasVar(filename1, var1)
  has2 = ncHasVar(filename2, var2)
  if has1 and not has2: return False
  if not has1 and has2: return False
  if (not has1) and (not has2): return True
  
  nc1 = Dataset(filename1)
  nc2 = Dataset(filename2)
  v1 = nc1.variables[var1]
  v2 = nc2.variables[var2]
  
  s1 = str(v1.shape)
  s2 = str(v2.shape)
  res = ( s1 == s2 )
  nc1.close()
  nc2.close()

  if (not res) and msg:
    msg.msg = 'Variable shapes %s%s <> %s%s in files "%s" and "%s".' % (var1, s1, var2, s2, filename1, filename2)

  return res

##################################################################
def ncValueEquals(filename, var, pos, value, msg=None, debug=False):
  """ 
  Test value at an index in variable such as '0,1,2'.
  pos: If string, then indices like '1,2,3' for 3d var, or a tuple or list [1,2,3].
  """
  if not ncHasVar(filename, var): return False
  
  nc = Dataset(filename)
  if type(pos) in [type([]), type(())]:
    slice = str(pos)[1:-1]
  else: # Single value or string already.
    slice = str(pos)
    
  v = eval('nc.variables[var][%s]' % (slice))
  if debug: print var,pos,'=',v
  res = (np.array(value).astype(v.dtype) == v)
  nc.close()

  if (not res) and msg:
    msg.msg = '%s%s = %s <> %s in file "%s".' % (var, str(pos), str(v), str(value), filename)

  return res

##################################################################
def ncValuesAre(filename, var, slice, any, value, msg=None):
  """ 
  Test slice if 'any' or 'all' are a digit value or 'min', 'max', 'nan', 'pinf', 'ninf', '-0', 'denorm', 'finite'. 
  """
  if not ncHasVar(filename, var, msg=msg): return False

  res = True
  msg.msg = ''
  nc = Dataset(filename)
  data = eval('nc.variables[var][%s]' % (slice))
  nc.close()
  dtype = data.dtype
  real = False
  
  if dtype == np.int8 or dtype == np.int16 or dtype == np.int32 or dtype == np.int64:
    tinfo = np.iinfo(dtype)
  else:
    tinfo = np.finfo(dtype)
    real = True
  
  if value == 'min':
    v = tinfo.min
    res = eval('(data == v).%s()' % (any))
  elif value == 'max':
    v = tinfo.max
    res = eval('(data == v).%s()' % (any))
  elif value == 'finite':
    res = eval('(np.isfinite(data)).%s()' % (any))
    if not res:
      msg.msg = '%s[%s].%s <> finite in "%s".' % (var, slice, any, filename)
  elif real:
    if value == 'nan':
      res = eval('(np.isnan(data)).%s()' % (any))
    elif value == 'pinf':
      res = eval('(np.isposinf(data)).%s()' % (any))
    elif value == 'ninf':
      res = eval('(np.isneginf(data)).%s()' % (any))
    elif value == '-0':
      res = eval('( (np.copysign(1,data)==-1) * (data==0) ).%s()' % (any))
    elif value == 'denorm':
      v = tinfo.tiny * tinfo.eps
      res = eval('(data == v).%s()' % (any))
  else: 
    try:
      v = float(value)
    except:
      msg.msg = 'ncValuesAre received invalid value "%s".' % (value)
      return False
            
  return res
  
##################################################################
def ncValuesEqual(filename1, var1, filename2, var2, slice1, slice2=None, msg=None):
  """ Test slice of values between 2 files and variables. """
  if not ncHasVar(filename1, var1, msg=msg): return False
  if not ncHasVar(filename2, var2, msg=msg): return False
  
  if slice2==None: slice2 = slice1
    
  nc1 = Dataset(filename1)
  nc2 = Dataset(filename2)
  res = eval('(nc1.variables[var1][%s] == nc2.variables[var2][%s]).all()' % (slice1, slice2))
  nc1.close()
  nc2.close()

  if (not res) and msg:
    msg.msg = '%s[%s] <> %s[%s] in files "%s" and "%s".' % (var1, slice1, var2, slice2, filename1, filename2)

  return res

##################################################################
def ncVarsEqual(filename1, var1, filename2, var2, atts=False, data=False, msg=None):
  """ Compare metadata and/or data between two variables in two files. """
  if not ncShapesEqual(filename1, var1, filename2, var2, msg=msg): return False
    
  res = True
  
  if atts:
    nc1 = Dataset(filename1)
    nc2 = Dataset(filename2)
    v1 = nc1.variables[var1]
    v2 = nc2.variables[var2]
    atts1 = v1.ncattrs()
    atts2 = v2.ncattrs()
    n1 = len(atts1)
    n2 = len(atts2)
    if n1 <> n2:
      if msg:
        msg.msg = 'Number of attributes of %s %d <> %d for %s in files "%s" and "%s".' % (var1, n1, n2, var2, filename1, filename2)
      res = False
    
    if res and (atts1 <> atts2):
      if msg:
        msg.msg = 'Attribute names of %s %s <> %s for %s in files "%s" and "%s".' % (var1, str(atts1), var2, str(atts2), filename1, filename2)
      res = False
    
    if res:
      names = []
      for name in atts1: names.append(name)
      
      for att in names:
        attval1 = getattr(v1,att)
        attval2 = getattr(v2,att)
        if attval1 <> attval2:
          if msg:
            msg.msg = 'Attribute %s.%s=%s <> %s.%s=%s in files "%s" and "%s".' % (var1, att, str(attval1), var2, att, str(attval2), filename1, filename2)
          res = False
          break

    nc1.close()
    nc2.close()

    if res:
      res = ncVarTypesEqual(filename1, var1, filename2, var2, msg=msg)
    
  if res and data:
    res = ncAllValuesEqual(filename1, var1, filename2, var2, msg=msg)
    
  if (not res) and msg:
    msg.msg = 'Variable %s <> %s in files "%s" and "%s".' % (var1, var2, filename1, filename2)

  return res
    
##################################################################
def ncVarTypeEquals(filename, name, typecode, msg=None):
  if not ncHasVar(filename, name, msg=msg): return False

  nc = Dataset(filename)
  dtype = nc.variables[name].dtype
  nc.close()
  
  res = (dtype == typecode)
  if (not res) and msg:
    msg.msg = 'Variable %s datatype %s <> %s in file "%s".' % (name, str(dtype), str(typecode), filename)
    
  return res
  
##################################################################
def ncVarTypesEqual(filename1, name1, filename2, name2, msg=None):
  if not ncHasVar(filename1, name1, msg=msg): return False
  if not ncHasVar(filename2, name2, msg=msg): return False

  nc1 = Dataset(filename1)
  nc2 = Dataset(filename2)
  dtype1 = nc1.variables[name1].dtype
  dtype2 = nc2.variables[name2].dtype
  nc1.close()
  nc2.close()
  
  res = (dtype1 == dtype2)
  if (not res) and msg:
    msg.msg = 'Variable %s datatype %s <> %s for variable %s in files "%s" and "%s".' % (name1, str(dtype1), str(dtype2), name2, filename1, filename2)
    
  return res

##################################################################
def ncAllVarsEqual(filename1, filename2, atts=False, data=False, msg=None):
  """ Compare all variables between two files. """
  res = True
  nc1 = Dataset(filename1)
  nc2 = Dataset(filename2)
  vars1 = nc1.variables
  vars2 = nc2.variables
  nv1 = len(vars1)
  nv2 = len(vars2)
  
  if nv1 <> nv2:
    if msg:
      msg.msg = 'Number of variables %d <> %d in files "%s" and "%s".' % (nv1, nv2, filename1, filename2)
      
    res = False
    
  if res:
    names1, names2 = [],[]
    for name in vars1: names1.append(name)
    for name in vars2: names2.append(name)
    if names1 <> names2:
      if msg:
        msg.msg = 'Variable names %s <> %s in files "%s" and "%s".' % (str(vars1), str(vars2), filename1, filename2)
      res = False
    
  if res:
    for var in vars1:
      if not ncVarsEqual(filename1, var, filename2, var, atts=atts, data=data, msg=msg):
        res = False
        break
  
  nc1.close()
  nc2.close()

  return res

