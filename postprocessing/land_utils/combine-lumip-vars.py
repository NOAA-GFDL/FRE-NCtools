#!/usr/bin/env python

import argparse
import sys
import numpy as np
import netCDF4 as nc
import time

def cloneAttrs(src,dst):
    '''
    Copies all attributes from src to dst. Returns dst.
    '''
    # it may be more efficient to use dst.setncatts(src.__dict__),
    # except I am not absolutely sure __dict__ doesn't contain something
    # else besides netcdf attributes
    for att in src.ncattrs() :
        dst.setncattr(att,src.getncattr(att))
    return dst

def cloneDim(src,dimname,dst,name=None):
    '''
    Clones dimension, dimension variable, and its data. Returns new 
    dimension.
    src     -- netcdf dataset
    dimname -- name of the dimension
    dst     -- netcdf data set
    name    -- new name for the dimension and variable.
    
    If the src dimension is unlimited, then the dst dimension is also
    unlimited, of zero length. The dimension variable data are not
    copied in this case.
    '''
    if name is None:
        name=dimname
    srcDim=src.dimensions[dimname]
    dst.createDimension(name,size=None if srcDim.isunlimited() else len(srcDim))
    if dimname in src.variables:
       srcVar=src.variables[dimname]
       # what about fill value, etc?
       dstVar=dst.createVariable(name,srcVar.dtype,(name,))
       # the line below doesn't work for MFDataset for unlimited dimension
       dstVar.setncatts(srcVar.__dict__)
       if not srcDim.isunlimited():
           dstVar[:] = srcVar[:]
    return dst.dimensions[name]

def cloneVar(src,varname,dst,name=None):
    '''
    Clones a variable definition from src to dst with all dimensions and 
    attributes, possibly renaming it. Variable data are not copied.
    Returns new variable.
    src     -- netcdf dataset
    varname -- name of the variable
    dst     -- destination netcdf dataset
    name    -- new name of the variable
    Note that if dimensions in dst are of different size, the size of the
    variable will match them.
    '''
    if name is None:
        name=varname
    srcVar=src.variables[varname]
    # what about fill_value, etc?
    dstVar=dst.createVariable(name,srcVar.dtype,srcVar.dimensions)    
    cloneAttrs(srcVar,dstVar)
    return dstVar

def fill_value(src):
    if '_FillValue' in src.ncattrs():
        return src.getncattr('_FillValue')
    elif 'missing_value' in src.ncattrs():
        return src.getncattr('mising_value')
    else:
        return None

parser=argparse.ArgumentParser(description='Combine LUT variables')
parser.add_argument('-o', '--output', required=True, help='output netcdf file')
parser.add_argument('--variables', metavar='VAR,VAR,...',
    help='''comma-separated list of variables to process. By default all variables from 
            input netcdf files are processesed.''')
parser.add_argument('-v','--verbose', dest='verb', action='count', help='increase verbosity')
parser.add_argument('input', nargs='+', 
    help='''input netcdf files for each land use tile, in the following order: psl, crp, 
            pst, urb. Urban file may be missing, in which case all urban data are filled 
            with missing values.''')
parser.add_argument('-H','--no-history', dest='history', action='store_false', 
    help='Do not append to "history" global attribute.')
args = parser.parse_args()

# names of the land use types, saved in the discrete coordinate variable in output netcdf
lu_names=['psl','crp','pst','urb']

if args.verb:
    print 'input files:'
    for lu,filename in zip(lu_names, args.input):
       print '    %s: "%s"'%(lu,filename)
    print 'output file:\n    "%s"'%args.output

# open input files
src = [nc.Dataset(arg,'r') for arg in args.input]

# create the output file using the first file in the list of input files as a template.
template = src[0]
dst = nc.Dataset(args.output,'w',format=template.data_model)

# clone all global attributes
# TODO: use correct "_NCProperties" attribute
# TODO: use correct "filename" attribute, or delete it altogether
for att in template.ncattrs() :
    dst.setncattr(att,template.getncattr(att))
# add to history attribute
if args.history :
    dst.history = (template.history+"\x0A" if 'history' in template.ncattrs() else '') + \
                  time.strftime("%c",time.localtime())+': '+' '.join(sys.argv)


# clone all dimensions from source file to destination
for d in template.dimensions.iterkeys() :
    if args.verb > 1:
        print 'cloning dimension "%s"'%d
    cloneDim(template,d,dst)

# copy unlimited dimension (time)
for name,dim in template.dimensions.iteritems():
   if dim.isunlimited():
       if args.verb > 1 :
           print 'copying unlimited dimension "%s" data'%name
       for i in xrange(dim.size):
           dst.variables[name][i] = template.variables[name][i]

# create land use type dimension
dst.createDimension('landusetype4',size = len(lu_names))

# create variable for land use type names and write the names into it
lu_name_len = max(len(x) for x in lu_names)
dst.createDimension('landusenameidx',size=lu_name_len)
namevar = dst.createVariable('landusename','c',('landusetype4','landusenameidx'))
namevar.long_name = 'names of land use types'
for i,name in enumerate(lu_names):
    namevar[i,:] = nc.stringtoarr(name,lu_name_len)


# find auxiliary (coordinate-related variables) in the netcdf, such as boundaries and 
# time averaging information variables: they are not going to be combined by land use,
# but are going to be copied to the output file
auxVars = set()
auxAttrs = {'bounds','edges','time_avg_info'} # attributes that may list auxiliary variables 
for var in template.variables.itervalues():
    for attr in auxAttrs:
        if attr not in var.ncattrs(): continue
        for v in var.getncattr(attr).split(','):
            if v not in template.dimensions: auxVars.add(v)

# find variables: if the list of variables is not provided on the command line, process
# all variables in the netcdf files, except dimension variables and averaging information 
# (average_DT, average_T1, average_T2)
if args.variables:
    vars = args.variables.split(',')
else:
    vars = template.variables.keys()
    vars = [ v for v in vars if v not in template.dimensions ]
    vars = [ v for v in vars if v not in auxVars ]

# copy all auxiliary (coordinate-related) variables
# exclude variables requested on command line from auxiliary variables
auxVars = [v for v in auxVars if v not in vars]
for varname in auxVars:
    if args.verb > 1:
        print 'processing auxiliary variable "%s"'%varname
    cloneVar(template, varname, dst)[:] = template.variables[varname][:]

# combine variables from different land use files
for varname in vars:
    if args.verb:
        print 'processing variable "%s"'%varname

    # created output variable
    invar  = template.variables[varname] # variable used as a template for destination variable

    if (template.dimensions[invar.dimensions[0]].isunlimited()) :
        outdim = (invar.dimensions[0],'landusetype4')+invar.dimensions[1:]
    else:
        outdim = ('landusetype4',)+invar[0].dimensions
    outvar=dst.createVariable(varname,invar.dtype,dimensions=outdim,fill_value=fill_value(invar))
    # clone relevant attributes from the template
    for att in invar.ncattrs() :
        outvar.setncattr(att,invar.getncattr(att))

    # add coordinates attribute (per CF convention, Section 6.1)
    outvar.coordinates = 'landusename'

    # copy data
    invars = [f.variables[varname] for f in src]
    if (template.dimensions[invar.dimensions[0]].isunlimited()) :
        for i in xrange(template.dimensions[invar.dimensions[0]].size) :
            for j,var in enumerate(invars) :
                outvar[i,j] = invars[j][i]
    else:
        for j,var in enumerate(invars) :
            outvar[j] = invars[j]

for f in src: f.close()
dst.close()