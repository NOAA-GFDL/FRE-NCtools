#!/usr/bin/env python2

import os
import argparse

import netCDF4
import numpy as np

def get_tiles_nc(fname):
    nc = netCDF4.Dataset(fname, 'w')

    Time = nc.createDimension('Time', None)
    var_Time = nc.createVariable('Time', 'double', 'Time')
    var_Time.long_name = 'Time'
    var_Time.units = 'time level'
    var_Time.cartesian_axis = 'T'
  
    zaxis_1 = nc.createDimension('zaxis_1', 49)
    var_zaxis_1 = nc.createVariable('zaxis_1', 'double', 'zaxis_1')
    var_zaxis_1.long_name = 'zaxis_1'
    var_zaxis_1.units = 'none'
    var_zaxis_1.cartesian_axis = 'Z'

    yaxis_1 = nc.createDimension('yaxis_1', 256)
    var_yaxis_1 = nc.createVariable('yaxis_1', 'double', 'yaxis_1')
    var_yaxis_1.long_name = 'yaxis_1'
    var_yaxis_1.units = 'none'
    var_yaxis_1.cartesian_axis = 'Y'

    xaxis_1 = nc.createDimension('xaxis_1', 256)
    var_xaxis_1 = nc.createVariable('xaxis_1', 'double', 'xaxis_1')
    var_xaxis_1.long_name = 'xaxis_1'
    var_xaxis_1.units = 'none'
    var_xaxis_1.cartesian_axis = 'X'

    var_o3 = nc.createVariable('o3', 'double', ('Time', 'zaxis_1', 'yaxis_1', 'xaxis_1'))
    var_o3.long_name = 'o3'
    var_o3.units = 'none'
    var_o3.checksum = 'A32650720153EC98'
    var_o3.interp_method = 'conserve_order1'

    return nc

def create_tiles_nc(fname):
    nc = get_tiles_nc(fname)

    return nc

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="NetCDF generator")
    parser.add_argument('fv_tracer', help="Destination file to create, will be overwritten if exists")
    args = parser.parse_args()

    
create_tiles_nc(args.fv_tracer)
