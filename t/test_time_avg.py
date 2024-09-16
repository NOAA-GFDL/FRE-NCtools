#!/usr/bin/env python3

import xarray as xr
import numpy as np
import os

def init_fake_data(a, count=0., factor=1.):
    with np.nditer(a, op_flags=['readwrite']) as it:
        for x in it:
            count = count + 1
            x[...] = count * factor

class TestTimeAvg():
    def create_input(self):
        print("Creating input files to test with")
        nx = 96
        ny = 96
        nt = 12

        a = np.zeros([nt, ny, nx])
        for k in range(0,nt):
          for j in range(0,ny):
            for i in range(0,nx):
                a[k, j, i] = (i+1.)*1000 + (j+1.)*10 + (k+1.)/100.
        self.a = a

        Time_attrs = {"units": "days since 2000-01-01", "axis": "T", "calendar_type": "JULIAN", "calendar": "julian"} 
        ds1 = xr.Dataset(
        data_vars={
            ## This is backwards #FORTRAN4LYFE
            "a": (("Time", "y", "x"), a),
            "Time": (("Time", np.arange(nt)*1.0+1, Time_attrs))
        },
        coords={
            "x": np.arange(nx)*1.0,
            "y": np.arange(ny)*1.0,
        },
        )
        ds1.to_netcdf("20120101.ice_shelf.nc", unlimited_dims="Time")


    def check_output(self):
      out_file = xr.open_dataset("wut.nc")
      if not (out_file.a.values == np.average(self.a, axis=0)).all():
          raise Exception("Test failed")


test_class = TestTimeAvg()
test_class.create_input()
try:
    os.system("TAVG.exe < input.nml")
except:
    raise Exception("Failed to run timavg")
test_class.check_output()