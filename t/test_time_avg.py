#!/usr/bin/env python3

import xarray as xr
import numpy as np
import os

class TestTimeAvg:
    def __init__(self, Test_case):
        self.Test_case = Test_case
        if self.Test_case == 1:
            self.In_file = "timavg_instantaneous_in.nc"
            self.Out_file = "timavg_instantaneous_out.nc"
            self.Namelist_file = "instantaneous.nml"
            self.Error_msg = "timavg with instantaneous input"
        elif self.Test_case == 2:
            self.In_file = "timavg_standard_in.nc"
            self.Out_file = "timavg_standard_out.nc"
            self.Namelist_file = "standard.nml"
            self.Error_msg = "timavg with standard input"


    def create_instantaneous(self):
        Time_attrs = {"units": "days since 2000-01-01", "axis": "T", "calendar_type": "JULIAN", "calendar": "julian"}
        ds1 = xr.Dataset(
        data_vars={
            ## This is backwards #FORTRAN4LYFE
            "a": (("Time", "y", "x"), self.a),
            "Time": (("Time", self.time_data, Time_attrs))
        },
        coords={
            "x": np.arange(self.nx)*1.0,
            "y": np.arange(self.ny)*1.0,
        },
        )
        ds1.to_netcdf(self.In_file, unlimited_dims="Time")

    def create_standard(self):
        Time_attrs = {"units": "days since 2000-01-01", "axis": "T", "calendar_type": "JULIAN", "calendar": "julian",
                      "bounds": "time_bnds"}
        ds1 = xr.Dataset(
        data_vars={
            ## This is backwards #FORTRAN4LYFE
            "a": (("Time", "y", "x"), self.a),
            "Time": (("Time", self.time_data, Time_attrs)),
            "time_bnds": (("Time", "nv"), self.time_bnds_data)
        },
        coords={
            "x": np.arange(self.nx)*1.0,
            "y": np.arange(self.ny)*1.0,
            "nv": np.arange(2)*1.0+1,
        },
        )
        ds1.to_netcdf(self.In_file, unlimited_dims="Time")


    def create_input(self):
        self.nx = 96
        self.ny = 96
        self.nt = 12
        self.time_data = np.arange(self.nt)*1.0+1

        self.a = np.zeros([self.nt, self.ny, self.nx])
        for k in range(0,self.nt):
            for j in range(0,self.ny):
              for i in range(0,self.nx):
                  self.a[k, j, i] = (i+1.)*1000 + (j+1.)*10 + (k+1.)/100.

        if self.Test_case == 1:
            self.create_instantaneous()
        elif self.Test_case == 2:
            self.time_bnds_data = np.zeros([self.nt, 2])
            self.time_bnds_data[:,1] = self.time_data[0:]
            self.time_bnds_data[1:,0] = self.time_data[0:-1]
            self.create_standard()


    def run_test(self):
      exit_code = os.system("TAVG.exe < " + self.Namelist_file)
      if exit_code != 0 :
          raise Exception(self.Error_msg + ":: Failed to complete")


    def check_output(self):
      out_file = xr.open_dataset(self.Out_file)
      if not (out_file.a.values == np.average(self.a, axis=0)).all():
          raise Exception(self.Error_msg + ":: answers were not the expected result!")


# Test time_average when the input is instantaneous (so no time bounds)
Test_Instantaneous = 1
test_class = TestTimeAvg(Test_Instantaneous)
test_class.create_input()
test_class.run_test()
test_class.check_output()

# Test time_average when the input is standard (with time_bounds, no average_* variables)
Test_Standard = 2
test_class = TestTimeAvg(Test_Standard)
test_class.create_input()
test_class.run_test()
test_class.check_output()
