#!/usr/bin/env python3

import netCDF4 as nc
import numpy as np

np.set_printoptions(formatter={'float': '{: 0.13f}'.format})
# Open the file
hg = nc.Dataset('horizontal_grid.tile2.nc', 'r')
lr = nc.Dataset('land_remap.tile2.nc', 'r')

# Read the data
lon_hg = hg.variables['x'][:][:]
lat_hg = hg.variables['y'][:][:]

lon_lr = lr.variables['lon'][:]
lat_lr = lr.variables['lat'][:]

# Close the file
hg.close()
lr.close()


#print(lon_hg[-2][1::2])
#print(lon_lr[:])

### Tile 1, 2
# for a, b in zip(lon_hg[-2][1::2], lon_lr):
#      diff = 0 if abs(a - b) < 1e-12 else a - b
#      print(f"{a:0.13f} {b:0.13f} {diff}")

### Tile 3, 4, 5
for a, b in zip(lon_hg[1][1::2], lon_lr):
     diff = 0 if abs(a - b) < 1e-12 else a - b
     print(f"{a:0.13f} {b:0.13f} {diff}")

### Tile 1, 2, 3, 4, 5
# for a, b in zip(lat_hg[1::2], lat_lr):
#   diff = 0 if abs(a[1] - b) < 1e-12 else a[1] - b
#   print(f"{a[1]:0.13f} {b:0.13f} {diff}")

