netcdf C48_mosaic {
dimensions:
	ntiles = 6 ;
	ncontact = 12 ;
	string = 255 ;
variables:
	char mosaic(string) ;
		mosaic:standard_name = "grid_mosaic_spec" ;
		mosaic:children = "gridtiles" ;
		mosaic:contact_regions = "contacts" ;
		mosaic:grid_descriptor = "" ;
	char gridlocation(string) ;
		gridlocation:standard_name = "grid_file_location" ;
	char gridfiles(ntiles, string) ;
	char gridtiles(ntiles, string) ;
	char contacts(ncontact, string) ;
		contacts:standard_name = "grid_contact_spec" ;
		contacts:contact_type = "boundary" ;
		contacts:alignment = "true" ;
		contacts:contact_index = "contact_index" ;
		contacts:orientation = "orient" ;
	char contact_index(ncontact, string) ;
		contact_index:standard_name = "starting_ending_point_index_of_contact" ;

// global attributes:
		:grid_version = "0.2" ;
		:code_version = "$Name: omsk_20071218_z1l $" ;
		:history = "/home/z1l/bin/tools_32bit/make_solo_mosaic --num_tiles 6 --dir . --mosaic_name atmos_mosaic" ;
data:

 mosaic = "atmos_mosaic" ;

 gridlocation = "./" ;

 gridfiles =
  "horizontal_grid.tile1.nc",
  "horizontal_grid.tile2.nc",
  "horizontal_grid.tile3.nc",
  "horizontal_grid.tile4.nc",
  "horizontal_grid.tile5.nc",
  "horizontal_grid.tile6.nc" ;

 gridtiles =
  "tile1",
  "tile2",
  "tile3",
  "tile4",
  "tile5",
  "tile6" ;

 contacts =
  "atmos_mosaic:tile1::atmos_mosaic:tile2",
  "atmos_mosaic:tile1::atmos_mosaic:tile3",
  "atmos_mosaic:tile1::atmos_mosaic:tile5",
  "atmos_mosaic:tile1::atmos_mosaic:tile6",
  "atmos_mosaic:tile2::atmos_mosaic:tile3",
  "atmos_mosaic:tile2::atmos_mosaic:tile4",
  "atmos_mosaic:tile2::atmos_mosaic:tile6",
  "atmos_mosaic:tile3::atmos_mosaic:tile4",
  "atmos_mosaic:tile3::atmos_mosaic:tile5",
  "atmos_mosaic:tile4::atmos_mosaic:tile5",
  "atmos_mosaic:tile4::atmos_mosaic:tile6",
  "atmos_mosaic:tile5::atmos_mosaic:tile6" ;

 contact_index =
  "96:96,1:96::1:1,1:96",
  "1:96,96:96::1:1,96:1",
  "1:1,1:96::96:1,96:96",
  "1:96,1:1::1:96,96:96",
  "1:96,96:96::1:96,1:1",
  "96:96,1:96::96:1,1:1",
  "1:96,1:1::96:96,96:1",
  "96:96,1:96::1:1,1:96",
  "1:96,96:96::1:1,96:1",
  "1:96,96:96::1:96,1:1",
  "96:96,1:96::96:1,1:1",
  "96:96,1:96::1:1,1:96" ;
}
