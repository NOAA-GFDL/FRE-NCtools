netcdf C48_mosaic {
dimensions:
	ntiles = 6 ;
	ncontact = 12 ;
	str = 255 ;
variables:
	char mosaic(str) ;
		mosaic:standard_name = "grid_mosaic_spec" ;
		mosaic:children = "gridtiles" ;
		mosaic:contact_regions = "contacts" ;
		mosaic:grid_descriptor = "" ;
	char gridlocation(str) ;
		gridlocation:standard_name = "grid_file_location" ;
	char gridfiles(ntiles, str) ;
	char gridtiles(ntiles, str) ;
	char contacts(ncontact, str) ;
		contacts:standard_name = "grid_contact_spec" ;
		contacts:contact_type = "boundary" ;
		contacts:alignment = "true" ;
		contacts:contact_index = "contact_index" ;
		contacts:orientation = "orient" ;
	char contact_index(ncontact, str) ;
		contact_index:standard_name = "starting_ending_point_index_of_contact" ;

// global attributes:
		:grid_version = "0.2" ;
		:code_version = "$Name: tikal_prerelease $" ;
		:history = "make_solo_mosaic --num_tiles 6 --dir ./ --mosaic C48_mosaic --tile_file C48_grid.tile1.nc,C48_grid.tile2.nc,C48_grid.tile3.nc,C48_grid.tile4.nc,C48_grid.tile5.nc,C48_grid.tile6.nc" ;
data:

 mosaic = "C48_mosaic" ;

 gridlocation = "./" ;

 gridfiles =
  "C48_grid.tile1.nc",
  "C48_grid.tile2.nc",
  "C48_grid.tile3.nc",
  "C48_grid.tile4.nc",
  "C48_grid.tile5.nc",
  "C48_grid.tile6.nc" ;

 gridtiles =
  "tile1",
  "tile2",
  "tile3",
  "tile4",
  "tile5",
  "tile6" ;

 contacts =
  "C48_mosaic:tile1::C48_mosaic:tile2",
  "C48_mosaic:tile1::C48_mosaic:tile3",
  "C48_mosaic:tile1::C48_mosaic:tile5",
  "C48_mosaic:tile1::C48_mosaic:tile6",
  "C48_mosaic:tile2::C48_mosaic:tile3",
  "C48_mosaic:tile2::C48_mosaic:tile4",
  "C48_mosaic:tile2::C48_mosaic:tile6",
  "C48_mosaic:tile3::C48_mosaic:tile4",
  "C48_mosaic:tile3::C48_mosaic:tile5",
  "C48_mosaic:tile4::C48_mosaic:tile5",
  "C48_mosaic:tile4::C48_mosaic:tile6",
  "C48_mosaic:tile5::C48_mosaic:tile6" ;

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
