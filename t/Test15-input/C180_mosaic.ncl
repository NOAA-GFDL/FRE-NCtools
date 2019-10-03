netcdf C180_mosaic {
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
		:code_version = "$Name: riga_201012 $" ;
		:history = "make_solo_mosaic --num_tiles 6 --dir ./ --mosaic_name C180_mosaic --tile_file C180_grid.tile1.nc,C180_grid.tile2.nc,C180_grid.tile3.nc,C180_grid.tile4.nc,C180_grid.tile5.nc,C180_grid.tile6.nc" ;
data:

 mosaic = "C180_mosaic" ;

 gridlocation = "./" ;

 gridfiles =
  "C180_grid.tile1.nc",
  "C180_grid.tile2.nc",
  "C180_grid.tile3.nc",
  "C180_grid.tile4.nc",
  "C180_grid.tile5.nc",
  "C180_grid.tile6.nc" ;

 gridtiles =
  "tile1",
  "tile2",
  "tile3",
  "tile4",
  "tile5",
  "tile6" ;

 contacts =
  "C180_mosaic:tile1::C180_mosaic:tile2",
  "C180_mosaic:tile1::C180_mosaic:tile3",
  "C180_mosaic:tile1::C180_mosaic:tile5",
  "C180_mosaic:tile1::C180_mosaic:tile6",
  "C180_mosaic:tile2::C180_mosaic:tile3",
  "C180_mosaic:tile2::C180_mosaic:tile4",
  "C180_mosaic:tile2::C180_mosaic:tile6",
  "C180_mosaic:tile3::C180_mosaic:tile4",
  "C180_mosaic:tile3::C180_mosaic:tile5",
  "C180_mosaic:tile4::C180_mosaic:tile5",
  "C180_mosaic:tile4::C180_mosaic:tile6",
  "C180_mosaic:tile5::C180_mosaic:tile6" ;

 contact_index =
  "360:360,1:360::1:1,1:360",
  "1:360,360:360::1:1,360:1",
  "1:1,1:360::360:1,360:360",
  "1:360,1:1::1:360,360:360",
  "1:360,360:360::1:360,1:1",
  "360:360,1:360::360:1,1:1",
  "1:360,1:1::360:360,360:1",
  "360:360,1:360::1:1,1:360",
  "1:360,360:360::1:1,360:1",
  "1:360,360:360::1:360,1:1",
  "360:360,1:360::360:1,1:1",
  "360:360,1:360::1:1,1:360" ;
}
