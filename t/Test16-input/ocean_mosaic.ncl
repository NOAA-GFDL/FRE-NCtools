netcdf ocean_mosaic {
dimensions:
	ntiles = 1 ;
	ncontact = 2 ;
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
		:history = "make_solo_mosaic --num_tiles 1 --dir ./ --mosaic_name ocean_mosaic --tile_file ocean_hgrid.nc --periodx 360" ;
data:

 mosaic = "ocean_mosaic" ;

 gridlocation = "./" ;

 gridfiles =
  "ocean_hgrid.nc" ;

 gridtiles =
  "tile1" ;

 contacts =
  "ocean_mosaic:tile1::ocean_mosaic:tile1",
  "ocean_mosaic:tile1::ocean_mosaic:tile1" ;

 contact_index =
  "720:720,1:400::1:1,1:400",
  "1:360,400:400::720:361,400:400" ;
}
