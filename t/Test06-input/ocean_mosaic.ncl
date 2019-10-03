netcdf ocean_mosaic {
dimensions:
	ntiles = 1 ;
	ncontact = 2 ;
	string = 255 ;
variables:
	char mosaic(string) ;
		mosaic:standard_name = "grid_mosaic_spec" ;
		mosaic:mosaic_spec_version = "0.2" ;
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
  "360:360,1:180::1:1,1:180",
  "1:180,180:180::360:181,180:180" ;
}
