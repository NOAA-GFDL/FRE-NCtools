netcdf grid_spec {
dimensions:
	string = 255 ;
	nfile_aXo = 6 ;
	nfile_aXl = 6 ;
	nfile_lXo = 6 ;
variables:
	char atm_mosaic_dir(string) ;
		atm_mosaic_dir:standard_name = "directory_storing_atmosphere_mosaic" ;
	char atm_mosaic_file(string) ;
		atm_mosaic_file:standard_name = "atmosphere_mosaic_file_name" ;
	char atm_mosaic(string) ;
		atm_mosaic:standard_name = "atmosphere_mosaic_name" ;
	char lnd_mosaic_dir(string) ;
		lnd_mosaic_dir:standard_name = "directory_storing_land_mosaic" ;
	char lnd_mosaic_file(string) ;
		lnd_mosaic_file:standard_name = "land_mosaic_file_name" ;
	char lnd_mosaic(string) ;
		lnd_mosaic:standard_name = "land_mosaic_name" ;
	char ocn_mosaic_dir(string) ;
		ocn_mosaic_dir:standard_name = "directory_storing_ocean_mosaic" ;
	char ocn_mosaic_file(string) ;
		ocn_mosaic_file:standard_name = "ocean_mosaic_file_name" ;
	char ocn_mosaic(string) ;
		ocn_mosaic:standard_name = "ocean_mosaic_name" ;
	char ocn_topog_dir(string) ;
		ocn_topog_dir:standard_name = "directory_storing_ocean_topog" ;
	char ocn_topog_file(string) ;
		ocn_topog_file:standard_name = "ocean_topog_file_name" ;
	char aXo_file(nfile_aXo, string) ;
		aXo_file:standard_name = "atmXocn_exchange_grid_file" ;
	char aXl_file(nfile_aXl, string) ;
		aXl_file:standard_name = "atmXlnd_exchange_grid_file" ;
	char lXo_file(nfile_lXo, string) ;
		lXo_file:standard_name = "lndXocn_exchange_grid_file" ;

// global attributes:
		:grid_version = "0.2" ;
		:code_version = "$Name: tikal_prerelease $" ;
		:history = "make_coupler_mosaic --atmos_mosaic C48_mosaic.nc --ocean_mosaic ocean_mosaic.nc --ocean_topog topog.nc --mosaic_name grid_spec" ;
data:

 atm_mosaic_dir = "./" ;

 atm_mosaic_file = "C48_mosaic.nc" ;

 atm_mosaic = "C48_mosaic" ;

 lnd_mosaic_dir = "./" ;

 lnd_mosaic_file = "C48_mosaic.nc" ;

 lnd_mosaic = "C48_mosaic" ;

 ocn_mosaic_dir = "./" ;

 ocn_mosaic_file = "ocean_mosaic.nc" ;

 ocn_mosaic = "ocean_mosaic" ;

 ocn_topog_dir = "./" ;

 ocn_topog_file = "topog.nc" ;

 aXo_file =
  "C48_mosaic_tile1Xocean_mosaic_tile1.nc",
  "C48_mosaic_tile2Xocean_mosaic_tile1.nc",
  "C48_mosaic_tile3Xocean_mosaic_tile1.nc",
  "C48_mosaic_tile4Xocean_mosaic_tile1.nc",
  "C48_mosaic_tile5Xocean_mosaic_tile1.nc",
  "C48_mosaic_tile6Xocean_mosaic_tile1.nc" ;

 aXl_file =
  "C48_mosaic_tile1XC48_mosaic_tile1.nc",
  "C48_mosaic_tile2XC48_mosaic_tile2.nc",
  "C48_mosaic_tile3XC48_mosaic_tile3.nc",
  "C48_mosaic_tile4XC48_mosaic_tile4.nc",
  "C48_mosaic_tile5XC48_mosaic_tile5.nc",
  "C48_mosaic_tile6XC48_mosaic_tile6.nc" ;

 lXo_file =
  "C48_mosaic_tile1Xocean_mosaic_tile1.nc",
  "C48_mosaic_tile2Xocean_mosaic_tile1.nc",
  "C48_mosaic_tile3Xocean_mosaic_tile1.nc",
  "C48_mosaic_tile4Xocean_mosaic_tile1.nc",
  "C48_mosaic_tile5Xocean_mosaic_tile1.nc",
  "C48_mosaic_tile6Xocean_mosaic_tile1.nc" ;
}
