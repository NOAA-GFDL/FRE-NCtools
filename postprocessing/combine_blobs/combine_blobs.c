#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include <string.h>
#include <fcntl.h>

// holds metadata information related to netcdf files
struct fileinfo {
  int ncid;  // id for the netcdf file
  int ndims;  // number of dimensions 
  int nvars;  // number of variables
  int ngatts;  // number of global attributes
  int unlimdimid;  // id of the unlimited dimension
  char varname[MAX_NC_VARS][MAX_NC_NAME];  // names of the variables
  long varsize[MAX_NC_VARS]; // sizes of the variables
  nc_type datatype[MAX_NC_VARS]; // data types of the variables
  int varndims[MAX_NC_VARS];  // number of dimensions for each variable
  int vardim[MAX_NC_VARS][MAX_NC_DIMS];  // dimension ids for each variable
  int natts[MAX_NC_VARS];  // number of attributes for each variable
  char dimname[MAX_NC_DIMS][MAX_NC_NAME];  // names of the dimensions
  size_t dimsize[MAX_NC_DIMS];  //sizes of the dimensions
};

// holds information necessary to extract slices of arrays
typedef struct {
  long start; // first index of slice
  long stop; // last index of slice
} slice;

// prints error message when netcdf functions fail
void handle_error(int status, char *message) {
  fprintf(stderr, "Oops: %s\nError: %s\n", message, nc_strerror(status));
  exit(-1);
}

void usage() {
  printf("combine_blobs: missing files\n");
  printf("combine_blobs input1.nc input2.nc ... output.nc\n");
  exit(2);
}

// calculate the size of all variables based on their dimensions
void inq_var_size(struct fileinfo *ncfile) {
  int v, d; // loop variables for netcdf variables and dimensions
  int size;
  int dimid;
  for (v = 0; v < ncfile->nvars; v++) {
    size = 1;
    for (d = 0; d < ncfile->varndims[v]; d++) {
      dimid = ncfile->vardim[v][d];
      size = size * ncfile->dimsize[dimid];
    } 
    ncfile->varsize[v] = size;
  }
}

void get_dimension_sizes(struct fileinfo *ncfile) {
  int d; 
  int status;
  for (d = 0; d < ncfile->ndims; d++) {
    status = nc_inq_dim(ncfile->ncid,d,ncfile->dimname[d],
                        &(ncfile->dimsize[d]));
    if (status != NC_NOERR) handle_error(status, "Problem making dimension inquiry.");
  }
}

// as the outputfile is written to and re-read, it's necessary to update our
// internal understanding of what is in the file (sizes of dims and vars)
void refresh_metadata(struct fileinfo *ncfile) {
  get_dimension_sizes(ncfile);
  inq_var_size(ncfile);
}

void open_file(char *filename, struct fileinfo *ncfile) {
  int v, d; // loop variables for netcdf variables and dimensions
  int status; // return status code for netcdf functions

  // open the netcdf file to read the metadata/data
  status = nc_open(filename,NC_NOWRITE,&ncfile->ncid);
  if (status != NC_NOERR) handle_error(status, "Problem opening netcdf file.");

  // collect some information about the file from its metadata
  nc_inq(ncfile->ncid, &ncfile->ndims, &ncfile->nvars, 
	 &ncfile->ngatts, &ncfile->unlimdimid);
  if (status != NC_NOERR) handle_error(status, "Problem making netcdf inquiry.");

  // collect some information about each variable from the metadata
  for (v = 0; v < ncfile->nvars; v++) {
    status = nc_inq_var(ncfile->ncid, v, ncfile->varname[v], &(ncfile->datatype[v]), 
			&(ncfile->varndims[v]), ncfile->vardim[v], &(ncfile->natts[v]));
    if (status != NC_NOERR) handle_error(status, "Problem making variable inquiry.");
  }

  // collect some information about each dimension from the metadata
  get_dimension_sizes(ncfile);

  // calculate the size of each variable
  inq_var_size(ncfile);
}

// given a variable name, retrieve the size of that netcdf variable
int get_varsize(char *varname, struct fileinfo *ncfile) {
  int v;
  for (v = 0; v < ncfile->nvars; v++) {
    if (strcmp(ncfile->varname[v], varname)) {
      break;
    }
  }
  return ncfile->varsize[v];
}

void close_file(struct fileinfo *ncfile) {
  int status;
  status = nc_close(ncfile->ncid);
  if (status != NC_NOERR) handle_error(status, "Problem closing file.");
  free(ncfile);
}

void append_slice(slice slices[], long start, long stop, int position) {
  slice varslice;
  varslice.start = start;
  varslice.stop = stop;
  slices[position] = varslice;
}

void write_values(slice slices[], int varid, struct fileinfo *ncinfile, struct fileinfo *ncoutfile) {
  int status;
  void *invarval; // buffer to hold variable data read from input file
  size_t type_size;

  status = nc_inq_type(ncinfile->ncid, ncinfile->datatype[varid], NULL, &type_size);
  if (status != NC_NOERR) handle_error(status, "Problem inquiring about types.");

  invarval = (void *) calloc(ncinfile->varsize[varid], type_size);

  int numdims = ncinfile->varndims[varid];

  // the nc_get_var command stated that only netdf4 files were legal, so
  // these vectors define the start and length of the array of data requested
  // from the get_vara function
  size_t put_start[numdims];
  size_t put_count[numdims];
  size_t get_start[numdims];
  size_t get_count[numdims];

  int i;
  long span;
  for (i = 0; i < numdims; i++) {
    get_start[i] = 0;
    get_count[i] = ncinfile->dimsize[i];
    put_start[i] = slices[i].start;
    span = slices[i].stop - slices[i].start;
    put_count[i] = span;
  }

  status = nc_get_vara(ncinfile->ncid, varid, get_start, get_count, invarval);
  if (status != NC_NOERR) handle_error(status, "Problem reading values from input file.");

  status = nc_put_vara(ncoutfile->ncid, varid, put_start, put_count, invarval);
  if (status != NC_NOERR) handle_error(status, "Problem writing values to output file.");
}

void create_outputfile(struct fileinfo *ncseedfile, struct fileinfo *ncoutfile, char *outfilename) {
  int status;
  status = nc_create(outfilename, 0, &(ncoutfile->ncid));
  if (status != NC_NOERR) handle_error(status, "Problem creating output file.");

  status = nc_inq(ncseedfile->ncid, &ncoutfile->ndims, &ncoutfile->nvars,
		  &ncoutfile->ngatts, &ncoutfile->unlimdimid);
  if (status != NC_NOERR) handle_error(status, "Problem inquring into seedfile.");

  char att_name[NC_MAX_NAME];
  size_t att_len;
  nc_type att_type;
  void *att_val;
  size_t type_size;

  int a;
  for (a = 0; a < ncseedfile->ngatts; a++) {
    status = nc_inq_attname(ncseedfile->ncid, NC_GLOBAL, a, att_name);
    if (status != NC_NOERR) handle_error(status, "Problem retrieving attribute name.");
    status = nc_inq_att(ncseedfile->ncid, NC_GLOBAL, att_name, &att_type, &att_len);
    if (status != NC_NOERR) handle_error(status, "Problem inquiring about global attribute from seed file.");
    status = nc_inq_type(ncseedfile->ncid, att_type, NULL, &type_size);
    if (status != NC_NOERR) handle_error(status, "Problem inquiring about attribute type.");
    att_val = (void *) calloc(att_len, type_size);
    status = nc_get_att(ncseedfile->ncid, NC_GLOBAL, att_name, att_val);
    if (status != NC_NOERR) handle_error(status, "Problem reading global attribute from seed file.");
    status = nc_put_att(ncoutfile->ncid, NC_GLOBAL, att_name, att_type, att_len, att_val);
    if (status != NC_NOERR) handle_error(status, "Problem writing global attribute to output file.");
    free(att_val);
  }
  int dimid;

  int d;
  for (d = 0; d < ncseedfile->ndims; d++) {
    status = nc_inq_dim(ncseedfile->ncid, d, ncoutfile->dimname[d], &(ncoutfile->dimsize[d]));
    if (status != NC_NOERR) handle_error(status, "Problem reading dimension from seed file.");
    if (d == ncseedfile->unlimdimid) {
      status = nc_def_dim(ncoutfile->ncid, ncoutfile->dimname[d], ncoutfile->dimsize[d], &dimid);
    } else {
      status = nc_def_dim(ncoutfile->ncid, ncoutfile->dimname[d], NC_UNLIMITED, &dimid);
    }
    if (status != NC_NOERR) handle_error(status, "Problem writing dimension to seed file.");
  }

  int varid;
  int numdims;

  int v;
  for (v = 0; v < ncseedfile->nvars; v++) {
    status = nc_inq_var(ncseedfile->ncid, v, ncoutfile->varname[v], &(ncoutfile->datatype[v]), 
			&(ncoutfile->varndims[v]), ncoutfile->vardim[v], &(ncoutfile->natts[v]));
    if (status != NC_NOERR) handle_error(status, "Problem inquiring into seed variable.");
    
    numdims = ncoutfile->varndims[v];

    int dimids[numdims];
    
    int d;
    for (d = 0; d < numdims; d++) {
      dimids[d] = d;
    }
    
    status = nc_def_var(ncoutfile->ncid, ncoutfile->varname[v], ncoutfile->datatype[v], 
			ncoutfile->varndims[v], dimids, &varid);
    if (status != NC_NOERR) handle_error(status, "Problem creating new variable in output file.");

    int va;
    for (va = 0; va < ncseedfile->natts[v]; va++) {
      status = nc_inq_attname(ncseedfile->ncid, v, va, att_name);
      if (status != NC_NOERR) handle_error(status, "Problem retrieving var attribute name.");
      status = nc_inq_att(ncseedfile->ncid, v, att_name, &att_type, &att_len);
      if (status != NC_NOERR) handle_error(status, "Problem inquiring about var attribute from seed file.");
      status = nc_inq_type(ncseedfile->ncid, att_type, NULL, &type_size);
      if (status != NC_NOERR) handle_error(status, "Problem inquiring about var attribute type.");
      att_val = (void *) calloc(att_len, type_size);
      status = nc_get_att(ncseedfile->ncid, v, att_name, att_val);
      if (status != NC_NOERR) handle_error(status, "Problem reading var attribute from seed file.");
      status = nc_put_att(ncoutfile->ncid, v, att_name, att_type, att_len, att_val);
      if (status != NC_NOERR) handle_error(status, "Problem writing var attribute to output file.");      
      free(att_val);
    }
  }
  status = nc_enddef(ncoutfile->ncid);
  if (status != NC_NOERR) handle_error(status, "Problem leaving define mode.");

  inq_var_size(ncoutfile);

}

int main(int argc, char *argv[]) {
  if (argc < 3) usage();
  char outfilename[2048];
  char infilename[2048];
  sprintf(outfilename, "%s", argv[argc - 1]);
  char seedfilename[2048];
  sprintf(seedfilename, "%s", argv[1]);
  struct fileinfo *ncoutfile;
  struct fileinfo *ncinfile;
  struct fileinfo *ncseedfile;
  int status; // return status for various netcdf functions
  int dimid;
  char *dimname;
  int varsize;

  // allocate a chunk of memory for the netcdf file info
  if ((ncoutfile=(struct fileinfo *)calloc(1, sizeof(struct fileinfo)))==NULL) {
      fprintf(stderr,"Error: cannot allocate enough memory for a file!\n"); 
      exit(1);
  }

  if ((ncseedfile=(struct fileinfo *)calloc(1, sizeof(struct fileinfo)))==NULL) {
      fprintf(stderr,"Error: cannot allocate enough memory for a file!\n"); 
      exit(1);
  }

  open_file(seedfilename, ncseedfile);
  create_outputfile(ncseedfile, ncoutfile, outfilename);
  close_file(ncseedfile);

  int i, v, d; // loop variables for input files, netcdf variables, and dimensions
  
  for (i = 1; i < argc - 1; i++) {
    sprintf(infilename, "%s", argv[i]);

    // allocate a chunk of memory for the netcdf file info
    if ((ncinfile=(struct fileinfo *)calloc(1, sizeof(struct fileinfo)))==NULL) {
      fprintf(stderr,"Error: cannot allocate enough memory for a file!\n");
      exit(1);
    }

    open_file(infilename, ncinfile);

    unsigned char flag = 0;    
    long start = 0;
    long stop = 0;    
    int status;
    long invarsize;
    long outvarsize;

    for (v = 0; v < ncinfile->nvars; v++) {
      slice slices[MAX_NC_DIMS];
      for (d = 0; d < ncinfile->varndims[v]; d++) {
	dimid = ncinfile->vardim[v][d];
	dimname = ncinfile->dimname[d];
	outvarsize = get_varsize(dimname, ncoutfile);
	invarsize = get_varsize(dimname, ncinfile);
	if (ncinfile->unlimdimid == dimid) {
	  if (flag) {
	    stop = outvarsize;
	    start = outvarsize - invarsize;
	    append_slice(slices, start, stop, d);
	  }
	  else {
	    start = outvarsize;
	    stop = outvarsize + invarsize;
	    flag = 1;
	    append_slice(slices, start, stop, d);
	  }
	}
	else {
	  start = 0;
	  stop = invarsize;
	  append_slice(slices, start, stop, d);	 
	}
      }
      if (stop != start) {
	write_values(slices, v, ncinfile, ncoutfile);
	refresh_metadata(ncoutfile); // need to refresh dim and var metadata
      }
    }	
    close_file(ncinfile);
  }
  close_file(ncoutfile);
  return 0;
}  
