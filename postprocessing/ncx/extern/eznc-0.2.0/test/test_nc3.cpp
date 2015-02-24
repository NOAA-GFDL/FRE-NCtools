/*
20111118 rsz Created. valgrind shows "Invalid free", so netcdf library not perfect memory-wise.
*/

#include "ezTest.hpp"
#include "netcdf.h"

#define handle_status(status) fprintf(stderr, "%s\n", nc_strerror(status));

std::string getFilename(ez::ezTestRunner& runner) {
  std::string f;
  
  if (runner.tmpdir.size()) {
    f.append(runner.tmpdir);
    
    if (runner.tmpdir.at(runner.tmpdir.size()-1) != '/')
      f.append("/");
  }
  
  f.append("tmp.nc");
  return f;
}

// This alone shows "Invalid free" issues with valgrind.
bool test_empty(ez::ezTestRunner& runner) {
  int ncid, status;
  //std::string filename = getFilename(runner);
    
  status = nc_create("/tmp/tmp.nc", NC_CLOBBER, &ncid);
  if (status) handle_status(status);
  ezassert( status == NC_NOERR );
  
  status = nc_enddef(ncid);
  if (status) handle_status(status);
  ezassert( status == NC_NOERR );
  
  status = nc_close(ncid);
  if (status) handle_status(status);
  ezassert( status == NC_NOERR );

  return true;
}

bool test_write3drec(ez::ezTestRunner& runner) {
  int ncid, status;
  int varid, dimids[4];
  std::string filename = getFilename(runner);
    
  status = nc_create(filename.c_str(), NC_CLOBBER, &ncid);
  if (status) handle_status(status);
  ezassert( status == NC_NOERR );
  
  status = nc_def_dim(ncid, "t", NC_UNLIMITED, &dimids[0]);
  if (status) handle_status(status);
  ezassert( status == NC_NOERR );
  
  status = nc_def_dim(ncid, "z", 10, &dimids[1]);
  if (status) handle_status(status);
  ezassert( status == NC_NOERR );

  status = nc_def_dim(ncid, "y", 20, &dimids[2]);
  if (status) handle_status(status);
  ezassert( status == NC_NOERR );

  status = nc_def_dim(ncid, "x", 30, &dimids[3]);
  if (status) handle_status(status);
  ezassert( status == NC_NOERR );

  status = nc_def_var(ncid, "v", NC_INT, 4, dimids, &varid);
  if (status) handle_status(status);
  ezassert( status == NC_NOERR );
  
  status = nc_set_fill(ncid, NC_FILL, 0);
  if (status) handle_status(status);
  ezassert( status == NC_NOERR );
  
  status = nc_enddef(ncid);
  if (status) handle_status(status);
  ezassert( status == NC_NOERR );

  size_t start[4], count[4];
  start[0] = 1;
  start[1] = start[2] = start[3] = 0;
  count[0] = 1;
  count[1] = 10;
  count[2] = 20;
  count[3] = 30;
  int data[10*20*30];
  std::fill(data, data+10*20*30, 0);
  status = nc_put_vara_int(ncid, varid, start, count, data);
  if (status) handle_status(status);
  ezassert( status == NC_NOERR );
    
  status = nc_close(ncid);
  if (status) handle_status(status);
  ezassert( status == NC_NOERR );

  return true;
}
  
int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  #define TEST(N) runner.add(&test_##N, "test_"#N);
  TEST(empty);
  TEST(write3drec);
  
  return runner.run();
}