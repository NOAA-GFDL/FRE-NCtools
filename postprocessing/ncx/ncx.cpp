/*
20111029 rsz Created.
*/
#include "ncx.hpp"

int main(int argc, const char* argv[]) {
  int status;
  NcX ncx;
  status = ncx.SetOptions(argc, argv);
  if (status) return status;
  
  status = ncx.Run();
  if (status) return status;
  
  return 0;
}