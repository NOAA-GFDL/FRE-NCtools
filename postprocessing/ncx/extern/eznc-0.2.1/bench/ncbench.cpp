#include "ezNc.hpp"

int main(int argc, const char* argv[]) {
  double ticTotal = ez::osQueryPerformance();

  if (argc < (1+15+1)) {
    printf("Reports time and memory used to write a configurable netcdf file.\n\n");
    printf("USAGE: ncbench [OPTIONS] out.nc\n\n");
    printf("OPTIONS:\n\n");
    printf("nt           Number of records.\n");
    printf("nx           Number of x points.\n");
    printf("ny           Number of y points.\n");
    printf("nz           Number of levels.\n");
    printf("format       File format choice: 32,64,hdf\n");
    printf("initialsz    Initial file size in bytes before writing data.\n");
    printf("bufrsizehint Request I/O blocksize. 0 will use system default.\n");
    printf("h_minfree    Header pad.\n");
    printf("v_align      Static data section alignment.\n");
    printf("v_minfree    Pad after static data section.\n");
    printf("r_align      Record section alignment.\n");
    printf("prefill      Enable prefill.\n");
    printf("chunksize    HDF feature. Total size of chunk cache, in bytes.\n");
    printf("chunknelems  HDF feature. Number of chunk slots in chunk table.\n");
    printf("chunkpreemt  HDF feature. Between 0 and 1 inclusive.\n");
    
    return 1;
  }
  
  size_t nt = (size_t)atoi(argv[1]);
  size_t nx = (size_t)atoi(argv[2]);
  size_t ny = (size_t)atoi(argv[3]);
  size_t nz = (size_t)atoi(argv[4]);
  const char* format = argv[5];
  size_t initialsz = (size_t)strtoul(argv[6],0,0);
  size_t bufrsizehint = (size_t)atoi(argv[7]);
  size_t h_minfree = (size_t)atoi(argv[8]);
  size_t v_align = (size_t)atoi(argv[9]);
  size_t v_minfree = (size_t)atoi(argv[10]);
  size_t r_align = (size_t)atoi(argv[11]);
  bool prefill = (bool)atoi(argv[12]);
  size_t chunksize = (size_t)atoi(argv[13]);
  size_t chunknelems = (size_t)atoi(argv[14]);
  float chunkpreemption = (float)atof(argv[15]);
  const char * filename = argv[16];
  
  double vmMeta, rssMeta, vmStatic, rssStatic, vmRec, rssRec, vmMax=0, rssMax=0;
  double tic, toc, tMeta, tStatic, tRec, tTotal;
  
  ez::ezNc nc;
  nc.filename = filename;
  nc.nrec = nt;
  nc.initialsz = initialsz;
  nc.bufrsizehint = bufrsizehint;
  nc.h_minfree = h_minfree;
  nc.v_align = v_align;
  nc.v_minfree = v_minfree;
  nc.r_align = r_align;
  nc.fill = prefill;
  nc.chunksize = chunksize;
  nc.chunknelems = chunknelems;
  nc.chunkpreemption = chunkpreemption;

  tic = ez::osQueryPerformance();
  
  switch(format[0]) {
    case '6': nc.create64Clobber(); break;
    case 'h': nc.createHDFClobber(); break;
    default: nc.createClobber(); break;
  }
  
  // Define dims, vars, fillvalues.
  int xd = nc.createDim("x", nx);
  int yd = nc.createDim("y", ny);
  int zd = nc.createDim("z", nz);
  int td = nc.createDim("t", NC_UNLIMITED);
  int dimids3d[3] = {zd,yd,xd};
  int dimids3dt[4] = {td,zd,yd,xd};

  int v1 = nc.createVar("var1", NC_BYTE, 3, dimids3d);
  int v2 = nc.createVar("var2", NC_SHORT, 3, dimids3d);
  int v3 = nc.createVar("var3", NC_INT, 3, dimids3d);
  int v4 = nc.createVar("var4", NC_FLOAT, 3, dimids3d);
  int v5 = nc.createVar("var5", NC_DOUBLE, 3, dimids3d);
  int v6 = nc.createVar("var6", NC_BYTE, 4, dimids3dt);
  int v7 = nc.createVar("var7", NC_SHORT, 4, dimids3dt);
  int v8 = nc.createVar("var8", NC_INT, 4, dimids3dt);
  int v9 = nc.createVar("var9", NC_FLOAT, 4, dimids3dt);
  int v10 = nc.createVar("var10", NC_DOUBLE, 4, dimids3dt);
    
  unsigned char cfill = 1;
  short sfill = 1;
  int ifill = 1;
  float ffill = 1.0;
  double dfill = 1.0;
  nc_put_att_uchar(nc.ncid, v1, "_FillValue", NC_BYTE, 1, &cfill);
  nc_put_att_short(nc.ncid, v2, "_FillValue", NC_SHORT, 1, &sfill);
  nc_put_att_int(nc.ncid, v3, "_FillValue", NC_INT, 1, &ifill);
  nc_put_att_float(nc.ncid, v4, "_FillValue", NC_FLOAT, 1, &ffill);
  nc_put_att_double(nc.ncid, v5, "_FillValue", NC_DOUBLE, 1, &dfill);
  nc_put_att_uchar(nc.ncid, v6, "_FillValue", NC_BYTE, 1, &cfill);
  nc_put_att_short(nc.ncid, v7, "_FillValue", NC_SHORT, 1, &sfill);
  nc_put_att_int(nc.ncid, v8, "_FillValue", NC_INT, 1, &ifill);
  nc_put_att_float(nc.ncid, v9, "_FillValue", NC_FLOAT, 1, &ffill);
  nc_put_att_double(nc.ncid, v10, "_FillValue", NC_DOUBLE, 1, &dfill);
  
  nc.enddef();
  
  size_t bytesStatic = nx*ny*nz*(sizeof(unsigned char)+sizeof(short)+sizeof(int)+sizeof(float)+sizeof(double));
  size_t bytesRec = nt*bytesStatic;
  
  toc = ez::osQueryPerformance();
  tMeta = toc - tic;
  tic = toc;
  
  std::vector<std::string> strings;
  strings.push_back("var1");
  strings.push_back("var2");
  strings.push_back("var3");
  strings.push_back("var4");
  strings.push_back("var5");
  strings.push_back("var6");
  strings.push_back("var7");
  strings.push_back("var8");
  strings.push_back("var9");
  strings.push_back("var10");
  ez::ezNcBuffers buf;
  nc.init(buf, strings);
  buf.allocate();
  buf.zeros();

  // Write static var data.
  size_t start[4], count[4];
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  start[3] = 0;
  count[0] = nz;
  count[1] = ny;
  count[2] = nx;

  nc.write(v1, &buf, start, count);
  nc.write(v2, &buf, start, count);
  nc.write(v3, &buf, start, count);
  nc.write(v4, &buf, start, count);
  nc.write(v5, &buf, start, count);

  toc = ez::osQueryPerformance();
  tStatic = toc - tic;
  tic = toc;
  ez::process_mem_usage(vmMeta, rssMeta);
  vmMax = std::max(vmMeta, vmMax);
  rssMax = std::max(rssMeta, rssMax);
  
  // Write rec var data.
  count[0] = 1;
  count[1] = nz;
  count[2] = ny;
  count[3] = nx;

  int i=0;
  for(; i < nt; ++i) {
    start[0] = i;
    nc.write(v6, &buf, start, count);
    nc.write(v7, &buf, start, count);
    nc.write(v8, &buf, start, count);
    nc.write(v9, &buf, start, count);
    nc.write(v10, &buf, start, count);
  }
  
  ez::process_mem_usage(vmRec, rssRec);  
  vmMax = std::max(vmRec, vmMax);
  rssMax = std::max(rssRec, rssMax);
  
  nc.close();
  toc = ez::osQueryPerformance();
  tRec = toc - tic;
  tic = toc;
  tTotal = toc - ticTotal;
  
  printf("Total(s),Meta(s),Static(s),Rec(s),StaticGB,RecGB,MBps,VMMax,RssMax\n");
  printf("%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%g,%g,%g\n", tTotal/1e6, tMeta/1e6, tStatic/1e6, tRec/1e6, bytesStatic/1.0e9, bytesRec/1.0e9, bytesRec/((tRec/1.0e6)*1.0e6), vmMax, rssMax);
  
  
  return 0;
}