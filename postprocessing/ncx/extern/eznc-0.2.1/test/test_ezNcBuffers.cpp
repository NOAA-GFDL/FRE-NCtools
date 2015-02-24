/*
20111111 rsz Created.
*/

#include "ezTest.hpp"
#include "ezNcBuffers.hpp"
#include "ezNcUtil.hpp"

bool test_allocate(ez::ezTestRunner& runner) {
  ez::ezNcBuffers buf;
  
  ezassert( buf.c == 0 );
  ezassert( buf.s == 0 );
  ezassert( buf.i == 0 );
  ezassert( buf.l == 0 );
  ezassert( buf.f == 0 );
  ezassert( buf.d == 0 );
  ezassert( buf.uc == 0 );
  ezassert( buf.us == 0 );
  ezassert( buf.ui == 0 );
  ezassert( buf.ul == 0 );
  
  buf.nc = 10;
  buf.ns = 11;
  buf.ni = 12;
  buf.nl = 13;
  buf.nf = 14;
  buf.nd = 15;
  buf.nuc = 16;
  buf.nus = 17;
  buf.nui = 18;
  buf.nul = 19;
  
  buf.allocate();
  
  ezassert( buf.c != 0 );
  ezassert( buf.s != 0 );
  ezassert( buf.i != 0 );
  ezassert( buf.l != 0 );
  ezassert( buf.f != 0 );
  ezassert( buf.d != 0 );
  ezassert( buf.uc != 0 );
  ezassert( buf.us != 0 );
  ezassert( buf.ui != 0 );
  ezassert( buf.ul != 0 );
  
  return true;
}

bool test_fill(ez::ezTestRunner& runner) {
  ez::ezNcBuffers buf;
  
  buf.nc = 10;
  buf.ns = 11;
  buf.ni = 12;
  buf.nl = 13;
  buf.nf = 14;
  buf.nd = 15;
  buf.nuc = 16;
  buf.nus = 17;
  buf.nui = 18;
  buf.nul = 19;
  
  buf.allocate();
  buf.fill(100);
  
  ezassert( buf.c[buf.nc-1] == 100 );
  ezassert( buf.s[buf.ns-1] == 100 );
  ezassert( buf.i[buf.ni-1] == 100 );
  ezassert( buf.l[buf.nl-1] == 100 );
  ezassert( buf.f[buf.nf-1] == 100 );
  ezassert( buf.d[buf.nd-1] == 100 );
  ezassert( buf.uc[buf.nuc-1] == 100 );
  ezassert( buf.us[buf.nus-1] == 100 );
  ezassert( buf.ui[buf.nui-1] == 100 );
  ezassert( buf.ul[buf.nul-1] == 100 );
  
  return true;
}

bool test_equals(ez::ezTestRunner& runner) {
  ez::ezNcBuffers buf1, buf2;

  ezassert( buf1.equals(buf2) );
    
  buf1.nc = 10;
  buf1.ns = 11;
  buf1.ni = 12;
  buf1.nl = 13;
  buf1.nf = 14;
  buf1.nd = 15;
  buf1.nuc = 16;
  buf1.nus = 17;
  buf1.nui = 18;
  buf1.nul = 19;
  
  buf1.allocate();
  buf1.fill(100);

  ezassert( buf1.equals(buf2) == false );
  
  buf2.nc = 110;
  buf2.ns = 111;
  buf2.ni = 112;
  buf2.nl = 113;
  buf2.nf = 114;
  buf2.nd = 115;
  buf2.nuc = 116;
  buf2.nus = 117;
  buf2.nui = 118;
  buf2.nul = 119;
  
  buf2.allocate();
  buf2.fill(100);

  ezassert( buf1.equals(buf2, NC_BYTE) == false );
  ezassert( buf1.equals(buf2, NC_CHAR) == false );
  ezassert( buf1.equals(buf2, NC_SHORT) == false );
  ezassert( buf1.equals(buf2, NC_INT) == false );
  ezassert( buf1.equals(buf2, NC_INT64) == false );
  ezassert( buf1.equals(buf2, NC_UINT) == false );  
  ezassert( buf1.equals(buf2, NC_UINT64) == false );
  ezassert( buf1.equals(buf2, NC_USHORT) == false );
  ezassert( buf1.equals(buf2, NC_FLOAT) == false );
  ezassert( buf1.equals(buf2, NC_DOUBLE) == false );
  ezassert( buf1.equals(buf2) == false );
  
  buf2.reset();
  buf2.nc = 10;
  buf2.ns = 11;
  buf2.ni = 12;
  buf2.nl = 13;
  buf2.nf = 14;
  buf2.nd = 15;
  buf2.nuc = 16;
  buf2.nus = 17;
  buf2.nui = 18;
  buf2.nul = 19;
  
  buf2.allocate();
  buf2.fill(0);
  
  ezassert( buf1.equals(buf2, NC_BYTE) == false );
  ezassert( buf1.equals(buf2, NC_CHAR) == false );
  ezassert( buf1.equals(buf2, NC_SHORT) == false );
  ezassert( buf1.equals(buf2, NC_INT) == false );
  ezassert( buf1.equals(buf2, NC_INT64) == false );
  ezassert( buf1.equals(buf2, NC_UINT) == false );  
  ezassert( buf1.equals(buf2, NC_UINT64) == false );
  ezassert( buf1.equals(buf2, NC_USHORT) == false );
  ezassert( buf1.equals(buf2, NC_FLOAT) == false );
  ezassert( buf1.equals(buf2, NC_DOUBLE) == false );
  ezassert( buf1.equals(buf2) == false );

  buf2.reset();
  buf2.nc = 10;
  buf2.ns = 11;
  buf2.ni = 12;
  buf2.nl = 13;
  buf2.nf = 14;
  buf2.nd = 15;
  buf2.nuc = 16;
  buf2.nus = 17;
  buf2.nui = 18;
  buf2.nul = 19;
  
  buf2.allocate();
  buf2.fill(100);
  
  ezassert( buf1.equals(buf2, NC_BYTE) );
  ezassert( buf1.equals(buf2, NC_CHAR) );
  ezassert( buf1.equals(buf2, NC_SHORT) );
  ezassert( buf1.equals(buf2, NC_INT) );
  ezassert( buf1.equals(buf2, NC_INT64) );
  ezassert( buf1.equals(buf2, NC_UINT) );  
  ezassert( buf1.equals(buf2, NC_UINT64) );
  ezassert( buf1.equals(buf2, NC_USHORT) );
  ezassert( buf1.equals(buf2, NC_FLOAT) );
  ezassert( buf1.equals(buf2, NC_DOUBLE) );
  ezassert( buf1.equals(buf2) );

  return true;  
}

bool test_find(ez::ezTestRunner& runner) {
  ez::ezNcBuffers buf;
  
  buf.nc = 10;
  buf.ns = 11;
  buf.ni = 12;
  buf.nl = 13;
  buf.nf = 14;
  buf.nd = 15;
  buf.nuc = 16;
  buf.nus = 17;
  buf.nui = 18;
  buf.nul = 19;
  
  buf.allocate();
  buf.fill(100);

  buf.c[5] = 199;
  buf.c[buf.nc-1] = 99;
  buf.uc[5] = 199;
  buf.uc[buf.nuc-1] = 99;
  buf.s[5] = 199;
  buf.s[buf.ns-1] = 99;
  buf.us[5] = 199;
  buf.us[buf.nus-1] = 99;
  buf.i[5] = 199;
  buf.i[buf.ni-1] = 99;
  buf.ui[5] = 199;
  buf.ui[buf.nui-1] = 99;
  buf.l[5] = 199;
  buf.l[buf.nl-1] = 99;
  buf.ul[5] = 199;
  buf.ul[buf.nul-1] = 99;
  buf.f[5] = 199;
  buf.f[buf.nf-1] = 99;
  buf.d[5] = 199;
  buf.d[buf.nd-1] = 99;

  std::vector<char> vc;
  std::vector<unsigned char> vuc;
  std::vector<short> vs;
  std::vector<unsigned short> vus;
  std::vector<int> vi;
  std::vector<unsigned int> vui;
  std::vector<long long> vl;
  std::vector<unsigned long long> vul;
  std::vector<float> vf;
  std::vector<double> vd;
  
  vc.push_back(99);
  vc.push_back(199);
  vuc.push_back(99);
  vuc.push_back(199);
  vs.push_back(99);
  vs.push_back(199);
  vus.push_back(99);
  vus.push_back(199);
  vi.push_back(99);
  vi.push_back(199);
  vui.push_back(99);
  vui.push_back(199);
  vl.push_back(99);
  vl.push_back(199);
  vul.push_back(99);
  vul.push_back(199);
  vf.push_back(99);
  vf.push_back(199);
  vd.push_back(99);
  vd.push_back(199);

  std::vector<int> idx;
  
  buf.find<unsigned char, int>(NC_BYTE, vuc, idx);
  ezassert( idx.size() == 2 );
  ezassert( idx[0] == (buf.nuc-1) );
  ezassert( idx[1] == 5 );

  buf.find<char, int>(NC_CHAR, vc, idx);
  ezassert( idx.size() == 2 );
  ezassert( idx[0] == (buf.nc-1) );
  ezassert( idx[1] == 5 );

  buf.find<short, int>(NC_SHORT, vs, idx);
  ezassert( idx.size() == 2 );
  ezassert( idx[0] == (buf.ns-1) );
  ezassert( idx[1] == 5 );

  buf.find<unsigned short, int>(NC_USHORT, vus, idx);
  ezassert( idx.size() == 2 );
  ezassert( idx[0] == (buf.nus-1) );
  ezassert( idx[1] == 5 );

  buf.find<int, int>(NC_INT, vi, idx);
  ezassert( idx.size() == 2 );
  ezassert( idx[0] == (buf.ni-1) );
  ezassert( idx[1] == 5 );

  buf.find<unsigned int, int>(NC_UINT, vui, idx);
  ezassert( idx.size() == 2 );
  ezassert( idx[0] == (buf.nui-1) );
  ezassert( idx[1] == 5 );

  buf.find<long long, int>(NC_INT64, vl, idx);
  ezassert( idx.size() == 2 );
  ezassert( idx[0] == (buf.nl-1) );
  ezassert( idx[1] == 5 );

  buf.find<unsigned long long, int>(NC_UINT64, vul, idx);
  ezassert( idx.size() == 2 );
  ezassert( idx[0] == (buf.nul-1) );
  ezassert( idx[1] == 5 );

  buf.find<float, int>(NC_FLOAT, vf, idx);
  ezassert( idx.size() == 2 );
  ezassert( idx[0] == (buf.nf-1) );
  ezassert( idx[1] == 5 );

  buf.find<double, int>(NC_DOUBLE, vd, idx);
  ezassert( idx.size() == 2 );
  ezassert( idx[0] == (buf.nd-1) );
  ezassert( idx[1] == 5 );
  
  return true;
}

bool test_findFirstGreaterEqual(ez::ezTestRunner& runner) {
  ez::ezNcBuffers buf;
  int i;
  
  buf.nc = 10;
  buf.ns = 11;
  buf.ni = 12;
  buf.nl = 13;
  buf.nf = 14;
  buf.nd = 15;
  buf.nuc = 16;
  buf.nus = 17;
  buf.nui = 18;
  buf.nul = 19;
  
  buf.allocate();
  ez::ezNcVariant v;
  
  v.type = NC_CHAR;
  v.set<char>(2);
  ez::seq<char>(buf.c, 10);
  i = buf.findFirstGreaterEqual(v);
  ezassert( i == 2 );

  v.type = NC_SHORT;
  v.set<short>(3);
  ez::seq<short>(buf.s, 11);
  i = buf.findFirstGreaterEqual(v);
  ezassert( i == 3 );

  v.type = NC_INT;
  v.set<int>(4);
  ez::seq<int>(buf.i, 12);
  i = buf.findFirstGreaterEqual(v);
  ezassert( i == 4 );

  v.type = NC_INT64;
  v.set<long long>(5);
  ez::seq<long long>(buf.l, 13);
  i = buf.findFirstGreaterEqual(v);
  ezassert( i == 5 );

  v.type = NC_FLOAT;
  v.set<float>(6);
  ez::seq<float>(buf.f, 14);
  i = buf.findFirstGreaterEqual(v);
  ezassert( i == 6 );

  v.type = NC_DOUBLE;
  v.set<double>(7);
  ez::seq<double>(buf.d, 15);
  i = buf.findFirstGreaterEqual(v);
  ezassert( i == 7 );

  v.type = NC_UBYTE;
  v.set<unsigned char>(8);
  ez::seq<unsigned char>(buf.uc, 16);
  i = buf.findFirstGreaterEqual(v);
  ezassert( i == 8 );

  v.type = NC_USHORT;
  v.set<unsigned short>(9);
  ez::seq<unsigned short>(buf.us, 17);
  i = buf.findFirstGreaterEqual(v);
  ezassert( i == 9 );

  v.type = NC_UINT;
  v.set<unsigned int>(10);
  ez::seq<unsigned int>(buf.ui, 18);
  i = buf.findFirstGreaterEqual(v);
  ezassert( i == 10 );

  v.type = NC_UINT64;
  v.set<unsigned long long>(11);
  ez::seq<unsigned long long>(buf.ul, 19);
  i = buf.findFirstGreaterEqual(v);
  ezassert( i == 11 );

  return true;
}

bool test_findFirstLessEqual(ez::ezTestRunner& runner) {
  ez::ezNcBuffers buf;
  int i;
  
  buf.nc = 10;
  buf.ns = 11;
  buf.ni = 12;
  buf.nl = 13;
  buf.nf = 14;
  buf.nd = 15;
  buf.nuc = 16;
  buf.nus = 17;
  buf.nui = 18;
  buf.nul = 19;
  
  buf.allocate();
  ez::ezNcVariant v;
  
  v.type = NC_CHAR;
  v.set<char>(2);
  ez::seq<char>(buf.c, 10, 1);
  i = buf.findFirstLessEqual(v);
  ezassert( i == 8 );

  v.type = NC_SHORT;
  v.set<short>(3);
  ez::seq<short>(buf.s, 11, 1);
  i = buf.findFirstLessEqual(v);
  ezassert( i == 8 );

  v.type = NC_INT;
  v.set<int>(4);
  ez::seq<int>(buf.i, 12, 1);
  i = buf.findFirstLessEqual(v);
  ezassert( i == 8 );

  v.type = NC_INT64;
  v.set<long long>(5);
  ez::seq<long long>(buf.l, 13, 1);
  i = buf.findFirstLessEqual(v);
  ezassert( i == 8 );

  v.type = NC_FLOAT;
  v.set<float>(6);
  ez::seq<float>(buf.f, 14, 1);
  i = buf.findFirstLessEqual(v);
  ezassert( i == 8 );

  v.type = NC_DOUBLE;
  v.set<double>(7);
  ez::seq<double>(buf.d, 15, 1);
  i = buf.findFirstLessEqual(v);
  ezassert( i == 8 );

  v.type = NC_UBYTE;
  v.set<unsigned char>(8);
  ez::seq<unsigned char>(buf.uc, 16, 1);
  i = buf.findFirstLessEqual(v);
  ezassert( i == 8 );

  v.type = NC_USHORT;
  v.set<unsigned short>(9);
  ez::seq<unsigned short>(buf.us, 17, 1);
  i = buf.findFirstLessEqual(v);
  ezassert( i == 8 );

  v.type = NC_UINT;
  v.set<unsigned int>(10);
  ez::seq<unsigned int>(buf.ui, 18, 1);
  i = buf.findFirstLessEqual(v);
  ezassert( i == 8 );

  v.type = NC_UINT64;
  v.set<unsigned long long>(11);
  ez::seq<unsigned long long>(buf.ul, 19, 1);
  i = buf.findFirstLessEqual(v);
  ezassert( i == 8 );

  return true;
}

bool test_findLastGreaterEqual(ez::ezTestRunner& runner) {
  ez::ezNcBuffers buf;
  int i;
  
  buf.nc = 10;
  buf.ns = 11;
  buf.ni = 12;
  buf.nl = 13;
  buf.nf = 14;
  buf.nd = 15;
  buf.nuc = 16;
  buf.nus = 17;
  buf.nui = 18;
  buf.nul = 19;
  
  buf.allocate();
  ez::ezNcVariant v;
  
  v.type = NC_CHAR;
  v.set<char>(2);
  ez::seq<char>(buf.c, 10, 1);
  i = buf.findLastGreaterEqual(v);
  ezassert( i == 8 );

  v.type = NC_SHORT;
  v.set<short>(3);
  ez::seq<short>(buf.s, 11, 1);
  i = buf.findLastGreaterEqual(v);
  ezassert( i == 8 );

  v.type = NC_INT;
  v.set<int>(4);
  ez::seq<int>(buf.i, 12, 1);
  i = buf.findLastGreaterEqual(v);
  ezassert( i == 8 );

  v.type = NC_INT64;
  v.set<long long>(5);
  ez::seq<long long>(buf.l, 13, 1);
  i = buf.findLastGreaterEqual(v);
  ezassert( i == 8 );

  v.type = NC_FLOAT;
  v.set<float>(6);
  ez::seq<float>(buf.f, 14, 1);
  i = buf.findLastGreaterEqual(v);
  ezassert( i == 8 );

  v.type = NC_DOUBLE;
  v.set<double>(7);
  ez::seq<double>(buf.d, 15, 1);
  i = buf.findLastGreaterEqual(v);
  ezassert( i == 8 );

  v.type = NC_UBYTE;
  v.set<unsigned char>(8);
  ez::seq<unsigned char>(buf.uc, 16, 1);
  i = buf.findLastGreaterEqual(v);
  ezassert( i == 8 );

  v.type = NC_USHORT;
  v.set<unsigned short>(9);
  ez::seq<unsigned short>(buf.us, 17, 1);
  i = buf.findLastGreaterEqual(v);
  ezassert( i == 8 );

  v.type = NC_UINT;
  v.set<unsigned int>(10);
  ez::seq<unsigned int>(buf.ui, 18, 1);
  i = buf.findLastGreaterEqual(v);
  ezassert( i == 8 );

  v.type = NC_UINT64;
  v.set<unsigned long long>(11);
  ez::seq<unsigned long long>(buf.ul, 19, 1);
  i = buf.findLastGreaterEqual(v);
  ezassert( i == 8 );

  return true;
}

bool test_findLastLessEqual(ez::ezTestRunner& runner) {
  ez::ezNcBuffers buf;
  int i;
  
  buf.nc = 10;
  buf.ns = 11;
  buf.ni = 12;
  buf.nl = 13;
  buf.nf = 14;
  buf.nd = 15;
  buf.nuc = 16;
  buf.nus = 17;
  buf.nui = 18;
  buf.nul = 19;
  
  buf.allocate();
  ez::ezNcVariant v;
  
  v.type = NC_CHAR;
  v.set<char>(2);
  ez::seq<char>(buf.c, 10);
  i = buf.findLastLessEqual(v);
  ezassert( i == 2 );

  v.type = NC_SHORT;
  v.set<short>(3);
  ez::seq<short>(buf.s, 11);
  i = buf.findLastLessEqual(v);
  ezassert( i == 3 );

  v.type = NC_INT;
  v.set<int>(4);
  ez::seq<int>(buf.i, 12);
  i = buf.findLastLessEqual(v);
  ezassert( i == 4 );

  v.type = NC_INT64;
  v.set<long long>(5);
  ez::seq<long long>(buf.l, 13);
  i = buf.findLastLessEqual(v);
  ezassert( i == 5 );

  v.type = NC_FLOAT;
  v.set<float>(6);
  ez::seq<float>(buf.f, 14);
  i = buf.findLastLessEqual(v);
  ezassert( i == 6 );

  v.type = NC_DOUBLE;
  v.set<double>(7);
  ez::seq<double>(buf.d, 15);
  i = buf.findLastLessEqual(v);
  ezassert( i == 7 );

  v.type = NC_UBYTE;
  v.set<unsigned char>(8);
  ez::seq<unsigned char>(buf.uc, 16);
  i = buf.findLastLessEqual(v);
  ezassert( i == 8 );

  v.type = NC_USHORT;
  v.set<unsigned short>(9);
  ez::seq<unsigned short>(buf.us, 17);
  i = buf.findLastLessEqual(v);
  ezassert( i == 9 );

  v.type = NC_UINT;
  v.set<unsigned int>(10);
  ez::seq<unsigned int>(buf.ui, 18);
  i = buf.findLastLessEqual(v);
  ezassert( i == 10 );

  v.type = NC_UINT64;
  v.set<unsigned long long>(11);
  ez::seq<unsigned long long>(buf.ul, 19);
  i = buf.findLastLessEqual(v);
  ezassert( i == 11 );

  return true;
}

bool test_increasing(ez::ezTestRunner& runner) {
  ez::ezNcBuffers buf;
  
  buf.nc = 10;
  buf.ns = 11;
  buf.ni = 12;
  buf.nl = 13;
  buf.nf = 14;
  buf.nd = 15;
  buf.nuc = 16;
  buf.nus = 17;
  buf.nui = 18;
  buf.nul = 19;
  
  buf.allocate();
  
  ez::seq<char>(buf.c, buf.nc);
  ezassert( buf.increasing(NC_CHAR) );
  buf.c[5] = 0;
  ezassert( buf.increasing(NC_CHAR, false) );
  ezassert( !buf.increasing(NC_CHAR) );
  
  ez::seq<short>(buf.s, buf.ns);
  ezassert( buf.increasing(NC_SHORT) );
  buf.s[5] = 0;
  ezassert( buf.increasing(NC_SHORT, false) );
  ezassert( !buf.increasing(NC_SHORT) );

  ez::seq<int>(buf.i, buf.ni);
  ezassert( buf.increasing(NC_INT) );
  buf.i[5] = 0;
  ezassert( buf.increasing(NC_INT, false) );
  ezassert( !buf.increasing(NC_INT) );

  ez::seq<long long>(buf.l, buf.nl);
  ezassert( buf.increasing(NC_INT64) );
  buf.l[5] = 0;
  ezassert( buf.increasing(NC_INT64, false) );
  ezassert( !buf.increasing(NC_INT64) );

  ez::seq<float>(buf.f, buf.nf);
  ezassert( buf.increasing(NC_FLOAT) );
  buf.f[5] = 0;
  ezassert( buf.increasing(NC_FLOAT, false) );
  ezassert( !buf.increasing(NC_FLOAT) );

  ez::seq<double>(buf.d, buf.nd);
  ezassert( buf.increasing(NC_DOUBLE) );
  buf.d[5] = 0;
  ezassert( buf.increasing(NC_DOUBLE, false) );
  ezassert( !buf.increasing(NC_DOUBLE) );

  ez::seq<unsigned char>(buf.uc, buf.nuc);
  ezassert( buf.increasing(NC_UBYTE) );
  buf.uc[5] = 0;
  ezassert( buf.increasing(NC_UBYTE, false) );
  ezassert( !buf.increasing(NC_UBYTE) );

  ez::seq<unsigned short>(buf.us, buf.nus);
  ezassert( buf.increasing(NC_USHORT) );
  buf.us[5] = 0;
  ezassert( buf.increasing(NC_USHORT, false) );
  ezassert( !buf.increasing(NC_USHORT) );

  ez::seq<unsigned int>(buf.ui, buf.nui);
  ezassert( buf.increasing(NC_UINT) );
  buf.ui[5] = 0;
  ezassert( buf.increasing(NC_UINT, false) );
  ezassert( !buf.increasing(NC_UINT) );

  ez::seq<unsigned long long>(buf.ul, buf.nul);
  ezassert( buf.increasing(NC_UINT64) );
  buf.ul[5] = 0;
  ezassert( buf.increasing(NC_UINT64, false) );
  ezassert( !buf.increasing(NC_UINT64) );

  return true;
}

bool test_reset(ez::ezTestRunner& runner) {
  ez::ezNcBuffers buf;
  
  buf.nc = 10;
  buf.ns = 11;
  buf.ni = 12;
  buf.nl = 13;
  buf.nf = 14;
  buf.nd = 15;
  buf.nuc = 16;
  buf.nus = 17;
  buf.nui = 18;
  buf.nul = 19;
  
  buf.allocate();
  buf.reset();
  
  ezassert( buf.c == 0 );
  ezassert( buf.nc == 0 );
  ezassert( buf.uc == 0 );
  ezassert( buf.nuc == 0 );
  ezassert( buf.s == 0 );
  ezassert( buf.ns == 0 );
  ezassert( buf.us == 0 );
  ezassert( buf.nus == 0 );
  ezassert( buf.i == 0 );
  ezassert( buf.ni == 0 );
  ezassert( buf.ui == 0 );
  ezassert( buf.nui == 0 );
  ezassert( buf.l == 0 );
  ezassert( buf.nl == 0 );
  ezassert( buf.ul == 0 );
  ezassert( buf.nul == 0 );
  ezassert( buf.f == 0 );
  ezassert( buf.nf == 0 );
  ezassert( buf.d == 0 );
  ezassert( buf.nd == 0 );
  
  return true;
}

bool test_setMaxSizePerRec(ez::ezTestRunner& runner) {
  ez::ezNcBuffers buf;
  std::vector<int> lens;
  std::vector<int> ids;
  int recid = 0;
  
  lens.push_back(10);
  lens.push_back(20);
  lens.push_back(30);
  lens.push_back(40);
  ids.push_back(1);
  ids.push_back(2);
  ids.push_back(3);
  ids.push_back(4);
  
  buf.setMaxSizePerRec(NC_CHAR, lens, ids, recid);
  ezassert( buf.nc == 10*20*30*40 );
  
  lens[0] = 1;
  buf.setMaxSizePerRec(NC_CHAR, lens, ids, recid);
  ezassert( buf.nc == 10*20*30*40 ); // Should be unchanged from last set.
  ezassert( buf.nuc == 0 );
  ezassert( buf.ns == 0 );
  ezassert( buf.nus == 0 );
  ezassert( buf.ni == 0 );
  ezassert( buf.nui == 0 );
  ezassert( buf.nl == 0 );
  ezassert( buf.nul == 0 );
  ezassert( buf.nf == 0 );
  ezassert( buf.nd == 0 );

  ids[0] = recid;
  lens[0] = 10;
  buf.setMaxSizePerRec(NC_INT, lens, ids, recid);
  ezassert( buf.ni == 20*30*40 );
  
  return true;
}

bool test_setSize(ez::ezTestRunner& runner) {
  ez::ezNcBuffers buf;
  
  ezassert( buf.nc == 0 );
  ezassert( buf.nuc == 0 );
  ezassert( buf.ns == 0 );
  ezassert( buf.nus == 0 );
  ezassert( buf.ni == 0 );
  ezassert( buf.nui == 0 );
  ezassert( buf.nl == 0 );
  ezassert( buf.nul == 0 );
  ezassert( buf.nf == 0 );
  ezassert( buf.nd == 0 );
  
  int nuc=10;
  int nc=11;
  int ns=12;
  int nus=13;
  int ni=14;
  int nui=15;
  int nl=16;
  int nul=17;
  int nf=18;
  int nd=19;
  
  buf.setSize(NC_BYTE, nuc);
  buf.setSize(NC_CHAR, nc);
  buf.setSize(NC_SHORT, ns);
  buf.setSize(NC_USHORT, nus);
  buf.setSize(NC_INT, ni);
  buf.setSize(NC_UINT, nui);
  buf.setSize(NC_INT64, nl);
  buf.setSize(NC_UINT64, nul);
  buf.setSize(NC_FLOAT, nf);
  buf.setSize(NC_DOUBLE, nd);
  
  ezassert( buf.nc == nc );
  ezassert( buf.nuc == nuc );
  ezassert( buf.ns == ns );
  ezassert( buf.nus == nus );
  ezassert( buf.ni == ni );
  ezassert( buf.nui == nui );
  ezassert( buf.nl == nl );
  ezassert( buf.nul == nul );
  ezassert( buf.nf == nf );
  ezassert( buf.nd == nd );

  return true;
}

bool test_zeros(ez::ezTestRunner& runner) {
  ez::ezNcBuffers buf;
  
  buf.nc = 10;
  buf.ns = 11;
  buf.ni = 12;
  buf.nl = 13;
  buf.nf = 14;
  buf.nd = 15;
  buf.nuc = 16;
  buf.nus = 17;
  buf.nui = 18;
  buf.nul = 19;
  
  buf.allocate();
  buf.fill(100);
  
  ezassert( buf.c[buf.nc-1] == 100 );
  ezassert( buf.s[buf.ns-1] == 100 );
  ezassert( buf.i[buf.ni-1] == 100 );
  ezassert( buf.l[buf.nl-1] == 100 );
  ezassert( buf.f[buf.nf-1] == 100 );
  ezassert( buf.d[buf.nd-1] == 100 );
  ezassert( buf.uc[buf.nuc-1] == 100 );
  ezassert( buf.us[buf.nus-1] == 100 );
  ezassert( buf.ui[buf.nui-1] == 100 );
  ezassert( buf.ul[buf.nul-1] == 100 );
  
  buf.zeros();
  
  ezassert( buf.c[buf.nc-1] == 0 );
  ezassert( buf.s[buf.ns-1] == 0 );
  ezassert( buf.i[buf.ni-1] == 0 );
  ezassert( buf.l[buf.nl-1] == 0 );
  ezassert( buf.f[buf.nf-1] == 0 );
  ezassert( buf.d[buf.nd-1] == 0 );
  ezassert( buf.uc[buf.nuc-1] == 0 );
  ezassert( buf.us[buf.nus-1] == 0 );
  ezassert( buf.ui[buf.nui-1] == 0 );
  ezassert( buf.ul[buf.nul-1] == 0 );

  return true;
}

int main(int argc, const char* argv[]) {
  ez::ezTestRunner runner;
  runner.getopt(argc, argv);
  #define TEST(N) runner.add(&test_##N, "test_"#N);
  TEST(allocate);
  TEST(equals);
  TEST(fill);
  TEST(find);
  TEST(findFirstGreaterEqual);
  TEST(findFirstLessEqual);
  TEST(findLastGreaterEqual);
  TEST(findLastLessEqual);
  TEST(increasing);
  TEST(reset);
  TEST(setMaxSizePerRec);
  TEST(setSize);
  TEST(zeros);
  
  return runner.run();
}