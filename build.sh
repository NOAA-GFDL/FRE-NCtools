autoreconf -i
./configure --prefix=$PREFIX
make -j
make install
