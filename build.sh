autoreconf -i
./configure --prefix=$PREFIX --disable-ocean-model-grid-generator
make -j 2
make install
