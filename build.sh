autoreconf -i
./configure --prefix=$PREFIX --disable-ocean-model-grid-generator
make -j
make install
