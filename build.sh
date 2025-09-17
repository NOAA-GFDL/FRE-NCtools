autoreconf -i
./configure --prefix=$PREFIX
make
make install
#cp configure.ac $SRC_DIR
#make check RUN_EXPENSIVE_TESTS=no
