echo 'builindg FRE-NCtools conda package...'
echo "Build directory:" && pwd
echo "Contents of Build directory:" && ls
autoreconf -i
./configure --prefix=$PREFIX
make
make install
#cp configure.ac $SRC_DIR
#make check RUN_EXPENSIVE_TESTS=no
