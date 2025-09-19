echo 'builindg FRE-NCtools conda package...'
echo "SRC_DIR / Build directory is:"
pwd
echo "Contents of SRC_DIR / Build directory are:"
ls

### this is sufficient
#autoreconf -i
#./configure --prefix=$PREFIX || cat config.log
#./configure --prefix=$PREFIX --with-mpi || cat config.log
#./configure --prefix=$PREFIX --enable-quad-precision || cat config.log
#./configure --prefix=$PREFIX --enable-quad-precision --with-mpi || cat config.log
#echo "compiling/building"
#make
#echo "installing into $PREFIX"
#make install

## to test, build-dir option, ala README
autoreconf -i
mkdir build && cd build
../configure --prefix=$PREFIX || cat config.log
#../configure --prefix=$PREFIX --with-mpi || cat config.log
#../configure --prefix=$PREFIX --enable-quad-precision || cat config.log
#../configure --prefix=$PREFIX --enable-quad-precision --with-mpi || cat config.logecho "compiling/building"
echo "compiling/building at $PWD"
make
echo "installing into $PREFIX"
make install
