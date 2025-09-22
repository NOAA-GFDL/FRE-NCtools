echo 'builindg FRE-NCtools conda package...'
echo "SRC_DIR / Build directory is:"
pwd
echo "Contents of SRC_DIR / Build directory are:"
ls

export LD_LIBRARY_PATH=${PREFIX}/lib
echo "PRE CONFIGURATION::"
echo ""
echo ""
echo "PATH is:"
echo $PATH
echo "LD_LIBRARY_PATH is:"
echo $LD_LIBRARY_PATH
echo ""

## this is sufficient
autoreconf -i
./configure --prefix=$PREFIX || cat config.log
#./configure --prefix=$PREFIX --with-mpi || cat config.log
#./configure --prefix=$PREFIX --enable-quad-precision || cat config.log
#./configure --prefix=$PREFIX --enable-quad-precision --with-mpi || cat config.log

export LD_LIBRARY_PATH=${PREFIX}/lib
echo "POST CONFIGURATION::"
echo ""
echo ""
echo "PATH is:"
echo $PATH
echo "LD_LIBRARY_PATH is:"
echo $LD_LIBRARY_PATH
echo ""

echo "compiling/building"
make
echo "installing into $PREFIX"
make install

### to test, build-dir option, ala README
#autoreconf -i
#mkdir build && cd build
#../configure --prefix=$PREFIX || cat config.log
##../configure --prefix=$PREFIX --with-mpi || cat config.log
##../configure --prefix=$PREFIX --enable-quad-precision || cat config.log
##../configure --prefix=$PREFIX --enable-quad-precision --with-mpi || cat config.logecho "compiling/building"
#make
#echo "installing into $PREFIX"
#make install
