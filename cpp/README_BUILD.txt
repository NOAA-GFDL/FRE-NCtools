#   I) load desired compilers and libs

#   IIA) #Compiling with CMake presets:
 # inspect the preset files under dir  FRE-NCtools/cpp
 # if you dont find one to you liking, modify one or make a
 # new one
 # List the available presets:
 cd FRE-NCtools/cpp
 cmake --list-presets
 #choose one of the presets, e.g. gpu-nvhpc-pr :
 cmake . --preset=gpu-nvhpc-pr -Dwith-gpu=ON
 #cd to the directory listed in the last output line - its
 #  where the build files are written.
 cd build_gpu_nvhpc_pr
 # run cmake in verbose mode
 cmake --build . -v
 # note the executables are currently copied to dir FRE-NCtools
-----------------
 #  IIB) Compiling using env vars and more cmake command line options.
 #set the compiler
 export CXX=/opt/nvidia/hpc_sdk/Linux_x86_64/23.7/compilers/bin/nvc++
 #pass to cmake the flags you want.
 cmake .. -Dwith-gpu="ON" -DCMAKE_CXX_FLAGS="-O2 -std=c++20 -stdpar=gpu"
 cmake  --build . -v
