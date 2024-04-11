#!/bin/bash

export PATH=$PATH:/home/Mikyung.Lee/FRE-NCTools/3_acc_read_map/fre-nctools/bin

nvc -g -gpu=debug,lineinfo -acc -std=c11 -I/home/Mikyung.Lee/FRE-NCTools/3_acc_read_map/tools/fregrid_acc \
    -I/home/Mikyung.Lee/FRE-NCTools/3_acc_read_map/tools/libfrencutils \
    /home/Mikyung.Lee/FRE-NCTools/3_acc_read_map/build/tools/fregrid_acc/conserve_interp_acc.o \
    /home/Mikyung.Lee/FRE-NCTools/3_acc_read_map/build/tools/libfrencutils/libfrencutils.a \
    /home/Mikyung.Lee/FRE-NCTools/3_acc_read_map/build/tools/libfrencutils_acc/libfrencutils_acc.a \
    `nc-config --cflags` `nc-config --libs` -Minfo=all \
    Test_read_remap_file.c
#./a.out 1
#./a.out 2
