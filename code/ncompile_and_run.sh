#!/bin/sh
exeName=$1
if [ "${exeName##*.}"x = "cu"x ];then
    exeName=${exeName%.*}
fi
#clang-format -i $exeName.cu -style=LLVM
nvcc -std=c++11 -o $exeName $exeName.cu -Xcompiler -fopenmp
echo compile [ $exeName.cu ] finished!
echo -----------------------------------------------------begin ↓
./$exeName $2 $3 $4 $5 $6 $7 $8 $9
echo -----------------------------------------------------eng   ↑