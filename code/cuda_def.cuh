#ifndef CUDA_DEF_CUH
#define CUDA_DEF_CUH

#include <stdio.h>

#define CUERR                                                                                                \
    {                                                                                                        \
        cudaError_t err;                                                                                     \
        if ((err = cudaGetLastError()) != cudaSuccess) {                                                     \
            printf("CUDA error: %s : %s, line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);           \
            exit(1);                                                                                         \
        }                                                                                                    \
    }

#define STARTTIME                                                                                            \
    cudaEvent_t Mstart, Mstop;                                                                               \
    float Mtime;                                                                                             \
    cudaEventCreate(&Mstart);                                                                                \
    cudaEventCreate(&Mstop);                                                                                 \
    cudaEventRecord(Mstart, 0);

#define STOPTIME                                                                                             \
    cudaEventRecord(Mstop, 0);                                                                               \
    cudaEventSynchronize(Mstop);                                                                             \
    cudaEventElapsedTime(&Mtime, Mstart, Mstop);                                                             \
    printf("gpu time=%f sec\n", Mtime / 1000.0);

#define CLOCKSTART time_t startClock = clock();

#define CLOCKSTOP                                                                                            \
    time_t stopClock = clock();                                                                              \
    printf("total time use=%f sec\n", (stopClock - startClock) / double(CLOCKS_PER_SEC));

#endif
