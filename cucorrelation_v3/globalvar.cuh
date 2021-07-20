#ifndef __GLOBALVAR_CUH__
#define __GLOBALVAR_CUH__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "spec.cuh"

#define HD cudaMemcpyHostToDevice
#define DH cudaMemcpyDeviceToHost
#define DD cudaMemcpyDeviceToDevice
#define soi 4
#define sor 4
#define soc 8
#define sop 8
#define sosph 88
#define sospi 344
#define cuds cudaDeviceSynchronize()

typedef struct{
    int x;
    int y;
}intx2;

////////////////io////////////////
char specpath[256],sacpath[256],tpath[256];
int sacl;
specinf *speclist;
int specfilenum;
char sactail[10]=".sacum\0";
char spectail[10]=".spec\0";

////////////////显存管理////////////////
float gpuuse=0.5;
int segmax,specmax,chimax,chinum;
intx2 *chin;
int a,n_a;
size_t gpusize,totgpusize;

////////////////内存/显存地址////////////////
cufftComplex *c1,*d_c1,*d_c2;
cufftReal *r1,*d_r1;
intx2 *mlist,*d_mlist;
specinf *d_speclist;

////////////////cuda////////////////
cudaDeviceProp gpuProp;
int maxth,blocknum;
int len[1];
int inembed[2],onembed[2];
cufftHandle plan;

////////////////cc////////////////
float cchalf;

////////////////setchn////////////////
int setchn;
char ncfchnname[8];

////////////////对称////////////////
int symmetric;

////////////////对称////////////////
int debug=0;

////////////////时间////////////////
clock_t s,t;

#endif