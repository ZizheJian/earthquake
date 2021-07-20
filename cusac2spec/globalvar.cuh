#ifndef __GLOBALVAR_CUH__
#define __GLOBALVAR_CUH__

#include <cuda_runtime.h>
#include <cufft.h>
#include "sac.cuh"
#include "spec.cuh"

#define HD cudaMemcpyHostToDevice
#define DH cudaMemcpyDeviceToHost
#define DD cudaMemcpyDeviceToDevice
#define PI 3.14159265358979323846
#define soi sizeof(int)
#define sor sizeof(cufftReal)
#define soc sizeof(cufftComplex)
#define cuds cudaDeviceSynchronize()
#define all (nfft*step+maxth-1)/maxth,maxth

//开始和结束时间
clock_t s,t;
//统计处理的文件数和数据量(float)
int filenum,datanum;

//n一个文件内float数,winn=win/dt,lagn=lag/dt,step分段数
int maxn,maxstep,maxnspec,planflag;
int n,winn,lagn,step,nfft,nspec;
//dt时间精度,win窗口长度(s),lag滑动距离(s)
float dt,win,lag,lagrate,df;
float *sacdata;
float *specdata;
int *segsize;//一个分段中原始数据量(float)

//-m3
//m3wn=m3w/dt,m3in=m3i/df
int m3wn,m31n,m32n,m33n,m34n;
//m3w第一次滤波后,时域下平滑半径
//m3top第一次滤波后,时域下削平时的最大值
//m3i第一次滤波时关键的频率
float m3w,m3top,m31,m32,m33,m34;

int wf,wf1n,wf2n,wf3n,wf4n;
float wf1,wf2,wf3,wf4;

int cmp;//是否使用指定cmpaz,cmpinc的开关
float cmpaz,cmpinc;

int ns=30;//第二次平滑半径

SACHEAD sachd;
SPECHEAD spechd;//输入的sac路径
char sacrootpath[256];//输入的spec路径
char specrootpath[256];
char sactail[10]=".sacum\0";
char spectail[10]=".spec\0";
int specl;//specrootpath的长度

int maxth;//最大线程数
cudaDeviceProp gpuProp;
int *d_segsize;//显存中的segsize
cufftReal *d_sac_r1,*d_sac_r2,*d_sac_r3;//显存中大小step*nfft*sor的空间
cufftComplex *d_sac_c1;//显存中step*(nspec+1)*soc
cufftHandle planR2C,planC2R;
int rank;
int len[1];
int inembed[2];
int istride;
int idist;
int onembed[2];
int ostride;
int odist;
int batch;

#endif