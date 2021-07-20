#ifndef __FFT_CU__
#define __FFT_CU__

#include "globalvar.cuh"

__global__ void cuMuli(float *data,int n)
{
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    if (id>=n)
        return;
    data[id]=data[id]*id;
}

__global__ void cuSum(float *data,int n,int gap)
{
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    if (id>=n)
        return;
    if (id%gap==0)
    {
        if (id+gap/2<n)
        {
            data[id]+=data[id+gap/2];
        }
    }
}

__global__ void currr(float *data,int n,float a,float b)
{
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    if (id>=n)
        return;
    data[id]-=(a*id+b);
}

void sumcontrol(float *data,int n,int muli)//muli为是否要乘上i的开关
{
    if (muli)
    {
        cuMuli<<<(n+maxth-1)/maxth,maxth>>>(data,n);
        cuds;
    }
    for (int i=2;i<=2*n;i*=2)
    {
        cuSum<<<(n+maxth-1)/maxth,maxth>>>(data,n,i);
        cuds;
    }
}

void rrrcontrol(float *data,int n,float* buffer)
{
    //buffer为给这个函数分配的一个寄存器,输入输出都是data
    float sum,summuli;//sum=sigma(data[i]),summuli=sigma(data[i]*i)
    cudaMemcpy(buffer,data,n*sor,DD);
    sumcontrol(buffer,n,0);
    cudaMemcpy(&sum,buffer,sor,DH);
    cudaMemcpy(buffer,data,n*sor,DD);
    sumcontrol(buffer,n,1);
    cudaMemcpy(&summuli,buffer,sor,DH);
    float n2=(n-1)*n*(2*n-1)/6;
    float n1=(n-1)*n/2;
    float n0=n;
    float a=(summuli*n0-sum*n1)/(n0*n2-n1*n1);
    float b=(summuli*n1-sum*n2)/(n1*n1-n0*n2);
    currr<<<(n+maxth-1)/maxth,maxth>>>(data,n,a,b);
    cuds;
}

__global__ void cuGenW1(cufftReal* d_w,int n,int f1,int f2,int f3,int f4)
{
    //生成第一次滤波用的权重数组
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    if (id>=n)
        return;
    if (id<f1)
        d_w[id]=0;
    else
        if (id<f2)
            d_w[id]=(1-cos(PI*(id-f1)/(f2-f1)))/2;
        else
            if (id<f3)
                d_w[id]=1;
            else
                if (id<f4)
                    d_w[id]=(1+cos(PI*(id-f3)/(f4-f3)))/2;
                else
                    d_w[id]=0;
}

__global__ void cuFilter(int step,int n,cufftComplex* c,cufftReal* w)//c乘上w,结果存在c里
{
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    if (id>=step*n)
        return;
    int ida=(id/n)*(n+1)+id%n;//频域数据的下标
    //因为之前将每个频域分段从nspec+1个complex压缩到nspec个
    //所以这里要跳过那个分段里第nspec个complex
    int idb=id%n;//权重的下标
    c[ida].x*=w[idb];
    c[ida].y*=w[idb];
}

__global__ void cuDivn(cufftReal *r,int step,int n)//总计step*n个float,全部除以n
{
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    if (id>=step*n)
        return;
    r[id]=r[id]/n;
}

__global__ void cuPos(cufftReal *a,int n)
{
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    if (id>=n)
        return;
    if (a[id]<0)
        a[id]*=-1;
}

__global__ void cuSmooth(cufftReal *b,cufftReal *a,int step,int n,int* segsize,int w)
{
    //b输出,a输入,w平滑半径
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    if (id>=step*n)
        return;
    b[id]=0;
    int n1=id/n;//第几个分段
    int n2=id%n;//分段内位置
    int s,t;//平滑范围的头尾
    if (segsize==NULL)//没有提供具体的分段大小,则认为每一段长n
    {
        s=(n2-w>=0)?(n2-w):0;
        t=(n2+w<n)?(n2+w):(n-1);
    }
    else
    {
        if (n2>=segsize[n1])
            return;
        s=(n2-w>=0)?(n2-w):0;
        t=(n2+w<segsize[n1])?(n2+w):(segsize[n1]-1);
    }
    for (int i=s;i<=t;i++)
        b[id]+=a[n1*n+i];
    b[id]/=(2*w+1);
}

__global__ void cuNormalr(cufftReal *b,cufftReal *a,int step,int n,int *segsize)
{
    //b原始数据,a平滑数据
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    if (id>=step*n)
        return;
    int n1=id/n;
    int n2=id%n;
    if (n2<segsize[n1])
        if (a[id]<1e-30)
            b[id]*=1e30;
        else
            b[id]/=a[id];
}

__global__ void cuNormalc(cufftComplex *b,cufftReal *a,int step,int n)
{
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    if (id>=step*n)
        return;
    if (a[id]<1e-30)
    {
        b[id].x*=1e30;
        b[id].y*=1e30;
    }
    else
    {
        b[id].x/=a[id];
        b[id].y/=a[id];
    }
}

__global__ void cuCutTop(cufftReal *a,int n,float m)
{
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    if (id>=n)
        return;
    if (a[id]>m)
        a[id]=m;
    if (a[id]<-m)
        a[id]=-m;
}

__global__ void cuGenW2(cufftReal* b,cufftComplex* a,int step,int n)
{
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    if (id>=step*n)
        return;
    int n1=id/n;
    int n2=id%n;
    if (n2==0)
        b[id]=abs(a[n1*n].x);
    if (n2==n-1)
        b[id]=abs(a[n1*n].y);
    if ((n2!=0) && (n2!=n-1))
        b[id]=sqrt(a[id].x*a[id].x+a[id].y*a[id].y);
}

#endif