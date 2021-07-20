#ifndef __CC_CU__
#define __CC_CU__

#include "fft.cu"

__global__ void mul(cufftComplex *c1,cufftComplex *c2,intx2 *mlist,int winw,int segmax,int spec)
{
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    int mi=id/(segmax*spec);
    int segi=(id/spec)%segmax;
    int speci=id%spec;
    int x=mlist[mi].x;
    int y=mlist[mi].y;
    if (mi>=winw)
        return;
    int xx=(x*segmax+segi)*spec+speci;
    int yy=(y*segmax+segi)*spec+speci;
    int zz=(mi*segmax+segi)*(spec+1)+speci;
    if (speci==0)
    {
        c2[zz].x=c1[xx].x*c1[yy].x;
        c2[zz].y=c1[xx].y*c1[yy].y;
    }
    else
    {
        c2[zz].x=(1-(speci%2)*2)*(c1[xx].x*c1[yy].x+c1[xx].y*c1[yy].y);
        c2[zz].y=(1-(speci%2)*2)*(c1[xx].x*c1[yy].y-c1[xx].y*c1[yy].x);
    }
}

__global__ void add(cufftComplex *c2,specinf *speclist,intx2 *mlist,int winw,int segmax,int spec,int gap,int gapnum)
{
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    int mi=id/(gapnum*spec);
    int segi=((id/spec)%gapnum)*gap*2;
    int speci=id%spec;
    int x=mlist[mi].x;
    int y=mlist[mi].y;
    if (mi>=winw)
        return;
    if ((segi+gap>=speclist[x].head.seg) || (segi+gap>=speclist[y].head.seg))
        return;
    int xx=(mi*segmax+segi)*(spec+1)+speci;
    int yy=(mi*segmax+segi+gap)*(spec+1)+speci;
    c2[xx].x+=c2[yy].x;
    c2[xx].y+=c2[yy].y;
}

__global__ void shift(cufftComplex *c2,int winw,int segmax,int spec)
{
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    int mi=id;
    if (mi>=winw)
        return;
    int xx=mi*segmax*(spec+1);
    int yy=mi*segmax*(spec+1)+spec;
    c2[yy].x=c2[xx].y;
    c2[yy].y=0;
    c2[xx].y=0;
}

__global__ void div(cufftReal *r1,specinf *speclist,intx2 *mlist,int winw,int spec)
{
    int id=blockIdx.x*blockDim.x+threadIdx.x;
    int mi=id/(spec*2);
    if (mi>=winw)
        return;
    int x=mlist[mi].x;
    int y=mlist[mi].y;
    int seg=speclist[x].head.seg;
    if (seg>speclist[y].head.seg)
        seg=speclist[y].head.seg;
    int a=2*spec*seg;
    r1[id]/=a;
}

void cc(cufftComplex *source,cufftComplex *station,cufftComplex *swap,int spec,float dt,float *cc,int cchalfn,int stepidx)
{
    memcpy(swap,station,spec*soc);
    swap[0]=cmplx(swap[0].x*source[0].x,swap[0].y*source[0].y);
    float aa=-1;
    for (int j=1;j<spec;j++)
    {
        swap[j]=cmltp(swap[j],conjg(source[j]));
        swap[j]=dmltp(aa,swap[j]);
        aa=-aa;
    }
    fftr(swap,spec,-dt);
    int ncc=2*cchalfn+1;
    memcpy(cc,(float*)swap+spec-cchalfn,ncc*sizeof(float));
}

#endif