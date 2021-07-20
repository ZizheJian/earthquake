#ifndef __SAC_CU__
#define __SAC_CU__

#include "sac.cuh"
#include "globalvar.cuh"

void readsac(FILE *sacfile)//读如sac的头和数据
{
    fread(&sachd,sizeof(SACHEAD),1,sacfile);
    filenum++;
    datanum+=sachd.npts;
    n=sachd.npts;
    dt=sachd.delta;
    winn=win/dt;
    lagn=lag/dt;
    lagrate=(float)(n-winn+lagn-1)/lagn;//和原代码不同,考虑到最后一小段
    step=lagrate+1;
    m3wn=m3w/dt;
    nfft=2;
    while (nfft<winn)
        nfft*=2;
    nspec=nfft/2;
    df=1/(dt*nfft);
    if (n>maxn)
    {
        maxn=n;
        if (sacdata!=NULL)
        {
            free(sacdata);
            sacdata=NULL;
        }
    }
    if (step>maxstep)
    {
        maxstep=step;
        if (segsize!=NULL)
        {
            free(segsize);
            segsize=NULL;
        }
        if (d_segsize!=NULL)
        {
            cudaFree(d_segsize);
            d_segsize=NULL;
        }
        if (d_sac_r1!=NULL)
        {
            cudaFree(d_sac_r1);
            d_sac_r1=NULL;
        }
        if (d_sac_r2!=NULL)
        {
            cudaFree(d_sac_r2);
            d_sac_r2=NULL;
        }
        if (d_sac_r3!=NULL)
        {
            cudaFree(d_sac_r3);
            d_sac_r3=NULL;
        }
        if (d_sac_c1!=NULL)
        {
            cudaFree(d_sac_c1);
            d_sac_c1=NULL;
        }
        if (planflag)
        {
            cufftDestroy(planR2C);
            cufftDestroy(planC2R);
            planflag=0;
        }
    }
    if (nspec>maxnspec)
    {
        maxnspec=nspec;
        if (specdata!=NULL)
        {
            free(specdata);
            specdata=NULL;
        }
        if (d_sac_r1!=NULL)
        {
            cudaFree(d_sac_r1);
            d_sac_r1=NULL;
        }
        if (d_sac_r2!=NULL)
        {
            cudaFree(d_sac_r2);
            d_sac_r2=NULL;
        }
        if (d_sac_r3!=NULL)
        {
            cudaFree(d_sac_r3);
            d_sac_r3=NULL;
        }
        if (d_sac_c1!=NULL)
        {
            cudaFree(d_sac_c1);
            d_sac_c1=NULL;
        }
        if (planflag)
        {
            cufftDestroy(planR2C);
            cufftDestroy(planC2R);
            planflag=0;
        }
    }
    if (sacdata==NULL)
        sacdata=(float*)malloc(maxn*sor);
    if (segsize==NULL)
        segsize=(int*)malloc(maxstep*soi);
    if (specdata==NULL)
        specdata=(float*)malloc(maxnspec*soc);
    if (d_segsize==NULL)
        cudaMalloc((void**)&d_segsize,maxstep*soi);
    if (d_sac_r1==NULL)
        cudaMalloc((void**)&d_sac_r1,maxstep*maxnspec*2*sor);
    if (d_sac_r2==NULL)
        cudaMalloc((void**)&d_sac_r2,maxstep*maxnspec*2*sor);
    if (d_sac_r3==NULL)
        cudaMalloc((void**)&d_sac_r3,maxstep*maxnspec*2*sor);
    if (d_sac_c1==NULL)
        cudaMalloc((void**)&d_sac_c1,maxstep*(maxnspec+1)*soc);
    if (!planflag)
    {
        rank=1;
        len[0]=nfft;
        inembed[0]=nfft;    inembed[1]=step;
        istride=1;
        idist=nfft;
        onembed[0]=nspec+1; onembed[1]=step;
        ostride=1;
        odist=nspec+1;
        batch=step;
        cufftPlanMany(&planR2C,rank,len,inembed,istride,idist,onembed,ostride,odist,CUFFT_R2C,batch);
        cufftPlanMany(&planC2R,rank,len,onembed,ostride,odist,inembed,istride,idist,CUFFT_C2R,batch);
        planflag=1;
    }
    fread(sacdata,n*sor,1,sacfile);
}

#endif