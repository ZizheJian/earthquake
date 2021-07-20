#ifndef __FFT_CU__
#define __FFT_CU__

#include "complex.cuh"

void fft(cufftComplex *a,int n,float dt)
{
    double pi=-PI;
    if (dt<0) 
        pi=PI;
    for (int m=n/2,j=0,i=1;i<n-1;i++)
    {
        int k;
        for (k=m;k<=j;k/=2)
            j-=k;
        j+=k;
        if(i<j) 
        {
            cufftComplex t=a[i];
            a[i]=a[j];
            a[j]=t;
        }
    }
    for (int m=1,step=2;m<n;m=step,step*=2)
    {
        cufftComplex u=One;
        cufftComplex w=cmplx(cos(pi/m),sin(pi/m));
        for (int j=0;j<m;j++)
        {
            for (int i=j;i<n;i+=step) 
            {
	            int k=i+m;
	            cufftComplex t=cmltp(a[k], u);
	            a[k]=cplus(a[i],cngtv(t));
	            a[i]=cplus(a[i],t);
            }
            u=cmltp(u,w);
        }
    }
    if (dt<0) 
        dt=-1/(n*dt);
    for (int i=0;i<n;i++) 
        a[i]=dmltp(dt,a[i]);
}

void fftr(cufftComplex *x,int n,float dt)
{
    int n2=n/2;
    float delw=PI/n;
    cufftComplex isg=IMAGE;
    if (dt>0)
    {
        delw=-delw;
        isg=cngtv(isg);
        fft(x,n,dt);
    }
    x[0]=cmplx(x[0].x+x[0].y,x[0].x-x[0].y);
    float w=delw;
    for (int i=1;i<n2;i++)
    {
        int j=n-i;
        cufftComplex t=conjg(x[j]);
        cufftComplex g=cplus(x[i],t);
        cufftComplex h=cplus(x[i],cngtv(t));
        h=cmltp(cmplx(cos(w),sin(w)),h);
        x[i]=dmltp(0.5,cplus(g,cmltp(isg,h)));
        x[j]=dmltp(0.5,cplus(conjg(g),cmltp(isg,conjg(h))));
        w+=delw;
    }
    x[n2]=conjg(x[n2]);
    if (dt<0) 
    {
        x[0]=dmltp(0.5,x[0]);
        fft(x,n,dt);
    }
}

#endif