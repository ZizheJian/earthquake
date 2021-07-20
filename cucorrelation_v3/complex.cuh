#ifndef __COMPLEX_CUH__
#define __COMPLEX_CUH__

#define PI 3.1415926
#define IMAGE cmplx(0,1)
#define One cmplx(1,0)

cufftComplex cplus(cufftComplex a,cufftComplex b)
{
    a.x += b.x;
    a.y += b.y;
    return(a);
}

cufftComplex cmltp(cufftComplex a, cufftComplex b)
{
    cufftComplex c;
    c.x = a.x*b.x - a.y*b.y;
    c.y = a.x*b.y + a.y*b.x;
    return(c);
}

cufftComplex cngtv(cufftComplex a)
{
    a.x = -a.x;
    a.y = -a.y;
    return(a);
}

cufftComplex conjg(cufftComplex a)
{
    a.y = -a.y;
    return(a);
}

cufftComplex dmltp(float a, cufftComplex b)
{
    b.x *= a;
    b.y *= a;
    return(b);
}

cufftComplex cmplx(float x, float y)
{
    cufftComplex a;
    a.x = x;
    a.y = y;
    return(a);
}

#endif