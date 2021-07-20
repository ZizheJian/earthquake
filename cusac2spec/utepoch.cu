#ifndef __UTEPOCH_CU__
#define __UTEPOCH_CU__

#include <stdio.h>
#include "sac.cuh"

long long int time2utepoch(int year,int jday,int hour,int min,int sec,int msec,int usec)
{
    long long int kk=1000000;
    long long int k=1000;
    long long int utepoch;
    int a4=year/4-!(year & 3);
    int a100=a4/25;
    int a400=a100/4;
    int intervening_leap_days=a4-a100+a400-477;
    int days=(365*(year-1970)+intervening_leap_days+jday-1);
    utepoch=(60*(60*((long long int)24*days+hour)+min)+sec)*kk+msec*k+usec;
    return utepoch;
}

long long int sactime2utepoch(SACHEAD hd,char mark)
{
    long long int kk=1000000;
    long long int utepoch;
    float marktime=0;
    int year,jday,hour,min,sec,msec;
    year=hd.nzyear;
    jday=hd.nzjday;
    hour=hd.nzhour;
    min=hd.nzmin;
    sec=hd.nzsec;
    msec=hd.nzmsec;
    if( mark=='b')
        marktime = hd.b;
    else 
        if(mark=='e')
            marktime=hd.e;
    utepoch=time2utepoch(year,jday,hour,min,sec,msec,0)+(long long int)((double)marktime*kk+0.5);
    return(utepoch);
}

#endif