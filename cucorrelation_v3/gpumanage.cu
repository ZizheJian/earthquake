#ifndef __GPUMANAGE_CU__
#define __GPUMANAGE_CU__

int seta()
{
    float avasize=gpuuse*gpusize;
    float a=segmax*(specmax+1)*soc+specmax*sor*2;
    float b=chimax*segmax*specmax*soc+chimax*sospi-avasize;
    return -b/a;
}

#endif