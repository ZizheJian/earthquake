#ifndef __VALID_CU__
#define __VALID_CU__

#include <stdio.h>

int validf(float in,float min,float max)
{
    if ((in>=min) && (in<=max))
        return 1;
    return 0;
}

#endif