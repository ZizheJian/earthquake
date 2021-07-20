#ifndef __SETNM_CU__
#define __SETNM_CU__

#include <stdio.h>

char nonet[8]="NONET\0";
char nosta[8]="NOSTA\0";
char voidnum[8]="-12345\0";
char chn[8]="CHN\0";

void setnm(char *dst,char *src,int n,char *dft)
{
    int i=0;
    int iuse=0;
    int isblank=1;
    char *ssrc=(char*)malloc(256*sizeof(char));
    while ((src[i]!='\0') && (i<n))
    {
        if(src[i]!=' ')
        {
            isblank=0;
            break;
        }
        i++;
    }
    if (isblank)
        strcpy(ssrc,dft);
    else
        strcpy(ssrc,src);
    for(i=0;i<n;i++)
        if(ssrc[i]!='\0')
            if(ssrc[i]!=' ')
            {
                dst[iuse]=ssrc[i];
                iuse++;
            }
        else
            break;
    if (iuse<n)
        dst[iuse]='\0';
    else
        dst[n-1]='\0';
    free(ssrc);
}

#endif