#ifndef __INI_CU__
#define __INI_CU__

#include "globalvar.cuh"

void usage()
{
    printf("sac2spec -i sac -o spec -w win -l lag -m3 halfwidth/wf1/wf2/wf3/wf4/cutN [-wf bf1/bf2/bf3/bf4] [-caz cmpaz -cinc cmpinc] [-ns nsmooth]\n");
    printf("-w  窗口宽度(7200s)\n");
    printf("-l  滑动距离(3600s)\n");
    printf("-m3 halfwidth/wf1/wf2/wf3/wf4/cutN: -m3 60/0.016/0.020/0.067/0.083/20\n");
    printf("-wf bf1/bf2/bf3/bf4:    -wf 0.0055/0.0066/0.400/0.454\n");
    printf("-caz -cinc  ???\n");
    printf("-ns 30\n");
}

void ini(int argc,char** argv)//获取maxth和输入指令
{
    filenum=0;
    datanum=0;
    cudaGetDeviceProperties(&gpuProp,0);
    maxth=gpuProp.maxThreadsPerBlock;
    maxn=0;
    maxstep=0;
    maxnspec=0;
    sacdata=NULL;
    specdata=NULL;
    segsize=NULL;
    d_segsize=NULL;
    d_sac_r1=NULL;
    d_sac_r2=NULL;
    d_sac_r3=NULL;
    d_sac_c1=NULL;
    planflag=0;
    if (argc <= 1) 
    {
        usage(); 
        exit(-1); 
    }
    for (int i=1;i<argc;i++) 
    {
        if (strcmp(argv[i],"-i")==0)
        {
            strcpy(sacrootpath,argv[i+1]);
            i++;
            continue;
        }
        if (strcmp(argv[i],"-o")==0)
        {
            strcpy(specrootpath,argv[i+1]);
            i++;
            continue;
        }
        if (strcmp(argv[i],"-w")==0)
        {
            win=atof(argv[i+1]);
            i++;
            continue;
        } 
        if (strcmp(argv[i],"-l")==0)
        {
            lag=atof(argv[i+1]);
            i++;
            continue;
        }
        if (strcmp(argv[i],"-m3")==0) 
        {  
            sscanf(argv[i+1],"%f/%f/%f/%f/%f/%f",&m3w,&m31,&m32,&m33,&m34,&m3top);
            i++;
            continue;
        }
        if (strcmp(argv[i],"-wf")==0) 
        {
            wf=1;
            sscanf(argv[i+1],"%f/%f/%f/%f",&wf1,&wf2,&wf3,&wf4);
            i++;
            continue;
        }
        if (strcmp(argv[i],"-caz")==0) 
        {
            cmpaz=atof(argv[i+1]);
            cmp++;
            i++;
            continue;
        }
        if (strcmp(argv[i],"-cinc")==0) 
        {
            cmpinc=atof(argv[i+1]);
            cmp++;
            i++;
            continue;
        }
        if (strcmp(argv[i],"-ns")==0) 
        {
            ns=atoi(argv[i+1]);
            i++;
            continue;
        }
    }
}

#endif