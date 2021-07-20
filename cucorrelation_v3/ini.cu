#ifndef __INI_CU__
#define __INI_CU__

#include "globalvar.cuh"
#include "file.cu"

void usage()
{
    printf("correlation -i specpath -o sacpath -c halfCCLength [-n KCMPNM] [-x] [-m gpumemuse] [-d debug]\n");
    printf("-i specpath\n");
    printf("-o sacpath\n");
    printf("-c half of cclenth: -C 3600\n");
    printf("-n kcmpnm\n");
    printf("-x output symmetric NCF or not\n");
    printf("-m proportion of GPU memory available, default=0.5\n");
    printf("-d use the single-thread algorithm to verify the result, larger number = more info\n");
}

void ini(int argc,char** argv)
{

    cudaGetDeviceProperties(&gpuProp,0);
    maxth=gpuProp.maxThreadsPerBlock;
    printf("maxthread=%d\n",maxth);
	cudaMemGetInfo(&gpusize,&totgpusize);
	printf("gpusize=%dMB\n",(int)(gpusize/1048576));
    if (argc<=1) 
    {
        usage(); 
        exit(-1); 
    }
    for (int i=1;i<argc;i++) 
    {
        if (strcmp(argv[i],"-i")==0)
        {
            strcpy(specpath,argv[i+1]);
            addslash(specpath);
            i++;
            continue;
        }
        if (strcmp(argv[i],"-o")==0)
        {
            strcpy(sacpath,argv[i+1]);
            addslash(sacpath);
            sacl=strlen(sacpath);
            i++;
            continue;
        }
        if (strcmp(argv[i],"-c")==0)
        {
            cchalf=atof(argv[i+1]);
            i++;
            continue;
        } 
        if (strcmp(argv[i],"-n")==0)
        {
            setchn=1;
            strcpy(ncfchnname,argv[i+1]);
            i++;
            continue;
        }
        if (strcmp(argv[i],"-x")==0)
        {
            symmetric=1;
            continue;
        }
        if (strcmp(argv[i],"-m")==0)
        {
            gpuuse=atof(argv[i+1]);
            i++;
            continue;
        }
        if (strcmp(argv[i],"-d")==0)
        {
            debug=atoi(argv[i+1]);
            i++;
            continue;
        }
    }
}

#endif