#ifndef __FILE_CU__
#define __FILE_CU__

#include "globalvar.cuh"

int isfolder(char *p)
{
    struct stat st;
    stat(p,&st);
    if (S_ISDIR(st.st_mode))
        return 1;
    return 0;
}

void addslash(char *p)
{
    if (p[strlen(p)-1]!='/')
        strcat(p,"/");
}

int speccp(specinf a,specinf b)
{
    int TOL=(a.head.dt+b.head.dt)*1000000+1;
    if (a.head.utepoch>b.head.utepoch+TOL)
        return 1;
    if (a.head.utepoch+TOL<b.head.utepoch)
        return -1;
    if (a.head.spec>b.head.spec)
        return 1;
    if (a.head.spec<b.head.spec)
        return -1;
    if (a.head.dt>b.head.dt)
        return 1;
    if (a.head.dt<b.head.dt)
        return -1;
    if (a.head.win>b.head.win)
        return 1;
    if (a.head.win<b.head.win)
        return -1;
    if (a.head.lag>b.head.lag)
        return 1;
    if (a.head.lag<b.head.lag)
        return -1;
    return 0;
}

void prescan(char *specpath,int t)
{
    DIR *currdir=opendir(specpath);
    int specl=strlen(specpath);
    struct dirent *ent;
    while ((ent=readdir(currdir))!=NULL)
    {
        if (ent->d_name[strlen(ent->d_name)-1]=='.')
            continue;
        strcat(specpath,ent->d_name);
        if (isfolder(specpath))
            prescan(specpath,t);
        else
        {
            if (t==1)
                specfilenum++;
            if (t==2)
            {
                strcpy(speclist[specfilenum].name,specpath);
                FILE *specfile=fopen(specpath,"rb");
                fread(&speclist[specfilenum].head,sizeof(SPECHEAD),1,specfile);
                fclose(specfile);
                specfilenum++;
            }
        }
        specpath[specl]='\0';
    }
    closedir(currdir);
}

#endif