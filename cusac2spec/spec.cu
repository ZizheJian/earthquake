#ifndef __SPEC_CU__
#define __SPEC_CU__

#include "spec.cuh"
#include "setnm.cu"
#include "utepoch.cu"
#include "valid.cu"
#include "file.cu"

void writespechead(FILE *specfile)//生成spec的头并写出
{
    spechd.type=1;
    if(strncmp(sachd.knetwk,"-12345  ",8)==0)
        sachd.knetwk[0]='\0';
    setnm(spechd.netwk,sachd.knetwk,8,nonet);
    if(strncmp(sachd.kstnm,"-12345  ",8)==0)
        sachd.kstnm[0]='\0';
    setnm(spechd.stanm,sachd.kstnm,8,nosta);
    if(strncmp(sachd.khole,"-12345  ",8)==0)
        sachd.khole[0]='\0';
    setnm(spechd.locid,sachd.khole,8,voidnum);
    spechd.nseg=step;
    spechd.winlen=win;
    spechd.laglen=lag;
    spechd.nspec=nspec;
    spechd.df=df;
    spechd.dt=dt;
    spechd.stlo=sachd.stlo;
    spechd.stla=sachd.stla;
    spechd.utepoch=sactime2utepoch(sachd,'b');
    setnm(spechd.chnnm,sachd.kcmpnm,8,chn);
    spechd.cmpaz=-12345.0;
    spechd.cmpinc=-12345.0;
    if(validf(sachd.cmpaz,-360.0,360.0))  
        spechd.cmpaz=sachd.cmpaz;
    if(validf(sachd.cmpinc,-360.0,360.0)==1)
        spechd.cmpinc=sachd.cmpinc;
    if(cmp==2) 
    { 
        spechd.cmpaz=cmpaz;
        spechd.cmpinc=cmpinc;
    }
    fwrite(&spechd,sizeof(SPECHEAD),1,specfile);
}

void spec_name_change(char *specrootpath,char *sacpath)//在specrootpath后面加上sac的文件名和".spec"
{
    if (isfolder(specrootpath))
    {
        sacpath[strlen(sacpath)-strlen(sactail)]='\0';//sac路径去掉后缀名
        int sacl;
        for (sacl=strlen(sacpath)-1;sacl>=0;sacl--)//求出sac文件开头的位置
            if (sacpath[sacl]=='/')
                break;
        sacl++;
        strcat(specrootpath,sacpath+sacl);
        strcat(specrootpath,spectail);
    }
}

#endif