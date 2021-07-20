/*
./cusac2spec -i ../result/sac.sacum -o ../result/spec.spec -w 8 -l 4 -m3 1/0.1/0.3/0.4/0.5/5 -wf 0.1/0.3/0.4/0.5 -ns 5
./cusac2spec -i ../result/smallsac -o ../result/smallspec -w 7200 -l 3600 -m3 60/0.016/0.020/0.067/0.083/20 -wf 0.0055/0.0066/0.400/0.454 -ns 30
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "globalvar.cuh"
#include "ini.cu"
#include "sac.cu"
#include "spec.cu"
#include "fft.cu"

void f(char *sacpath)
{
    FILE *sacfile=fopen(sacpath,"rb");
    readsac(sacfile);//读入sac的头和数据
    fclose(sacfile);
    spec_name_change(specrootpath,sacpath);//生成spec的完整文件名
    FILE *specfile=fopen(specrootpath,"wb");
    specrootpath[specl]='\0';//把spec路径改回去
    writespechead(specfile);//写spec的头

    cudaMemcpy(d_sac_r1,sacdata,n*sor,HD);//d_sac_r1保存原始数据
    ////////////////分段存储+rrr////////////////
    for (int i=step-1;i>=0;i--)
    {
        if (winn+lagn*i<=n)
            segsize[i]=winn;
        else
            segsize[i]=n-lagn*i;
        cudaMemset(d_sac_r2+i*nfft,0,nfft*sor);
        cudaMemcpy(d_sac_r2+i*nfft,d_sac_r1+i*lagn,segsize[i]*sor,DD);
        rrrcontrol(d_sac_r2+i*nfft,segsize[i],d_sac_r1+i*nfft);///rrr=rdc+rtr
        //这里和原代码不同
        //原代码中先对所有数据作rrr,再对分段数据作rrr
        //这里直接对分段数据作rrr
    }
    cudaMemcpy(d_segsize,segsize,step*soi,HD);
    cudaMemcpy(d_sac_r1,d_sac_r2,step*nfft*sor,DD);//d_sac_r1,d_sac_r2保存分段数据
    ////////////////第一次FFT////////////////
    cufftExecR2C(planR2C,d_sac_r2,d_sac_c1);
    for (int i=0;i<step;i++)//每段数据压缩到nspec位
    {
        cufftComplex t;
        cudaMemcpy(&t,d_sac_c1+i*(nspec+1)+nspec,soc,DH);
        t.y=t.x;
        t.x=0;
        cudaMemcpy(d_sac_c1+i*(nspec+1),&t,soc,HD);
    }
    //d_sac_r2保存分段频域数据
    ////////////////第一次滤波////////////////
    m31n=m31/df;
    m32n=m32/df;
    m33n=m33/df;
    m34n=m34/df;
    if (m34n>=nspec)
        m34n=nspec-1;
    cuGenW1<<<(nspec+maxth-1)/maxth,maxth>>>(d_sac_r3,nspec,m31n,m32n,m33n,m34n);
    //d_sac_r3保存权重
    cuds;
    cuFilter<<<(step*nspec+maxth-1)/maxth,maxth>>>(step,nspec,d_sac_c1,d_sac_r3);
    cuds;
    ////////////////第二次FFT////////////////
    for (int i=0;i<step;i++)//每段数据扩展到nspec+1位
    {
        cufftComplex t;
        cudaMemcpy(&t,d_sac_c1+i*(nspec+1),soc,DH);
        t.x=t.y;
        t.y=0;
        cudaMemcpy(d_sac_c1+i*(nspec+1)+nspec,&t,soc,HD);
        t.x=0;
        cudaMemcpy(d_sac_c1+i*(nspec+1),&t,soc,HD);
    }
    cufftExecC2R(planC2R,d_sac_c1,d_sac_r2);
    cuDivn<<<all>>>(d_sac_r2,step,nfft);
    cuds;
    //d_sac_r2保存滤波后的分段时域数据
    ////////////////临域取平均////////////////
    cuPos<<<all>>>(d_sac_r2,step*nfft);//全部取绝对值
    cuds;
    cudaMemcpy(d_sac_r3,d_sac_r2,step*nfft*sor,DD);
    cuSmooth<<<all>>>(d_sac_r2,d_sac_r3,step,nfft,d_segsize,m3wn);
    ////////////////第一次归一化+削平////////////////
    cuNormalr<<<all>>>(d_sac_r1,d_sac_r2,step,nfft,d_segsize);
    //这里d_sac_r1为刚分好段时的时域数据,d_sac_r2为平滑后的数据
    cuds;
    cuCutTop<<<all>>>(d_sac_r1,step*nfft,m3top);
    cuds;
    ////////////////频谱白化////////////////
    if (wf)
    {
        ////////////////第三次FFT////////////////
        cufftExecR2C(planR2C,d_sac_r1,d_sac_c1);
        for (int i=0;i<step;i++)//这里没有执行rrr,频域中一些之前为0的点现在不为0,压缩方法不一样
        {
            cufftComplex t1,t2;
            cudaMemcpy(&t1,d_sac_c1+i*(nspec+1),soc,DH);
            cudaMemcpy(&t2,d_sac_c1+i*(nspec+1)+nspec,soc,DH);
            t1.y=t2.x;
            cudaMemcpy(d_sac_c1+i*(nspec+1),&t1,soc,HD);
        }
        ////////////////第二次归一化////////////////
        cuGenW2<<<(step*(nspec+1)+maxth-1)/maxth,maxth>>>(d_sac_r2,d_sac_c1,step,nspec+1);
        cuSmooth<<<(step*(nspec+1)+maxth-1)/maxth,maxth>>>(d_sac_r3,d_sac_r2,step,nspec+1,NULL,ns);
        cuNormalc<<<all>>>(d_sac_c1,d_sac_r3,step,nspec+1);
        //此时频域数据在d_sac_c1里
        ////////////////第二次滤波////////////////
        wf1n=wf1/df;
        wf2n=wf2/df;
        wf3n=wf3/df;
        wf4n=wf4/df;
        if (wf4n>=nspec)
            wf4n=nspec-1;
        cuGenW1<<<(nspec+maxth-1)/maxth,maxth>>>(d_sac_r2,nspec,wf1n,wf2n,wf3n,wf4n);
        cuds;
        cuFilter<<<(step*nspec+maxth-1)/maxth,maxth>>>(step,nspec,d_sac_c1,d_sac_r2);
        cuds;
    }
    for (int i=0;i<step;i++)
    {
        cudaMemcpy(specdata,d_sac_c1+i*(nspec+1),nspec*soc,DH);
        fwrite(specdata,nspec*soc,1,specfile);
    }

    fclose(specfile);
}

void brousesac(char *sacpath)
{
    if (strcmp(sacpath+strlen(sacpath)-strlen(sactail),sactail)==0)//输入路径是个sacum文件
    {
        f(sacpath);
        return;
    }
    if (!isfolder(sacpath))//不是sacum文件,且不是文件夹
        return;
    DIR *currdir=opendir(sacpath);
    addslash(sacpath);
    struct dirent *ent;
    int sacl=strlen(sacpath);//当输入路径是文件夹时有效,记录当前文件夹的路径的长度
    while ((ent=readdir(currdir))!=NULL)//遍历当前文件夹里所有文件(夹)
    {
        if (ent->d_name[strlen(ent->d_name)-1]=='.')
            continue;
        strcat(sacpath,ent->d_name);//现在sacpath为准备访问的文件(夹)名
        if (isfolder(sacpath))
            brousesac(sacpath);//如果是文件夹,则递归
        else
            f(sacpath);
        sacpath[sacl]='\0';//最后把加上的"当前访问的文件(夹)名"去掉
    }
    closedir(currdir);
}

int main(int argc,char **argv)
{
    s=clock();
    ini(argc,argv);
    if (isfolder(specrootpath))
        addslash(specrootpath);
    //如果输入路径是一个文件,则只处理这一个文件
    //如果输入路径是一个文件夹,则处理包含的所有文件
    //如果输出路径是一个文件,则要求输入路径也是文件,并直接写到文件内(不考虑文件名)
    //如果输出路径是一个文件夹,则在文件夹内新建和输入文件同名的文件(后缀名不同)并输出
    if ((isfolder(sacrootpath)) && (!isfolder(specrootpath)))
        return 0;
    specl=strlen(specrootpath);
    brousesac(sacrootpath);
    t=clock();
    printf("time=%fs   filenum=%d   datanum=%d\n",(float)(t-s)/1000000,filenum,datanum);
}