/*
./cucorrelation -i ../result/smallspec -o ../result/smallcc -c 3600 -x -m 0.5 -d 2
*/

#include "globalvar.cuh"
#include "ini.cu"
#include "file.cu"
#include "cc.cu"
#include "gpumanage.cu"

int main(int argc,char **argv)
{
    s=time(NULL);
    ini(argc,argv);
    ////////////////读入文件名和头////////////////
    specfilenum=0;
    prescan(specpath,1);
    speclist=(specinf*)malloc(specfilenum*sizeof(specinf));
    specfilenum=0;
    prescan(specpath,2);
    printf("specfilenum=%d\n",specfilenum);
    ////////////////文件分为团////////////////
    for (int i=0;i<specfilenum;i++)//排序
        for (int j=i+1;j<specfilenum;j++)
            if (speccp(speclist[i],speclist[j])==1)
            {
                specinf t=speclist[i];
                speclist[i]=speclist[j];
                speclist[j]=t;
            }
    chimax=1;
    chinum=1;
    for (int i=1,j=1;i<=specfilenum;i++)//统计最大团和团数
        if (i!=specfilenum)
            if (speccp(speclist[i],speclist[i-1])==0)
                j++;
            else
            {
                if (chimax<j)
                    chimax=j;
                j=1;
                chinum++;
            }
        else
            if (chimax<j)
                chimax=j;
    chin=(intx2*)malloc(chinum*soi*2);
    chin[0].x=0;
    for (int i=1,j=0;i<specfilenum;i++)//计算团的两端位置
        if (speccp(speclist[i],speclist[i-1])!=0)
        {
            chin[j].y=i-1;
            j++;
            chin[j].x=i;
        }
    chin[chinum-1].y=specfilenum-1;
    printf("chinum=%d, chimax=%d\n",chinum,chimax);
    ////////////////确定单个文件大小////////////////
    specmax=0;
    segmax=0;
    for (int i=0;i<specfilenum;i++)
    {
        if (specmax<speclist[i].head.spec)
            specmax=speclist[i].head.spec;
        if (segmax<speclist[i].head.seg)
            segmax=speclist[i].head.seg;
    }
    printf("segmax=%d, specmax=%d\n",segmax,specmax);
    ////////////////确定窗口大小并分配////////////////
    a=seta();
    printf("a=%d\n",a);
    c1=(cufftComplex*)malloc(chimax*segmax*specmax*soc);
    r1=(cufftReal*)malloc(a*specmax*sor*2);
    mlist=(intx2*)malloc(a*soi*2);
    cudaMalloc((void**)&d_c1,chimax*segmax*specmax*soc);
    cudaMalloc((void**)&d_c2,a*segmax*(specmax+1)*soc);
    cudaMalloc((void**)&d_r1,a*specmax*sor*2);
    cudaMalloc((void**)&d_speclist,chimax*sospi);
    cudaMalloc((void**)&d_mlist,a*soi*2);
    ////////////////cc////////////////
    for (int chii=0;chii<chinum;chii++)
    {
        int mx=1;
        int my=0;
        printf("chii=%d\n",chii);
        ////////////////读入////////////////
        int spec=speclist[chin[chii].x].head.spec;
        float dt=speclist[chin[chii].x].head.dt;
        int cchalfn=cchalf/dt;
        int chincur=chin[chii].y-chin[chii].x+1;
        for (int filei_a=chin[chii].x;filei_a<=chin[chii].y;filei_a++)
        {
            int filei_r=filei_a-chin[chii].x;
            FILE *specfile=fopen(speclist[filei_a].name,"rb");
            SPECHEAD trashhead;
            fread(&trashhead,sosph,1,specfile);
            for (int segi=0;segi<speclist[filei_a].head.seg;segi++)
                fread(c1+segmax*spec*filei_r+segi*spec,spec*soc,1,specfile);
            fclose(specfile);
        }
        cudaMemcpy(d_c1,c1,chincur*segmax*spec*soc,HD);
        cudaMemcpy(d_speclist,speclist+chin[chii].x,chincur*sospi,HD);
        for (int wini=0;wini<(chincur*(chincur-1)/2+a-1)/a;wini++)
        {
            ////////////////分配任务////////////////
            int wins=wini*a;
            int wint=min(wins+a-1,chincur*(chincur-1)/2);
            int winw=wint-wins+1;
            for (int mi=0;mi<winw;mi++)
            {
                mlist[mi].x=mx;
                mlist[mi].y=my;
                if (my==mx-1)
                {
                    my=0;
                    mx++;
                }
                else
                    my++;
            }
            printf("\twini=%d\n",wini);
            cudaMemcpy(d_mlist,mlist,a*soi*2,HD);
            ////////////////乘////////////////
            blocknum=(winw*segmax*spec+maxth-1)/maxth;
            mul<<<blocknum,maxth>>>(d_c1,d_c2,d_mlist,winw,segmax,spec);
            cuds;
            ////////////////加////////////////
            for (int gap=1;gap<segmax;gap*=2)
            {
                int gapnum=(segmax+gap-1)/(2*gap);
                blocknum=(winw*gapnum*spec+maxth-1)/maxth;
                add<<<blocknum,maxth>>>(d_c2,d_speclist,d_mlist,winw,segmax,spec,gap,gapnum);
                cuds;
            }
            printf("\t\tadd complete\n");
            for (int conj=0;conj<1;conj++)
            {
                ////////////////共轭////////////////
                ////////////////切换复数存储格式////////////////
                blocknum=(winw+maxth-1)/maxth;
                shift<<<blocknum,maxth>>>(d_c2,winw,segmax,spec);
                cuds;
                printf("\t\tshift complete\n");
                ////////////////fft////////////////
                len[0]=spec*2;
                inembed[0]=spec*2;
                inembed[1]=winw;
                onembed[0]=spec+1;
                onembed[1]=winw;
                cufftPlanMany(&plan,1,len,onembed,1,segmax*(spec+1),inembed,1,spec*2,CUFFT_C2R,winw);
                cufftExecC2R(plan,d_c2,d_r1);
                cufftDestroy(plan);
                printf("\t\tfft complete\n");
                ////////////////除////////////////
                blocknum=(winw*spec*2+maxth-1)/maxth;
                div<<<blocknum,maxth>>>(d_r1,d_speclist,d_mlist,winw,spec);
                cuds;
                printf("\t\tdiv complete\n");
                ////////////////输出////////////////
                cudaMemcpy(r1,d_r1,winw*spec*2*sor,DH);
                for (int mi=0;mi<winw;mi++)
                {   
                    int srci=chin[chii].x+mlist[mi].x;
                    int stai=chin[chii].x+mlist[mi].y;
                    for (int k=strlen(speclist[srci].name)-1;k>=0;k--)
                        if (speclist[srci].name[k]=='/')
                        {
                            strcpy(tpath,speclist[srci].name+k+1);
                            break;
                        }
                    tpath[strlen(tpath)-strlen(spectail)]='\0';
                    strcat(sacpath,tpath);
                    strcat(sacpath,"+");
                    for (int k=strlen(speclist[stai].name)-1;k>=0;k--)
                        if (speclist[stai].name[k]=='/')
                        {
                            strcpy(tpath,speclist[stai].name+k+1);
                            break;
                        }
                    tpath[strlen(tpath)-strlen(spectail)]='\0';
                    strcat(sacpath,tpath);
                    strcat(sacpath,sactail);
                    FILE *sacfile=fopen(sacpath,"wb");
                    if (sacfile==NULL)
                        printf("open failed, sac=%s\n",sacpath);
                    else
                    {
                        fwrite(r1+mi*spec*2+spec-cchalfn,(2*cchalfn+1)*sor,1,sacfile);
                        fclose(sacfile);
                    }
                    if (debug>=1)
                    {
                        cufftComplex* swap=(cufftComplex*)malloc(spec*soc);
                        float* finalcc=(float*)malloc((2*cchalfn+1)*sor);
                        float* ingcc=(float*)malloc((2*cchalfn+1)*sor);
                        memset(finalcc,0,(2*cchalfn+1)*sor);
                        for (int step=0;step<min(speclist[srci].head.seg,speclist[stai].head.seg);step++)
                        {
                            cc(c1+mlist[mi].x*segmax*spec+step*spec,c1+mlist[mi].y*segmax*spec+step*spec,swap,spec,dt,ingcc,cchalfn,step);
                            for (int p=0;p<2*cchalfn+1;p++)
                                finalcc[p]+=ingcc[p];
                        }
                        float diff=0;
                        float avg=0;
                        for (int p=0;p<2*cchalfn+1;p++)
                        {
                            finalcc[p]/=min(speclist[srci].head.seg,speclist[stai].head.seg);
                            diff=max(diff,abs(finalcc[p]-r1[mi*spec*2+spec-cchalfn+p]));
                            avg+=abs(r1[mi*spec*2+spec-cchalfn+p]);
                        }
                        if (debug>=2)
                        {
                            printf("multi thread:\n");
                            for (int p=0;p<2*cchalfn+1;p++)
                                printf("%f ",r1[mi*spec*2+spec-cchalfn+p]);
                            printf("\nsingle thread:\n");
                            for (int p=0;p<2*cchalfn+1;p++)
                                printf("%f ",finalcc[p]);
                            getchar();
                        }
                        avg/=(2*cchalfn+1);
                        printf("src=%d,  \tsta=%d,  \t%s: diff=%f, avg=%f\n",srci,stai,sacpath,diff,avg);
                        if (diff>=0.00001)
                            getchar();
                        free(swap);
                        free(finalcc);
                        free(ingcc);
                    }
                    sacpath[sacl]='\0';
                }
            }
        }
    }
    t=time(NULL);
    printf("total time=%ds\n",(int)(t-s));
}