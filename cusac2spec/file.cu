#ifndef __FILE_CU__
#define __FILE_CU__

int isfolder(char *p)//判断路径是否是文件
{
    struct stat st;
    stat(p,&st);
    if (S_ISDIR(st.st_mode))
        return 1;
    return 0;
}

void addslash(char *p)//确保输入路径末尾有一个'/'
{
    if (p[strlen(p)-1]!='/')
        strcat(p,"/");
}

#endif