#ifndef __SPEC_CUH__
#define __SPEC_CUH__

typedef struct{
    int type;
    char netwk[8];
    char stanm[8];
    char locid[8];
    char chnnm[8];
    float stlo;
    float stla;
    float cmpaz;
    float cmpinc;
    float win;
    float lag;
    int seg;
    int spec;
    float df;
    float dt;
    long long int utepoch;
}SPECHEAD;

typedef struct{
    char name[256];
    SPECHEAD head;
}specinf;

#endif