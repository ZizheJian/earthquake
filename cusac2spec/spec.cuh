#ifndef __SPEC_CUH__
#define __SPEC_CUH__

typedef struct spec_head{
    int type;
    char netwk[8];
    char stanm[8];
    char locid[8];
    char chnnm[8];
    float stlo;
    float stla;
    float cmpaz;
    float cmpinc;
    float winlen;
    float laglen;
    int nseg;
    int nspec;
    float df;
    float dt;
    long long int utepoch;
}SPECHEAD;

#endif