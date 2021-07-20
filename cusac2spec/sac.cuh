#ifndef __SAC_CUH__
#define __SAC_CUH__

typedef struct sac_head
{
    float delta;
    float depmin;
    float depmax;
    float scale;
    float odelta;
    float b;
    float e;
    float o;
    float a;
    float internal1;
    float t0;
    float t1;
    float t2;
    float t3;
    float t4;
    float t5;
    float t6;
    float t7;
    float t8;
    float t9;
    float f;
    float resp0;
    float resp1;
    float resp2;
    float resp3;
    float resp4;
    float resp5;
    float resp6;
    float resp7;
    float resp8;
    float resp9;
    float stla;
    float stlo;
    float stel;
    float stdp;
    float evla;
    float evlo;
    float evel;
    float evdp;
    float unused1;
    float user0;
    float user1;
    float user2;
    float user3;
    float user4;
    float user5;
    float user6;
    float user7;
    float user8;
    float user9;
    float dist;
    float az;
    float baz;
    float gcarc;
    float internal2;
    float internal3;
    float depmen;
    float cmpaz;
    float cmpinc;
    float unused2;
    float unused3;
    float unused4;
    float unused5;
    float unused6;
    float unused7;
    float unused8;
    float unused9;
    float unused10;
    float unused11;
    float unused12;
    int nzyear;
    int nzjday;
    int nzhour;
    int nzmin;
    int nzsec;
    int nzmsec;
    int internal4;
    int internal5;
    int	internal6;
    int	npts;
    int	internal7;
    int	internal8;
    int	unused13;
    int	unused14;
    int	unused15;
    int	iftype;
    int	idep;
    int	iztype;
    int	unused16;
    int	iinst;
    int	istreg;
    int	ievreg;
    int	ievtyp;
    int	iqual;
    int	isynth;
    int	unused17;
    int	unused18;
    int	unused19;
    int	unused20;
    int	unused21;
    int	unused22;
    int	unused23;
    int	unused24;
    int	unused25;
    int	unused26;
    int	leven;
    int	lpspol;
    int	lovrok;
    int	lcalda;
    int	unused27;
    char kstnm[8];
    char kevnm[16];
    char khole[8];
    char ko[8];
    char ka[8];
    char kt0[8];
    char kt1[8];
    char kt2[8];
    char kt3[8];
    char kt4[8];
    char kt5[8];
    char kt6[8];
    char kt7[8];
    char kt8[8];
    char kt9[8];
    char kf[8];
    char kuser0[8];
    char kuser1[8];
    char kuser2[8];
    char kcmpnm[8];
    char knetwk[8];
    char kdatrd[8];
    char kinst[8];
}SACHEAD;

static SACHEAD sacnull={
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,     6,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    -12345,-12345,-12345,-12345,-12345,
    {'-','1','2','3','4','5',' ',' '},
    {'-','1','2','3','4','5',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '},
    {'-','1','2','3','4','5',' ',' '},{'-','1','2','3','4','5',' ',' '},
    {'-','1','2','3','4','5',' ',' '},{'-','1','2','3','4','5',' ',' '},
    {'-','1','2','3','4','5',' ',' '},{'-','1','2','3','4','5',' ',' '},
    {'-','1','2','3','4','5',' ',' '},{'-','1','2','3','4','5',' ',' '},
    {'-','1','2','3','4','5',' ',' '},{'-','1','2','3','4','5',' ',' '},
    {'-','1','2','3','4','5',' ',' '},{'-','1','2','3','4','5',' ',' '},
    {'-','1','2','3','4','5',' ',' '},{'-','1','2','3','4','5',' ',' '},
    {'-','1','2','3','4','5',' ',' '},{'-','1','2','3','4','5',' ',' '},
    {'-','1','2','3','4','5',' ',' '},{'-','1','2','3','4','5',' ',' '},
    {'-','1','2','3','4','5',' ',' '},{'-','1','2','3','4','5',' ',' '},
    {'-','1','2','3','4','5',' ',' '}
};

#endif