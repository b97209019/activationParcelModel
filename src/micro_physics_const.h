#ifndef __MICRO_PHYSICS_CONST_H__
#define __MICRO_PHYSICS_CONST_H__

typedef long double real8;

#define M_1D3  0.33333333333333333333L
#define M_2D3  0.66666666666666666667L
#define M_4D3  1.33333333333333333333L

//#define M_PI    3.14159265358979324L
#define M_4PI   12.56637061435917295L
#define M_4PID3  4.18879020478639098L
#define M_PID6   0.52359877559829887L
//#define M_SQRT2   1.41421356237309505L
//#define M_3DSQRT2  2.1213203435596424L
#define M_SQRT3    1.73205080756887729L
#define M_CBRT2    1.259921049894582L
#define M_SQRTPI   1.77245385090551603L
#define M_SQRTPID2 1.25331413731550025L
#define M_SQRT2PI  2.50662827463100050L
#define M_4D27     0.14814814814814814L

#define SS2ETDEW 15.39414112281038
#define POWERLAWSIZELIMIT 0

// [kg/mol]
#define CMWDRY .028964
#define CMWH2O .018015

// [J/(kg.K)]
#define CR    8.3143
#define CRDRY 287.05
#define CRH2O 461.51
#define CCpDRY 1005.
#define CCpH2O 1850.

#define CE       .6220 // CRDRY/CRH2O
#define CrhoH2O0C 999.8396

#if defined __AS__ // for SO4(HN3)2
#define CalHygro
#define CrhoSALTDRY 1769.
#define CrhoSALTWET 2093.
#define CMWSALT  .13214
#define Cnu 3 // number of desovled ions

#elif defined __NACL__ //for NaCl
#define CalHygro
#define CrhoSALTDRY 2165.
#define CrhoSALTWET 2726.
#define CMWSALT  .05844
#define Cnu 2 // number of desovled ions

#elif defined __S__//for sulfate
#define CalHygro
#define CrhoSALTDRY 1841.
#define CrhoSALTWET 2510.
#define CMWSALT  .09606
#define Cnu 3 // number of desovled ions

#elif defined __N__ //for nitrate
#define CalHygro
#define CrhoSALTDRY 1503.
#define CrhoSALTWET 1909.
#define CMWSALT  .0620049
#define Cnu 2 // number of desovled ions

#elif defined __HCL__ //for hydrochloric
#define CalHygro
#define CrhoSALTDRY 1187.
#define CrhoSALTWET 1782.
#define CMWSALT  .03646
#define Cnu 2 // number of desovled ions

#elif defined __BC__ //for hydrochloric
#define CalHygro
#define CrhoSALTDRY 1187.
#define CrhoSALTWET 1782.
#define CMWSALT  .05844
#define Cnu 2 // number of desovled ions
#endif

#define CK0C  273.15
#define CP0 1.01325e5

#define CVSDYN0 1.718e-5
#define CHMFP0  6.62e-8
#define CG    9.8

#define CRCUTMAX 3.e-3L
#define CRCUTmin 1.0e-9L
#define CRCUT1   1.5e-9L

#ifdef __TCK__
#undef CRCUTMAX
#undef CBIN_Q
//#define CRCUTMAX 1.5e-6L
//#define CBIN_Q  1.148698354997035L
//#define CRCUTMAX 4e-7L
#define CBIN_Q  1.122462048309373L
#endif

#endif
