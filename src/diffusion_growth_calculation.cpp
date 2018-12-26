#include"diffusion_growth_header_file.h"
void set_aerosol(real8 n_aerosol_mode, real8 *Nts,
                 real8 *SIGMAs, real8 *MUs, real8 *kaps,
                 real8 *dry_r3s, real8 *ns, real8 *kapxrhow0Cs){
// creat bin limit
    real8 *bin_rcut = new real8[g_nbin+1];
    bin_rcut[0] = CRCUT1/g_binq;
    for (int i = 0; i != g_nbin; ++i){
        dry_r3s[i] = ns[i] = kapxrhow0Cs[i] = 0;
        bin_rcut[i+1] = bin_rcut[i] * g_binq;
    }
    bin_rcut[g_nbin] = CRCUTMAX;   // 3.e-3;
    bin_rcut[0] = CRCUTmin;        // 2.e-10;

//  for (int j = 0; j != g_nbin+1; ++j)
//    printf("%d, %Le\n", j, bin_rcut[j]);
    for (int i = 0; i != n_aerosol_mode; ++i){
        real8 hNi = Nts[i]*2.45e4/CMWDRY*.5,
              sqrt2sigma = (M_SQRT2*SIGMAs[i]),
              logmu = logl(MUs[i]),
              r3erfckernel = - 1.5*sqrt2sigma,
              r3kernel = expl(3.*logmu+4.5*SIGMAs[i]*SIGMAs[i]),
              nerfckernel1, nerfckernel2, n1, n2, r31, r32,
              kapxrho0C = kaps[i]*CrhoH2O0C;
/*
          nerfckernel = (logl(bin_rcut[0])-logmu)/(sqrt2sigma),
          n1 = erfcl( nerfckernel ), n2,
          r31 = erfcl( nerfckernel + r3erfckernel ), r32;
*/
        for (int j = 0; j != g_nbin; ++j){
            nerfckernel1 = (logl(bin_rcut[j  ])-logmu)/(sqrt2sigma);
            nerfckernel2 = (logl(bin_rcut[j+1])-logmu)/(sqrt2sigma);
            if ( nerfckernel1 < 0 && nerfckernel2 < 0 ){
                n2  = erfcl( -nerfckernel1 );                    // 1-erf()
                n1  = erfcl( -nerfckernel2 );                    // 1-erf()
                r32 = erfcl( -nerfckernel1 - r3erfckernel );     // 1-erf()
                r31 = erfcl( -nerfckernel2 - r3erfckernel );     // 1-erf()
            }
            else {
                n1  = erfcl( nerfckernel1 );                    // 1-erf()
                n2  = erfcl( nerfckernel2 );                    // 1-erf()
                r31 = erfcl( nerfckernel1 + r3erfckernel );     // 1-erf()
                r32 = erfcl( nerfckernel2 + r3erfckernel );     // 1-erf()
            }
            ns[j] += hNi*(n1 - n2);
            real8 tmp1 = hNi*r3kernel*(r31 - r32);
            dry_r3s[j] += tmp1;
            kapxrhow0Cs[j] += tmp1*kapxrho0C;
        }
    }
// The particle_ns[g_nbin] stores the total number of aerosol
    ns[g_nbin] = 0;
    for (int i = 0; i != g_nbin; ++i){
        if (ns[i] != 0 and dry_r3s[i] != 0){
            kapxrhow0Cs[i] /= dry_r3s[i];
            dry_r3s[i] /= ns[i];
        } else
            dry_r3s[i] = ns[i] = kapxrhow0Cs[i] = 0;
        ns[g_nbin] += ns[i];
    }
    return;
}
void set_aerosol(char *aerosol_spectrum_file_name,
                 real8 *dry_r3s, real8 *ns, real8 *kapxrhow0Cs){
    ns[g_nbin] = 0;
    char aerosol_spectrum_full_file_name[64];
    strcpy(aerosol_spectrum_full_file_name, "../input/");
    strcat(aerosol_spectrum_full_file_name, aerosol_spectrum_file_name);
    FILE *aerosol_file_ID = fopen(aerosol_spectrum_full_file_name, "r");
    for (int i = 0; i != g_nbin; ++i){
        fscanf(aerosol_file_ID, "%Le, %Le, %Le\n", ns+i, dry_r3s+i, kapxrhow0Cs+i);
        ns[g_nbin] += ns[i];
    }
}
void set_aerosol(real8 ra1, real8 ra2, TPecoef tpecoef,
                 real8 *kaps, real8 *dry_r3s, real8 *ns,
                 real8 *kapxrhow0Cs){
// creat bin limit
    real8 *bin_rcut = new real8[g_nbin+1];
    bin_rcut[0] = CRCUT1/g_binq;
    for (int i = 0; i != g_nbin; ++i){
        dry_r3s[i] = ns[i] = 0;
        bin_rcut[i+1] = bin_rcut[i] * g_binq;
    }
    bin_rcut[g_nbin] = CRCUTMAX;  // 3.e-3;
    bin_rcut[0] = CRCUTmin;        // 2.e-10;
  
    real8 N2 = 0, r32 = 0,
          N1 =  ra1/(ra2-1)*powl(bin_rcut[g_nbin],(1-ra2)),
          r31 = ra1/(ra2-4)*powl(bin_rcut[g_nbin],(4-ra2));
// The [g_nbin] store the total number of aerosol
    ns[g_nbin] = 0;
    for (int j = g_nbin-1; j != -1; --j){
        N2 =       ra1/(ra2-1)*powl(bin_rcut[j],(1-ra2));
        r32 = ra1/(ra2-4)*powl(bin_rcut[j],(4-ra2));
        real8 tmp_n = N2 - N1, tmp_r3 = r32 - r31, tmp_r3bar = tmp_r3/tmp_n;
        N1 = N2; r31 = r32;
        if (tmp_n != 0 and tmp_r3!= 0){
            ns[j] = tmp_n;
            dry_r3s[j] = tmp_r3 / tmp_n;
        }
        else
            dry_r3s[j] = ns[j] = 0;
        kapxrhow0Cs[j] = kaps[0]*CrhoH2O0C;
        ns[g_nbin] += ns[j];
    }
    return;
}
void set_aerosol(real8 c, real8 k, real8 rLimit, TPecoef tpecoef,
                 real8 *kaps, real8 *dry_r3s, real8 *ns,
                 real8 *kapxrhow0Cs){
// creat bin limit
    real8 *bin_rcut = new real8[g_nbin+1];
    bin_rcut[0] = CRCUT1/g_binq;
    for (int i = 0; i != g_nbin; ++i){
        dry_r3s[i] = ns[i] = 0;
        bin_rcut[i+1] = bin_rcut[i] * g_binq;
    }
    bin_rcut[g_nbin] = CRCUTMAX;  // 3.e-3;
    bin_rcut[0] = CRCUTmin;        // 2.e-10;

// N = c ssCritical^k
    c *= 2.45e4/CMWDRY;
    const real8 tmp2 = M_4D27*pow(tpecoef.alpha,3)/kaps[0],
                b = c*pow(SS2ETDEW*SS2ETDEW*tmp2,k*.5),
                bm = b*k/(k-2), tmp3 = 1.5*(2-k),
                tmp4 = bm;
    real8 N1 = b*pow(bin_rcut[g_nbin],-1.5*k),
          r31 = tmp4*pow(bin_rcut[g_nbin],tmp3),
          N2 = 0, r32 = 0;
// The [g_nbin] store the total number of aerosol
    ns[g_nbin] = 0;
    for (int j = g_nbin-1; j != -1; --j){
        N2 = b*pow(bin_rcut[j],-1.5*k);
        r32 = tmp4*pow(bin_rcut[j],tmp3);
        real8 tmp_n = N2 - N1, tmp_r3 = r32 - r31;
//        printf("N1 %Le M1 %Le\n", N1, M1);
        N1 = N2; r31 = r32;
        if (tmp_n != 0 && tmp_r3 != 0 && bin_rcut[j] < rLimit ){
            ns[j] = tmp_n;
            dry_r3s[j] = tmp_r3 / tmp_n;
            kapxrhow0Cs[j] = kaps[0]*CrhoH2O0C;
        }
        else
            dry_r3s[j] = ns[j] = kapxrhow0Cs[j] = 0;
        ns[g_nbin] += ns[j];
//        printf("rc %Le rs %Le ns %Le\n", bin_rcut[j], powl(dry_r3s[j], M_1D3), ns[j]);
    }
//  ns[g_nbin] -= ns[g_nbin-1];
//  dry_r3s[g_nbin-1] = ns[g_nbin-1] = 0;
    return;
}
void set_aerosol_to_eq(real8 *dry_r3s, real8 *kapxrhow0Cs,
                       real8 *h2om, TPecoef tpecoef){
    if (tpecoef.S > 1){
        printf("error: S > 1\n");
        exit(-1);
    }
    for (int i = 0; i != g_nbin; ++i){
        if (dry_r3s[i] <= 0) continue;
        real8 alf = tpecoef.alpha, r3 = dry_r3s[i],
              kap = kapxrhow0Cs[i]/tpecoef.rho_water,
              lnS = logl(tpecoef.S), j = 0,
              x1=sqrt(3*kap*r3/alf), x2=x1*1.001, f1, f2, x3,
              r = pow(r3, M_1D3);
        f1 = r*pow(kap/(alf/x1 - lnS) +1, M_1D3)-x1;
        do {
            f2 = r*powl(kap/(alf/x2 - lnS) +1, M_1D3)-x2;
            x3 = x2 - f2*(x2-x1)/(f2-f1);
            x1 = x2; x2 = x3; f1 = f2; ++j;
        } while (fabs(x1/x2 -1) > 1e-10 and j < 30);
        if (j>20) printf("warning (to_eq) iters:%d\n", j);
        h2om[i] = M_4PID3*(pow(x3,3)-r3)*tpecoef.rho_water;
        continue;

        real8 Vs = M_4PID3*dry_r3s[i], logS = logl(tpecoef.S),
              tmp1 = Vs*kap*tpecoef.rho_water, preh2om, rwet;

        h2om[i] = Vs*1e3;
//        printf("S:%Le :", tpecoef.S);
        do {
            preh2om = h2om[i];
            rwet = pow((Vs+preh2om/tpecoef.rho_water)/M_4PID3, M_1D3);
            h2om[i] = tmp1/(tpecoef.alpha/rwet-logS);
        } while (fabs(1-preh2om/h2om[i]) > 1e-6);
/*        Particle par(&tpecoef, Vs, h2om[i]);
        printf("S_eq:%Le, h2om:%Le", par.S_eq, h2om[i]);
        printf("\n"); */
    }
/*
    printf("h2om: ");
    for (int i = 0; i < g_nbin; ++i)
        printf("%Le, ", h2om[i]);
    printf("\n"); */
}


real8 dmi(Particle &particle, real8 dt){
    TPecoef *tpecoefp = particle.tpecoef_pointer;

    real8 tmp1 = tpecoefp->latent_heat_vapor/tpecoefp->K,
          tmp2 = CRH2O*tpecoefp->K/
        (tpecoefp->saturated_vapor_pressure*particle.Dvstar),
          tmp3 = tmp1/particle.kastar*(tmp1/CRH2O-1),
          initial_dS = tpecoefp->S - particle.S_eq,
          dmib = M_4PI*particle.wet_radius*initial_dS/
        (tmp2+tmp3)*dt,
          initial_h2om = particle.h2om;
    if (initial_h2om + dmib <= 0){
        dmib = -initial_h2om*(1-1e-14);
//      printf("dmib makes water mass negtive and has been corrected.\n");
    }
    particle.update_h2om( initial_h2om + dmib);
    real8 dSb = tpecoefp->S - particle.S_eq;
    if (dSb*initial_dS >= 0){
        if (isnan(dmib))
            printf("dmib: %Le\n", dmib);
        return dmib;
    }else{
        real8 alf = tpecoefp->alpha, r3 = particle.saltr3,
              kap = particle.kap,
              lnS = logl(tpecoefp->S), h2om,
              x1=sqrt(3*kap*r3/alf), x2=x1*1.0001, f1, f2, x3,
              r = pow(r3, M_1D3);
        int j = 0;
        f1 = r*pow(kap/(alf/x1 - lnS) +1, M_1D3)-x1;
        do {
            f2 = r*pow(kap/(alf/x2 - lnS) +1, M_1D3)-x2;
            x3 = x2 - f2*(x2-x1)/(f2-f1);
            x1 = x2; x2 = x3; f1 = f2; ++j;
        } while (fabs(x1/x2 -1) > 1e-10);
        h2om = M_4PID3*(pow(x3,3)-r3)*tpecoefp->rho_water;
        particle.update_h2om( h2om );
        return h2om - initial_h2om;
    }
}

real8 fdwv(real8 *dry_r3s, real8 *ns, real8 *kapxrhow0Cs,
          real8 *h2om, real8 *delta_h2om, 
          TPecoef &tpecoef, real8 dt){
    real8 dwv = 0;
#ifdef PARALLEL
#pragma omp parallel for reduction( -:dwv)
#endif
    for (int i = 0; i < g_nbin; ++i){
        if (dry_r3s[i] == 0)
            continue;
        Particle particle(&tpecoef, dry_r3s[i],
            kapxrhow0Cs[i], h2om[i]);
        delta_h2om[i] = dmi(particle, dt);
        dwv -= ns[i] * delta_h2om[i];
        h2om[i] = particle.h2om;
        if (isnan(delta_h2om[i]))
            printf("delta_h2om:%Le\n", delta_h2om[i]);
    }
    return dwv;
}

real8 cal_qw(real8* ns, real8* h2om){
    real8 qw = 0;
    return qw;
#ifdef PARALLEL
#pragma omp parallel for reduction( +:qw)
#endif
    for (int i = 0; i < g_nbin; ++i){
        if (isnan(h2om[i]))
            printf("h2om[%i]:%Le\n", i, h2om[i]);
        if (isnan(ns[i]))
            printf("ns[%i]:%Le\n", i, ns[i]);
        qw += ns[i]*h2om[i];
    }
    if (isnan(qw))
        printf("qw:%Le\n", qw);
    return qw;
}


