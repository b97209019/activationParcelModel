#ifndef __DIFFUSION_GROWTH_HEADER_FILE_H__
#define __DIFFUSION_GROWTH_HEADER_FILE_H__
#define ONE_DIGIT   "%17.10Le, " 
#define ONE_DIGIT_N "%17.10Le\n" 
#include"micro_physics_const.h"
#include<string>
#include<iostream>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<omp.h>
using namespace std;

extern int g_nbin;
extern real8 g_binq;
extern real8 g_ac;
typedef struct TPECOEF{
    public:
        real8 K, p, S, rho_water, cw, saturated_vapor_pressure,
              latent_heat_vapor, normaled_temperature,
              normaled_pressure, dynamic_viscosity,
              mean_free_path, Dv, ka, surface_tension, alpha, 
              thermal_velocity_vapor, e, wv, xv,
              thermal_velocity_air, rho_air,
              latent_heat_sublim,
              A1, A2, A3, A;
        TPECOEF(){}
        TPECOEF(real8, real8, real8);
        void update_TPe(const string &);
        void update_TPe(real8, real8 *, real8 *,
                        real8, real8, real8 qw = 0);
        void update_withS();
        void update_withwv();
        void parcel_descend(real8, real8, real8 &);
} TPecoef;

typedef struct PARTICLE{
    public:
        TPecoef *tpecoef_pointer;
        real8 Dvstar, kastar, h2om, saltr3, kap,
              salt_apparent_volume, S_eq, wet_radius,
              Kohler_critical_h2om, Kohler_critical_S_eq;
        PARTICLE(TPecoef *tpecoef, real8 saltr3, real8 kapxrhow0C,
                 real8 h2om, bool calculate_Kohler_critical = false);
        void update_h2om(real8 h2om);
} Particle;
 
void set_aerosol(real8 n_aerosol_mode, real8 *Nts,
                 real8 *SIGMAs, real8 *MUs, real8 *kaps,
                 real8 *dry_r3s, real8 *ns, real8 *kapxrhow0Cs);
void set_aerosol(char *aerosol_spectrum_file_name,
                 real8 *dry_r3s, real8 *ns, real8 *kapxrhow0Cs);
void set_aerosol(real8 ra1, real8 ra2, TPecoef tpecoef,
                 real8 *kaps, real8 *dry_r3s, real8 *ns,
                 real8 *kapxrhow0Cs);
void set_aerosol(real8 c, real8 k, real8 rLimit, TPecoef tpecoef,
                 real8 *kaps, real8 *dry_r3s, real8 *ns,
                 real8 *kapxrhow0Cs);


void set_aerosol_to_eq(real8 *, real8 *, real8 *, TPecoef);
real8 dmi(Particle &, real8);
real8 fdwv(real8 *, real8 *, real8 *, real8 *, real8 *, 
           TPecoef &, real8);
real8 cal_qw(real8 *, real8 *);
void read_input_data(real8 &, real8 &, real8 &, real8 &,
                     real8 &, real8 &, real8 &, real8 &,
                     real8 &, real8 &, char *, char *, char *);
void read_lognormal_distribution(char *aerosol_file_name,
     short &n_aerosol_mode, real8 *Ns, real8 *SIGMAs, real8 *MUs);
void read_kappa(char *kappa_file_name, real8 *kap0Cs);
void read_powerlaw_distribution(char *aerosol_file_name, real8 &c, real8 &k, real8 &rLimit);
void create_files(real8 *, real8 *, real8 *, char *, char *);
void output_data(TPecoef, real8 *, real8 *, real8 *, real8 *,
                 real8, real8, real8, real8, real8, real8 &,
                 char *, char *);
void output_data_1variable(FILE *particle_file_ID,
                           const char variable_name[],
                           real8 *particle_variable,
                           real8 model_time);
real8 S_crit(real8 Vdry, real8 alf, real8 kap, real8 &a);
real8 rcut(real8 lnS, real8 alf, real8 kap, real8 &a);
#endif
