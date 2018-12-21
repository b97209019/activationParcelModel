#include"diffusion_growth_header_file.h"
real8 g_ac, g_binq, g_aT;
int g_nbin;
void read_input_data(real8 &W0, real8 &aW, real8 &dt,
     real8 &t_start, real8 &t_stop, real8 &print_dt,
     real8 &t_snd, real8 &p_snd, real8 &rh_snd,
     real8 &z_snd, char *aerosol_type, char *aerosol_file_name,
     char *kappa_file_name){
    char input_file_name[] = "../input/input.data",
         line[64];
    FILE *input_file_ID = fopen(input_file_name, "r");
    fgets(line, 64, input_file_ID);
    short io_get_vars =
        fscanf(input_file_ID, "%Lf, %Lf, %Lg, %Lf, %Lf, %Lf\n",
        &W0, &aW, &dt, &t_start, &t_stop, &print_dt);
    fgets(line, 64, input_file_ID);
    io_get_vars +=
        fscanf(input_file_ID, "%Lf %Lf %Lf %Lf %Lf %Lf\n",
        &p_snd, &t_snd, &rh_snd, &z_snd, &g_ac, &g_aT);
    fgets(line, 64, input_file_ID);
    io_get_vars +=
        fscanf(input_file_ID, "%s %s %i\n", aerosol_type, 
               aerosol_file_name, &g_nbin);
    fgets(line, 64, input_file_ID);
    io_get_vars +=
        fscanf(input_file_ID, "%s\n", kappa_file_name);
    fclose(input_file_ID);
    if( io_get_vars != 16 ){
        printf( "Fail to read file: %s\n",
            input_file_name );
        exit(-1);
    }
    g_binq = pow(2,50./(3.*g_nbin));
    return;
}
void read_lognormal_distribution(char *aerosol_file_name,
     short &n_aerosol_mode, real8 *Ns, real8 *SIGMAs,
     real8 *MUs){
    char line[128], aerosol_full_file_name[64];
    strcpy(aerosol_full_file_name, "../input/");
    strcat(aerosol_full_file_name, aerosol_file_name);
    FILE *aerosol_file_ID = fopen(aerosol_full_file_name, "r");
    fgets(line, 128, aerosol_file_ID);
    fgets(line, 128, aerosol_file_ID);
    short io_get_vars = sscanf(line, "%i", &n_aerosol_mode);
    if (n_aerosol_mode > 3){
        printf( "Error: more than trimode\n");
        exit(-1);
    }
    for (int i = 0; i != n_aerosol_mode; ++i){
        fgets(line, 128, aerosol_file_ID);
        io_get_vars +=
            sscanf(line, "%Lg %Lg %Lg %*s\n",
    	Ns+i, SIGMAs+i, MUs+i);
    }
    fclose(aerosol_file_ID);
  
    if( io_get_vars != 1+3*n_aerosol_mode ){
        printf( "Fail to read file: %s\n",
            aerosol_file_name );
        exit(-1);
    }
    return;
}
void read_kappa(char *kappa_file_name, real8 *kaps){
    char line[128], kappa_full_file_name[64];
    strcpy(kappa_full_file_name, "../input/");
    strcat(kappa_full_file_name, kappa_file_name);
    FILE *kappa_file_ID = fopen(kappa_full_file_name, "r");
    fgets(line, 128, kappa_file_ID);
    fgets(line, 128, kappa_file_ID);
    char given_kappa[10];
    short io_get_vars = sscanf(line, "%s", given_kappa);
    if (strcmp(given_kappa, "False") == 0){
        for (int i = 0; i != 3; ++i)
            kaps[i] = Cnu*CrhoSALTDRY*CMWH2O/(CrhoH2O0C*CMWSALT);
        return;
    } else if (strcmp(given_kappa, "True") != 0){
        printf("error found in kappa.txt: unknown keyword for option is_prescribed.\n");
        exit(-1);
    }
    fgets(line, 128, kappa_file_ID);
    int n_aerosol_mode;
    io_get_vars += sscanf(line, "%i", &n_aerosol_mode);
    if (n_aerosol_mode > 3){
        printf( "Error: more than trimode\n");
        exit(-1);
    }
    for (int i = 0; i != n_aerosol_mode; ++i){
        fgets(line, 128, kappa_file_ID);
        io_get_vars += sscanf(line, "%Lg\n", kaps+i);
    }
    fclose(kappa_file_ID);
  
    if( io_get_vars != 2+n_aerosol_mode ){
        printf( "Fail to read file: %s\n",
            kappa_file_name );
        exit(-1);
    }
    return;
}
void read_powerlaw_distribution(
    char *aerosol_file_name, real8 &c, real8 &k, real8 &rLimit){
    char line[128], aerosol_full_file_name[64];
    strcpy(aerosol_full_file_name, "../input/");
    strcat(aerosol_full_file_name, aerosol_file_name);
    FILE *aerosol_file_ID = fopen(aerosol_full_file_name, "r");
    fgets(line, 128, aerosol_file_ID);
    fgets(line, 128, aerosol_file_ID);
    short io_get_vars = sscanf(line, "%Lg %Lg %Lg\n", &c, &k, &rLimit);
    if( io_get_vars != 3){
        printf( "Fail to read file: %s\n",
            aerosol_file_name );
        exit(-1);
    }
    return;
}

void create_files(real8 *dry_r3s, real8 *ns,
                  real8 *kapxrhow0Cs,
                  char *ambient_file_name,
                  char *particle_file_name){
    FILE *ambient_file_ID = fopen(ambient_file_name, "w");
    fprintf(ambient_file_ID, "%-17s, %-17s, %-17s, %-17s, "
        "%-17s, %-17s, %-17s, %-17s, %-17s, %-17s, %-17s, "
        "%-17s, %-17s\n", "time", "K", "p", "S", "height", "wv",
        "dwvdt", "wc", "A1", "W", "A2", "activated_n",
        "fdwv_calc_t");
    fclose(ambient_file_ID);
  
    FILE *particle_file_ID = fopen(particle_file_name,"w");
    fprintf(particle_file_ID, "%-10s, %-17s, ", "variable",
        "time\\bin");
    for (int i = 0; i != g_nbin-1; ++i)
        fprintf(particle_file_ID, "%17i, ", i+1);
    fprintf(particle_file_ID, "%17i\n", g_nbin);
    output_data_1variable(particle_file_ID,
        "saltr3", dry_r3s, 0L);
    output_data_1variable(particle_file_ID,
        "n", ns, 0L);
    output_data_1variable(particle_file_ID,
        "kapxrhow0C", kapxrhow0Cs, 0L);
  
    fclose(particle_file_ID);
    return;
}
void output_data(TPecoef tpecoef, real8 *dry_r3s,
    real8 *kapxrhow0Cs,
    real8 *h2om, real8 *ns, real8 dwvdt,
    real8 height, real8 W, real8 calculate_fdwv_time,
    real8 model_time, real8 &outNc,
    char *ambient_file_name, char *particle_file_name){
// create some variables about S* r*
    real8 *Kohler_critical_S_eqs = new real8[g_nbin],
          *Kohler_critical_h2oms = new real8[g_nbin],
          *wet_radiuss = new real8[g_nbin], *S_eqs = new real8[g_nbin],
          wc = 0;
#ifdef PARALLEL
#pragma omp parallel for
#endif
    for (int i = 0; i < g_nbin; ++i){
        Particle particle(&tpecoef, dry_r3s[i], kapxrhow0Cs[i],
                          h2om[i], true);
        Kohler_critical_S_eqs[i] = particle.Kohler_critical_S_eq;
        Kohler_critical_h2oms[i] = particle.Kohler_critical_h2om;
        wet_radiuss[i] =  particle.wet_radius;
        S_eqs[i] = particle.S_eq;
        if (tpecoef.S/S_eqs[i] -1 > 1e-1)
            wc += h2om[i]*ns[i];
    }
    FILE *particle_file_ID = fopen(particle_file_name,
        "a");
    output_data_1variable(particle_file_ID,
        "h2om", h2om, model_time);
    output_data_1variable(particle_file_ID,
        "r_wet", wet_radiuss, model_time);
    output_data_1variable(particle_file_ID,
        "S_eq", S_eqs, model_time);
    output_data_1variable(particle_file_ID,
        "S_eq*", Kohler_critical_S_eqs, model_time);
    output_data_1variable(particle_file_ID,
        "h2om*", Kohler_critical_h2oms, model_time);
    fclose(particle_file_ID);
    real8 not_activated_number = 0;
    for (int i = 0; i != g_nbin; ++i){
        if (h2om[i] < Kohler_critical_h2oms[i])
            not_activated_number += ns[i];
        else if (ns[i] < 1e-30)
            continue;
//    else break;
    }
    FILE *ambient_file_ID = fopen(ambient_file_name, "a");
    real8 cp_i = CCpDRY+CCpH2O*tpecoef.wv,
          A1 = (CE*tpecoef.latent_heat_vapor/
        (cp_i*tpecoef.K) - (1+tpecoef.wv)/(1+tpecoef.wv/CE))
        *CG/(CRDRY*tpecoef.K),
          pme = tpecoef.p - tpecoef.e,
          A2 = (pme)*(pme)/(CE*tpecoef.p*tpecoef.e) +
        tpecoef.latent_heat_vapor*
        tpecoef.latent_heat_vapor/
        (cp_i*CRH2O*tpecoef.K*tpecoef.K);
    outNc = ns[g_nbin] - not_activated_number;
    fprintf(ambient_file_ID, ONE_DIGIT ONE_DIGIT
        ONE_DIGIT ONE_DIGIT ONE_DIGIT ONE_DIGIT
        ONE_DIGIT ONE_DIGIT ONE_DIGIT ONE_DIGIT
        ONE_DIGIT ONE_DIGIT ONE_DIGIT_N , model_time,
        tpecoef.K, tpecoef.p, tpecoef.S, height,
        tpecoef.wv, dwvdt, wc, A1, W, A2, outNc,
        calculate_fdwv_time);
    fclose(ambient_file_ID);
    return;
}
void output_data_1variable(FILE *particle_file_ID,
    const char variable_name[], real8 *particle_variable,
    real8 model_time){
    fprintf(particle_file_ID, "%-10s, " ONE_DIGIT,
        variable_name, model_time);
    for (int i = 0; i != g_nbin-1; ++i){
        fprintf(particle_file_ID, ONE_DIGIT,
            particle_variable[i]);
    }
    fprintf(particle_file_ID, ONE_DIGIT_N,
        particle_variable[g_nbin-1]);
    return;
}


