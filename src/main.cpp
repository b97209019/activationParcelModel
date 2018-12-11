#include"diffusion_growth_header_file.h"
using namespace std;
int main(int argc, char* argv[]){
    real8 W, W0, aW, dt, t_start, t_stop, print_dt,
          t_snd, p_snd, rh_snd, z_snd, initial_height;
  
    short n_aerosol_mode, io_get_vars;
    char aerosol_file_name[32], aerosol_type[32], kappa_file_name[32];
    read_input_data(W0, aW, dt, t_start, t_stop, print_dt, t_snd,
        p_snd, rh_snd, z_snd, aerosol_type, aerosol_file_name, kappa_file_name);
// setup TPecoef
    real8 K = t_snd +CK0C, p = p_snd*100, S = rh_snd/100.;
    TPecoef tpecoef(K,p,S);
    tpecoef.parcel_descend(W0, aW, initial_height);
// setup aerosol spectrum
    real8 *ns = new real8[g_nbin+1], *dry_r3s = new real8[g_nbin],
          *h2om = new real8[g_nbin], *kapxrhow0Cs = new real8[g_nbin],
          Ns[3], SIGMAs[3], MUs[3], c, k, rLimit, kaps[3];
    read_kappa(kappa_file_name, kaps);
    if (strcmp(aerosol_type, "lognormal") == 0){
        read_lognormal_distribution(aerosol_file_name, n_aerosol_mode, Ns, SIGMAs, MUs);
        set_aerosol(n_aerosol_mode, Ns, SIGMAs, MUs, kaps, dry_r3s, ns, kapxrhow0Cs);
    } else if (strcmp(aerosol_type, "powerlaw") == 0){
        read_powerlaw_distribution(aerosol_file_name, c, k, rLimit);
        if (rLimit > 0)
            set_aerosol(c, k, rLimit, tpecoef, kaps, dry_r3s, ns, kapxrhow0Cs);
        else{
            real8 ra1 = c, ra2 = k;
            set_aerosol(ra1, ra2, tpecoef, kaps, dry_r3s, ns, kapxrhow0Cs);
        }
    } else if (strcmp(aerosol_type, "fileDefined") == 0){
        set_aerosol(aerosol_file_name, dry_r3s, ns, kapxrhow0Cs);
    } else{
        printf("Error: unknown aerosol_type!\n");
        exit(-3);
    }
    set_aerosol_to_eq(dry_r3s, kapxrhow0Cs, h2om, tpecoef);
  
    real8 preS = -1;
    real8 print_indicator = 0, calculate_time = 0,
          calculate_fdwv_time = 0, total_wtime, fdwv_wtime;
    char ambient_file_name[128], particle_file_name[128],
         postfix[37];
    sprintf(postfix, "%s_T%02iP%02iW%010.7Lfa%010.7Lf",
            aerosol_file_name, int(t_snd),
            int(p_snd/100), W0, aW);
    sprintf(ambient_file_name, "../output/ambient_"
            "%15s.csv\0", postfix);
    sprintf(particle_file_name, "../output/particle_"
            "%15s.csv\0", postfix);
    create_files(dry_r3s, ns, kapxrhow0Cs, ambient_file_name,
                 particle_file_name);
    bool Smax_reached = false;
    int fine_output_indicator = 0;
    real8 dwv_k[4], half_dt = dt*.5, T_k[4], lnp_k[4],
          qw, *h2om_k1 = new real8[g_nbin],
          *h2om_k2 = new real8[g_nbin],
          *h2om_k3 = new real8[g_nbin], *h2om_k4 = new real8[g_nbin],
          *h2om_tmp = new real8[g_nbin], time_Smax, lastPNc = 0;
    TPecoef
          tpecoef_tmp(tpecoef.K, tpecoef.p, tpecoef.S);
    for (real8 model_time = t_start; model_time < t_stop;
          model_time += dt){
        if (aW < 0)
            W = W*(1- expl(aW*model_time));
        else if (aW >= 0)
            W = W0 + aW*model_time;
        total_wtime = omp_get_wtime();
        tpecoef_tmp.K = tpecoef.K;
        tpecoef_tmp.p = tpecoef.p;
        tpecoef_tmp.wv = tpecoef.wv;
        tpecoef_tmp.update_TPe("wv");
#ifdef PARALLEL
#pragma omp parallel for
#endif
        for (int j = 0; j < g_nbin; ++j){
            h2om_tmp[j] = h2om[j];
        }
        fdwv_wtime = omp_get_wtime();
        dwv_k[0] =
            fdwv(dry_r3s, ns, kapxrhow0Cs, h2om_tmp,
                 h2om_k1, tpecoef_tmp, half_dt);
        calculate_fdwv_time += omp_get_wtime() - fdwv_wtime;
        qw = cal_qw(ns, h2om_tmp);
        tpecoef_tmp.update_TPe(
            dwv_k[0], T_k, lnp_k, W, half_dt, qw);
    
        fdwv_wtime = omp_get_wtime();
        dwv_k[1] =
            fdwv(dry_r3s, ns, kapxrhow0Cs, h2om_tmp,
                 h2om_k2, tpecoef_tmp, half_dt);
        calculate_fdwv_time += omp_get_wtime() - fdwv_wtime;
        qw = cal_qw(ns, h2om_tmp);
        tpecoef_tmp.update_TPe(
            dwv_k[1], T_k+1, lnp_k+1, W, half_dt, qw);
     
        tpecoef_tmp.K = tpecoef.K+T_k[1];
        tpecoef_tmp.p = tpecoef.p*expl(lnp_k[1]);
        tpecoef_tmp.wv = tpecoef.wv+dwv_k[1];
        tpecoef_tmp.update_TPe("wv");
#ifdef PARALLEL
#pragma omp parallel for
#endif
        for (int j = 0; j < g_nbin; ++j)
          h2om_tmp[j] = h2om[j] + h2om_k2[j];
        fdwv_wtime = omp_get_wtime();
        dwv_k[2] =
            fdwv(dry_r3s, ns, kapxrhow0Cs, h2om_tmp,
                 h2om_k3, tpecoef_tmp, dt);
        calculate_fdwv_time += omp_get_wtime() - fdwv_wtime;
        qw = cal_qw(ns, h2om_tmp);
        tpecoef_tmp.update_TPe(
            dwv_k[2], T_k+2, lnp_k+2, W, dt, qw);
    
        tpecoef_tmp.K = tpecoef.K+T_k[2];
        tpecoef_tmp.p = tpecoef.p*expl(lnp_k[2]);
        tpecoef_tmp.wv = tpecoef.wv+dwv_k[2];
        tpecoef_tmp.update_TPe("wv");
#ifdef PARALLEL
#pragma omp parallel for
#endif
        for (int j = 0; j < g_nbin; ++j)
          h2om_tmp[j] = h2om[j] + h2om_k3[j];
        fdwv_wtime = omp_get_wtime();
        dwv_k[3] =
            fdwv(dry_r3s, ns, kapxrhow0Cs, h2om_tmp,
                 h2om_k4, tpecoef_tmp, half_dt);
        calculate_fdwv_time += omp_get_wtime() - fdwv_wtime;
        qw = cal_qw(ns, h2om_tmp);
        tpecoef_tmp.update_TPe(
            dwv_k[3], T_k+3, lnp_k+3, W, half_dt, qw);
        calculate_time += omp_get_wtime() - total_wtime;
// pre_ for previous time step variable
        real8 *&pre_h2om = h2om_tmp;
        TPecoef &pre_tpecoef = tpecoef_tmp;
#ifdef PARALLEL
#pragma omp parallel for
#endif
        for (int j = 0; j < g_nbin; ++j){
            pre_h2om[j] = h2om[j];
            real8 dh2om = (h2om_k1[j] + 2*h2om_k2[j] +
                h2om_k3[j] + h2om_k4[j])*M_1D3;
            if (h2om[j]+dh2om < 0){
                printf("h2omm< 0 %15Le h2omm%15Le salt %15Le\n",
                    dh2om, h2om[j], dry_r3s[j]);
            }
            h2om[j] += dh2om;
        }
        real8 dwv = (dwv_k[0] + 2*dwv_k[1] +
              dwv_k[2] + dwv_k[3])*M_1D3;
        pre_tpecoef.K = tpecoef.K;
        pre_tpecoef.p = tpecoef.p;
        pre_tpecoef.wv = tpecoef.wv;
        pre_tpecoef.update_TPe("wv");
        real8 dwvdt = dwv/dt;
        tpecoef.K += (T_k[0]+2*T_k[1]+T_k[2]+T_k[3])*M_1D3;
        tpecoef.p *= expl((lnp_k[0]+2*lnp_k[1]+
                      lnp_k[2]+lnp_k[3])*M_1D3);
        tpecoef.wv += dwv;
        tpecoef.update_TPe("wv");
/*
    if (Smax_reached && pre_tpecoef.S > 1 && pre_tpecoef.S < tpecoef.S){
        printf("S climb again\n");
        exit(0);
    }
*/
        bool S_first_descend= false;
// check if had reached Smax at previous time
        if (!Smax_reached && pre_tpecoef.S > 1 && pre_tpecoef.S > tpecoef.S){
            S_first_descend = true;
            Smax_reached = true;
//            fine_output_indicator = 120+print_indicator-int(model_time);
            print_indicator = int(model_time) +1;
        }
// output at the time reached Smax
// caution! calculated_fdwv_time is not the
// previous time consumed and it does
// not cover the whole print-out interval.
        if (S_first_descend){
            time_Smax = model_time-dt;
            real8 parcel_height = initial_height + (aW < 0? 
                  W0*(time_Smax-exp(aW*time_Smax)/aW):
                  (W0 + .5*aW*time_Smax)*time_Smax), outNc;
            output_data(
                pre_tpecoef, dry_r3s, kapxrhow0Cs, pre_h2om,
                ns, dwvdt, parcel_height, W, 
                calculate_fdwv_time, time_Smax, outNc,
                ambient_file_name, particle_file_name);
// output to other file after Smax
            sprintf(ambient_file_name, "../output/ambient_"
                "%15s_postSmax.csv\0", postfix);
            sprintf(particle_file_name, "../output/particle_"
                "%15s_postSmax.csv\0", postfix);
            create_files(dry_r3s, ns, kapxrhow0Cs,
                         ambient_file_name,
                         particle_file_name);
// change t_stop to model_time+200s
            t_stop = model_time + 200;
//        t_stop = model_time + 3600;
        }
        real8 parcel_height = initial_height + (aW < 0? 
              W0*(model_time-exp(aW*model_time)/aW):
              (W0 + .5*aW*model_time)*model_time);
//    if (Smax_reached && parcel_height > 2000)
//      exit(0);
// failsafe
        if (tpecoef.K  < 200 or isnan(tpecoef.K)){
            printf("error: temperature abnormal\n");
            exit(-1);
        }
// normal output every print_dt
        if (model_time >= print_indicator){
            real8 outNc;
            output_data(
                  tpecoef, dry_r3s, kapxrhow0Cs, h2om, ns,
                  dwvdt, parcel_height, W, calculate_fdwv_time,
                  model_time, outNc, ambient_file_name,
                  particle_file_name);
            if (outNc > lastPNc || outNc == 0 || !Smax_reached ){
                t_stop = model_time + 600;
//              t_stop = model_time + 1600;
                lastPNc = outNc;
            }
            calculate_time = 0; 
            calculate_fdwv_time = 0; 
            if (fine_output_indicator > 0){
                print_indicator += 1;
                --fine_output_indicator;
            } else print_indicator += print_dt;
        }
    }
}
