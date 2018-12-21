#include"diffusion_growth_header_file.h"
TPecoef::TPECOEF(real8 K, real8 p, real8 S){
    this->K = K; this->p = p; this->S = S;
    this->update_TPe("S");
}
void TPecoef::update_TPe(const string &s){
    real8 K = this->K, p = this->p, C = K - CK0C,
          X = (C-35)*(C-35);
    this->rho_water = C >= 0.?
        (999.84+C*(18.225+C*(-7.9222E-3+C*(-5.5448E-5+
        C*(1.4976E-7-C*3.933E-10)))))/(1.+0.01816*C)
        : 999.84+C*(0.086-C*0.0108);
    this->cw = C > 0? 4178 + ( 1.298e-2 + 1.591e-5*X )*X :
          4218. + 0.3471 * C*C;
  
    this->saturated_vapor_pressure =  K > 223.15?
        610.78+C*(44.36518521+C*(1.428945805+C*(2.65064847E-2+
    C*(3.031240396E-4+C*(2.034080948E-6+C*6.136820929E-9)))))
        : 6.3374*expl(25.584*(K-223.15)/K);
    if ( s.compare("S") == 0)
        update_withS();
    else if ( s.compare("wv") == 0)
        update_withwv();
    real8 S = this->S;
    this->rho_air = (p+(CE-1)*this->e)/(CRDRY*K);
    this->latent_heat_vapor = 
        2.501E6 * powl( CK0C/K , 0.167+3.67E-4*K );
    this->latent_heat_sublim = 
        2.6332E+6 + K*(1727. - 3.625*K);
    this->normaled_temperature = K/CK0C;
    this->normaled_pressure = 101325/p;
    this->dynamic_viscosity = 1.718e-5 + 
        (C > 0 ?  4.9e-8*C : C*(4.9e-8-1.2e-10*C));
    this->mean_free_path = 6.62e-8*normaled_pressure*
        sqrtl(this->normaled_temperature)*
        dynamic_viscosity/1.718e-5;
    this->Dv = 2.11e-5*powl(this->normaled_temperature,
               1.94)* normaled_pressure;
    this->ka = 2.3823e-2+7.11756e-5*C;
    this->surface_tension = 0.0761 - 1.55E-4 * C;
    this->alpha = 2*this->surface_tension/(CRH2O*this->rho_water*K);
    this->thermal_velocity_vapor = sqrtl(CRH2O*K/(M_PI+M_PI));
    this->thermal_velocity_air= sqrtl(CRDRY*K/(M_PI+M_PI));
    real8 cp_i = (CCpDRY + CCpH2O* this->wv),
          tmp1 = 1/(this->wv*(1+this->wv/CE)),
          tmp2 = this->latent_heat_vapor/(cp_i*CRH2O*K*K);
    this->A1 = (CE*this->latent_heat_vapor/(cp_i*K)
               -(1+this->wv)/(1+this->wv/CE))*CG/(CRDRY*K);
    this->A2 = tmp1 - tmp2*this->latent_heat_vapor;
    this->A3 = tmp1 - tmp2*this->latent_heat_sublim;
}
void TPecoef::update_withS(){
    this->e = S*this->saturated_vapor_pressure;
    this->xv = this->e / (p - this->e);
    this->wv = this->xv *CE;
    return;
}
void TPecoef::update_withwv(){
    this->xv = this->wv/CE;
    this->e = this->p*this->xv/(1 + this->xv);
    this->S = this->e/this->saturated_vapor_pressure;
    return;
}
void TPecoef::update_TPe(real8 dwv, real8 *T_k,
    real8 *lnp_k, real8 W, real8 dt, real8 qw/* = 0*/){
    real8 mcgwdt = -CG*W*dt;
    *lnp_k = CE*(1+this->wv)/
        ((CE+this->wv)*CRDRY*this->K)*mcgwdt;
    if (qw == 0) {
        real8 cp_i = (CCpDRY + CCpH2O* this->wv);
        *T_k = mcgwdt/cp_i - this->latent_heat_vapor
            /CCpH2O*logl(1+CCpH2O*dwv/cp_i);
    }else{
        real8 cp_i = (CCpDRY + CCpH2O* this->wv+this->cw*qw);
        *T_k = (mcgwdt - this->latent_heat_vapor
            *dwv)/cp_i;
    }
    if (isnan(*T_k)){
        printf("dT:%Le, dwv:%Le\n", *T_k, dwv);
        printf("lv:%Le, cw:%Le\n", this->latent_heat_vapor, this->cw);
        printf("qw:%Le, wv:%Le\n", qw, this->wv);
    }
    if (this->K + *T_k < 0)
        printf("dT:%Le\n", T_k);
    this->K += *T_k;
    this->p *= expl(*lnp_k);
    this->wv += dwv;
    this->update_TPe("wv");
}

void TPecoef::parcel_descend(real8 W0, real8 aW, real8 &height){
  real8 K0 = this->K, cp = (CCpDRY+CCpH2O*this->wv),
//          a = cp*(CE+CE*this->wv)/((CE+this->wv)*CRDRY),
          a = cp*(1+this->wv)/((CE+this->wv)*CRH2O),
          W = W0 + 30*aW;
//    height = 100./ W < 120. ? -120.*W : -100.;
    if ( aW <= 0 )
        height = -100.*W0;
    else if ( aW > 0 )
        height = -100.*W;
      
    this->K -= CG*height/cp;
    this->p *= powl((this->K/K0),a);
    this->update_TPe("wv");
    return;
}

Particle::PARTICLE( TPecoef *tpecoefp, real8 saltr3, real8 kapxrhow0C, real8 h2om,
    bool calculate_Kohler_critical/* = false*/){
    this->tpecoef_pointer = tpecoefp;
    this->saltr3 = saltr3;
    this->kap = kapxrhow0C/tpecoefp->rho_water;
    this->salt_apparent_volume = M_4PID3*saltr3;
    if (h2om < 0) printf("warning h2om less than zero\n");
    if (calculate_Kohler_critical){
        real8 alf = tpecoefp->alpha, rwet;
        this->Kohler_critical_S_eq =
        S_crit(this->salt_apparent_volume, alf, this->kap, rwet);
        this->Kohler_critical_h2om =
            (M_4PID3*pow(rwet, 3) - this->salt_apparent_volume)*
            tpecoefp->rho_water;
    }
    this->update_h2om(h2om);
}
void Particle::update_h2om(real8 h2om){
    TPecoef *tpecoefp = this->tpecoef_pointer;
    this->h2om = h2om > 0? h2om : 0;
    this->wet_radius = cbrtl(
        (this->salt_apparent_volume
        + this->h2om / tpecoefp->rho_water)
        /M_4PID3 ); 
    real8 r_lambda = this->wet_radius+tpecoefp->mean_free_path,
    q1 = this->wet_radius/r_lambda ,
    q2 = this->wet_radius*tpecoefp->thermal_velocity_vapor,
    q3 = this->wet_radius*tpecoefp->thermal_velocity_air;
    this->Dvstar = tpecoefp->Dv/(q1+tpecoefp->Dv/(g_ac*q2));
    this->kastar = tpecoefp->ka/(q1+tpecoefp->ka/
           (tpecoefp->rho_air*CCpDRY*g_aT*q3));
    this->S_eq = expl(
           tpecoefp->alpha/this->wet_radius
           - this->kap*M_4PID3*saltr3*tpecoefp->rho_water/(this->h2om));
}

real8 S_crit(real8 Vdry, real8 alf, real8 kap, real8 &a){
    real8 d = -Vdry/M_4PID3, p = d*kap/alf, b = -pow(-3*p, .5),
	  q = 2*pow(b,3)/27 + d, tmp1 = pow(-p/3, .5),
	  cos3th = 3*q/(2*p*tmp1);
    a = 2*tmp1*cosh( acosh(cos3th)/3 ) - b/3.;
    return exp(alf/a - kap/( pow(a, 3)/(-d) - 1 ));
}

real8 rcut(real8 lnS, real8 alf, real8 kap, real8 &a){
    real8 tmp1 = pow(1+3*lnS/kap, .5);
    a = alf*(3*lnS-kap*(1+tmp1))/(3*lnS*(lnS-kap));
    real8 tmp2 = alf-lnS*a;
    return a*pow( tmp2/(a*kap+tmp2), M_1D3);
}
