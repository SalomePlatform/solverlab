//////////////////////////////////////////////////////////////////////////////
//
// File:        fluid.h
// Directory:   
// Version:     
//
//////////////////////////////////////////////////////////////////////////////

#ifndef fluid_included
#define fluid_included

#include <math.h>
#include <Noms.h>
#include <SFichier.h>
#include <Process.h>

class fluid
{
public :
    fluid() { fluide = 0; };
    void set_fluid(Nom &nom)
    {
        Noms liste(2); liste(0) = "Na"; liste(1) = "LBe"; 
        fluide = liste.rang(nom);
        if (fluide == -1)
        {
            if (!Process::me()) Cerr << "fluide " << nom << " non reconnu parmi : " << liste << finl;
            Process::exit();
        }
    }

    //proprietes physiques
    /* Lois physiques du sodium issues de fits sur DEN/DANS/DM2S/STMF/LMES/RT/12-018/A */
    int fluide;
    #define Tk ((T) + 273.15)
    #define Ksil (fluide ? 3.022e-11 : 1e-9)

    /* inverse de la densite du liquide (hors pression) */
    inline double IRhoL(double T)
    {
      return fluide ? (9.03e-5 + Tk * (1.003e-8 + Tk * 2.01e-12))
                    : (9.829728e-4 + Tk * (2.641186e-7 + Tk * (-3.340743e-11 + Tk * 6.680973e-14)));
    }

    /* sa derivee */
    inline double DTIRhoL(double T)
    {
      return fluide ? (1.003e-8 + 2.01e-12 * 2)
                    : (2.641186e-7 + Tk * (-3.340743e-11 * 2 + Tk * 6.680973e-14 * 3));
    }

    /* densite du liquide */
    inline double RhoL(double T, double P)
    {
      return exp(Ksil * (P - 1e5)) / IRhoL(T);
    }

    /* enthalpie du liquide */
    inline double HL(double T, double P)
    {
      return (fluide ? 97.98e3 + Tk * (159 + Tk * (-2.72e-2 / 2 + Tk * 7.12e-6 / 3))
                     : 2992600 / Tk - 365487.2 + Tk * (1658.2 + Tk * (-.42395 + Tk * 1.4847e-4)))
          + (IRhoL(T) - Tk * DTIRhoL(T)) * (1 - exp(Ksil*(1e5 - P))) / Ksil;
    }
    
    /* derivees */
    inline double DTHL(double T, double P)
    {
      return (fluide ? 159 + Tk * (-2.72e-2 + Tk * 7.12e-6)
                     : -2992600 / (Tk*Tk) + 1658.2 + Tk * (-.42395 * 2 + Tk * 1.4847e-4 * 3))
             - Tk * (fluide ? 0 : -3.340743e-11 * 2 + Tk * 6.680973e-14 * 6) * (1 - exp(Ksil*(1e5 - P))) / Ksil;
    }
    
    /* inverses par methode de Newton */
    inline double TLh(double h, double T0, double P)
    {
      double T = T0, sec; int i;
      for (i = 0; i < 100 && dabs( sec = HL(T, P) - h ) > 1e-8; i++)
        T -= sec / DTHL(T, P);
      return i < 100 ? T : -1e10;
    }
    /* derivees */
    inline double DTRhoL(double T, double P)
    {
      return - exp(Ksil * (P - 1e5)) * DTIRhoL(T) / (IRhoL(T) * IRhoL(T));
    }
    inline double DPRhoL(double T, double P)
    {
        return Ksil * RhoL(T, P);
    }
    
    /* derivees */
    inline double DPHL(double T, double P)
    {
        return (IRhoL(T) - Tk * DTIRhoL(T)) * exp(Ksil * (1e5 - P));
    }
    
    /* energie volumique du liquide */
    inline double EL(double T, double P)
    {
      return HL(T, P) * RhoL(T, P);
    }
    /* derivees */
    inline double DTEL(double T, double P)
    {
      return DTHL(T, P) * RhoL(T, P) + HL(T, P) * DTRhoL(T, P);
    }
    inline double DPEL(double T, double P)
    {
      return DPHL(T, P) * RhoL(T, P) + HL(T, P) * DPRhoL(T, P);
    }
    
    inline double TLe(double e, double T0, double P)
    {
      double T = T0, sec; int i;
      for (i = 0; i < 100 && dabs( sec = EL(T, P) - e ) > 1e-3; i++)
        T -= sec / DTEL(T, P);
      return i < 100 ? T : -1e10;
    }
    
    /* viscosite du liquide*/
    inline double MuL( double T )
    {
      return fluide ? 4.94e-4 * exp(754.1 / Tk)
                    : exp( -6.4406 - 0.3958 * log(Tk) + 556.835 / Tk );
    }
    /* conductivite du liquide */
    inline double LambdaL( double T )
    {
      return fluide ? 3.61 + Tk * (1.517e-2 - Tk * 1.741e-6)
                    : max(124.67 + Tk * (-.11381 + Tk * (5.5226e-5 - Tk * 1.1848e-8)), 0.045);
    }
    /* tension superficielle */
    inline double SigmaL(double T)
    {
        return 0.2405 * pow(1 - Tk / 2503.7, 1.126);
    }
    /* Nusselt du liquide -> Skupinski */
    inline double NuL(double T, double P, double w, double Dh)
    {
      double Pe = RhoL(T, P) * dabs(w) * Dh * DTHL(T, P) / LambdaL(T);
      return 4.82 + 0.0185 * pow(Pe, 0.827);
    }
    /* temperature de saturation */
    inline double Tsat( double P )
    {
      double A = 7.8130, B = 11209, C = 5.249e5;
      return 2 * C / ( -B + sqrt(B*B + 4 * A * C - 4 * C * log(P / 1e6))) - 273.15;
    }
    /* sa derivee */
    inline double DTsat( double P )
    {
      double A = 7.8130, B = 11209, C = 5.249e5, Tsk = Tsat(P) + 273.15;
      return Tsk * Tsk / ( P * sqrt(B*B + 4 * A * C - 4 * C * log(P / 1e6)));
    }
    /* pression de vapeur saturante */
    inline double Psat( double T )
    {
      double A = 7.8130, B = 11209, C = 5.249e5;
      return 1e6 * exp(A - B / Tk - C / (Tk * Tk));
    }
    /* sa derivee */
    inline double DPsat( double T )
    {
      double B = 11209, C = 5.249e5;
      return (B / (Tk * Tk) + 2 * C / (Tk * Tk * Tk)) * Psat(T);
    }
    /* enthalpie massique de saturation */
    inline double Hsat( double P )
    {
      return HL(Tsat(P), P);
    }
    /* sa derivee */
    inline double DHsat( double P )
    {
      return DTsat(P) * DTHL(Tsat(P), P) + DPHL(Tsat(P), P);
    }
    /* chaleur latente, a prendre a Tsat */
    inline double Lvap( double P )
    {
      double Tc = 2503.7, Tsk = Tsat(P) + 273.15;
      return 3.9337e5 * ( 1 - Tsk / Tc) + 4.3986e6 * pow( 1 - Tsk / Tc, .29302);
    }
    
    /* sa derivee */
    inline double DLvap( double P )
    {
      double Tc = 2503.7, Tsk = Tsat(P) + 273.15;
      return DTsat(P) * (-3.9337e5 / Tc - 4.3986e6  * .29302 * pow( 1 - Tsk / Tc, .29302 - 1) / Tc);
    }
    
    /* densite de la vapeur : (gaz parfait) * f1(Tsat) * f2(DTsat)*/
    #define  f1(Ts) (2.49121 + Ts * (-5.53796e-3 + Ts * (7.5465e-6 + Ts * (-4.20217e-9 + Ts * 8.59212e-13))))
    #define Df1(Ts) (-5.53796e-3 + Ts * (2 * 7.5465e-6 + Ts * (-3 * 4.20217e-9 + Ts * 4 * 8.59212e-13)))
    #define  f2(dT) (1 + dT * (-5e-4 + dT * (6.25e-7 - dT * 4.1359e-25)))
    #define Df2(dT) (-5e-4 + dT * (2 * 6.25e-7 - dT * 3 * 4.1359e-25))
    inline double RhoV( double T, double P )
    {
      double Tsk = Tsat(P) + 273.15, dT = Tk - Tsk;
      return P * 2.7650313e-3 * f1(Tsk) * f2(dT) / Tk;
    }
    /* ses derivees */
    inline double DTRhoV( double T, double P )
    {
      double Tsk = Tsat(P) + 273.15, dT = Tk - Tsk;
      return P * 2.7650313e-3 * f1(Tsk) * (Df2(dT) - f2(dT) / Tk) / Tk;
    }
    inline double DPRhoV( double T, double P )
    {
      double Tsk = Tsat(P) + 273.15, dT = Tk - Tsk;
      return 2.7650313e-3 * (f1(Tsk) * f2(dT) + P * DTsat(P) * (Df1(Tsk) * f2(dT) - f1(Tsk) * Df2(dT))) / Tk;
    }
    #undef f1
    #undef Df1
    #undef f2
    #undef Df2
    /* inverse de la densite de la vapeur */
    inline double IRhoV( double T, double P )
    {
      return 1 / RhoV(T, P);
    }
    /* ses derivees */
    inline double DTIRhoV( double T, double P )
    {
      return -DTRhoV(T, P) / (RhoV(T, P) * RhoV(T, P));
    }
    inline double DPIRhoV( double T, double P )
    {
      return -DPRhoV(T, P) / (RhoV(T, P) * RhoV(T, P));
    }
    /* enthalpie de la vapeur */
    #define  CpVs(Ts) (-1223.89 + Ts * (14.1191 + Ts * (-1.62025e-2 + Ts * 5.76923e-6)))
    #define DCpVs(Ts) (14.1191 + Ts * (- 2 * 1.62025e-2 + Ts * 3 * 5.76923e-6))
    inline double HV( double T, double P)
    {
      double Ts = Tsat(P), dT = T - Ts, Cp0 = 910, k = -4.6e-3;
      return HL(Ts, P) + Lvap(P) + Cp0 * dT + (Cp0 - CpVs(Ts)) * (1 - exp(k * dT)) / k;
    }
    /* ses derivees */
    inline double DTHV( double T, double P)
    {
        double Ts = Tsat(P), dT = T - Ts, Cp0 = 910, k = -4.6e-3;
        return Cp0 - (Cp0 - CpVs(Ts)) * exp(k * dT);
    }
    inline double DPHV( double T, double P)
    {
        double Ts = Tsat(P), dT = T - Ts, Cp0 = 910, k = -4.6e-3;
        return DPHL(Ts, P) + DLvap(P) + DTsat(P) * (DTHL(Ts, P) - Cp0 - DCpVs(Ts) * (1 - exp(k * dT)) / k + (Cp0 - CpVs(Ts)) * exp(k * dT));
    }
    #undef CpVs
    #undef DCpVs
    /* energie volumique de la vapeur */
    inline double EV( double T, double P)
    {
      return RhoV(T, P) * HV(T, P);
    }
    /* ses derivees */
    inline double DTEV(double T, double P)
    {
      return DTRhoV(T, P) * HV(T, P) + RhoV(T, P) * DTHV(T, P);
    }
    inline double DPEV(double T, double P)
    {
      return DPRhoV(T, P) * HV(T, P) + RhoV(T, P) * DPHV(T, P);
    }
    
    /* conductivite de la vapeur */
    inline double LambdaV( double T )
    {
      return 0.045;
    }
    
    /* viscosite de la vapeur */
    inline double MuV( double T )
    {
      return 2.5e-6 + 1.5e-8 * Tk;
    }
    
    /* Nusselt de la vapeur -> Dittus/ Boetler */
    inline double NuV(double T, double P, double w, double Dh)
    {
      double Re = RhoV(T, P) * dabs(w) * Dh / MuV(T), Pr = MuV(T) * DTHV(T, P) / LambdaV(T);
      return 0.023 * pow(Re, 0.8) * pow(Pr, 0.4);
    }
    
    inline double Rho( double T, double x, double P )
    {
      return 1 / ( (1-x)/RhoL(T, P) + x / RhoV(T, P) );
    }
    
    #undef Tk
    #undef Ksil
};
#endif
