/*--------------------------------------------------------------------------- 
 * Alternative2_estimation.mod
 * 
 * This is the dynare code that estimates the model in "Quantifying Stock
 * and Flow Effects of QE."
 * 
 * ...........................................................................
 * Copyright 2018-2020 by Nao Sudo and Masaki Tanaka (Bank of Japan)
 * 
 *---------------------------------------------------------------------------*/

  load Calibrated_Parameters;

//---------------------------------------------------------------------------- 
// Declare the variables
//---------------------------------------------------------------------------- 
//Endogenous variables
  var y K K_bar u I q rk N mc w w0_u w0_r pi Xpn Xpd c c_u c_r MUC MUC_u MUC_r 
      r rL rd bL bLCB bLP PL m m_u m_r g a_u a_r z zl lmbdp_t lmbdw_t 
      mu zeta_ex zeta_st zeta_fl zeta RP_ex RP_st RP_fl RP rLEH dy dI dc dw
      EXPr1 EXPr2 EXPr3 EXPr4 er0 er1 er2 er3 er4 rr rrLEH rrL mv mvP mvCB
      QE eQE;

//Observation variables
  var dy_obs dI_obs dc_obs dN_obs dw_obs pi_obs r_obs rL_obs EXPr1_obs EXPr2_obs
      EXPr3_obs EXPr4_obs dbLCB_obs eps_QE_obs eps_QEs_obs;

//Exogenous variables
  varexo eps_a_u eps_a_r eps_QE eps_QEs eps_z eps_zl eps_zeta eps_r eps_bL eps_g
         eps_mu eps_lmbdp eps_lmbdw eps_r1 eps_r2 eps_r3 eps_r4;

//Exogenous variables
  varexo me_dw me_pi;

//Calibrated values
  parameters alpha beta delta D y_ss g_ss bLratio mvCB_ss Rd_ss;
    set_param_value('alpha', alpha)    ; set_param_value('beta', beta)      ;
    set_param_value('delta', delta)    ; set_param_value('D', D)            ;
    set_param_value('y_ss', y_ss)      ; set_param_value('g_ss', g_ss)      ;
    set_param_value('bLratio', bLratio); set_param_value('mvCB_ss', mvCB_ss);
    set_param_value('Rd_ss', Rd_ss)    ;
 
//Estimated values
//Structual parameters
  parameters sigmau sigmar h nu omega tau1_p tau2_p tau3_p tau4_p tau5_p psi eta
             Cratio zetap zetaw lmbdp lmbdw delta0 delta1;

//Policy parameters
  parameters rho_r rho_pi;

//Shock parameters
  parameters rho_bL rho_QE rho_a_u rho_a_r rho_z rho_zeta rho_g rho_lmbdp
             rho_lmbdw rho_mu rho_eQE;

//Steady-state values
  parameters pi_ss_p gamma_p zeta_ss_p;

//---------------------------------------------------------------------------- 
// Declare the model
//---------------------------------------------------------------------------- 
model(linear);
 // Define parameters and steady-state values as local values
   # tau1    = tau1_p/1000;
   # tau2    = tau2_p/1000;
   # tau3    = tau3_p/1000;
   # tau4    = tau4_p/1000;
   # tau5    = tau5_p/1000;
   # gamma   = gamma_p/400;
   # pi_ss   = pi_ss_p/400;
   # zeta_ss = zeta_ss_p/400;
   # Pi_ss   = exp(pi_ss) ;
   # R_ss    = exp(gamma) * Pi_ss / beta;
   # RL_ss   = (1 + zeta_ss) * R_ss;
   # kappa   = RL_ss * (1 - 1 / D);
   # PL_ss   = 1 / (RL_ss - kappa);
   # mc_ss   = 1 / (1 + lmbdp);
   # rk_ss   = beta^(-1) * exp(gamma) - (1 - delta);
   # w_ss    = alpha^(alpha/(1 - alpha)) * (1 - alpha) *
               rk_ss^(-alpha/(1 - alpha)) * mc_ss^(1/(1 - alpha));
   # K_ss    = alpha / (1 + lmbdp) / rk_ss;
   # K_bar_ss= exp(gamma) * K_ss;
   # I_ss    = (1 - (1 - delta) * exp(-gamma)) * K_bar_ss;
   # N_ss    = (alpha/(1-alpha) * w_ss/rk_ss)^(- alpha);
   # c_ss    = y_ss - g_ss - I_ss;
   # cr_ss   = c_ss / (Cratio * omega + (1 - omega));
   # cu_ss   = cr_ss * Cratio;
   # MUCu_ss = (1 - beta * h)*(1 - h)^(-sigmau)*cu_ss^(-sigmau);
   # MUCr_ss = (1 - beta * h)*(1 - h)^(-sigmar)*cr_ss^(-sigmar);
   # MUC_ss  = omega * MUCu_ss + (1 - omega) * MUCr_ss;
   # r_ss    = log(R_ss);
   # rL_ss   = log(RL_ss);
   # m_ss    = mvCB_ss;
   # phiMinv = m_ss/(omega * (MUCu_ss * (R_ss - Rd_ss)/R_ss)^(-1/delta1) +
               (1 - omega)*(MUCr_ss * (RL_ss - Rd_ss)/RL_ss)^(-1/delta1));
   # mu_ss   = phiMinv * (MUCu_ss * (R_ss  - Rd_ss)/R_ss )^(-1/delta1);
   # mr_ss   = phiMinv * (MUCr_ss * (RL_ss - Rd_ss)/RL_ss)^(-1/delta1);
   # bLCB_ss = mvCB_ss /PL_ss;
   # bLP_ss  = bLCB_ss / bLratio;
   # bL_ss   = bLCB_ss + bLP_ss;
   # w_ratio = (MUCu_ss / MUCr_ss)^(-1/ (1 + (1 + lmbdw)/lmbdw*nu));
   # chiw    = omega / (omega + (1 - omega) * w_ratio^(1 / lmbdw));

 // Define some local values for the real money demand functions
   # mu1  = delta0 / (1 + delta0 * (1 + beta));
   # mu2  = beta * mu1;
   # mu3  = - (1 / delta1) * ( 1 / (1 + delta0 * (1 + beta)) );
   # mu4u = - (Rd_ss / (delta1 * (R_ss  - Rd_ss)) ) * (1 / (1 + delta0 * (1 + beta)) );
   # mu4r = - (Rd_ss / (delta1 * (RL_ss - Rd_ss)) ) * (1 / (1 + delta0 * (1 + beta)) );

//********************************************//
//              State Equations               //
//********************************************//

// 1. Production function
   y = alpha * K + (1 -alpha) * N;

 // 2. Real marginal cost
   mc = alpha * rk + (1 - alpha) * w;

 // 3. FOC for cost minimization 
   K = w - rk + N; 

 // 4. Effective capital
   K = - z + u + K_bar(-1);

 // 5. Law of motion of capital accumulation
   K_bar = (1 - delta)/exp(gamma) * (K_bar(-1) - z) 
           + (1 - (1 - delta)/exp(gamma)) * (I + mu);

 // 6. FOC for optimal capital utilization
   rk = psi * u ;

 // 7. Investment decision
   0 = q + mu - exp(2*gamma) * eta * (I - I(-1) + z) 
      + beta * exp(2*gamma) * eta * (I(+1) - I + z(+1));

 // 8. Real Tobin's Q
   q = beta / exp(gamma) * (rk_ss * (rk(+1) + u(+1)) + (1 - delta) * q(+1))
      + MUC(+1) - MUC - z(+1);

 // 9-11. Phillips curve
   Xpn = (1 - beta * zetap) * (MUC + y + mc + lmbdp_t) 
            + beta * zetap  * ((1 + lmbdp)/lmbdp * pi(+1) + Xpn(+1));   
   Xpd = (1 - beta * zetap) * (MUC + y) 
            + beta * zetap  * ( 1/lmbdp * pi(+1) + Xpd(+1)); 
   pi = (1 - zetap) / zetap * (Xpn - Xpd);

 //12. Marginal utility of consumption: Unrestricted
   MUC_u = 1/(1 - beta * h) * ( a_u - beta * h * a_u(+1) - sigmau/(1 - h) *
           ( (1 + beta * h^2) * c_u - beta * h * c_u(+1) - h * c_u(-1) ) );

 //13. Marginal utility of consumption: Restricted
   MUC_r = 1/(1 - beta * h) * ( a_r - beta * h * a_r(+1) - sigmar/(1 - h) *
           ( (1 + beta * h^2) * c_r - beta * h * c_r(+1) - h * c_r(-1) ) );

 //14. Euler equation: Unrestrected, Long-term bonds
   MUC_u = rL(+1) + PL(+1) - PL - pi(+1) + MUC_u(+1) - z(+1) - zeta;

 //15. Euler equation: Unrestrected, Short-term bonds
   MUC_u = r - pi(+1) + MUC_u(+1) - z(+1);

 //16. Euler equation: Restrected, Long-term bonds
   MUC_r = rL(+1) + PL(+1) - PL - pi(+1) + MUC_r(+1) - z(+1);

 //17. Aggregation of marginal utility of consumption
   MUC_ss * MUC = omega * MUCu_ss * MUC_u + (1 - omega) * MUCr_ss * MUC_r;

 //18. Real wage setting: Unrestricted
   w0_u = beta*zetaw * (w0_u(+1) + pi(+1) + z(+1)) + (1 - beta*zetaw) * w
        + (1 - beta*zetaw)*lmbdw/(lmbdw + nu*(1 + lmbdw)) * (lmbdw_t + a_u -MUC_u + nu * N - w);

 //19. Real wage setting: Restricted
   w0_r = beta*zetaw * (w0_r(+1) + pi(+1) + z(+1)) + (1 - beta*zetaw) * w
        + (1 - beta*zetaw)*lmbdw/(lmbdw + nu*(1 + lmbdw)) * (lmbdw_t + a_r -MUC_r + nu * N - w);

 //20. Law of motion of real wages
   w =  (1 - zetaw)* (chiw * w0_u + (1 - chiw) * w0_r )
        + zetaw * (w(-1) - pi - z);

 //21. Real money demand: Unrestricted
   m_u = mu1 * m_u(-1) + mu2 * m_u(+1) + mu3 * (MUC_u - a_u) + mu4u *(r - rd);

 //22. Real money demand: Restricted
   m_r = mu1 * m_r(-1) + mu2 * m_r(+1) + mu3 * (MUC_r - a_r) + mu4r * (PL(+1) - PL + rL(+1) - rd);

 //23. Aggregate real money demand
   m_ss * m = omega * mu_ss * m_u + (1 - omega) * mr_ss * m_r;

 //24. Monetary policy
   r = rho_r * r(-1) + (1 - rho_r)* (rho_pi * pi) - er0;

 //25. Transaction cost
   zeta = zeta_st + zeta_fl +zeta_ex;

 //26. "Stock effect" component
   zeta_st = tau1 * mvP - tau2 * m_u - tau3 * m_r;

 //27. "Flow effect" component
   zeta_fl = tau4 *(- eps_QEs) + tau5 *(- eps_QE);

 //28. Aggregate real consumption
   c_ss * c = omega * cu_ss * c_u + (1 - omega) * cr_ss * c_r;

 //29. Resource constraint for the final good market
   y = c_ss * c + I_ss * I + g_ss * g + rk_ss * K_bar_ss / exp(gamma) * u ;

 //30. Long-term bond policy
   mv = rho_bL * mv(-1) + eps_bL;

 //31. Resource constraint for the long-term bond market
   bL_ss * bL = bLCB_ss * bLCB + bLP_ss * bLP;

 //32. QE Operation
   mvCB = rho_QE * mvCB(-1) + eQE + eps_QEs;

 //33. The central bank's B/S
   m = mvCB;

 //Exogenous processes
   //34. Preference shock process: Unrestricted
     a_u = rho_a_u * a_u(-1) + eps_a_u;

   //35. Preference shock process: Restricted
     a_r = rho_a_r * a_r(-1) + eps_a_r;

   //36. Price markup shock process
     lmbdp_t = rho_lmbdp * lmbdp_t(-1) + eps_lmbdp;

   //37. Price markup shock process
     lmbdw_t = rho_lmbdw * lmbdw_t(-1) + eps_lmbdw;

   //38. Productivity shock process: persistent + temporary
     z = zl + eps_z;
   
   //39. Productivity shock process: persistent
     zl = rho_z * zl(-1) + eps_zl;

   //40. Exogenous transaction cost shock process
     zeta_ex = rho_zeta * zeta_ex(-1) - eps_zeta;

   //41. Gov't expenditure process
     g = rho_g * g(-1) + eps_g;

   //42. Investment shock process
     mu   = rho_mu * mu(-1) + eps_mu;

   //43-47. News shocks on monetary policy
     er0 = er1(-1) + eps_r;
     er1 = er2(-1) + eps_r1;
     er2 = er3(-1) + eps_r2;
     er3 = er4(-1) + eps_r3;
     er4 = eps_r4;

   //48. QE shock Process
     eQE = rho_eQE * eQE(-1) + eps_QE;

//Define some auxiliary variables
   //49. Output growth
     dy = y - y(-1) + z;

   //50. Real wage growth
     dw = w - w(-1) + z;

   //51. Real investment growth
     dI = I - I(-1) + z;

   //52. Real consumption growth
     dc = c - c(-1) + z;

   //53. Price of long-term bonds
     PL = - D *rL;

   //54. Market value of outstanding long-term bonds
     mv = PL + bL;

   //55. Market value of long-term bonds held by the private sector
     mvP = PL + bLP;

   //56. Market value of long-term bonds held by the central bank
     mvCB = PL + bLCB;

   //57. Risk premium
     RP(+1) = D/(D - 1) * RP - 1/(D - 1) * zeta;

   //58. Stock effect component of risk premium
     RP_st(+1) = D/(D - 1) * RP_st - 1/(D - 1) * zeta_st;

   //59. Flow effect component of risk premium
     RP_fl(+1) = D/(D - 1) * RP_fl - 1/(D - 1) * zeta_fl;

   //60. Exogenous component of risk premium
     RP_ex(+1) = D/(D - 1) * RP_ex - 1/(D - 1) * zeta_ex;

   //61. Risk free long-term rate (EH: Expectation Hypothesis)
     rLEH(+1) = D/(D - 1) * rLEH - 1/(D - 1) * r;

   //62. Short-term real interest rate
     rr    = r - pi(+1);

   //63. Long-term real yield, based on the expectations hypothesis
     rrLEH = D/(D - 1) * rrLEH - 1/(D - 1) * rr;

   //64. Long-term real yield
     rrL   = rrLEH + RP;

   //65-68. Expected short-term interst rates
     EXPr1 = r(+1);
     EXPr2 = EXPr1(+1);
     EXPr3 = EXPr2(+1);
     EXPr4 = EXPr3(+1);

//********************************************//
//            Observation Equations           //
//********************************************//
 //69. Observed output growth
   dy_obs = 100 * gamma + dy;

 //70. Observed consumption growth
   dc_obs = 100 * gamma + dc;

 //71. Observed investment growth
   dI_obs = 100 * gamma + dI;

 //72. Observed inflation
   pi_obs = 100 * pi_ss + pi + me_pi;

 //73. Observed labor supply
   dN_obs  = N - N(-1);
 
 //74. Observed real wage growth
   dw_obs = 100 * gamma + dw + me_dw;

 //75. Observed nominal short-term interest rate
   r_obs  = 100 * r_ss + r;

 //76. Observed nominal long-term interest rate
   rL_obs = 100 * rL_ss + rL;

 //77-80. Observed expected nominal short-term interest rate
   EXPr1_obs  = 100 * r_ss + EXPr1;
   EXPr2_obs  = 100 * r_ss + EXPr2;
   EXPr3_obs  = 100 * r_ss + EXPr3;
   EXPr4_obs  = 100 * r_ss + EXPr4;

 //81. Observed real long-term bonds held by the central bank
   dbLCB_obs = 100 * gamma + mvCB - mvCB(-1) + z;

 //82. Observed QE shock
   eps_QE_obs  = eps_QE;
   eps_QEs_obs = eps_QEs;

// Define an additional auxiliary variable
// (Sum of stock and flow effects)
   QE = RP_st + RP_fl;
end;

//---------------------------------------------------------------------------- 
// Define the exogenous shocks
//---------------------------------------------------------------------------- 
shocks;
    var eps_a_u    ; stderr 0.5;
    var eps_a_r    ; stderr 0.5;
    var eps_z      ; stderr 0.5;
    var eps_zl     ; stderr 0.5;
    var eps_zeta   ; stderr 0.5;
    var eps_r      ; stderr 0.5;
    var eps_bL     ; stderr 0.5;
    var eps_QE     ; stderr 0.5;
    var eps_QEs    ; stderr 0.5;
    var eps_g      ; stderr 0.5;
    var eps_lmbdp  ; stderr 0.5;
    var eps_lmbdw  ; stderr 0.5;
    var eps_mu     ; stderr 0.5;
    var eps_r1     ; stderr 0.5;
    var eps_r2     ; stderr 0.5;
    var eps_r3     ; stderr 0.5;
    var eps_r4     ; stderr 0.5;
    var me_dw      ; stderr 0.5;
    var me_pi      ; stderr 0.5;
end;

-/---------------------------------------------------------------------------- 
// Estimate Parameters
//---------------------------------------------------------------------------- 
//********************************************//
//              Declare Priors                //
//********************************************//
estimated_params;
  // structual parameters
    sigmau  ,  gamma_pdf, 2.0000, 0.5000; // Chen et al.(2012)
    sigmar  ,  gamma_pdf, 2.0000, 0.5000; // Chen et al.(2012)
    h       ,   beta_pdf, 0.7000, 0.0500; // Fueki et al.(2010)
    nu      ,  gamma_pdf, 2.0000, 0.5000; // Fueki et al.(2010)
    omega   ,   beta_pdf, 0.5000, 0.2000; // 
    psi     ,  gamma_pdf, 1.0000, 0.5000; // Fueki et al.(2010)
    eta     ,  gamma_pdf, 4.0000, 1.0000; // Smets and Wouters (2007)
    zetap   ,   beta_pdf, 0.5000, 0.1000; // Chen et al.(2012)
    zetaw   ,   beta_pdf, 0.5000, 0.1000; // Chen et al.(2012)
    lmbdp   ,  gamma_pdf, 0.1500, 0.0200; // Justiniano et al.(2010)
    lmbdw   ,  gamma_pdf, 0.1500, 0.0200; // Justiniano et al.(2010)
    Cratio  ,  gamma_pdf, 1.0000, 0.0500; // Chen et al.(2012)
    tau1_p  , normal_pdf, 1.5000, 0.5000; // Chen et al.(2012)
    tau2_p  , normal_pdf, 1.5000, 0.5000; // Chen et al.(2012)
    tau3_p  , normal_pdf, 1.5000, 0.5000; // Chen et al.(2012)
    tau4_p  , normal_pdf, 2.5000, 0.5000; // Chen et al.(2012)
    tau5_p  , normal_pdf, 2.5000, 0.5000; // Chen et al.(2012)
    delta1  ,  gamma_pdf, 1.8200, 1.0000; // ALSN(2004)
    delta0  ,  gamma_pdf, 4.3600, 1.0000; // ALSN(2004)

  // policy parameters
    rho_r     ,   beta_pdf, 0.7500, 0.0500; //
    rho_pi    ,  gamma_pdf, 1.5000, 0.2000; // 

  // steady state values
    gamma_p   , gamma_pdf, 1.3000, 0.2000; // Data Mean
    pi_ss_p   , gamma_pdf, 0.3000, 0.2000; // Data Mean
    zeta_ss_p , gamma_pdf, 1.0000, 0.2000; // Data Mean

  // standard deviation of shocks
    stderr eps_bL      , inv_gamma_pdf, 0.50, inf; // 
    stderr eps_QE      , inv_gamma_pdf, 5.00, inf; // 
    stderr eps_QEs     , inv_gamma_pdf, 5.00, inf; // 
    stderr eps_r       , inv_gamma_pdf, 0.10, inf; // Fueki et al.(2010)
    stderr eps_a_u     , inv_gamma_pdf, 5.00, inf; // Fueki et al.(2010)
    stderr eps_a_r     , inv_gamma_pdf, 5.00, inf; // Fueki et al.(2010)
    stderr eps_z       , inv_gamma_pdf, 5.00, inf; // Fueki et al.(2010)
    stderr eps_zl      , inv_gamma_pdf, 0.10, inf; // Fueki et al.(2010)
    stderr eps_g       , inv_gamma_pdf, 1.00, inf; // Fueki et al.(2010)
    stderr eps_zeta    , inv_gamma_pdf, 0.50, inf; //
    stderr eps_lmbdp   , inv_gamma_pdf, 1.00, inf; // Fueki et al.(2010)
    stderr eps_lmbdw   , inv_gamma_pdf, 1.00, inf; // Fueki et al.(2010)
    stderr eps_mu      , inv_gamma_pdf, 3.00, inf; // Fueki et al.(2010)
    stderr eps_r1      , inv_gamma_pdf, 0.10, inf; // same to that of eps_r
    stderr eps_r2      , inv_gamma_pdf, 0.10, inf; // same to that of eps_r
    stderr eps_r3      , inv_gamma_pdf, 0.10, inf; // same to that of eps_r
    stderr eps_r4      , inv_gamma_pdf, 0.10, inf; // same to that of eps_r
    stderr me_dw       , inv_gamma_pdf, 0.30, inf; //
    stderr me_pi       , inv_gamma_pdf, 0.30, inf; //

  // persistence of shocks
    rho_bL    , beta_pdf, 0.7000, 0.1500;
    rho_QE    , beta_pdf, 0.8000, 0.1500;
    rho_a_u   , beta_pdf, 0.8000, 0.1500; // Fueki et al.(2010)
    rho_a_r   , beta_pdf, 0.8000, 0.1500; // Fueki et al.(2010)
    rho_z     , beta_pdf, 0.9800, 0.0100; // Fueki et al.(2010)
    rho_zeta  , beta_pdf, 0.7000, 0.1500;
    rho_g     , beta_pdf, 0.5000, 0.1500; // Fueki et al.(2010)
    rho_lmbdp , beta_pdf, 0.5000, 0.1500; // Fueki et al.(2010)
    rho_lmbdw , beta_pdf, 0.5000, 0.1500; // Fueki et al.(2010)
    rho_mu    , beta_pdf, 0.7000, 0.1500; // Fueki et al.(2010)
    rho_eQE   , beta_pdf, 0.7000, 0.1500; 
end;

estimated_params_init;
stderr eps_bL,0.197904757;
stderr eps_QE,3.203406321;
stderr eps_QEs,0.590879388;
stderr eps_r,0.022429579;
stderr eps_a_u,4.039351158;
stderr eps_a_r,9.658459644;
stderr eps_z,1.761657307;
stderr eps_zl,0.094868742;
stderr eps_g,3.034648309;
stderr eps_zeta,0.191253893;
stderr eps_lmbdp,1.980692459;
stderr eps_lmbdw,1.250132546;
stderr eps_mu,6.997915539;
stderr eps_r1,0.018106849;
stderr eps_r2,0.012978015;
stderr eps_r3,0.016283754;
stderr eps_r4,0.01238161;
stderr me_dw,1.149174402;
stderr me_pi,0.23756498;
sigmau,2.228793776;
sigmar,1.146142964;
h,0.65346347;
nu,1.007726654;
omega,0.63809498;
psi,1.460671684;
eta,4.814769761;
zetap,0.900940135;
zetaw,0.712095446;
lmbdp,0.188527043;
lmbdw,0.136243652;
Cratio,1.006828756;
tau1_p,1.567997371;
tau2_p,1.846202638;
tau3_p,0.31709969;
delta1,1.321654207;
delta0,4.573742817;
rho_r,0.820835659;
rho_pi,2.146723105;
gamma_p,1.017093263;
pi_ss_p,1.079779941;
zeta_ss_p,1.050763799;
rho_bL,0.743838185;
rho_QE,0.797904141;
rho_a_u,0.994098185;
rho_a_r,0.114055691;
rho_z,0.987521198;
rho_zeta,0.881307988;
rho_g,0.878880147;
rho_lmbdp,0.818000249;
rho_lmbdw,0.546106646;
rho_mu,0.481129371;
rho_eQE,0.990627772;
end;

//********************************************//
//       Declare Observed Variables           //
//********************************************//
varobs dy_obs dc_obs dI_obs dN_obs dw_obs pi_obs r_obs rL_obs EXPr1_obs
       EXPr2_obs EXPr3_obs EXPr4_obs dbLCB_obs eps_QE_obs eps_QEs_obs;

//********************************************//
//             Run Estimation                 //
//********************************************//
estimation(
datafile=data_86Q3to17Q4, 
mode_compute=0, mode_file=Alternative2_estimation_mode,
mh_replic=200000, mh_nblocks=2,
mh_jscale=0.18,
smoother
)y pi r rL RP RP_st RP_fl QE eQE mvCB mv;
