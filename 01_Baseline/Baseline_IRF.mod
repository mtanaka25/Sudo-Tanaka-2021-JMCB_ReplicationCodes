/*--------------------------------------------------------------------------- 
 * Baseline_IRF.mod
 * 
 * This is the dynare code that calculates IRF of the model in "Quantifying 
 * Stock and Flow Effects of QE."
 *
 * ...........................................................................
 * Copyright 2018-2020 by Nao Sudo and Masaki Tanaka (Bank of Japan)
 *
 *---------------------------------------------------------------------------*/

//---------------------------------------------------------------------------- 
// Declare the variables
//---------------------------------------------------------------------------- 
//Endogenous variables
  var y K K_bar u I q rk N mc w w0_u w0_r pi Xpn Xpd c c_u c_r MUC MUC_u MUC_r 
      r rL rd bL bLCB bLP PL m m_u m_r g a_u a_r z zl lmbdp_t lmbdw_t 
      mu zeta_ex zeta_st zeta_fl zeta RP_ex RP_st RP_fl RP rLEH dy dI dc dw
      EXPr1 EXPr2 EXPr3 EXPr4 er0 er1 er2 er3 er4 rr rrLEH rrL mv mvP mvCB
      QE eQE;

//Exogenous variables
  varexo eps_a_u eps_a_r eps_QE eps_QEs eps_z eps_zl eps_zeta eps_r eps_bL eps_g
         eps_mu eps_lmbdp eps_lmbdw eps_r1 eps_r2 eps_r3 eps_r4;

//Calibrated values
  parameters alpha beta delta D;
    set_param_value('alpha', alpha)    ; set_param_value('beta', beta)   ;
    set_param_value('delta', delta)    ; set_param_value('D', D)         ;
 

//Structual parameters
  parameters sigmau sigmar h nu tau1 tau2 tau3 tau4 tau5 tau6 eta psi zetap
             zetaw lmbdp lmbdw delta0 delta1;
    set_param_value('sigmau', sigmau) ; set_param_value('sigmar', sigmar);
    set_param_value('h', h)           ; set_param_value('nu', nu)        ;
    set_param_value('tau1', tau1)     ; set_param_value('tau2', tau2)    ; 
    set_param_value('tau3',tau3)      ; set_param_value('tau4', tau4)    ; 
    set_param_value('tau5',tau5)      ; set_param_value('eta', eta)      ;
    set_param_value('psi', psi)       ; set_param_value('zetap', zetap)  ;
    set_param_value('lmbdp', lmbdp)   ; set_param_value('lmbdw', lmbdw)  ;
    set_param_value('zetaw', zetaw)   ; set_param_value('delta0', delta0);
    set_param_value('delta1', delta1) ; set_param_value('tau6', tau6)    ;

//Policy parameters
  parameters rho_r rho_pi;
    set_param_value('rho_r', rho_r)   ; set_param_value('rho_pi', rho_pi);

//Shock parameters
  parameters rho_bL rho_QE rho_a_u rho_a_r rho_z rho_zeta rho_g rho_lmbdp
             rho_lmbdw rho_mu rho_eQE
             sig_bL sig_QE sig_QEs sig_a_u sig_a_r sig_z sig_zl sig_r sig_zeta sig_g
             sig_lmbdp sig_lmbdw sig_mu sig_r1 sig_r2 sig_r3 sig_r4;
    set_param_value('rho_bL', rho_bL)      ; set_param_value('rho_QE', rho_QE)      ;
    set_param_value('rho_a_u', rho_a_u)    ; set_param_value('rho_a_r', rho_a_r)    ;
    set_param_value('rho_z', rho_z)        ; set_param_value('rho_zeta', rho_zeta)  ;
    set_param_value('rho_g', rho_g)        ; set_param_value('rho_lmbdp', rho_lmbdp);
    set_param_value('rho_lmbdw', rho_lmbdw); set_param_value('rho_mu', rho_mu)      ;
    set_param_value('rho_eQE', rho_eQE)    ; set_param_value('sig_bL', sig_bL)      ;
    set_param_value('sig_QE', sig_QE)      ; set_param_value('sig_QEs', sig_QEs)    ;
    set_param_value('sig_a_u', sig_a_u)    ; set_param_value('sig_a_r', sig_a_r)    ;
    set_param_value('sig_z', sig_z)        ; set_param_value('sig_zl', sig_zl)      ;
    set_param_value('sig_r', sig_r)        ; set_param_value('sig_zeta', sig_zeta)  ;
    set_param_value('sig_g', sig_g)        ; set_param_value('sig_lmbdp', sig_lmbdp); 
    set_param_value('sig_lmbdw', sig_lmbdw); set_param_value('sig_mu', sig_mu)      ;
    set_param_value('sig_r1', sig_r1)      ; set_param_value('sig_r2', sig_r2)      ;
    set_param_value('sig_r3', sig_r3)      ; set_param_value('sig_r4', sig_r4)      ;


//Steady-state values
  parameters  pi_ss Pi_ss gamma R_ss RL_ss Rd_ss PL_ss m_ss mu_ss mr_ss I_ss 
              c_ss cu_ss cr_ss rk_ss K_bar_ss MUC_ss MUCu_ss MUCr_ss bLP_ss
              bLCB_ss omega g_ss bL_ss;
    set_param_value('pi_ss', pi_ss)    ; set_param_value('Pi_ss', Pi_ss)      ;
    set_param_value('gamma', gamma)    ; set_param_value('Rd_ss', Rd_ss)      ;
    set_param_value('R_ss', R_ss)      ; set_param_value('RL_ss', RL_ss)      ;
    set_param_value('PL_ss', PL_ss)    ; set_param_value('m_ss', m_ss)        ;
    set_param_value('mu_ss', mu_ss)    ; set_param_value('mr_ss', mr_ss)      ;
    set_param_value('I_ss', I_ss)      ; set_param_value('c_ss', c_ss)        ;
    set_param_value('cu_ss', cu_ss)    ; set_param_value('cr_ss', cr_ss)      ;
    set_param_value('rk_ss', rk_ss)    ; set_param_value('K_bar_ss', K_bar_ss);
    set_param_value('MUC_ss', MUC_ss)  ; set_param_value('MUCu_ss', MUCu_ss)  ;
    set_param_value('MUCr_ss', MUCr_ss); set_param_value('bLP_ss', bLP_ss)    ;
    set_param_value('bLCB_ss', bLCB_ss); set_param_value('omega', omega)      ;
    set_param_value('g_ss', g_ss)      ; set_param_value('bL_ss', bL_ss)      ;

//---------------------------------------------------------------------------- 
// Declare the model
//---------------------------------------------------------------------------- 
model(linear);
   # w_ratio = (MUCu_ss / MUCr_ss)^(-1/ (1 + (1 + lmbdw)/lmbdw*nu));
   # chiw    = omega / (omega + (1 - omega) * w_ratio^(1/lmbdw));

// Define some local values for the real money demand functions
   # mu1  = delta0 / (1 + delta0 * (1 + beta));
   # mu2  = beta * mu1;
   # mu3  = - (1 / delta1) * ( 1 / (1 + delta0 * (1 + beta)) );
   # mu4u = - (Rd_ss / (delta1 * (R_ss - Rd_ss)) ) * (1 / ( 1 + delta0 * (1 + beta)) );
   # mu4r = - (Rd_ss / (delta1 * (RL_ss - Rd_ss)) ) * (1 / ( 1 + delta0 * (1 + beta)) );

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
   MUC_u = 1/(1 - beta*h) * ( a_u - beta * h * a_u(+1) - sigmau/(1 - h) *
           ( (1 + beta * h^2) * c_u - beta * h * c_u(+1) - h * c_u(-1) ) );

 //13. Marginal utility of consumption: Restricted
   MUC_r = 1/(1 - beta*h) * ( a_r - beta * h * a_r(+1) - sigmar/(1 - h) *
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
   zeta_fl = tau4 * (mvP - mvP(-1))
             - tau5 * (m_u - m_u(-1)) - tau6 * (m_r - m_r(-1));

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

   //38. Total Productivity shock
     z = zl + eps_z;

   //39. Persistent Productivity Shock Process
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

   //48. Persistent QE Shock process
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

   //69. Sum of stock and flow effect
     QE = RP_st + RP_fl;
end;

//---------------------------------------------------------------------------- 
// Compute the steady state and check if Blanchard-Kahn Condition is verified
//---------------------------------------------------------------------------- 
steady;
check;

//---------------------------------------------------------------------------- 
// Define the exogenous shocks
//---------------------------------------------------------------------------- 
shocks;
    var eps_a_u   ; stderr sig_a_u   ;
    var eps_a_r   ; stderr sig_a_r   ;
    var eps_z     ; stderr sig_z     ;
    var eps_zl    ; stderr sig_zl    ;
    var eps_zeta  ; stderr sig_zeta  ;
    var eps_r     ; stderr sig_r     ;
    var eps_bL    ; stderr sig_bL    ;
    var eps_QE    ; stderr sig_QE    ;
    var eps_QEs   ; stderr sig_QEs   ;
    var eps_g     ; stderr sig_g     ;
    var eps_lmbdp ; stderr sig_lmbdp ;
    var eps_lmbdw ; stderr sig_lmbdw ;
    var eps_mu    ; stderr sig_mu    ;
    var eps_r1    ; stderr sig_r1    ;
    var eps_r2    ; stderr sig_r2    ;
    var eps_r3    ; stderr sig_r3    ;
    var eps_r4    ; stderr sig_r4    ;
end;

//---------------------------------------------------------------------------- 
// Run stochastic simulations
//---------------------------------------------------------------------------- 
stoch_simul(order=1, irf=30, irf_shocks=(eps_r))dy;