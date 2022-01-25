% MakeParamsFilePost.m
%
% Assign posterior mean to the estimated parameters.
%
% ...........................................................................
% Create: September 26, 2017
%
%%
%--------------------------------------------------------------------
% a. Define parameters
%--------------------------------------------------------------------
% a-1. Calibrated Parameters
    load('Calibrated_Parameters.mat');
    
% a-2. Estimated structual parameters 
    sigmau  = x2(nowLine, 20);
    sigmar  = x2(nowLine, 21);
    h       = x2(nowLine, 22);
    nu      = x2(nowLine, 23);
    omega   = x2(nowLine, 24);
    psi     = x2(nowLine, 25);
    eta     = x2(nowLine, 26);
    zetap   = x2(nowLine, 27);
    zetaw   = x2(nowLine, 28);
    lmbdp   = x2(nowLine, 29);
    lmbdw   = x2(nowLine, 30);
    Cratio  = x2(nowLine, 31);
    tau1    = x2(nowLine, 32)/1000;
    tau2    = x2(nowLine, 33)/1000;
    tau3    = x2(nowLine, 34)/1000;
    delta1  = x2(nowLine, 35);
    delta0  = x2(nowLine, 36);
    
% a-3. Estimated policy parameters
    rho_r     = x2(nowLine, 37);
    rho_pi    = x2(nowLine, 38);

% a-4. Estimated steady-state values
    gamma_p   = x2(nowLine, 39);
    pi_ss_p   = x2(nowLine, 40);
    zeta_ss_p = x2(nowLine, 41);
    
% a-5. Estimated shock parameters
    rho_bL    = x2(nowLine, 42);
    %rho_QE    = x2(nowLine, 43);
    rho_QE    = 0.896675567; % Posterior mean in the baseline case
    rho_a_u   = x2(nowLine, 44);
    rho_a_r   = x2(nowLine, 45);
    rho_z     = x2(nowLine, 46);
    rho_zeta  = x2(nowLine, 47);
    rho_g     = x2(nowLine, 48);
    rho_lmbdp = x2(nowLine, 49);
    rho_lmbdw = x2(nowLine, 50);
    rho_mu    = x2(nowLine, 51);
    rho_eQE   = x2(nowLine, 52);

    sig_bL      = x2(nowLine, 1);
    sig_QE      = x2(nowLine, 2);
    sig_QEs     = x2(nowLine, 3);
    sig_r       = x2(nowLine, 4);
    sig_a_u     = x2(nowLine, 5);
    sig_a_r     = x2(nowLine, 6);
    sig_z       = x2(nowLine, 7);
    sig_zl      = x2(nowLine, 8);
    sig_g       = x2(nowLine, 9);
    sig_zeta    = x2(nowLine, 10);
    sig_lmbdp   = x2(nowLine, 11);
    sig_lmbdw   = x2(nowLine, 12);
    sig_mu      = x2(nowLine, 13);
    sig_r1      = x2(nowLine, 14);
    sig_r2      = x2(nowLine, 15);
    sig_r3      = x2(nowLine, 16);
    sig_r4      = x2(nowLine, 17);
    sig_me_dw   = x2(nowLine, 18);
    sig_me_pi   = x2(nowLine, 19);

%%
%--------------------------------------------------------------------
% b. Compute steady state values
%--------------------------------------------------------------------
% Technology growth rate
    gamma = gamma_p/400;

% Net inflation rate (quarterly)
    pi_ss = pi_ss_p/400;

% Steady state transaction cost(quarterly)
    zeta_ss = zeta_ss_p/400;
    
% Gross inflation rate (quarterly)
    Pi_ss = exp(pi_ss);

% Short-term interst rate (<= Eular equation: Unrestrected, Short)
    R_ss  = exp(gamma) * Pi_ss / beta;

% Long-term interst rate (<= Eular equation: Restrected, Long)
    RL_ss  = (1 + zeta_ss) * R_ss;

% coupon  (æ D = rL_ss / (rL - kappa) )
    kappa = RL_ss * (1 - 1 / D);

% Market price of long-term bonds
    PL_ss = 1 / (RL_ss - kappa);

% Real marginal cost (<= Optimal relative prices)
    mc_ss = 1 / (1 + lmbdp);

% Real return on capital (<= FOC for optimal capital stock)
    rk_ss = beta^(-1) * exp(gamma) - (1 - delta);

% Real wages (<= Real marginal cost func.)
    w_ss = alpha^(alpha/(1 - alpha)) * (1 - alpha) *...
           rk_ss^(-alpha/(1 - alpha)) * mc_ss^(1/(1 - alpha));

% Effective capital (<= Capital demand func.)
    K_ss = alpha / (1 + lmbdp) / rk_ss;

% Capital stock
    K_bar_ss = exp(gamma) * K_ss;

% Investment (<= Capital accumulation)
    I_ss = (1 - (1 - delta) * exp(-gamma)) * K_bar_ss;

% Labor (<= Labor demand func.)
    N_ss = (alpha/(1-alpha) * w_ss/rk_ss)^(- alpha);

% Consumption (<= Resource Constraint)
    c_ss  = y_ss - g_ss - I_ss;
    cr_ss = c_ss / (Cratio * omega + (1 - omega));
    cu_ss = cr_ss * Cratio;

% Marginal utility of consumption (<= FOC for HH's utility maximization)
    MUCu_ss = (1 - beta * h)*(1 - h)^(-sigmau)*cu_ss^(-sigmau);
    MUCr_ss = (1 - beta * h)*(1 - h)^(-sigmar)*cr_ss^(-sigmar);
    MUC_ss  = omega * MUCu_ss + (1 - omega) * MUCr_ss;
    
% Real money demand
    m_ss    = mvCB_ss;
    phiMinv = m_ss/(omega * (MUCu_ss * (R_ss - Rd_ss)/R_ss)^(-1/delta1) +...
               (1 - omega)*(MUCr_ss * (RL_ss - Rd_ss)/RL_ss)^(-1/delta1));
    mu_ss   = phiMinv * (MUCu_ss * (R_ss - Rd_ss)/R_ss)^(-1/delta1);
    mr_ss   = phiMinv * (MUCr_ss * (RL_ss - Rd_ss)/RL_ss)^(-1/delta1);

% Long-term bonds held by the private sector
    bLCB_ss = mvCB_ss / PL_ss;
    bLP_ss  = bLCB_ss / bLratio;
    bL_ss  = bLCB_ss + bLP_ss;
