% Step2_IRF_without_commitment.m
%
% runs simulation and plots results.
% 
% ...................................................................
% Copyright 2018-2020 by Nao Sudo and Masaki Tanaka (Bank of Japan)
%


%%
clear all;
close all;

%%
%--------------------------------------------------------------------
% A. Options
%--------------------------------------------------------------------
    SimCode     = 'Alternative1_IRF'        ; % Name of the mod file: simulation

    LenZLB        = 0;    % length of ZLB comittment
    Shock2Plot    = {'eps_QEs'};

    LenIRFs       = 21;   % # of quarters on which IRFs are computed
    nDraws        = 1000; % # of draws to be used to compute credible bands
        
    nMP       = 24;  % id # of Monetary Policy equation
    nRS       = 22;  % id # of short-term interest rate
    nRL       = 23;  % id # of long-term interest rate
    nQE       = 67;  % id # of central bank's B/S size
    nEpsQE    =  4;  % id # of temporary QE shock
    Var2Plot  = [47, 45, 46];
    
%%
%--------------------------------------------------------------------
% B. Compute the policy functions and IRFs under the commitment to ZLB  
%--------------------------------------------------------------------
% B-1.
    load(sprintf('%s\\Alternative1_est_original\\Alt1_original_draws.mat',  pwd));   
    nowLine = 1;
    
    irfs     = zeros(nDraws, size(Var2Plot,2), LenIRFs, size(Shock2Plot,2));
    
    wtbar    = waitbar(0,'Preparing...','Name','Computing IRFs with credible bands...');

    % iteration: from B-2 to B-4 
    for j = 1 : nDraws
        % B-2. Simulate the model with parameters set as posterior mean
        if j == 1
         	Set_parameters_1;
          	eval(sprintf('dynare %s noclearall nograph', SimCode));
        else
            Set_parameters_2;
            [~,~,M_,~,oo_] = resol(0, M_, options_, oo_);
        end
    
        % B-3. Simulate the model with parameters set as posterior mean
        [PHI1, PHI2]    = zlbPolicyFunc(sprintf('%s_dynamic', SimCode),...
                          LenZLB, nMP, nRS);
 
        % B-4. Compute Impulse Response Functions
        % Adjust shock sizes (10% of GDP)
        M_.Sigma_e(nEpsQE, nEpsQE) = (119.65 / PHI2(nQE, nEpsQE))^2;
  
        irfj = zlbIRF(M_, PHI1, PHI2, 'shocks'   , Shock2Plot,...
                                      'LengthIRF', LenIRFs);
            
        for v = 1 : size(Var2Plot, 2)
            irfs(j, v , :, :) = irfj(Var2Plot(1,v),:,:);
        end
            
        nowLine = nowLine + 1;
        waitbar(j/nDraws, wtbar, sprintf('Now: %s/%s',num2str(j),num2str(nDraws)));
    end
    delete(wtbar)
    
%%
%--------------------------------------------------------------------
% C. Plot the IRFs
%--------------------------------------------------------------------
% C-0. Display
figure(1)
% TP
subplot(3,1,1)
Ydata = squeeze(irfs(:,1,:))*400;
PlotDistBands(Ydata);
axis tight
title('Term Premium')

% Stock
subplot(3,1,2)
Ydata = squeeze(irfs(:,2,:))*400;
PlotDistBands(Ydata);
axis tight
title('Stock Effect')

% Flow
subplot(3,1,3)
Ydata = squeeze(irfs(:,3,:))*400;
PlotDistBands(Ydata);
axis tight
title('Flow Effect')

