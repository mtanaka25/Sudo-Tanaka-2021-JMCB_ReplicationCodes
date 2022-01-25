% Step4_historical_decomp.m
%
% calculates historical decompositions.
% ...................................................................
% Copyright 2018-2020 by Nao Sudo and Masaki Tanaka (Bank of Japan)
%


%%
clear all;
close all;

%%
%-------------------------------------------------------------------------
% A. Setting
%-------------------------------------------------------------------------
    SimCode     = 'Baseline_HD'         ; % Name of the mod file: simulation

    % # of parameter draws
    nDraws     = 1000;

%%
%-------------------------------------------------------------------------
% B. Calclating Shock Decompositions
%-------------------------------------------------------------------------
    load(sprintf('%s\\Baseline_est_original\\Baseline_original_draws.mat',  pwd));   
    nowLine = 1;
    
    % Allocate memory 
    HD_Stock_3D = NaN(126, 21, nDraws);
    HD_Flow_3D  = NaN(126, 21, nDraws);
   
    wtbar    = waitbar(0,'Preparing...','Name','Calculating shock decompositions...');
    for j = 1 : nDraws

        Set_parameters_3;

        eval(sprintf('dynare %s noclearall nograph', SimCode));
        close all

        HD_Stock_3D(:,:,j) = squeeze(oo_.shock_decomposition(45, :, :))';
        HD_Flow_3D(:,:,j)  = squeeze(oo_.shock_decomposition(46, :, :))';
        
        nowLine = nowLine + 1;
        waitbar(j/nDraws, wtbar, sprintf('Now: %s/%s',num2str(j),num2str(nDraws)));
    end
    delete(wtbar)
 
%-------------------------------------------------------------------------
% C. Write out the result
%-------------------------------------------------------------------------
    fname  = 'Historical_Decomposition_Baseline.xlsx';

    warning off MATLAB:xlswrite:AddSheet

    HD_Stock = prctile(HD_Stock_3D, [92.5, 50, 2.5], 3);    
    HD_Flow  = prctile(HD_Flow_3D , [92.5, 50, 2.5], 3);   
    load('data_86Q3to17Q4.mat', 'rL_obs');
    
    Stock_Effect = sum(HD_Stock(:,[3, 4], 2), 2) * 4;
    Flow_Effect  = sum(HD_Flow(:,[3, 4], 2), 2) * 4;
    Others       = 4 * rL_obs - (Stock_Effect + Flow_Effect);   
 
    QE_Total_med  = Stock_Effect + Flow_Effect;
    QE_Total_high = (sum(HD_Stock(:,[3, 4], 1), 2) + sum(HD_Flow(:,[3, 4], 1), 2)) * 4;
    QE_Total_low  = (sum(HD_Stock(:,[3, 4], 3), 2) + sum(HD_Flow(:,[3, 4], 3), 2)) * 4;

    xlswrite(fname, [QE_Total_high, QE_Total_med, QE_Total_low]    , 'Figure5_1', 'A1');
    xlswrite(fname, [4 * rL_obs, Stock_Effect, Flow_Effect, Others], 'Figure5_2', 'A1');