% Step3_historical_decomp.m
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
    SimCode     = 'Alternative2_HD'         ; % Name of the mod file: simulation

    % # of parameter draws
    nDraws     = 1000;

%%
%-------------------------------------------------------------------------
% B. Calclating Shock Decompositions
%-------------------------------------------------------------------------
    load(sprintf('%s\\Alternative2_est_original\\Alt2_original_draws.mat',  pwd));   
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
    fname  = 'Historical_Decomposition_Alt2.xlsx';

    warning off MATLAB:xlswrite:AddSheet

    HD_Stock     = prctile(HD_Stock_3D, 50, 3);
    Stock_Effect = sum(HD_Stock(:,[3, 4]), 2) * 4;
    
    HD_Flow     = prctile(HD_Flow_3D, 50, 3);   
    Flow_Effect = sum(HD_Flow(:,[3, 4]), 2) * 4;
    
    load('data_86Q3to17Q4.mat', 'rL_obs');
    QE_Total = Stock_Effect + Flow_Effect;
    Others   = 4 * rL_obs - (Stock_Effect + Flow_Effect);
    
    xlswrite(fname, [4 * rL_obs, Stock_Effect, Flow_Effect, Others], 'Fig7_Alt2', 'A1');
 

