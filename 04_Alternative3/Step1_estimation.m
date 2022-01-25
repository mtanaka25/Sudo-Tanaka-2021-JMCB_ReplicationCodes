% Step1_estimation.m
%
% Run DYNARE and estimate the model. 
%
% ...................................................................
% Copyright 2018-2020 by Nao Sudo and Masaki Tanaka (Bank of Japan)
%
%%
clear all;
close all;

%%
%--------------------------------------------------------------------
% A. Option
%--------------------------------------------------------------------
EstCode     = 'Alternative3_estimation' ; % Name of the mod file: estimation

%%
%--------------------------------------------------------------------
% B. Estimate parameters with DYNARE
%--------------------------------------------------------------------
eval(sprintf('dynare %s noclearall nograph', EstCode));
