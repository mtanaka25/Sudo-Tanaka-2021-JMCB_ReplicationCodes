function irfs = zlbIRF_for_Alt3(M_, PHI1, PHI2, varargin)
%
% computes impulse responses under the central bank's commitment to
% zero lower bound (modified for Alternative model 3).
%
%
% [Required Inputs]
%  M_         : Model structure
%               -- Dynare automatically generates
%  PHI1, PHI2 : Coefficient matrixes of policy function
%               -- see also PolicyFuncUnderZLB
%  
% [Options]
%  shocks     : names of shocks for which IRFs are computed
%               -- Default: All shocks
%  LengthIRF  : Length of impulse response fuctions computed
%               -- Default: 25
%  cutoff     : Threshold under which values of impulse response are
%               considered as zero
%               -- Default: 1e-10
%
% [Output]
%  irfs : Impulse response functions
% ...................................................................
% Copyright 2018-2020 by Nao Sudo and Masaki Tanaka (Bank of Japan)
%

%% ------------------------------------------------------------------------

%% 1. Check Inputs
if nargin < 3
    msg =['Error: Not enough input arguments. \n',...
          'At least, model structure and two coefficient matrixes of policy \n',...
          'function must be specified. For more detail, see IrfUnderZLB''s help.'];
    error(msg, class(msg));
end

o = struct;
for jo = 1:(length(varargin)/2) % Load Options
    o.(varargin{(jo-1)*2+1}) = varargin{jo*2};
end

for Exoj = 1:M_.exo_nbr
    o.ExoVars(1, Exoj) = cellstr(deblank(M_.exo_names(Exoj,:)));
end


%% 2. Apply default options to unspecified ones
oDefault.shocks    = o.ExoVars;
oDefault.LengthIRF = 25;
oDefault.cutoff    = 1e-10;

oList = fieldnames(oDefault);
for jo = 1:length(oList)
    oName = oList{jo};
    if ~isfield(o,oName)
        o.(oName) = oDefault.(oName);
    end
end

%% 3. Comupte Impulse Response Functions
nshock    = size(o.shocks, 2);
irfs      = zeros(M_.endo_nbr + M_.nsfwrd, o.LengthIRF, nshock);
LengthZLB = size(PHI1, 3);

for j = 1 : nshock
    irfj = zeros(M_.endo_nbr + M_.nsfwrd, o.LengthIRF + 1);
    
    if isfield(o,'shockMat')
        e2 = o.shockMat;
    else
        sj   = o.shocks(1,j);
        e2   = zeros(M_.exo_nbr, o.LengthIRF + 1);
        e2(strcmp(sj, o.ExoVars), 1) = M_.Sigma_e(strcmp(sj, o.ExoVars), strcmp(sj, o.ExoVars))^(1/2);
    end

    for k = 1:o.LengthIRF
      if k < LengthZLB
         irfj(:, k + 1) = PHI1(:, : , k) * irfj(:, k) + PHI2(:, : , k) * e2(:, k);
      else
         irfj(:, k + 1) = PHI1(:, : , LengthZLB) * irfj(:, k) + PHI2(:, : , LengthZLB) * e2(:, k);
      end  
    end
    
    irfs(:, :, j) = irfj(:, 2:o.LengthIRF+1);
end
    irfs(abs(irfs) <= o.cutoff) = 0;

%% ------------------------------------------------------------------------
