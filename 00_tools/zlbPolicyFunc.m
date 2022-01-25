function [PHI1, PHI2] = zlbPolicyFunc(fun, zlb, TR, R)
%
% computes policy functions under the commitment to the ZLB
%
% Inputs :
%   fun        : 'FILENAME_dynamic.m'
%                 -- Dynare automatically generates this matlab file.
%   zlb        : # of quarters on which the monetary policy commits to the ZLB
%   TR         : id # of the monetary policy equation
%   R          : id # of short-term interest rate  
%
% Outputs :
%   PHI1, PHI2 : Coefficient matixes of policy functions
%
% ...................................................................
% Copyright 2018-2020 by Nao Sudo and Masaki Tanaka (Bank of Japan)
%

%% ------------------------------------------------------------------------

%% 0. Load global structure
 global M_ oo_
 
%% 1. Load the policy function without ZLB 
    zlb     = zlb + 1;
    nExo    = M_.exo_nbr;
    nEndo   = M_.endo_nbr;
    nBwd    = M_.nspred;
    nFwd    = M_.nsfwrd;
    nAll    = nEndo + nBwd + nFwd;
    Lag     = M_.lead_lag_incidence(1:2,:);
    Lead    = M_.lead_lag_incidence(2:3,:) - max(M_.lead_lag_incidence(1,:));
    order1  = oo_.dr.inv_order_var;
    order2  = oo_.dr.state_var;
    PHI1    = zeros(nEndo + nFwd, nEndo + nFwd, zlb);
    PHI2    = zeros(nEndo + nFwd, nExo, zlb);

    for i = 1 : nEndo + nFwd
        if i < nEndo
            for j = 1 : size(oo_.dr.ghx, 2)
                PHI1(i, order2(j), zlb) = oo_.dr.ghx(order1(i), j);
            end
        elseif i == nEndo
            for j = 1 : size(oo_.dr.ghx, 2)
                PHI1(i, order2(j), zlb) = oo_.dr.ghx(order1(i), j);
            end
            ghx = PHI1(1:nEndo, 1:nEndo, zlb)*PHI1(1:nEndo, 1:nEndo, zlb); 
        else
            [~, col] = find(Lead == i);
            i0 = Lead(1, col);  
            PHI1(i, 1:nEndo , zlb) = ghx(i0, :);          
        end
    end
    
    for i = 1 : size(oo_.dr.ghu, 1)
        for j = 1 : size(oo_.dr.ghu, 2)
            PHI2(i, j, zlb) = oo_.dr.ghu(order1(i), j);
        end  
    end
    
%% 2. Rewrite model as canonical form
%--------------------------------------------------------------------
%[Canonical form]
%  GAMMA4 * z = GAMMA1 * z(-1) + GAMMA2 * eps + GAMMA3 * eta
%   --   z    [nEndo + nFwd by 1]           : vector of state variables
%   --  eps   [nExo by 1]                   : vector of exo. innovations
%   --  eta   [nFwd by 1]                   : vector of forcast errors
%   -- GAMMA4 [nEndo + nFwd by nEndo + nFwd]: coef. matrix on z
%   -- GAMMA1 [nEndo + nFwd by nEndo + nFwd]: coef. matrix on z(-1)
%   -- GAMMA2 [nEndo + nFwd by nExo]        : coef. matrix on eps
%   -- GAMMA3 [nEndo + nFwr by nFwd]        : coef. matrix on eta
%--------------------------------------------------------------------
  % Running 'FILENAME_dynamic.m' to get a jacobian matrix
  [~, jcb] = feval(fun, zeros(1, nAll), zeros(1, nExo), M_.params, oo_.steady_state, 1);
  
  % Decomposing the jacobian into GAMMA1, ..., GAMMA4 
  GAMMA4   = jcb(: , nBwd+1 : nAll);
  GAMMA1   = zeros(nEndo + nFwd, nEndo + nFwd);
      for i = 1: nBwd
          [~, col] = find(Lag == i);
          GAMMA1(1:nEndo ,col) = -jcb(: , i);   
      end     
  GAMMA2   = -jcb(: , nAll+1 : nAll+nExo);
  GAMMA3   = zeros(nEndo + nFwd, nFwd);
  
  % Equations for forcast errors
  for i = 1 : nFwd
      raw      = size(jcb,1) + i;
      id1      = nEndo + i;
      [~, col] = find(Lead == id1);
      id2      = Lead(1,col); 
      
      GAMMA4(raw, id2) = 1;
      GAMMA1(raw, id1) = 1;
      GAMMA2(raw, :)   = 0;
      GAMMA3(raw, i)   = 1;
  end
  
%% 3. Compute policy functions under ZLB constraint (perfect foresight solution)

  % Allow short-term interest rate to remain at zero
    GAMMA4(TR,:)   = zeros(1, nEndo+nFwd);
    GAMMA1(TR,:)   = zeros(1, nEndo+nFwd);
    GAMMA4(TR, R) = 1; 
    GAMMA1(TR, R) = 1;

  for k = 1 : zlb - 1
    GAMMA1j  = GAMMA1(nEndo+1 : nEndo+nFwd, : ) ;
    GAMMA4j  = GAMMA4(nEndo+1 : nEndo+nFwd, : ) ;
    GAMMA4jt = GAMMA4j * PHI1( : , : , zlb-k+1 ) - GAMMA1j; 
    
    GAMMA1t  = vertcat(GAMMA1(1 : nEndo, : ), zeros(nFwd, nEndo+nFwd));
    GAMMA2t  = vertcat(GAMMA2(1 : nEndo, : ), zeros(nFwd, nExo));
    GAMMA4t  = vertcat(GAMMA4(1 : nEndo, : ), GAMMA4jt);
  
    PHI1( : , : , zlb-k)    = GAMMA4t \ GAMMA1t;
    PHI2( : , : , zlb-k)    = GAMMA4t \ GAMMA2t;
  end

%% ------------------------------------------------------------------------
