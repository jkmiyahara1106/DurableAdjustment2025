function param = SetParameters
%SetParameters Set parameters

param.Display     = 2;
param.MakePlots   = 1;
param.ComputeMPC  = 0;
param.IterateKFE  = 0;

% preferences
param.risk_aver   = 3.5;
param.rho         = 0.025/4; %0.01; %0.005275; %quarterly 
param.alpha       = 0.8; %share in non-durable
% zeta: 1 - exponent in homothetic formulation of utility
param.zeta        = 1-param.alpha*(1-param.risk_aver); 

%returns
param.r           = 0.01/4; %quarterly
param.r_risk      = 0.07/4; %quarterly
param.sigma2      = 0.1655^2; %quarterly
param.sigma       = sqrt(param.sigma2); %quarterly

% transaction costs
param.prop_cost = 0.001;
param.adj_arriv = 365/4; % set the hjb as initial guess
param.risky_share = 0.5;
param.upper_bound_cost_fun = 100000;

% asset grids
param.na          = 500; %100;
param.amax        = log(20); %400; 
param.borrow_lim  = log(0.001);
param.agrid_par   = 1; %1 for linear, 0 for L-shaped

% computation
param.maxiter_hjb = 300;
param.tol_hjb     = 1.0e-7;
param.delta_hjb   = 1.0e8;
param.mindV = 1.0e-8;

param.maxiter_kfe = 1000;
param.tol_kfe     = 1.0e-8;
param.delta_kfe   = 1.0e8;

%mpc options
param.cumconT     = 1; %duration to measure cumulative consumption
param.delta_mpc   = 0.01; %time steps for cumulative consumption
param.mpcamount1  = 1.0e-10; %approximate thoeretical mpc
param.mpcamount2  = 0.10; % ten percent of average income: approx $5000

%% UTILITY FUNCTION

if param.risk_aver==1
    param.u = @(c) param.alpha*log(c);
else    
    param.u = @(c)(c.^(1-param.zeta)-1)./(1-param.risk_aver);
end    
param.u1 = @(c) param.alpha*c.^(-param.zeta);
param.u1inv = @(u) (u/param.alpha).^(-1./param.zeta);
    
end