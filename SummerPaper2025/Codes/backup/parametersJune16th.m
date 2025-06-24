%% PARAMETERS

% preferences
risk_aver   = 3.5;
rho         = 0.025/4; %0.01; %0.005275; %quarterly 
alpha       = 0.8; %share in non-durable
zeta        = 1-alpha*(1-risk_aver);

%returns
r           = 0.01/4; %quarterly
r_risk      = 0.07/4; %quarterly
sigma2      = 0.1655^2; %quarterly
sigma       = sqrt(sigma2);

% transaction costs
prop_cost = 0.05;
adj_arriv = 365/4; % set the hjb as initial guess
risky_share = 0.3;

% asset grids
na          = 500; %100;
amax        = 20; %400; 
borrow_lim  = 0.01;
agrid_par   = 1; %1 for linear, 0 for L-shaped

% computation
maxiter_hjb = 200;
tol_hjb     = 1.0e-7;
delta_hjb   = 1.0e8;
mindV = 1.0e-8;

maxiter_kfe = 1000;
tol_kfe     = 1.0e-8;
delta_kfe   = 1.0e8;