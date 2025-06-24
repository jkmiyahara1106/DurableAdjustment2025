% Finite difference implicit updating with investment risk
% Ken Miyahara 2025
% Based on Greg Kaplan 2024

clear;
close all;

addpath('Results')
addpath('Codes')
%addpath('Figures\')

%% OPTIONS
Display     = 2;
MakePlots   = 1;
ComputeMPC  = 0;
IterateKFE  = 0;

%% PARAMETERS

% preferences
risk_aver   = 3.5;
rho         = 0.025/4; %0.01; %0.005275; %quarterly 
alpha       = 0.8; %share in non-durable
zeta        = 1-alpha*(1-risk_aver);

%returns
r           = 0.01/4; %quarterly
r_risk      = 0.07/4; %quarterly
sigma2      = 0.1655^2/4; %quarterly
sigma       = sqrt(sigma2);
risky_share = 0.7;
margin_prop = 3;

% transaction costs
prop_cost = 0.05;
adj_arriv_u = 1; % adjust up opportunities
adj_arriv_d = 1; % adjust down opportunities 
psi_val_u = [ 0; 5; 10; 20 ; 30 ; 50; 100];
pmf_psi_u = [ 0; 0.2; 0.2; 0.2; 0.2; 0.1; 0.1];
cdf_psi_u = cumsum(pmf_psi_u);
psi_cum_u = cumsum(psi_val_u.*pmf_psi_u);

dist_up.vals = psi_val_u*100;
dist_up.cdf = cdf_psi_u;
dist_up.costCum =  psi_cum_u*100;

dist_down.vals = psi_val_u*100;
dist_down.cdf = cdf_psi_u;
dist_down.costCum =  psi_cum_u*100;

% asset grids
nx          = 100; %100;
xmax        = 4;%log(20); %400; 
borrow_lim  = -2;%log(0.1);
agrid_par   = 1; %1 for linear, 0 for L-shaped

% computation
maxiter_hjb = 100;
tol_hjb     = 1.0e-7;
delta_hjb   = 1.0e8;
mindV = 1.0e-10;

maxiter_kfe = 1000;
tol_kfe     = 1.0e-8;
delta_kfe   = 1.0e8;

%mpc options
cumconT     = 1; %duration to measure cumulative consumption
delta_mpc   = 0.01; %time steps for cumulative consumption
mpcamount1  = 1.0e-10; %approximate thoeretical mpc
mpcamount2  = 0.10; % ten percent of average income: approx $5000
    

%% SET UP GRIDS

% assets
xgrid = linspace(0,1,nx)';
xgrid = xgrid.^(1./agrid_par);
xgrid = borrow_lim + (xmax-borrow_lim).*xgrid;

% asset grid spacing: for partial derivatives
dxgrid  =  diff(xgrid);
%dxgridf = [dxgrid; dxgrid(nx-1)];
%dxgridb = [dxgrid(1); dxgrid];
dx      = dxgrid(1);

% trapezoidal rule: for KFE and moments
xdelta = ones(nx,1);
xdelta = xdelta*dxgrid(1);

%% UTILITY FUNCTION

if risk_aver <= 1
    disp('Risk aversion must be larger than 1')
else    
    u = @(c,x) (c.*exp(x)).^(1-zeta)./(1-risk_aver);
end    

u1 = @(c,x) alpha*c.^(-zeta).*exp((1-zeta)*x); % marginal utility
u1inv = @(vp,x) (vp.*exp(-(1-zeta)*x)/alpha).^(-1./zeta); % inverse of marginal utility
opt_expos = @(dv,d2v) min((r_risk-r).*sigma.*dv./(sigma2.*max(dv-d2v,mindV)),margin_prop*sigma);
opt_expos_uncon = @(dv,d2v) (r_risk-r).*sigma.*dv./(sigma2.*max(dv-d2v,mindV));
int_income = @(expos) r + (r_risk-r)*expos/sigma - expos.^2/2;
bdgt_const = @(dv,d2v,x) int_income(opt_expos_uncon(max(dv,mindV),d2v)) - u1inv(max(dv,mindV),x) ;

%% INITIALIZE VALUE FUNCTION

Vguess = zeros(nx,1);
c0 = r;
Vguess(:) = u(c0,xgrid)./rho;

% ITERATE ON VALUE FUNCTION
load('Vguess.mat')
%V    = Vguess;
Vdiff = 1;
iter = 0;
dVf= 0*V;
dVb= 0*V;

while iter <= maxiter_hjb && Vdiff>tol_hjb
    iter = iter + 1;
    
    % forward difference
    dVf(1:nx-1) = max((V(2:nx)-V(1:nx-1))/dx,mindV);
    c0 = 5*r;
    dVf(nx) = max(u1(c0,xgrid(nx)),mindV); %state constraint
    % backward difference
    dVb(2:nx) = max((V(2:nx)-V(1:nx-1))./dx,mindV);
    c0 = r;
    dVb(1) = max(u1(c0,xgrid(1)),mindV); %state constraint

    % Central difference second derivative
    d2V = min((dVf - dVb) /dx,-mindV);
    
    % forward portfolio 
    exposure_f = opt_expos(dVf,d2V);
    
    % backward portfolio
    exposure_b = opt_expos(dVb,d2V);

    %consumption and savings with forward difference
    conf = u1inv(dVf,xgrid);
    driftf = int_income(exposure_f) - conf;
    Hf = u(conf,xgrid) + dVf.*driftf + d2V.*exposure_f.^2/2;
    
    %consumption and savings with backward difference
    conb = u1inv(dVb,xgrid);
    driftb = int_income(exposure_b) -conb;
    Hb = u(conb,xgrid) + dVb.*driftb + d2V.*exposure_b.^2/2;
    
    %consumption and derivative with adot = 0
    dV0 = dVf*0;
    %find dV consistent with zero drift
    for it = 1:nx
    dV0(it) = fzero( @(dv) bdgt_const(dv,d2V(it),xgrid(it)),dVf(it));
    end
    exposure_0 = opt_expos(dV0,d2V);
    con0 = int_income(exposure_0);
   %dV0 = u1(con0,xgrid);
    H0 = u(con0,xgrid) +  d2V.*exposure_0.^2/2;
    
    % choice of forward or backward differences based on sign of drift    
    Ineither = (1-(driftf>0)) .* (1-(driftb<0));
    Iunique = (driftb<0).*(1-(driftf>0)) + (1-(driftb<0)).*(driftf>0);
    Iboth = (driftb<0).*(driftf>0);
    Ib = Iunique.*(driftb<0).*(Hb>H0) + Iboth.*(Hb>Hf).*(Hb>H0);
    If = Iunique.*(driftf>0).*(Hf>H0) + Iboth.*(Hf>Hb).*(Hf>H0);
    I0 = 1-Ib-If;
    
    %consumption, savings and utility
    con   = conf.*If + conb.*Ib + con0.*I0;
    expos = exposure_f.*If + exposure_b.*Ib + exposure_0.*I0;
    dV = dVb.*Ib + dVf.*If + dV0.*I0;
    drift   = driftf.*If + driftb.*Ib; 

    util  = u(con,xgrid);

    % Impulse Hamiltonian
    % find optimum and maximized continuation value
    M = V./( (1+exp(xgrid)) .^ (1-risk_aver) ) ;
    [max_val , ind_max] = max(M);
    % define benefit of adjusting
    dval = (exp(xgrid)-prop_cost+1).^(1-risk_aver).*max_val-V;
    % compute impulse Hamiltonian
    
    % compute adjustment hazard
    up = xgrid<xgrid(ind_max); %adjust up
    down = xgrid>xgrid(ind_max); %adjust down
    
    ind_bar_d = find_ind(dist_down.vals,dval); %index that truncates dist
    ind_bar_u = find_ind(dist_up.vals,dval); %index that truncates dist
    
    adj_hazard = adj_arriv_d.*down.*dist_down.cdf(ind_bar_d) + ...
        adj_arriv_u.*up.*dist_up.cdf(ind_bar_u);
    
    util_cost = adj_arriv_d.*down.*dist_down.costCum(ind_bar_d) + ...
        adj_arriv_u.*up.*dist_up.costCum(ind_bar_u);

    %construct A matrix: tri-diagonal elements
    Alowdiag = -Ib.*driftb./dx + expos.^2/2/dx^2;
    Adiag = -If.*driftf./dx + Ib.*driftb./dx-expos.^2/dx^2;
        Adiag(1) = Adiag(1) + expos(1).^2/2/dx^2 - Ib(1).*driftb(1)./dx;
        Adiag(end) = Adiag(end) + expos(end).^2/2/dx^2 + If(end).*driftf(end)./dx;
    Aupdiag = If.*driftf./dx + expos.^2/2/dx^2;

    %use spdiags to create A matrix 
    Ahjb = spdiags(Adiag(:),0,nx,nx) + ...
                spdiags(Alowdiag(2:nx),-1,nx,nx) + ...
                spdiags([0;Aupdiag(1:nx-1)],1,nx,nx);
    
    Akfe = Ahjb - spdiags(adj_hazard,0,nx,nx);
    Akfe(:,ind_max) = Akfe(:,ind_max) + adj_hazard;

    Ahjb = Ahjb - spdiags(adj_hazard,0,nx,nx);
    Ahjb(:,ind_max) = Ahjb(:,ind_max) + adj_hazard.*(exp(xgrid)-prop_cost+1).^(1-risk_aver)/((1+exp(xgrid(ind_max))) ^ (1-risk_aver) );

    if max(abs(sum(Akfe,2)))>10^(-8)
        disp('Ill-posed Transition matrix')
        return
    end

    B = (rho + 1./delta_hjb)*speye(nx) - Ahjb;
        
    % solve linear system
    Vnew = B \ (util - util_cost + V./delta_hjb);
   
    Vdiff = max(abs(Vnew-V));
    if Display >=1
        disp(['HJB iteration ' int2str(iter), ' diff: ' num2str(Vdiff)]);
    end

    V = Vnew;
end 
M = V./( (1+exp(xgrid)) .^ (1-risk_aver) ) ;
[max_val , ind_max] = max(M);
% define benefit of adjusting
dval = (exp(xgrid)-prop_cost+1).^(1-risk_aver).*max_val-V;


%% SOLVE KFE

if IterateKFE==0
    gvecadj = [Akfe'; ones(1,nx)] \ [zeros(nx,1); 1];
   
elseif IterateKFE==1    

    %initialize at ergodic income distribution at a=0, adjusting for non-uniform grids
    
    gvecadj = ones(nx,1)/nx;

    gdiff = 1;
    iter = 0;
    %iterate to convergence
    while iter <= maxiter_kfe && gdiff>tol_kfe
        iter = iter + 1;
        gvecadjnew = (speye(nx) - delta_kfe.* Akfe') \ gvecadj;

        gdiff = max(abs(gvecadjnew-gvecadj));
        if Display >=1
            disp(['KFE iteration ' int2str(iter), ' diff: ' num2str(gdiff)]);
        end

        gvecadj = gvecadjnew;
    end    

end
gmat    = gvecadj./xdelta;
mean_a = sum(xgrid.*gvecadj);
freq_adj = sum(adj_hazard.*gvecadj);
dist_adj = adj_hazard.*gvecadj/freq_adj;
disp(mean_a)
disp(freq_adj)

%% MAKE PLOTS
if MakePlots ==1 
    figure(1);
    
    % consumption policy function
    subplot(2,4,1);
    plot(xgrid,drift(:),'b-','LineWidth',1);
    grid;
    xlim([borrow_lim xmax]);
    title('Drift Policy Function');
    %legend('Lowest income state','Highest income state');

    % savings policy function
    subplot(2,4,2);
    plot(xgrid,expos(:)/sigma,'b-','LineWidth',1);
    hold on;
    xline(xgrid(ind_max), 'r--') 
    %plot(xgrid,risky_share*ones(nx,1),'k','LineWidth',0.5);
    hold off;
    grid;
    xlim([borrow_lim xmax]);
    title('Optimal Portfolio Function');
    
    % consumption policy function: zoomed in
    subplot(2,4,3:4);
    plot(log(exp(xgrid)-prop_cost+1) - log(exp(xgrid(ind_max))+1),dist_adj,'-b','LineWidth',2);
    grid;
    %xlim(borrow_lim + [0 1]);
    title('Distribution of log Durable Adjustment');
    
    % % savings policy function: zoomed in
    % subplot(2,4,4);
    % plot(xgrid,drift(:,1),'o-b','LineWidth',2);
    % hold on;
    % plot(xgrid,zeros(nx,1),'k','LineWidth',0.5);
    % hold off;
    % grid;
    % xlim(borrow_lim + [0 1]);
    % title('Savings: Zoomed');
    
    subplot(2,4,5)
    plot(xgrid,adj_hazard,'b-','LineWidth',1);
    grid;
    xlim([borrow_lim xmax]);
    title('Hazard function');
    
    subplot(2,4,6:8)
    plot(xgrid,gmat,'b-','LineWidth',1);
    grid;
    xlim([borrow_lim xmax]);
    title('Distribution');

end

