% Finite difference implicit updating with investment risk
% Ken Miyahara 2025
% Based on Greg Kaplan 2024

clear;
close all;

addpath('Results\')
addpath('Codes\')
addpath('Figures\')

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
sigma2      = 0.1655^2; %quarterly
sigma       = sqrt(sigma2);

% transaction costs
prop_cost = 0.05;
adj_arriv = 365/4; % set the hjb as initial guess
risky_share = 0.25;

% asset grids
na          = 500; %100;
amax        = 10; %400; 
borrow_lim  = 0.00000001;
agrid_par   = 1; %1 for linear, 0 for L-shaped

% computation
maxiter_hjb = 200;
tol_hjb     = 1.0e-7;
delta_hjb   = 1.0e8;
mindV = 1.0e-8;

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
agrid = linspace(0,1,na)';
agrid = agrid.^(1./agrid_par);
agrid = borrow_lim + (amax-borrow_lim).*agrid;

% asset grid spacing: for partial derivatives
dagrid  =  diff(agrid);
dagridf = [dagrid; dagrid(na-1)];
dagridb = [dagrid(1); dagrid];
da      = dagrid(1);

% trapezoidal rule: for KFE and moments
adelta = ones(na,1);
adelta = adelta*dagrid(1);

%% UTILITY FUNCTION

if risk_aver==1
    u = @(c)alpha*log(c);
else    
    u = @(c)(c.^(1-zeta)-1)./(1-risk_aver);
end    
u1 = @(c) alpha*c.^(-zeta);
u1inv = @(u) (u/alpha).^(-1./zeta);

%% INITIALIZE VALUE FUNCTION

%Vguess = zeros(na,1);
%Vguess(:) = u(r.*agrid + (r_risk-r)*risky_share.*agrid)./rho;

% ITERATE ON VALUE FUNCTION
load('./Results/Vguess.mat')
%V    = Vguess;
Vdiff = 1;
iter = 0;
dVf= 0*V;
dVb= 0*V;

while iter <= maxiter_hjb && Vdiff>tol_hjb
    iter = iter + 1;
    
    % forward difference
    dVf(1:na-1) = max((V(2:na)-V(1:na-1))/da,mindV);
    exposure_f = risky_share*agrid*sigma;
    dVf(na) = max(u1(r.*agrid(na)+(r_risk-r)*exposure_f(na)/sigma),mindV); %state constraint

    % backward difference
    dVb(2:na) = max((V(2:na)-V(1:na-1))./da,mindV);
    exposure_b = risky_share*agrid*sigma;
    dVb(1) = max(u1(r.*agrid(1) + (r_risk-r)*exposure_b(1)/sigma),mindV); %state constraint
    
    % Central difference second derivative
    d2V = min((dVf - dVb) /da,-mindV);

    %consumption and savings with forward difference
    conf = u1inv(dVf);
    savf = r.*agrid  + (r_risk-r)*exposure_f/sigma - conf;
    Hf = u(conf) + dVf.*savf;
    
    %consumption and savings with backward difference
    conb = u1inv(dVb);
    savb = r.*agrid  + (r_risk-r)*exposure_b/sigma - conb;
    Hb = u(conb) + dVb.*savb;
    
    %consumption and derivative with adot = 0
    exposure_0 = risky_share*agrid*sigma;
    con0 = r.*agrid  + (r_risk-r)*exposure_0/sigma;
    dV0 = u1(con0);
    H0 = u(con0);
    
    % choice of forward or backward differences based on sign of drift    
    Ineither = (1-(savf>0)) .* (1-(savb<0));
    Iunique = (savb<0).*(1-(savf>0)) + (1-(savb<0)).*(savf>0);
    Iboth = (savb<0).*(savf>0);
    Ib = Iunique.*(savb<0).*(Hb>H0) + Iboth.*(Hb>Hf).*(Hb>H0);
    If = Iunique.*(savf>0).*(Hf>H0) + Iboth.*(Hf>Hb).*(Hf>H0);
    I0 = 1-Ib-If;
    
    %consumption, savings and utility
    con   = conf.*If + conb.*Ib + con0.*I0;
    expos = exposure_f.*If + exposure_b.*Ib + exposure_0.*I0;
    sav   = savf.*If + savb.*Ib;    
    util  = u(con);

    % Impulse Hamiltonian
    iH_fun = @(dv) adj_arriv*(dv>0).*dv;
    % find optimum and maximized continuation value
    M = V./( (1+agrid) .^ (1-risk_aver) ) ;
    [max_val , ind_max_rat] = max(M);
    % define benefit of adjusting
    dval = (agrid-prop_cost+1).^(1-risk_aver).*max_val-V;
    % compute impulse Hamiltonian
    iH = iH_fun(dval);
    % compute adjustment hazard
    adj_hazard = adj_arriv.*(dval>0);

    %construct A matrix: tri-diagonal elements
    Alowdiag = -Ib.*savb./da + expos.^2/2/da^2;
    Adiag = -If.*savf./da + Ib.*savb./da-expos.^2/da^2;
        Adiag(1) = Adiag(1) + expos(1).^2/2/da^2;
        Adiag(end) = Adiag(end) + expos(end).^2/2/da^2;
    Aupdiag = If.*savf./da + expos.^2/2/da^2;

    %use spdiags to create A matrix 
    Ahjb = spdiags(Adiag(:),0,na,na) + ...
                spdiags(Alowdiag(2:na),-1,na,na) + ...
                spdiags([0;Aupdiag(1:na-1)],1,na,na);
    
    Akfe = Ahjb - spdiags(adj_hazard,0,na,na);
    Akfe(:,ind_max_rat) = Akfe(:,ind_max_rat) + adj_hazard;

    Ahjb = Ahjb - spdiags(adj_hazard,0,na,na);
    Ahjb(:,ind_max_rat) = Ahjb(:,ind_max_rat) + adj_hazard.*(agrid-prop_cost+1).^(1-risk_aver)/((1+agrid(ind_max_rat)) ^ (1-risk_aver) );

    if max(abs(sum(Akfe,2)))>10^(-8)
        disp('Ill-posed Transition matrix')
        return
    end

    B = (rho + 1./delta_hjb)*speye(na) - Ahjb;
        
    % solve linear system
    Vnew = B \ (util + V./delta_hjb);
   
    Vdiff = max(abs(Vnew-V));
    if Display >=1
        disp(['HJB iteration ' int2str(iter), ' diff: ' num2str(Vdiff)]);
    end

    V = Vnew;
end 
M = V./( (1+agrid) .^ (1-risk_aver) ) ;
[max_val , ind_max_rat] = max(M);
% define benefit of adjusting
dval = (agrid-prop_cost+1).^(1-risk_aver).*max_val-V;


%% SOLVE KFE

if IterateKFE==0
    gvecadj = [Akfe'; ones(1,na)] \ [zeros(na,1); 1];
   
elseif IterateKFE==1    

    %initialize at ergodic income distribution at a=0, adjusting for non-uniform grids
    
    gvecadj = ones(na,1)/na;

    gdiff = 1;
    iter = 0;
    %iterate to convergence
    while iter <= maxiter_kfe && gdiff>tol_kfe
        iter = iter + 1;
        gvecadjnew = (speye(na) - delta_kfe.* Akfe') \ gvecadj;

        gdiff = max(abs(gvecadjnew-gvecadj));
        if Display >=1
            disp(['KFE iteration ' int2str(iter), ' diff: ' num2str(gdiff)]);
        end

        gvecadj = gvecadjnew;
    end    

end
gmat    = gvecadj./adelta;
mean_a = sum(agrid.*gvecadj);
freq_adj = sum(adj_hazard.*gvecadj);
disp(mean_a)
disp(freq_adj)

%% MAKE PLOTS
if MakePlots ==1 
    figure(1);
    
    % consumption policy function
    subplot(2,4,1);
    plot(agrid,con(:),'b-','LineWidth',1);
    grid;
    xlim([borrow_lim amax]);
    title('Consumption Policy Function');
    %legend('Lowest income state','Highest income state');

    % savings policy function
    subplot(2,4,2);
    plot(agrid,sav(:),'b-','LineWidth',1);
    hold on;
    plot(agrid,zeros(na,1),'k','LineWidth',0.5);
    hold off;
    grid;
    xlim([borrow_lim amax]);
    title('Savings Policy Function');
    
    % consumption policy function: zoomed in
    subplot(2,4,3);
    plot(agrid,con(:,1),'o-b','LineWidth',2);
    grid;
    xlim(borrow_lim + [0 1]);
    title('Consumption: Zoomed');
    
    % savings policy function: zoomed in
    subplot(2,4,4);
    plot(agrid,sav(:,1),'o-b','LineWidth',2);
    hold on;
    plot(agrid,zeros(na,1),'k','LineWidth',0.5);
    hold off;
    grid;
    xlim(borrow_lim + [0 1]);
    title('Savings: Zoomed');
    
    subplot(2,4,5)
    plot(agrid,V,'b-','LineWidth',1);
    grid;
    xlim([borrow_lim amax]);
    title('Value Function');
    
    subplot(2,4,6:8)
    plot(agrid,gmat,'b-','LineWidth',1);
    grid;
    xlim([borrow_lim amax]);
    title('Distribution');

end

