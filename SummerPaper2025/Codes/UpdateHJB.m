function [Vnew,Ahjb,Akfe,adj_hazard] = UpdateHJB(V,param,grids)
%UpdateHJB Update HJB


% forward difference
dVf(1:param.na-1,1) = max((V(2:param.na)-V(1:param.na-1))/grids.da,param.mindV);
exposure_f = param.risky_share*grids.agrid*param.sigma;
dVf(param.na,1) = max(param.u1(param.r.*grids.agrid(param.na)+(param.r_risk-param.r)*exposure_f(param.na)/param.sigma),param.mindV); %state constraint

% backward difference
dVb(2:param.na,1) = max((V(2:param.na)-V(1:param.na-1))./grids.da,param.mindV);
exposure_b = param.risky_share*grids.agrid*param.sigma;
dVb(1,1) = max(param.u1(param.r.*grids.agrid(1) + (param.r_risk-param.r)*exposure_b(1)/param.sigma),param.mindV); %state constraint

% Central difference second derivative
d2V = min((dVf - dVb) /grids.da,-param.mindV);

%consumption and savings with forward difference
conf = param.u1inv(dVf);
savf = param.r.*grids.agrid  + (param.r_risk-param.r)*exposure_f/param.sigma - conf;
Hf = param.u(conf) + dVf.*savf + d2V/2.*exposure_f.^2;

%consumption and savings with backward difference
conb = param.u1inv(dVb);
savb = param.r.*grids.agrid  + (param.r_risk-param.r)*exposure_b/param.sigma - conb;
Hb = param.u(conb) + dVb.*savb + d2V/2.*exposure_b.^2;

%consumption and derivative with adot = 0
exposure_0 = param.risky_share*grids.agrid*param.sigma;
con0 = param.r.*grids.agrid  + (param.r_risk-param.r)*exposure_0/param.sigma;
dV0 = param.u1(con0);
H0 = param.u(con0) + + d2V/2.*exposure_0.^2;

% choice of forward or backward differences based on sign of drift    
%Ineither = (1-(savf>0)) .* (1-(savb<0));
Iunique = (savb<0).*(1-(savf>0)) + (1-(savb<0)).*(savf>0);
Iboth = (savb<0).*(savf>0);
Ib = Iunique.*(savb<0).*(Hb>H0) + Iboth.*(Hb>Hf).*(Hb>H0);
If = Iunique.*(savf>0).*(Hf>H0) + Iboth.*(Hf>Hb).*(Hf>H0);
I0 = 1-Ib-If;

%consumption, savings and utility
con   = conf.*If + conb.*Ib + con0.*I0;
expos = exposure_f.*If + exposure_b.*Ib + exposure_0.*I0;
sav   = savf.*If + savb.*Ib;    
util  = param.u(con);

% Hazard function H'(u)
adj_hazard_fun = @(dv) param.adj_arriv * (dv <= 0) .* 0 + ...
         param.adj_arriv * (dv > 0 & dv <= param.upper_bound_cost_fun) .* dv/param.upper_bound_cost_fun  + ...
         param.adj_arriv * (dv > param.upper_bound_cost_fun) .* 1;

% find optimum and maximized continuation value
M = V./( (1+grids.agrid) .^ (1-param.risk_aver) ) ;
[max_val , ind_max_rat] = max(M);

% define benefit of adjusting
dval = (grids.agrid-param.prop_cost+1).^(1-param.risk_aver).*max_val-V;

% compute adjustment hazard
adj_hazard = adj_hazard_fun(dval);

% flow utility cost
util_cost = adj_hazard .*( (dval>0 & dval<=param.upper_bound_cost_fun).*dval/2  +...
    (dval>param.upper_bound_cost_fun).* param.upper_bound_cost_fun/2);

%construct A matrix: tri-diagoparam.nal elements
Alowdiag = -Ib.*sav./grids.da + expos.^2/2/grids.da^2;
Adiag = -If.*sav./grids.da + Ib.*sav./grids.da-expos.^2/grids.da^2;
    Adiag(1) = Adiag(1) + expos(1).^2/2/grids.da^2;
    Adiag(end) = Adiag(end) + expos(end).^2/2/grids.da^2;
Aupdiag = If.*sav./grids.da + expos.^2/2/grids.da^2;

%use spdiags to create A matrix 
Ahjb = spdiags(Adiag(:),0,param.na,param.na) + ...
            spdiags(Alowdiag(2:param.na),-1,param.na,param.na) + ...
            spdiags([0;Aupdiag(1:param.na-1)],1,param.na,param.na);

Akfe = Ahjb - spdiags(adj_hazard,0,param.na,param.na);
Akfe(:,ind_max_rat) = Akfe(:,ind_max_rat) + adj_hazard;

Ahjb = Ahjb - spdiags(adj_hazard,0,param.na,param.na);
Ahjb(:,ind_max_rat) = Ahjb(:,ind_max_rat) + adj_hazard.*(grids.agrid-param.prop_cost+1).^(1-param.risk_aver)/((1+grids.agrid(ind_max_rat)) ^ (1-param.risk_aver) );

if max(abs(sum(Akfe,2)))>10^(-8)
    disp('Ill-posed Transition matrix')
    return
end

B = (param.rho + 1./param.delta_hjb)*speye(param.na) - Ahjb;
    
% solve linear system
Vnew = B \ (util - util_cost + V./param.delta_hjb);
end