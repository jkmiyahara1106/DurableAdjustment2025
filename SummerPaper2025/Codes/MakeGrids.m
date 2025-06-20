function grids = MakeGrids(param)
%MakeGrids Make grids

% assets
agrid = linspace(0,1,param.na)';
agrid = agrid.^(1./param.agrid_par);
agrid = param.borrow_lim + (param.amax-param.borrow_lim).*agrid;

% asset grid spacing: for partial derivatives
dagrid  =  diff(agrid);
dagridf = [dagrid; dagrid(param.na-1)];
dagridb = [dagrid(1); dagrid];
da      = dagrid(1);

% trapezoidal rule: for KFE and moments
adelta = ones(param.na,1);
adelta = adelta*dagrid(1);

grids.agrid = agrid;
grids.dagridf = dagridf;
grids.dagridb = dagridb;
grids.dagrid = dagrid;
grids.da = da;
grids.adelta = adelta;

end