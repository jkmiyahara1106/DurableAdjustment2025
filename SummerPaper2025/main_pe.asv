% Finite difference implicit updating with investment risk
% Ken Miyahara 2025
% Based on Greg Kaplan 2024

clear;
close all;
addpath("Codes\")
addpath("Results\")
%%
% INITIALIZE VALUE FUNCTION

%Vguess = zeros(na,1);
%Vguess(:) = u(r.*agrid + (r_risk-r)*risky_share.*agrid)./rho;
param = SetParameters;
grids = MakeGrids(param);

% ITERATE ON VALUE FUNCTION
load('./Results/Vguess2.mat')
%V    = Vguess;
Vdiff = 1;
iter = 0;
dVf= 0*V;
dVb= 0*V;

while iter <= param.maxiter_hjb && Vdiff>param.tol_hjb
    iter = iter + 1;
    
    Vnew = UpdateHJB(V,param,grids);
   
    [Vdiff,ind_diff] = max(abs(Vnew-V));
    if param.Display >=1
        disp(['HJB iteration ' int2str(iter), ' diff: ' num2str(Vdiff)]);
    end

    V = Vnew;
end 

[Vnew,Ahjb,Akfe,adj_hazard] = UpdateHJB(V,param,grids);

% M = V./( (1+agrid) .^ (1-risk_aver) ) ;
% [max_val , ind_max_rat] = max(M);
% % define benefit of adjusting
% dval = (agrid-prop_cost+1).^(1-risk_aver).*max_val-V;


%% SOLVE KFE

if param.IterateKFE==0
    gvecadj = [Akfe'; ones(1,param.na)] \ [zeros(param.na,1); 1];
   
elseif param.IterateKFE==1    

    %initialize at ergodic income distribution at a=0, adjusting for non-uniform grids
    
    gvecadj = ones(param.na,1)/param.na;

    gdiff = 1;
    iter = 0;
    %iterate to convergence
    while iter <= param.maxiter_kfe && gdiff>param.tol_kfe
        iter = iter + 1;
        gvecadjnew = (speye(param.na) -param.delta_kfe.* Akfe') \ gvecadj;

        gdiff = max(abs(gvecadjnew-gvecadj));
        if Display >=1
            disp(['KFE iteration ' int2str(iter), ' diff: ' num2str(gdiff)]);
        end

        gvecadj = gvecadjnew;
    end    

end
gmat    = gvecadj./grids.adelta;
mean_a = sum(grids.agrid.*gvecadj);
freq_adj = sum(adj_hazard.*gvecadj);
dist_adj = gmat.*adj_hazard/freq_adj;
disp(mean_a)
disp(freq_adj)

% %% Find barriers
% wlow = agrid(find(dval<=0,1,'first'));
% what = agrid(ind_max_rat);
% whigh = agrid(find(dval<=0,1,'last'));
% 
% %% Define barrier region for focused plots
% barrier_region = [wlow, whigh];
% barrier_margin = 0.5 * (whigh - wlow); % 10% margin around barriers
% plot_xlim = [wlow - barrier_margin, whigh + barrier_margin];
% 
% %% MAKE PLOTS
% 
% %% Create a separate figure focusing only on the barrier region
% if MakePlots == 1
%     figure(1);
% 
%     % Find indices for barrier region
%     barrier_indices = 1:length(agrid);%agrid >= (wlow - barrier_margin) & agrid <= (whigh + barrier_margin);
%     agrid_barrier = agrid(barrier_indices);
% 
%     % Policy functions in barrier region
%     subplot(2,3,1);
%     plot(agrid_barrier, con(barrier_indices,1), 'b-', 'LineWidth', 2);
%     hold on;
%     yrange = ylim;
%     plot([wlow wlow], yrange, 'r--', 'LineWidth', 1.5);
%     plot([what what], yrange, 'g--', 'LineWidth', 1.5);
%     plot([whigh whigh], yrange, 'r--', 'LineWidth', 1.5);
%     xlim([borrow_lim 3])
%     ylim([con(1) con(find(agrid>3,1,'first'))])
%     hold off;
%     grid;
%     title('Consumption Policy');
%     xlabel('Wealth-to-Durable');
%     ylabel('Consumption');
% 
%     subplot(2,3,2);
%     plot(-log(agrid(ind_max_rat))+log(agrid_barrier), dist_adj(barrier_indices,1), 'b-', 'LineWidth', 2);
%     hold on;
%     %plot(agrid_barrier, zeros(sum(barrier_indices),1), 'k-', 'LineWidth', 0.5);
%     yrange = ylim;
%     plot(-log(agrid(ind_max_rat))+log([wlow wlow]), yrange, 'r--', 'LineWidth', 1.5);
%     plot(-log(agrid(ind_max_rat))+log([what what]), yrange, 'g--', 'LineWidth', 1.5);
%     plot(-log(agrid(ind_max_rat))+log([whigh whigh]), yrange, 'r--', 'LineWidth', 1.5);
%     xlim([-.8 1.2])
%     xticks([-0.8 -0.5 0 0.5 1.2])
%     hold off;
%     grid;
%     title('Distribution of $1/w$ adjustments',Interpreter='latex');
%     xlabel('Durable-to-wealth');
%     ylabel('Density');
% 
%     subplot(2,3,3);
%     plot(agrid_barrier, 1-exp(-adj_hazard(barrier_indices)/3), 'b-', 'LineWidth', 2);
%     hold on;
%     yrange = ylim;
%     plot([wlow wlow], yrange, 'r--', 'LineWidth', 1.5);
%     plot([what what], yrange, 'g--', 'LineWidth', 1.5);
%     plot([whigh whigh], yrange, 'r--', 'LineWidth', 1.5);
%     xlim([borrow_lim 3])
%     hold off;
%     grid;
%     title('Hazard of adjustment');
%     xlabel('Wealth-to-Durable');
%     ylabel('Probability (per month)');
% 
%     subplot(2,3,4);
%     plot(agrid_barrier, gmat(barrier_indices), 'b-', 'LineWidth', 2);
%     hold on;
%     yrange = ylim;
%     plot([wlow wlow], yrange, 'r--', 'LineWidth', 1.5);
%     plot([what what], yrange, 'g--', 'LineWidth', 1.5);
%     plot([whigh whigh], yrange, 'r--', 'LineWidth', 1.5);
%     xlim([borrow_lim 3])
%     hold off;
%     grid;
%     title('Wealth-to-durable Distribution');
%     xlabel('Wealth-to-Durable');
%     ylabel('Density');
% 
%     sgtitle('Durable Adjustment Model', ...
%         'FontSize', 14, 'FontWeight', 'bold');
% 
%     subplot(2,3,5);
%     plot(log(agrid-prop_cost+1) - log(1+agrid(ind_max_rat)), dist_adj(barrier_indices,1), 'b-', 'LineWidth', 2);
%     hold on;
%     yrange = ylim;
%     plot(log([wlow wlow]-prop_cost+1)-log(1+agrid(ind_max_rat)), yrange, 'r--', 'LineWidth', 1.5);
%     plot(log([what what]-prop_cost+1)-log(1+agrid(ind_max_rat)), yrange, 'g--', 'LineWidth', 1.5);
%     plot(log([whigh whigh]-prop_cost+1)-log(1+agrid(ind_max_rat)), yrange, 'r--', 'LineWidth', 1.5);
%     xlim([-0.5 0.5])
%     hold off;
%     grid;
%     title('Distribution of durable adjustment');
%     xlabel('Durable adjustment');
%     ylabel('Density');
% 
%     sgtitle('Durable Adjustment Model', ...
%         'FontSize', 14, 'FontWeight', 'bold');
% 
%     %print('DurableAdjModel.png', '-dpng', '-r300');
% end
% 
% sum((log(agrid-prop_cost+1) - log(1+agrid(ind_max_rat))).*dist_adj)