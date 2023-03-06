function [objective,cancercost,combatantcost] = cost_to_go(sol)
%solution is a Nx6 matrix of N time points, where
% sol(n,1) is the time at the n-th time point
% sol(n,2) is the helper cell concentration
% sol(n,3) is the combatant cell concentration
% sol(n,4) is the cancer cell concentration
% sol(n,5) is the dendritic cell concentration
% sol(n,6) is the IL-2 food concentration

%interpolate incoming data
% interpolationgrid = transpose(sol(1,1):0.01:sol(size(sol,1),1));
% helper_solution = interp1(sol(:,1),sol(:,2),interpolationgrid,"spline");
% combatant_solution = interp1(sol(:,1),sol(:,3),interpolationgrid,"spline");
% cancer_solution = interp1(sol(:,1),sol(:,4),interpolationgrid,"spline");
% dendritic_solution = interp1(sol(:,1),sol(:,5),interpolationgrid,"spline");
% food_solution = interp1(sol(:,1),sol(:,6),interpolationgrid,"spline");

%test

combatant_square_sum = trapz(sol(:,1),sol(:,3).^2);
cancer_square_sum = trapz(sol(:,1),sol(:,4).^2);


objective = combatant_square_sum + cancer_square_sum;
%objective = sol(end,4);
cancercost = cancer_square_sum;
combatantcost = combatant_square_sum;

%objective = rand;
end