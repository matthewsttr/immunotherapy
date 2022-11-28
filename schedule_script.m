clc;
clear all;
close all;



horizon=50000;

%create a dosing schedule for x fraction of the horizon, with n equally
%spaced doses with specified level

n=floor(0.838*horizon/100);
x=0.8;
initCond = [0.1 ,0.1 ,1 ,9.5 ,0.01];
level=0.5;
scheduletimes = linspace(round(horizon/n),round(horizon*x),n);
scheduledoses=ones(1,n)*level;
schedule2 = horzcat(transpose(scheduletimes), transpose(scheduledoses));
%%%%

for i=1:30
omega=0.1*i;
n=floor(omega*horizon/100);
x=0.99;
scheduletimes = linspace(round(horizon/n),round(horizon*x),n);
scheduledoses=ones(1,n)*level;
%schedules{i} = horzcat(transpose(scheduletimes), transpose(scheduledoses));
scheduledoses = horzcat(transpose(scheduletimes), transpose(scheduledoses));
ithsolution=therapy_dosed(horizon,initCond,scheduledoses);
ithcost(i)=cost_to_go(ithsolution);
end

ithcost

figure()
plot(ithcost)
title('ithcost')
xlabel('n')
ylabel('cost')
legend('cost')


% therapy_cancer(horizon,initCond);
% schedule2sol = therapy_dosed(horizon, initCond,schedule2);
% 
% cost_to_go(schedule2sol)
