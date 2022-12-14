function [outputArg1] = therapy_cancer(horizon,initCond)
%THERAPY Summary of this function goes here
%   Detailed explanation goes here

%H params a(1)
a0 = 10^(-4);
b0 = 0.005;
c0 = 10;
d0 = 10^(-2);
f0 = 1;

%C params a(2)
a1 = 10^(-4);
b1 = 0.005;
c1 = 10;
d1 = 10^(-2);
f1 = 1;

%M params a(3)
d2 = 0.02;
e2 = 0.1;
f2 = 1;

%D params a(4)
e3 = 0.1;

%I params a(5)
a4 = 10^(-2);
c4 = 10^(-7);
e4 = 10^(-2);

%initial conditions
% H0 = 0.01;
% C0 = 0.01;
% M0 = 1;
% D0 = 0;
% I0 = 0.01;
H0=initCond(1);
C0=initCond(2);
M0=initCond(3);
D0=initCond(4);
I0=initCond(5);

%time discretization
%tspan = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
tspan = linspace(1,horizon,horizon);

f = @(t,a) [
    a0 - b0*a(1) + c0*a(4)*d0*a(1)*(1-a(1)/f0);
    a1 - b1*a(2) + c1*a(5)*(a(3)+a(4))*d1*a(2)*(1-a(2)/f1);
    d2*a(3)*(1-a(3)/f2)-e2*a(3)*a(2);
    -e3*a(4)*a(2);
    a4*a(1)*a(4)-c4*a(2)*a(5)-e4*a(5);
    ];
g = @(t,a) [
    a0 - b0*a(1) + c0*a(4)*d0*a(1)*(1-a(1)/f0); %a(1) = H
    a1 - b1*a(2) + c1*a(5)*(a(3)+a(4))*d1*a(2)*(1-a(2)/f1); %a(2) = C
    d2*log(f2/a(3))*a(3)-e2*a(3)*a(2); %a(3) =M
    -e3*a(4)*a(2); % a(4) = D
    a4*a(1)*a(4)-c4*a(2)*a(5)-e4*a(5); % a(5) = I
    ];

[timelogistic,logistic] = ode45(f,tspan,[H0,C0,M0,D0,I0]);
[~,gompertz] = ode45(g,tspan,[H0,C0,M0,D0,I0]);

%[timelogistic,timegompertz,logistic(:,3),gompertz(:,3)]

% figure() 
% plot(timelogistic,logistic(:,1),'-o',timelogistic,gompertz(:,1),'.')
% title('H(t) via ODE45 - CD4+ (helper T-cell) concentration')
% xlabel('Time t (hours)')
% ylabel('Concentration')
% legend('Logistic Growth','Gompertz Growth')
% 
% figure()
% plot(timelogistic,logistic(:,2),'-o',timelogistic,gompertz(:,2),'.')
% title('C(t) via ODE45 - CD8+ (combatants) concentration')
% xlabel('Time t (hours)')
% ylabel('Concentration')
% legend('Logistic Growth','Gompertz Growth')


figure()
plot(timelogistic,log(logistic(:,3)),'-o',timelogistic,log(gompertz(:,3)),'.')
title('M(t) via ODE45 - Cancer concentration')
xlabel('Time t (hours)')
ylabel('Concentration')
legend('Logistic Growth','Gompertz Growth')
xlim([0,horizon])
ylim([-30,30])

% figure()
% plot(timelogistic,logistic(:,4),'-o',timelogistic,gompertz(:,4),'.')
% title('D(t) via ODE45 - dendritic cell concentration')
% xlabel('Time t (hours)')
% ylabel('Concentration')
% legend('Logistic Growth','Gompertz Growth')
% 
% figure()
% plot(timelogistic,logistic(:,5),'-o',timelogistic,gompertz(:,5),'.')
% title('I(t) via ODE45 - food concentration')
% xlabel('Time t (hours)')
% ylabel('Concentration')
% legend('Logistic Growth','Gompertz Growth')
 outputArg1 = [timelogistic,logistic]; 
end