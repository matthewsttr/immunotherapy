function [cancerConc] = therapy_dosed(horizon, initCond, schedule)
% horizon   -    float:      represents final time to simulate to
% initCond  -    1x5 float:  initial parameters for system
% schedule  -    Nx2 float:  matrix of dosetime i and dendritic cell
%                            concentration targets at time i:
%                [ t0, target0;
%                  t1, target1;
%                  ...          ]
% output/cancerConc -   float representing concentration of cancer at
%                       horizon time

%H params a(1)
a0 = 10^(-4); b0 = 0.005; c0 = 10; d0 = 10^(-2); f0 = 1;

%C params a(2)
a1 = 10^(-4); b1 = 0.005; c1 = 10; d1 = 10^(-2); f1 = 1;

%M params a(3)
d2 = 0.02; e2 = 0.1; f2 = 1;

%D params a(4)
e3 = 0.1;

%I params a(5)
a4 = 10^(-2); c4 = 10^(-7); e4 = 10^(-2);

%initial conditions
launchConditions(:) = initCond(:);
%disp(launchConditions);

tspan = linspace(1,horizon,horizon);

f = @(t,a) [
    a0 - b0*a(1) + c0*a(4)*d0*a(1)*(1-a(1)/f0);
    a1 - b1*a(2) + c1*a(5)*(a(3)+a(4))*d1*a(2)*(1-a(2)/f1);
    d2*a(3)*(1-a(3)/f2)-e2*a(3)*a(2);
    -e3*a(4)*a(2);
    a4*a(1)*a(4)-c4*a(2)*a(5)-e4*a(5);
    ];

%numDoses will be used to determine how many times to restart ode45
[numDoses,schedulewidth] = size(schedule);
if schedulewidth ~= 2
    fprintf('ERROR: schedule matrix must be Nx2 \n')
end

runningsolution = zeros(6);

for currentDose = 1:numDoses
    if currentDose==1
        tspan = [0 schedule(currentDose,1)]
        [temptime,tempsolution]=ode45(f,tspan,launchConditions);
        runningsolution=vertcat(runningsolution,[temptime,tempsolution]);
        [endtimeindex,~] = size(runningsolution);
        launchConditions = runningsolution(endtimeindex,2:6);
        launchConditions(4) = launchConditions(4) + schedule(currentDose,2);
        disp(launchConditions)
    else 
        tspan=[schedule(currentDose-1,1) schedule(currentDose,1)]
        [temptime,tempsolution]=ode45(f,tspan,launchConditions);
        runningsolution=vertcat(runningsolution,[temptime,tempsolution]);
        [endtimeindex,~] = size(runningsolution);
        launchConditions = runningsolution(endtimeindex,2:6); 
        launchConditions(4) = launchConditions(4) + schedule(currentDose,2);
        disp(launchConditions)
    end
end

tspan=[schedule(numDoses,1), horizon]
[temptime,tempsolution]=ode45(f,tspan,launchConditions);
runningsolution=vertcat(runningsolution,[temptime,tempsolution]);
 
% disp(runningsolution)

% figure()
% plot(runningsolution(:,1),log(runningsolution(:,4)),'-o',runningsolution(:,1),log(runningsolution(:,5)),'-.')
% title('M(t) and D(t) via ODE45 - Cancer concentration')
% xlabel('Time t (hours)')
% ylabel('Concentration')
% legend('Cancer Concentration','Dendritic Cell Concentration')
% xlim([0,horizon])
% ylim([-30,30])

end %function
