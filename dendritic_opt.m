% % step 0
% T = 1000;       % horizon time (hours)
% N = 3;          % number of injection
% V = 0.2;        % Vaccine quantity
% x0 = [0 ,0 ,1 ,0.1 ,0];    % initial cell population
% s0 = [300,600,900];   % initial schedule of doses (hours)
% 
% % step 1
% % solve the original system with x0
% therapy_dosed



%%%%%%%%%%%%%%%%%
%create a dosing schedule for x fraction of the horizon, with n equally
%spaced doses with specified level V
clear all;
close all;
clc;
format long;

horizon=4000;
N=80;
doselevel=0.5;
initCond = [0 ,0 ,1 ,0.1 ,0];
tolerance=10^-8;

dosecut = linspace(0,horizon,N+2);
scheduletimes= dosecut(2:(N+2)-1);

% dosecut = horizon*transpose(rand(N,1));
% dosecut = sort(dosecut);
% scheduletimes = dosecut;

%scheduletimes = [1000,1100,1500,1700]

scheduledoses=ones(1,N)*doselevel;
schedule = horzcat(transpose(scheduletimes), transpose(scheduledoses));

stepsize = 1;

%solve the system and calculate the cost
currentschedule = schedule;
currentsolution = therapy_dosed(horizon,initCond,currentschedule);
[currentcost,~,~] = cost_to_go(currentsolution);
running_optimization_data = [currentcost,transpose(currentschedule(:,1));0,transpose(currentschedule(:,2))];
dendriticplotter(currentsolution,'original schedule')
cumulativetimes = transpose(currentschedule(:,1));

costdiff(1)=currentcost;
optstep=1;
oldschedule=10^29*currentschedule;
oldoldschedule=10^30*currentschedule;
oldcost=10^29;
oldoldcost=10^30;

while costdiff(optstep) > tolerance
optstep=optstep+1;
    for i=1:N
        steppedschedule = currentschedule;
        steppedschedule(i,1) = steppedschedule(i,1) + stepsize;
        steppedsolution = therapy_dosed(horizon,initCond,steppedschedule);
        [steppedcost(i),b,c] = cost_to_go(steppedsolution);
        steppedcost(i) = steppedcost(i)-currentcost;
        [steppedcost(i),b,c];
        %dendriticplotter(steppedsolution)
    end
    
    gradient=transpose(steppedcost./norm(steppedcost));
    if N>1 
        oldoldschedule=oldschedule; 
    end
    oldschedule=currentschedule;
    

    % the gradient step
    currentschedule(:,1) = currentschedule(:,1) - stepsize*gradient 
    % if the schedule keeps going in circles, just finish
    if norm(currentschedule(:,1)-oldoldschedule(:,1)) < stepsize/N
        fprintf('stuck in a rut')
        break; 
    end

    % if the gradient step puts the current schedule out of bounds, knock
    % it back in bounds to the closest point. This is an example of the p
    % projective gradient method. this is simple since the feasible region
    % is a square box
    for i=1:N 
        if currentschedule(i,1) < 0 
            currentschedule(i,1) = i*stepsize; %sketchy to include the i, but it keeps two zero-plunging doses from piling on each-other
        end 
    end
        

    cumulativetimes = vertcat(cumulativetimes,transpose(currentschedule(:,1)));
    currentsolution = therapy_dosed(horizon,initCond,currentschedule); %solve the gradient-stepped problem
    oldoldcost=oldcost;
    oldcost=currentcost;
    [currentcost,~,~] = cost_to_go(currentsolution); % get the cost after the gradient step
    costdiff(optstep) = oldcost-currentcost;
    costdoublediff(optstep) = oldoldcost-currentcost;
    running_optimization_data=vertcat(running_optimization_data,[currentcost,transpose(currentschedule(:,1));0,transpose(currentschedule(:,2))]);
    if any(currentschedule(:,1)<0) 
           break; 
    end
    if any(currentschedule(:,1)>horizon) 
           break;
    end
end    

running_optimization_data(:,1);
running_optimization_data;
plotdata = running_optimization_data(1:2:end,:);
figure()
plot(plotdata(:,1))
title('Cost Function vs Optimization Steps')
dendriticplotter(currentsolution,'optimized schedule')
%opt_data = running_optimization_data(1:2:end,:)  % odd matrix
%B = running_optimization_data(2:2:end,:)  % even matrix
figure()
plot(1:size(cumulativetimes),cumulativetimes)
ylim([0 horizon])
%titlestr = ['Dose Times vs Optimization Steps' newline 'initial schedule: ' newline mat2str(round(dosecut))];
titlestr = 'Dose Times vs Optimization Steps';
title(titlestr)
xlabel("optimization steps")
ylabel("scheduled dosetimes")
format default


%dim = [-0.5 -0.5 0.3 0.3];
%annotation('textbox',dim,'String',str,'FitBoxToText','on');


optstep
currentcost
oldcost
dosecut