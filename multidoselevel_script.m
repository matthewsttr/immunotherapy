clc;
clear all;
close all;



horizon=4000;

%create a dosing schedule for x fraction of the horizon, with n equally
%spaced doses with specified level
n=2;

initCond = [0 ,0 ,1 ,0.1 ,0];
level=0.5;
dosecut = linspace(0,horizon,n+2);
scheduletimes= dosecut(2:(n+2)-1);

scheduledoses=ones(1,n)*level;
schedule2 = horzcat(transpose(scheduletimes), transpose(scheduledoses));
%%%%

max_numDoses=100
for currentlevelindex=1:5
    for n=1:max_numDoses
        initCond = [0 ,0 ,1 ,0.1 ,0];
        level=0.2*currentlevelindex;
        dosecut = linspace(0,horizon,n+2);
        scheduletimes= dosecut(2:(n+2)-1);
        scheduledoses=ones(1,n)*level;
        %schedules{i} = horzcat(transpose(scheduletimes), transpose(scheduledoses));
        scheduledoses = horzcat(transpose(scheduletimes), transpose(scheduledoses));
        nthsolution=therapy_dosed(horizon,initCond,scheduledoses);
        [nthcost(n,1,currentlevelindex),nthcost(n,2,currentlevelindex),nthcost(n,3,currentlevelindex)]=cost_to_go(nthsolution);
    end
end

nthcost
plotnvector = transpose(linspace(1,max_numDoses,max_numDoses));
figure()
hold on
for currentlevelindex=1:5
    plot(plotnvector,nthcost(:,1,currentlevelindex),'-')
    x=size(plotnvector,1)
    y=nthcost(x,1,currentlevelindex)
    labeltext=strcat("level=",num2str(currentlevelindex*0.2))
    text(x,y,labeltext);
    %plot(plotnvector,nthcost(:,1,currentlevelindex),'-','Color','#77AC30')
    %plot(plotnvector,nthcost(:,2,currentlevelindex),'-','Color','blue')
    %plot(plotnvector,nthcost(:,3,currentlevelindex),'-','Color','red')
end
hold off

title('Quadratic Total Cost')
xlabel('total number of doses')
ylabel('cost index')
%legend('total cost', 'cancer cost', 'combatant cost')

% therapy_cancer(horizon,initCond);
% schedule2sol = therapy_dosed(horizon, initCond,schedule2);
% 
% cost_to_go(schedule2sol)