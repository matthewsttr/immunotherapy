clc;
clear all;
close all;



horizon=25000

%create a dosing schedule for x fraction of the horizon, with n equally
%spaced doses with specified level

n=200;
x=0.8;
initCond = [0.1 ,0.1 ,1 ,9.5 ,0.01]
level=0.5;
scheduletimes = linspace(round(horizon/n),round(horizon*x),n)
scheduledoses=ones(1,n)*level
schedule2 = horzcat(transpose(scheduletimes), transpose(scheduledoses))
%%%%


therapy_cancer(horizon,initCond)
therapy_dosed(horizon, initCond,schedule2)