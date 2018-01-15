close all; clc; clear all;
%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script for numerical simulations of the GMP differentiation model 
% presented in "Dynamic Analysis of Cytokine-Induced Differentiation of
% Granulocyte-Monocyte Progenitor Cells"
% Dec. 7, 2017. Bronson R. Weston et al.
% This script requires 'StochModeling.m' and 'SystemODE.m' to run.
% Results of script include the following:
% 1) A plot of protein concentrations over the time of a numerical 
%    simulation of a default (average) cell under the specified conditions. 
%    Result is displayed as a figure.
% 2) A simulation report of the population compositions for a 
%    stochastically generated GMP population of size NumberOfCells. Result
%    is displayed in the command window. 
%  Simulation conditions should be specified in the area below.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%



%% DEFINE SIMULATION CONDITIONS HERE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
GMCSF=0.6; GCSF=0.6; MCSF=0.6;
timeofsim = 70;
NumberOfCells= 100;

%The following are conditions specified in figures from the paper
%Figure 3A: GMCSF=0; GCSF=0; MCSF=1; timeofsim = 70;
%Figure 3B: GMCSF=0; GCSF=1; MCSF=0; timeofsim = 70;
%Figure 7A: GMCSF=0.6; GCSF=0; MCSF=0; timeofsim = 60;
%Figure 7B: GMCSF=1.2; GCSF=0; MCSF=0; timeofsim = 60;
%Figure 8A: GMCSF=0.9; GCSF=0; MCSF=0; timeofsim = 110;

%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%

global S1 S2 S3;
S1=0; S2=0; S3=0;
K=7.5; sz1=2.5; sz2=1.25; sz3=1;

%simulation from initial conditions to ensure settlement to the GMP state
Ci=0.14242; Pi=0.13313; Gi= 0.088018; Ii=0.1494; Ei=0.091107; GMRi=0.0686; MRi=0.0449; GRi=0.07725;
[t,y] = ode45(@SystemODE, [0 20], [Ci, Pi, Gi, Ii, Ei, GMRi, MRi, GRi]);
Ci=y(end,1); Pi=y(end,2); Gi=y(end,3); Ii=y(end,4); Ei=y(end,5); GMRi=y(end,6);  MRi=y(end,7); GRi=y(end,8);

figure('Position', [150, 250, 950, 450]);
[t,y] = ode45(@SystemODE, [0 5], [Ci, Pi, Gi, Ii, Ei, GMRi, MRi, GRi ]);
hold on
plot(t, y(:,1), 'r', 'LineWidth',sz1);  
b=y(:,4)-y(:,1)+(1/K);
c=-1*y(:,1)/K;
CF=(-1*b+(((b.*b)-4*c).^(1/2)))/2;
plot(t,CF, 'k', 'LineWidth',sz1)
plot( t, y(:,3), 'r', 'LineWidth', sz2)
plot(t, y(:,8),'--r', 'LineWidth', sz2)
plot(t, y(:,2), 'b', 'LineWidth',sz1)
plot(t, y(:,5) , 'b', 'LineWidth', sz2) 
plot(t, y(:,4), 'c', 'LineWidth', sz2) 
plot(t, y(:,7), '--b', 'LineWidth', sz2)
plot(t, y(:,6) , '--m', 'LineWidth', sz2)
plot([5 5], [0, 1.05], 'k','LineWidth',sz3)

%Simulation of specified conditions for default cell
S1=GMCSF; S2=MCSF; S3=GCSF;
[t,y] = ode45(@SystemODE, [5 timeofsim+5], [y(end,1), y(end,2), y(end,3), y(end,4), y(end,5), y(end,6), y(end,7), y(end,8)]);
plot(t, y(:,1), 'r', 'LineWidth',sz1);
b=y(:,4)-y(:,1)+(1/K);
c=-1*y(:,1)/K;
CF=(-1*b+(((b.*b)-4*c).^(1/2)))/2;
plot(t,CF, 'k', 'LineWidth',sz1)
plot( t, y(:,3), 'r', 'LineWidth', sz2)
plot(t, y(:,8),'--r', 'LineWidth', sz2)
plot(t, y(:,2), 'b', 'LineWidth',sz1)
plot(t, y(:,5) , 'b', 'LineWidth', sz2) 
plot(t, y(:,4), 'c', 'LineWidth', sz2) 
plot(t, y(:,7), '--b', 'LineWidth', sz2)
plot(t, y(:,6) , '--m', 'LineWidth', sz2)
set(gca, 'fontsize',18)
xlabel('Time', 'fontsize', 18)
ylabel('Protein Activity', 'fontsize', 18)
ylim([0 1.05])
xlim([0 timeofsim+5])
title({['GM-CSF = ', num2str(GMCSF)], ['M-CSF = ', num2str(MCSF)], ['G-CSF = ', num2str(GCSF)]}, 'fontsize', 16)
leg=legend('[C/EBP\alpha]_T', '[C/EBP\alpha]_F', 'Gfi-1', 'G-CSFR', 'PU.1', 'Egr-2', 'IRF8/AP1', 'M-CSFR', 'GM-CSFR', 'Location', 'eastoutside');%, 'Orientation','horizontal');
box on

%Stochastic Population Simulation
StochModeling(NumberOfCells, timeofsim, GMCSF, MCSF, GCSF)