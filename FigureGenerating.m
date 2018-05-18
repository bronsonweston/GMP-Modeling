%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script for reproducing figures published in the paper,
% "Mathematical Analysis of Cytokine-Induced Differentiation of 
% Granulocyte-Monocyte Progenitor Cells"
% May 18, 2018. Bronson R. Weston et al.
% 
% This script requires a user to uncomment code corresponding to the
% figure they wish to generate.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
close all; clear all; clc;
global S1 S2 S3
S1=0; S2=0; S3=0;

%%%%%%%%%%%%%%%%%%%%% Figure 3A&B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TimeCourseSimulation(0, 1,0,70)
% TimeCourseSimulation(0, 0,1,70)

%%%%%%%%%%%%%%%%%%%%% Figure 3C&D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BarChart(0, 1,0)
% BarChart(0, 0,1)

%%%%%%%%%%%%%%%%%%%%% Figure 4A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NullclinesProjections(0, 0, 0, 'yes')

%%%%%%%%%%%%%%%%%%%%% Figure 4B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NullclinesProjections(0, 1, 0, 'yes')

%%%%%%%%%%%%%%%%%%%%% Figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bifurcation('M-CSF')

%%%%%%%%%%%%%%%%%%%%% Figure 6A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NullclinesProjections(0, 0, 1, 'yes')

%%%%%%%%%%%%%%%%%%%%% Figure 6B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nullclines3D(0, 0, 1, 'yes')

%%%%%%%%%%%%%%%%%%%%% Figure 6C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bifurcation('G-CSF')

%%%%%%%%%%%%%%%%%%%%% Figure 7A&B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TimeCourseSimulation(0.6, 0,0,60)
% TimeCourseSimulation(1.2, 0, 0,60)

%%%%%%%%%%%%%%%%%%%%% Figure 7C&D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BarChart(0.6, 0,0)
% BarChart(1.2, 0, 0)

%%%%%%%%%%%%%%%%%%%%% Figure 7E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure7E()

%%%%%%%%%%%%%%%%%%%%% Figure 8A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TimeCourseSimulation(0.9, 0, 0,110)

%%%%%%%%%%%%%%%%%%%%% Figure 8B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BarChart(0.9, 0, 0)

%%%%%%%%%%%%%%%%%%%%% Figure 8C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bifurcation('GM-CSF')


%%%%%%%%%%%%%%%%%%%%% Figure 9A-D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NullclinesProjections(0.5, 0, 0, 'yes')
% NullclinesProjections(0.8, 0, 0, 'yes')
% NullclinesProjections(0.9, 0, 0, 'yes')
% NullclinesProjections(1.2, 0, 0, 'yes')

%%%%%%%%%%%%%%%%%%%%% Figure 10A-E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Note: HeatMapping functions take a significant amount of time to
%%%%% compute. For quicker computation times, decrease the first input
%%%%% (the number of simulated cells per cytokine combination) and/or
%%%%% change the axes range/interval (currently [0 0.05 1.2] if changed
%%%%% to [0 0.1 1.2] for both axes would be 4x faster)

% HeatMapping(500, [0 0.05 1.4], 'M-CSF', [0 0.05 1.4], 'G-CSF', 'GM-CSFR', 'Gran+Mono', 0)
% HeatMapping(500, [0 0.05 1.2], 'GM-CSF', [0 0.05 1.2], 'G-CSF', 'GM-CSFR', 'Gran+Mono', 0)
% HeatMapping(500, [0 0.05 1.2], 'GM-CSF & M-CSF', [0 0.05 1.2], 'G-CSF', 'GM-CSFR', 'Gran+Mono', 0)
% HeatMapping(500, [0 0.05 1.2], 'GM-CSF', [0 0.05 1.2], 'M-CSF', 'GM-CSFR', 'Gran+Mono', 0)
% HeatMapping(500, [0 0.05 1.2], 'G-CSF & GM-CSF', [0 0.05 1.2], 'M-CSF', 'GM-CSFR', 'Gran+Mono', 0)




