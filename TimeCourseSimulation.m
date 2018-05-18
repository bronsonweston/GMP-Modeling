function [ output_args ] = TimeCourseSimulation(GMCSF, MCSF,GCSF,timeofsim)
%Returns a simulation of protein activity over the specified simulation
%time, timeofsim, and given the specified cytokine concentrations. Requires
%SystemODE.m
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
leg=legend('C/EBP_T', 'C/EBP_F', 'Gfi-1', 'G-CSFR', 'PU.1', 'Egr-2', 'IRF8', 'M-CSFR', 'GM-CSFR', 'Location', 'eastoutside');%, 'Orientation','horizontal');
box on

end

