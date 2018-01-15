function [ ] = Figure7C()
%This function reproduces the results of Figure 7C from Weston et al. 2018
global S1 S2 S3
sigma = 3.2;  
Wco=-0.79; Wcc=2.71; %effect on C 
Wgo=-0.75; Wcg=1.6; Weg=-1.27; %effect on G
Wpo=-0.80; Wpp=1.79; Wgp =-1.22; Wcp=0.99; %effect on P
Wio=-0.73; Wpi=1.40; %effect on I
Weo=-0.8; Wpe=1.45; Wge=-1.27; %Effect on E
Wgmro= -1.2; Wcgmr= 0.9; Wpgmr= 2.3; % effect on GM-CSFR
Wmro= -1.25; Wpmr= 2; Wcmr= 0.5; Wemr = 1; Wgmr = -1.2; % effect on M-CSFR
Wgro= -1.0; Wpgr= 0.4; Wcgr= 1.1; Wggr = 0.9; % effect on G-CSFR
Kg = 0.9; Kgm = 0.65; Km=0.45; K=7.5;%Binding Constants
Vgm=.75; Vm=0.85; Vgc=0.85; Vgg=0.75; Vgmi=-0.62; %Receptor Rates
sz1=2.5; sz2=1.25; sz3=1;

Ci=0.14242; Pi=0.13313; Gi= 0.088018; Ii=0.1494; Ei=0.091107; GMRi=0.0686; MRi=0.0449; GRi=0.07725;
[t,y] = ode45(@SystemODE, [0 50], [Ci, Pi, Gi, Ii, Ei, GMRi, MRi, GRi]);
Ci=y(end,1); Pi=y(end,2); Gi=y(end,3); Ii=y(end,4); Ei=y(end,5); GMRi=y(end,6);  MRi=y(end,7); GRi=y(end,8);

figure('Position', [25, 100, 760, 450]);
S1 = 0; S2 = 0; S3 = 0;
[t,y] = ode45(@SystemODE, [0 5], [Ci, Pi, Gi, Ii, Ei, GMRi, MRi, GRi ]);
GMSignal=S1*y(:,6)/(S1+Kgm);
plot(t,GMSignal, '--b', 'LineWidth',sz1)
hold on
[t,y] = ode45(@SystemODE, [1 5], [Ci, Pi, Gi, Ii, Ei, GMRi, MRi, GRi ]);

plot(t,GMSignal, '--r', 'LineWidth',sz1)

hold on

plot([5 5], [0, 0.55], 'k','LineWidth',sz3)

S1=0.6; [t,y] = ode45(@SystemODE, [5 60], [Ci, Pi, Gi, Ii, Ei, GMRi, MRi, GRi]);
GMSignal=S1*y(:,6)/(S1+Kgm);
plot(t,GMSignal, '--b', 'LineWidth',sz1)
plot([5 5], [0,GMSignal(1)] , '--b', 'LineWidth',sz1)

S1=1.2; [t,y] = ode45(@SystemODE, [5 60], [Ci, Pi, Gi, Ii, Ei, GMRi, MRi, GRi]);
GMSignal=S1*y(:,6)/(S1+Kgm);
plot(t,GMSignal, '--r', 'LineWidth',sz1)
plot([5 5], [0.015,GMSignal(1)] , '--r', 'LineWidth',sz1)

xlabel('Time')
%plot(t, y(:,3), 'g', t, y(:,4), 'c', t, y(:,5) , '--c')
set(gca, 'fontsize',18)
xlabel('Time', 'fontsize', 18)
ylabel('GM-CSF Signal Strength', 'fontsize', 18)
ylim([0 0.55])
box on
S1 = 0; S2 = 0; S3 = 0; 

end

