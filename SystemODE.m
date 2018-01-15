function [ dy ] = SystemODE( t, y)
%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This is the ODE function for the model presented in
% "Dynamic Analysis of Cytokine-Induced Differentiation of
% Granulocyte-Monocyte Progenitor Cells"
% Dec. 7, 2017. Bronson R. Weston et al.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
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

dy = zeros(8,1);
C = y(1); P=y(2); G=y(3); I=y(4); E=y(5); GMCSFR = y(6); MCSFR = y(7); GCSFR = y(8);

a=1; b=I-C+1/K; c=-1*C/K;
Cf = (-1*b+sqrt(b^2-4*a*c))/2;
dCdt= (1/(1+exp(-1*sigma*(Wco + Wcc*Cf  + Vgm*S1*GMCSFR/(S1+Kgm) + Vgc*S3*GCSFR/(Kg+S3)))))-C;
dPdt= ((1/(1+exp(-1*sigma*(Wpo + Wpp*P + Wgp*G + Wcp*Cf + Vm*S2*MCSFR/(Km+S2)))))-P);
dGdt= ((1/(1+exp(-1*sigma*(Wgo + Wcg*Cf + Weg*E + Vgg*S3*GCSFR/(Kg+S3)))))-G);
dIdt= ((1/(1+exp(-1*sigma*(Wio + Wpi*P + Vgmi*(S1*GMCSFR)/(S1+Kgm)+ Vgmi*(S3*GCSFR)/(S3+Kg)))))-I);
dEdt= ((1/(1+exp(-1*sigma*(Weo + Wpe*P + Wge*G))))-E);
dGMCSFRdt = (1/10)*((1/(1+exp(-1*sigma*(Wgmro + Wpgmr*P + Wcgmr*Cf))))-GMCSFR);
dMCSFRdt = (1/10)*((1/(1+exp(-1*sigma*(Wmro + Wpmr*P + Wcmr*Cf + Wemr*E + Wgmr*G))))-MCSFR);
dGCSFRdt = (1/10)*((1/(1+exp(-1*sigma*(Wgro + Wpgr*P + Wcgr*Cf + Wggr*G))))-GCSFR);


dy(1) = dCdt;
dy(2) = dPdt;
dy(3) = dGdt;
dy(4) = dIdt;
dy(5) = dEdt;
dy(6) = dGMCSFRdt;
dy(7) = dMCSFRdt;
dy(8) = dGCSFRdt;
end

