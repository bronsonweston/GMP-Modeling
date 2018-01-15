function [ oneOrZero, e ] = Stable(parm0)
%This function evaluates the stability of the steady state solution. It
%returns 1 if it is stable, and 0 if it is unstable. 
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
Int=0.01;
Cf=parm0(1);
P=parm0(2);
dCdt=[];
dPdt=[];
Gall= FindGCytokines(P, Cf, 0, 1.0, Int, S1, S2, S3);
for i=1:length(Gall)
    G=Gall(i);
    GMCSFR= (1/(1+exp(-1*sigma*(Wgmro + Wpgmr*P + Wcgmr*Cf))));
    GCSFR= (1/(1+exp(-1*sigma*(Wgro + Wpgr*P + Wcgr*Cf + Wggr*G))));
    I= (1/(1+exp(-1*sigma*(Wio + Wpi*P + Vgmi*(S1*GMCSFR)/(S1+Kgm)+ Vgmi*(S3*GCSFR)/(S3+Kg)))));
    E= 1/(1+exp(-1*sigma*(Weo + Wpe*P + Wge*G)));
    MCSFR= 1/(1+exp(-1*sigma*(Wmro + Wpmr*P + Wcmr*Cf + Wemr*E + Wgmr*G)));
    %         Ci= ((K*(I+Cf)+1)/(K+1/Cf));
    Ci= ((I+1/K)^2-(2*Cf+I+1/K)^2)/(-4*(Cf+1/K));
    dC= (1/(1+exp(-1*sigma*(Wco + Wcc*Cf  + Vgm*S1*GMCSFR/(S1+Kgm) + Vgc*S3*GCSFR/(Kg+S3)))))-Ci;
    dP=(1/(1+exp(-1*sigma*(Wpo + Wpp*P + Wgp*G + Wcp*Cf + Vm*S2*MCSFR/(Km+S2)))))-P;
    if isempty(dCdt)
        dCdt= dC;
        dPdt = dP;
        C=Ci;
        Gi=Gall(i);
        GMCSFRi= (1/(1+exp(-1*sigma*(Wgmro + Wpgmr*P + Wcgmr*Cf))));
        GCSFRi= (1/(1+exp(-1*sigma*(Wgro + Wpgr*P + Wcgr*Cf + Wggr*G))));
        Ii= (1/(1+exp(-1*sigma*(Wio + Wpi*P + Vgmi*(S1*GMCSFR)/(S1+Kgm)+ Vgmi*(S3*GCSFR)/(S3+Kg)))));
        Ei= 1/(1+exp(-1*sigma*(Weo + Wpe*P + Wge*G)));
        MCSFRi= 1/(1+exp(-1*sigma*(Wmro + Wpmr*P + Wcmr*Cf + Wemr*E + Wgmr*G)));
    elseif norm([dC,dP])<norm([dCdt,dPdt]) %&& Ci <= 1
        dCdt=dC;
        dPdt=dP;
        C=Ci;
        Gi=Gall(i);
        GMCSFRi= (1/(1+exp(-1*sigma*(Wgmro + Wpgmr*P + Wcgmr*Cf))));
        GCSFRi= (1/(1+exp(-1*sigma*(Wgro + Wpgr*P + Wcgr*Cf + Wggr*G))));
        Ii= (1/(1+exp(-1*sigma*(Wio + Wpi*P + Vgmi*(S1*GMCSFR)/(S1+Kgm)+ Vgmi*(S3*GCSFR)/(S3+Kg)))));
        Ei= 1/(1+exp(-1*sigma*(Weo + Wpe*P + Wge*G)));
        MCSFRi= 1/(1+exp(-1*sigma*(Wmro + Wpmr*P + Wcmr*Cf + Wemr*E + Wgmr*G)));
    end
end
InitialConditions= [C,P,Gi,Ii,Ei,GMCSFRi,MCSFRi,GCSFRi];
j=1;
stepsize=0.000001;
for i=1:length(InitialConditions(1,:))
    InitialConditions(j,i) = InitialConditions(j,i) + stepsize;
    J(:,i)= SystemODE(0,InitialConditions);
    InitialConditions(j,i) = InitialConditions(j,i) - 2*stepsize;
    J(:,i)= (J(:,i)+(-1*SystemODE(0,InitialConditions)))./2;
    InitialConditions(j,i) = InitialConditions(j,i) + stepsize; %returns step to normal
end
J;
e=eig(J);
oneOrZero=all(real(e)<0);
end
