function [ dCdt, dPdt ] = dcfdp( Cf, P,Int )
%This function takes the [C/EBP]f and [PU.1] concentrations and calculates
%the rate of change of those variabels, returned as dCdt and dPdt
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
    C= ((K*(I+Cf)+1)/(K+1/Cf));
    dC= (1/(1+exp(-1*sigma*(Wco + Wcc*Cf  + Vgm*S1*GMCSFR/(S1+Kgm) + Vgc*S3*GCSFR/(Kg+S3)))))-C;
    dP=(1/(1+exp(-1*sigma*(Wpo + Wpp*P + Wgp*G + Wcp*Cf + Vm*S2*MCSFR/(Km+S2)))))-P;
    if isempty(dCdt) 
        dCdt= dC;
        dPdt = dP;
    elseif norm([dC,dP])<norm([dCdt,dPdt])
        dCdt=dC;
        dPdt=dP;
    end
    end
end

