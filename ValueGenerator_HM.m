function [GMHcount, DefaultConc] = ValueGenerator_HM(NumberOfCells, StandDev, c1, c1value, c2, c2value,constantThirdStimulant)
%[GMHcount] = Stoch_Generator(NumberOfCells, StandDev, GMCSF, MCSF, GCSF, s4)
global S1 S2 S3
K=7.5;
S1=0; S2=0; S3=0;
GMCSF=constantThirdStimulant; MCSF=constantThirdStimulant; GCSF=constantThirdStimulant; 


GMHcount = zeros(1,4);
if strcmp(c1,c2)
    error('Only 1 cytokine type specified')
end

if strcmp(c1, 'GM-CSF')
    GMCSF=c1value;
elseif strcmp(c1, 'M-CSF')
    MCSF = c1value;
elseif strcmp(c1, 'G-CSF')
    GCSF= c1value;
elseif strcmp(c1, 'GM-CSF & M-CSF') || strcmp(c1, 'M-CSF & GM-CSF')
    GMCSF= c1value;
    MCSF= c1value;
elseif strcmp(c1, 'G-CSF & M-CSF') || strcmp(c1, 'M-CSF & G-CSF')
    GCSF= c1value;
    MCSF= c1value;
elseif strcmp(c1, 'GM-CSF & G-CSF') || strcmp(c1, 'G-CSF & GM-CSF')
    GCSF= c1value;
    GMCSF= c1value;
elseif strcmp(c1, 'rateconstant')
    Rateconstant=c1value;
else
    error('Cytokine 1 invalid')
end
   
if strcmp(c2, 'GM-CSF')
    GMCSF=c2value;
elseif strcmp(c2, 'M-CSF')
    MCSF = c2value;
elseif strcmp(c2, 'G-CSF')
    GCSF= c2value;
elseif strcmp(c2, 'rateconstant')
    Rateconstant=c2value;
elseif ~strcmp(c2, 'None')
    error('Cytokine 2 invalid')
end

% Simulates other cytokines in the agar solution
if GMCSF == 0
    GMCSF=0;
end

Ci=0.14242; Pi=0.13313; Gi= 0.088018; Ii=0.1494; Ei=0.091107; GMRi=0.0686; MRi=0.0449; GRi=0.07725;
[~,y] = ode45(@(t,y) SystemODEfast(t,y,S1,S2,S3), [0 50], [Ci, Pi, Gi, Ii, Ei, GMRi, MRi, GRi]);
Ci=y(end,1); Pi=y(end,2); Gi=y(end,3); Ii=y(end,4); Ei=y(end,5); GMRi=y(end,6);  MRi=y(end,7); GRi=y(end,8);
S1 = GMCSF; S2 = MCSF; S3 = GCSF;
[~,y] = ode45(@(t,y) SystemODEfast(t,y,S1,S2,S3), [0 150], [Ci, Pi, Gi, Ii, Ei, GMRi, MRi, GRi]);
DefaultConc = y(end,:);
S1=0; S2=0; S3=0;

for i=1:NumberOfCells
    S1 = 0; S2 = 0; S3 = 0;
    tCi = normrnd(Ci, StandDev*Ci);
    tPi = normrnd(Pi, StandDev*Pi);
    tGi = normrnd(Gi, StandDev*Gi);
    tIi = normrnd(Ii, StandDev*Ii);
    tEi = normrnd(Ei, StandDev*Ei);
    tGMRi = normrnd(GMRi, StandDev*GMRi);
    tMRi = normrnd(MRi, StandDev*MRi);
    tGRi = normrnd(GRi, StandDev*GRi);
    S1 = GMCSF; S2 = MCSF; S3 = GCSF;
    [~,y2] = ode45(@(t,y) SystemODEfast(t,y,S1,S2,S3), [0 150], [tCi, tPi, tGi, tIi, tEi, tGMRi, tMRi, tGRi]);
    S1 = 0; S2 = 0; S3 = 0; 
    b=y2(:,4)-y2(:,1)+(1/K);
    c=-1*y2(:,1)/K;
    CF=(-1*b+(((b.*b)-4*c).^(1/2)))/2;
    if CF(end,1)/y2(end,2) > 3
        GMHcount(1,1) = GMHcount(1,1) + 1;
    elseif CF(end,1)/y2(end,2) < 1/3
        GMHcount(1,2) = GMHcount(1,2) + 1;
    elseif y2(end,1) > 0.4
            GMHcount(1,3) = GMHcount(1,3) + 1; 
    else
        GMHcount(1,4) = GMHcount(1,4) + 1; %undifferentiated cells
    end
end
end

