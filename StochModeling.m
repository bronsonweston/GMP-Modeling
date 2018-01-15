function [] = StochModeling(NumberOfCells, timeofsim, GMCSF, MCSF, GCSF)
%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Function for stochastic simulations of the GMP differentiation model 
% presented in "Dynamic Analysis of Cytokine-Induced Differentiation of
% Granulocyte-Monocyte Progenitor Cells"
% Dec. 7, 2017. Bronson R. Weston et al.
% This function requires 'SystemODE.m' to run.
% This function generates a population of cells with stochastically
% generated initial conditions and simulates their differentiation over the
% specified time. The function prints out the final population composition
% in the command window.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
global S1 S2 S3
S1 = 0; S2 = 0; S3 = 0; K=7.5;
StandDev=0.2;
Ci=0.14242; Pi=0.13313; Gi= 0.088018; Ii=0.1494; Ei=0.091107; GMRi=0.0686; MRi=0.0449; GRi=0.07725;
[~,y] = ode45(@SystemODE, [0 50], [Ci, Pi, Gi, Ii, Ei, GMRi, MRi, GRi]);
Ci=y(end,1); Pi=y(end,2); Gi=y(end,3); Ii=y(end,4); Ei=y(end,5); GMRi=y(end,6);  MRi=y(end,7); GRi=y(end,8);

Data = zeros(NumberOfCells, 9);
GMHcount = zeros(1, 3);
for i=1:length(Data(:,1))
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
    [~,y2] = ode45(@SystemODE, [0 timeofsim], [tCi, tPi, tGi, tIi, tEi, tGMRi, tMRi, tGRi]);
    b=y2(:,4)-y2(:,1)+(1/K);
    c=-1*y2(:,1)/K;
    CF=(-1*b+(((b.*b)-4*c).^(1/2)))/2;
    if CF(end,1)/y2(end,2) > 3
        GMHcount(1,1) = GMHcount(1,1) + 1;
    elseif CF(end,1)/y2(end,2) < 1/3
        GMHcount(1,2) = GMHcount(1,2) + 1;
    else
        if y2(end,1) > 0.4 && y2(end,2) > 0.4
            GMHcount(1,3) = GMHcount(1,3) + 1; %hybrid cells
        end
    end
end
disp('Simulation Report:') 
disp(['     Cytokine Strength: GM-CSF=', num2str(S1), '  M-CSF=', num2str(S2), '  G-CSF=', num2str(S3)]) 
disp(['     GMPs:  ', num2str(100*(NumberOfCells-(GMHcount(1,1)+GMHcount(1,2)+GMHcount(1,3)))/NumberOfCells), '%'])
disp(['     Granulocyte Progenitors:  ', num2str(100*GMHcount(1,1)/NumberOfCells), '%'])
disp(['     Monocytes:  ', num2str(100*GMHcount(1,2)/NumberOfCells), '%'])
disp(['     M-MDSCs:  ', num2str(100*GMHcount(1,3)/NumberOfCells), '%'])

S1 = 0; S2 = 0; S3 = 0;

end

