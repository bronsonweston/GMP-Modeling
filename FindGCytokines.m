function [GSolution] = FindGCytokines( P, Cf, Ymin, Ymax, Int, S1, S2, S3, tracker )
%This function finds all [Gfi-1] solutions within the specified interval 
%given the cytokine conditions and a specified [PU.1] and {C/EBP]f 
%concentration
global Iteration IterationB IterationC
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
if ~exist('tracker','var')
    tracker=1;
end

Tolerance = 1e-9;
if isempty(Iteration)
    Iteration= 0; %if the solver cannot find the solution it will iterate over a set number of times. Starts counter if not already started
end
GSolution=[];
tempSol=[];

%Functions for solver to work on
E= @(G) 1/(1+exp(-1*sigma*(Weo + Wpe*P + Wge*G)));
GCSFR= @(G) 1/(1+exp(-1*sigma*(Wgro + Wpgr*P + Wcgr*Cf + Wggr*G)));
Function = @(G) (1/(1+exp(-1*sigma*(Wgo + Wcg*Cf + Weg*E(G) + Vgg*S3*GCSFR(G)/(Kg+S3)))))-G;

%declare initial 
stoppoint = 1;
interval = Int;
cmin = Ymin;
cmax = Ymax;
cValue = 1;
%Find inital estimates of solutions for Function
PotentialSS = zeros(1, 4); %first column x value, second colum y value, third is slope, and fourth is direction
counter=0;
for i = 1:1:((cmax-cmin)/interval+1)
    value = Function((i-1)*interval+cmin);
    if isreal(value) 
        counter=counter+1;
        PotentialSS(counter,1) = (i-1)*interval+cmin;
        PotentialSS(counter,2) = value;
    end
end

%Calculate slope and direction
%First calculate the inital and end slopes 
try 
    PotentialSS(1,3) = (PotentialSS(2,2)-PotentialSS(1,2));
    PotentialSS(1,4) = PotentialSS(1,3)*PotentialSS(1,2);
    PotentialSS(end,3) = (PotentialSS(end,2)-PotentialSS(end-1,2));
    PotentialSS(end,4) = PotentialSS(end,3)*PotentialSS(end,2);
catch %if there is only one PotentialSS
    [Gsol]=FindGCytokines(P, Cf,PotentialSS(1,1)-Int,PotentialSS(1,1)+Int,Int/10, S1, S2, S3, tracker+1);
    GSolution = [GSolution; Gsol];
    return
end

for i=2:length(PotentialSS(:,1))-1
    PotentialSS(i,3) = (PotentialSS(i+1,2)-PotentialSS(i-1,2));
    PotentialSS(i,4) = PotentialSS(i,3)*PotentialSS(i,2);
end

cPointer = zeros(length(PotentialSS(:,1))-1,1);
for i= length(PotentialSS(:,1))-1:-1:1
    P1 = PotentialSS(i,2);
    P2 = PotentialSS(i+1,2);

    if P1 == 0 
        cPointer(i) = PotentialSS(i,1);
    elseif (P1 > 0 && P2 < 0) || (P1 < 0 && P2 > 0)
        try
            if PotentialSS(i,4)<0
                %A=PotentialSS(i-2:i+2,:)
                cPointer(i) = PotentialSS(i,1);
            elseif PotentialSS(i,4)>0 && PotentialSS(i-1,4)<0
                cPointer(i) = PotentialSS(i,1);
            end
        catch
            disp('Local minima or maxima, no zero found');
        end
    elseif PotentialSS(i,4)<0 && PotentialSS(i+1,4)>0  %minimum, might contain zero
        try
            intersect = (PotentialSS(i+1,2)-PotentialSS(i-1,2)-PotentialSS(i+1,4)*2)/(PotentialSS(i-1,4)-PotentialSS(i+1,4));
            estInt= (PotentialSS(i-1,4)*intersect+PotentialSS(i-1,2))*PotentialSS(i-1,2);
            if  estInt <= 0%find intersection of tangent lines and multiply it by y to see if it should cross zero
                [Gsol]=FindGCytokines(P, Cf,PotentialSS(i-1,1),PotentialSS(i+1,1),Int/10, S1, S2, S3, tracker+1);
                GSolution = [GSolution; Gsol];
            end
        catch
            intersect = (PotentialSS(i+1,2)-PotentialSS(i,2)-PotentialSS(i+1,4)*2)/(PotentialSS(i,4)-PotentialSS(i+1,4));
            estInt= (PotentialSS(i,4)*intersect+PotentialSS(i,2))*PotentialSS(i,2);
            if estInt <= 0
                [Gsol]=FindGCytokines(P, Cf,PotentialSS(i,1),PotentialSS(i+1,1),Int/10, S1, S2, S3, tracker+1);
                GSolution = [GSolution; Gsol];
            end
        end
    else
        cPointer(i) = [];
    end
end
 
if tracker==1
    IterationB=0;
end

% Find solutions for Function
for x=1:length(cPointer)
    cmin = cPointer(x);
    cmax = cPointer(x)+Int;
    stoppoint = 1;
    minValue = 1000;
    interval = Int/10;
    while stoppoint == 1
        c = cmin:interval:cmax;
        for i = 1:length(c)
            value= Function(c(i));
            if isreal(value)
                pointer = value;
                if abs(pointer) < abs(minValue)
                    minValue = pointer;
                    cValue = c(i);
                end
            end
        end
        cmin = cValue - interval;
        cmax = cValue + interval;
        interval = interval/10;
        if interval < 10^-17
            break
        end
        if abs(minValue) < Tolerance 
            stoppoint = 0;
            tempSol = [tempSol;cValue];
        end
    end
end
%Remove all repeating Solutions
tempSol=[tempSol;GSolution];
if ~isempty(tempSol)
    tempSol2 = tempSol(1);
else
    tempSol2 = [];
end
for i=1:length(tempSol)
    add=1;
    for x=1:length(tempSol2)
        if abs(tempSol(i)-tempSol2(x))<=Tolerance
            add=0;
        end
    end
    if add==1
        tempSol2=[tempSol2; tempSol(i)];
    end
end
 % finished removing all repeating values

% solution is currently in G, convert to PU.1
GSolution=tempSol2;

GSolution = sort(GSolution);
GSolution1= zeros(length(GSolution),1);

if isempty(GSolution) && Iteration < 10 && abs(PotentialSS(end,3))-abs(PotentialSS(end-1,3))<0
    Iteration = Iteration + 1;
    IterationC= 1;
    PotentialSS(end,1);
    PotentialSS(end-1,1);
    [GSolution] = FindGCytokines(P,Cf,PotentialSS(end-1,1), PotentialSS(end,1), Int/10, S1, S2, S3, tracker+1 );
    return
end

for i=length(GSolution):-1:1
    if GSolution(i) < 0 
        GSolution(i)=[];
    end
end
% GSolution
if tracker==1
    Iteration = 0;
%     G= GSolution
    IterationC=0;
end

end
