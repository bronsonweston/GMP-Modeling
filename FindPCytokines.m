function [CSolution, GSolution] = FindPCytokines(Cf, Ymin, Ymax, Int, S1, S2, S3, tracker, tripointer )
%This function finds all [PU.1] solutions within the specified interval 
%given the cytokine conditions and a specified [C/EBP]f concentration
global Iteration
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

if ~exist('tripointer','var')
    tripointer=0;
end
Tolerance = 1e-7;
if isempty(Iteration)
    Iteration= 0; %if the solver cannot find the solution it will iterate over a set number of times. Starts counter if not already started
end
CSolution=[];
GSolution=[];
tempSol=[];
gtempSol=[];

%Functions for solver to work on
G= @(P) FindGCytokines(P, Cf, 0, 1.0, 0.01, S1, S2, S3);
E= @(P) 1/(1+exp(-1*sigma*(Weo + Wpe*P + Wge*G(P))));
MCSFR= @(P) 1/(1+exp(-1*sigma*(Wmro + Wpmr*P + Wcmr*Cf + Wemr*E(P) + Wgmr*G(P))));
Function = @(P) (1/(1+exp(-1*sigma*(Wpo + Wpp*P + Wgp*G(P) + Wcp*Cf + Vm*S2*MCSFR(P)/(Km+S2)))))-P;

%declare initial
stoppoint = 1;
interval = Int;
cmin = Ymin;
cmax = Ymax;
cValue = 1;
% PotentialSS1 stores all potential solutions
%first column = x value, second colum = y value, third is slope, and fourth is direction
PotentialSS1 = zeros(1, 4); %first column x value, second colum y value, third is slope, and fourth is direction
%In case function G produces more than 1 solution, PotentialSS2 and
%Potential SS3 are used to store the other solutions
PotentialSS2=[];
PotentialSS3=[];
counter=0;
alternatecounter=0;
PnTracker=1;
Extra=[];
if tripointer==0 || tripointer==1 %If in iteration, tripointer is used to determine which G solution is of interest (if there are multiple G solutions)
    %Initial scan for solutions
    for i = 1:1:((cmax-cmin)/interval+1)
        try
            value = Function((i-1)*interval+cmin);
            if PnTracker==1
                if isreal(value)
                    counter=counter+1;
                    PotentialSS1(counter,1) = (i-1)*interval+cmin;
                    PotentialSS1(counter,2) = value;
                else
                end
            else
                if tripointer~=0
                    break
                end
                alternatecounter= alternatecounter+1;
                PotentialSS2(alternatecounter,1) = (i-1)*interval+cmin;
                PotentialSS2(alternatecounter,2) = value;
            end
            
        catch % Catch error due to multiple G solutions 
            Pn=(i-1)*interval+cmin;
            Gn= FindGCytokines(Pn, Cf, 0, 1.0, 0.01, S1, S2, S3);
            En= (1+exp(-1*sigma*(Weo + Wpe*Pn + Wge*Gn))).^-1;
            MCSFRn= (1+exp(-1*sigma*(Wmro + Wpmr*Pn + Wcmr*Cf + Wemr*En + Wgmr*Gn))).^-1;
            Valuen = ((1+exp(-1*sigma*(Wpo + Wpp*Pn + Wgp*Gn + Wcp*Cf + Vm*S2*MCSFRn/(Km+S2)))).^-1)-Pn;
            try
                target=Valuen(1);
                targetslope=(Valuen(1)-PotentialSS1(counter,2));
                slope=PotentialSS1(counter,2)-PotentialSS1(counter-1,2);
                pointer=1;
                for j=2:length(Valuen)
                    if abs(slope-(Valuen(j)-PotentialSS1(counter,2))) < abs(slope-targetslope)
                        target=Valuen(j);
                        targetslope=(Valuen(j)-PotentialSS1(counter,2));
                        pointer=j;
                    end
                end
                Valuen(pointer)=[];
                value=target;
                counter=counter+1;
                alternatecounter= alternatecounter+1;
                PotentialSS1(counter,1) = Pn;
                PotentialSS1(counter,2) = value;
                PotentialSS2(alternatecounter,1)=Pn;
                PotentialSS3(alternatecounter,1)=Pn;
                PotentialSS2(alternatecounter,2)=Valuen(1);
                PotentialSS3(alternatecounter,2)=Valuen(2);
                PnTracker=3;
            catch
                counter=counter+1;
                alternatecounter= alternatecounter+1;
                PotentialSS1(counter,1) = Pn;
                PotentialSS1(counter,2) = Valuen(3);
                PotentialSS2(alternatecounter,1)=Pn;
                PotentialSS3(alternatecounter,1)=Pn;
                PotentialSS2(alternatecounter,2)=Valuen(1);
                PotentialSS3(alternatecounter,2)=Valuen(2);
            end
        end
    end
else %If iteration used and solution is converging for a specific G solution (must be more than 1 G solution if this is used)
    try
        disp('tripointer')
        Pn=(i-1)*interval+cmin;
        Gn= FindGCytokines(Pn, Cf, 0, 1.0, 0.01, S1, S2, S3);
        En= (1+exp(-1*sigma*(Weo + Wpe*Pn + Wge*Gn))).^-1;
        MCSFRn= (1+exp(-1*sigma*(Wmro + Wpmr*Pn + Wcmr*Cf + Wemr*En + Wgmr*Gn))).^-1;
        Valuen = ((1+exp(-1*sigma*(Wpo + Wpp*Pn + Wgp*Gn + Wcp*Cf + Vm*S2*MCSFRn/(Km+S2)))).^-1)-Pn;
        PotentialSS1(counter,1) = Pn;
        PotentialSS1(counter,2) = Valuen(tripointer);
    catch
    end
end

for l=1:PnTracker %tracks which G solution we are solving for
    if l==1
        PotentialSS = PotentialSS1;
    elseif l==2
        PotentialSS = PotentialSS2;
    else
        PotentialSS=PotentialSS3;
    end
        
    %First calculate the inital and end slopes
    try
        PotentialSS(1,3) = (PotentialSS(2,2)-PotentialSS(1,2));
        PotentialSS(1,4) = PotentialSS(1,3)*PotentialSS(1,2);
        PotentialSS(end,3) = (PotentialSS(end,2)-PotentialSS(end-1,2));
        PotentialSS(end,4) = PotentialSS(end,3)*PotentialSS(end,2);
    catch %if there is only one PotentialSS
        continue
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
                    [Csol, Gsol]=FindPCytokines(Cf,PotentialSS(i-1,1),PotentialSS(i+1,1),Int/10, S1, S2, S3, tracker+1);
                    CSolution = [CSolution; Csol];
                    GSolution = [GSolution; Gsol];
                end
            catch
                intersect = (PotentialSS(i+1,2)-PotentialSS(i,2)-PotentialSS(i+1,4)*2)/(PotentialSS(i,4)-PotentialSS(i+1,4));
                estInt= (PotentialSS(i,4)*intersect+PotentialSS(i,2))*PotentialSS(i,2);
                if estInt <= 0
                    [CSol, Gsol]=FindPCytokines(Cf,PotentialSS(i,1),PotentialSS(i+1,1),Int/10, S1, S2, S3, tracker+1);
                    CSolution = [CSolution; Csol];
                    GSolution = [GSolution; Gsol];
                end
            end
        else
            cPointer(i) = [];
        end
    end
    
    if PotentialSS(1,4) < 0%(PotentialSS(end,3)+PotentialSS(end,2)<0 && PotentialSS(end,2) > 0) || (PotentialSS(end,3)+PotentialSS(end,2)>0 && PotentialSS(end,2) < 0)
        %         l
        cPointer=[cPointer;PotentialSS(end,1)];
    end

% Find solutions for Function
    for x=1:length(cPointer)
        cmin = cPointer(x);
        cmax = cPointer(x)+Int;
        stoppoint = 1;
        minValue = 1000;
        interval = Int/10;
        tristable=0;
        while stoppoint == 1
            checker=0;
            c = cmin:interval:cmax;
            for i = 1:length(c)
                try
                    c(i);
                    value= Function(c(i));
                    if isreal(value)
                        pointer = value;
                        if abs(pointer) < abs(minValue)
                            minValue = pointer;
                            cValue = c(i);
                            gValue = G(c(i));
                        end
                    end
                catch
                    Pn=c(i);
                    Gn= FindGCytokines(Pn, Cf, 0, 1.0, 0.01, S1, S2, S3);
                    if length(Gn)==2
                        Gn= FindGCytokines(Pn, Cf, 0, 1.0, 0.001, S1, S2, S3);
                    end
                    En= (1+exp(-1*sigma*(Weo + Wpe*Pn + Wge*Gn))).^-1;
                    MCSFRn= (1+exp(-1*sigma*(Wmro + Wpmr*Pn + Wcmr*Cf + Wemr*En + Wgmr*Gn))).^-1;
                    Valuen = ((1+exp(-1*sigma*(Wpo + Wpp*Pn + Wgp*Gn + Wcp*Cf + Vm*S2*MCSFRn/(Km+S2)))).^-1)-Pn;
                    target=Valuen(1);
                    gvalue=Gn(1);
                    for j=2:length(Valuen)
                        if abs(Valuen(j))<abs(target)
                            target=Valuen(j);
                            gvalue=Gn(j);
                        end
                    end
                    value=target;
                    pointer = value;
                    if abs(pointer) < abs(minValue)
                        minValue = pointer;
                        cValue = c(i);
                        gValue = gvalue;
                    end
                end
            end
            if checker==0
                cmin = cValue - interval;
                cmax = cValue + interval;
                interval = interval/10;
            end
            if interval < 10^-17 %Solution cannot be found
                break
            end
            if abs(minValue) < Tolerance
                stoppoint = 0;
                tempSol = [tempSol;cValue];
                gtempSol= [gtempSol;gValue];
            end
        end
    end
end


%Remove all repeating Solutions
tempSol=[tempSol;CSolution];
gtempSol=[gtempSol;GSolution];

if ~isempty(tempSol)
    tempSol2 = tempSol(1);
    gtempSol2= gtempSol(1);
else
    tempSol2 = [];
    gtempSol2=[];
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
        gtempSol2=[gtempSol2; gtempSol(i)];
    end
end
% finished removing all repeating values

CSolution=tempSol2;
GSolution=[];
CSolution = sort(CSolution);
%re-sort GSolution according to CSolution
for x=1:length(CSolution)
    for i=1:length(tempSol2)
        if tempSol2(i)==CSolution(x)
            GSolution=[GSolution;gtempSol2(i)];
        end
    end
end

if isempty(CSolution) && Iteration < 10 && abs(PotentialSS(end,3))-abs(PotentialSS(end-1,3))<0
    Iteration = Iteration + 1;
    PotentialSS(end,1);
    PotentialSS(end-1,1);
    return
end

for i=length(CSolution):-1:1
    if CSolution(i) < 0
        CSolution(i)=[];
        GSolution(i)=[];
    end
end

if tracker==1
    Iteration = 0;
end
end