function [] = NullclinesProjections(s1, s2, s3, Load, Title)
%%
% This function generates nullclines in the PU.1:C/EBP phaseplane. It also
% simulates the differentiation of 100 stochastically generated cells and
% plots their trajectories over time within the phase plane.
% Code requires 'FindPCytokines.m', 'FindCfCytokines.m',
% 'GraphMatrixCompiler.m', 'CalculateCdC.m', 'CalculateCdP.m',
% 'SystemODE.m', 'FindGCytokines.m' and 'FindPCytokines.m'
%%
global S1 S2 S3 PVALUES CVALUES
NumberOfCells=100;
StandDev=0.2;
K=7.5;
%Initial step size for iterative solvers
Pint=0.01; 
Cint=0.01;


Cvalues = -0.02:0.01:1.11; %interval of solution
Pvalues = zeros(length(Cvalues),4);
Solutions=FindPCytokines(Cvalues(1), 0, 1, Pint, s1, s2, s3 );
S1=s1; S2=s2; S3=s3;
i=0;
backwords=0;

fig1 = figure; ax1=gca;
fig2 = figure; ax2=gca;
hold on

%If told to load file and file is available, it loads the data. Otherwise,
%the function finds the solution
if strcmp(Load,'yes')
    try 
        load([pwd '/NullclineData/GMCSF = ' num2str(s1) ' MCSF = ' num2str(s2) ' G-CSF = ' num2str(s3) '.mat'], 'PVALUES', 'CVALUES')
    catch
        Load='no';
    end
else
    Load='no';
end

%Find PU.1 nullcline solution if not loaded
if strcmp(Load,'no')
    while i<length(Cvalues)
            i=i+1;
            Cf= Cvalues(i)
            OldSolutions=Solutions;
            Solutions=FindPCytokines(Cf, 0, 1, Pint, s1, s2, s3 );
            if rem(length(Solutions),2)==0
                Solutions=FindPCytokines(Cf, 0, 1, Pint/10, s1, s2, s3 );
                if rem(length(Solutions),2)==0
                    Cvalues(i)=[];
                    try
                        Pvalues(i,:)=[];
                    catch
                    end
                    disp('deleted value, even number of solutions')
                    i=i-1;
                    Solutions=OldSolutions;
                    continue
                end
            end
            if length(Solutions)~=length(OldSolutions)
                if backwords==0
                    Cvalues=[Cvalues(1:i-2),Cvalues(i-1):0.001:Cvalues(i),Cvalues(i+1:end)];
                    Pvalues=[Pvalues;zeros(length(Cvalues(i-1):0.001:Cvalues(i))-2,length(Pvalues(1,:)))];
                    i=i-1;
                    backwords=1
                    Solutions=OldSolutions;
                    Cf=Cvalues(i);
                else
                    backwords=0
                end
            end
            Solutions=sort(Solutions);
            Pvalues(i,1)= Cf;
            for j=1:length(Solutions)
                Pvalues(i,j+1)=Solutions(j);
            end
    end
    PVALUES=Pvalues;
else
    Pvalues = PVALUES;
end

%Compiles all solutions into 1 line
[Pvalues1,Pvalues2]= GraphMatrixCompiler(Pvalues);

%Plots PU.1 solution
if isempty(Pvalues2)
    hold(ax1, 'on')
    plot(ax1, Pvalues1(:,1),Pvalues1(:,2), 'b', 'LineWidth',2);
    C=CalculateCdP(Pvalues1(:,1),Pvalues1(:,2),[s1,s2,s3]); %Calculates the results for C/EBPtotal
    hold(ax2, 'on')
    plot(ax2, C,Pvalues1(:,2), 'b', 'LineWidth',2);
else %in case a second, independent solution was found
    hold(ax1, 'on')
    plot(ax1, Pvalues1(:,1),Pvalues1(:,2), 'b', 'LineWidth',2);
    plot(ax1, Pvalues2(:,1),Pvalues2(:,2), 'b', 'LineWidth',2);
    C=CalculateCdP(Pvalues1(:,1),Pvalues1(:,2),[s1,s2,s3]); %Calculates the results for C/EBPtotal
    hold(ax2, 'on')
    plot(ax2, C,Pvalues1(:,2), 'b', 'LineWidth',2);
    C=CalculateCdP(Pvalues2(:,1),Pvalues2(:,2),[s1,s2,s3]); %Calculates the results for C/EBPtotal
    plot(ax2, C,Pvalues2(:,2), 'b', 'LineWidth',2);
end


%Find C/EBPf nullcline solution if not loaded
if strcmp(Load,'no')
    PvaluesC = -0.02:0.01:1.11;
    CvaluesC = zeros(length(PvaluesC),4);
    Solutions=FindCfCytokines(PvaluesC(1), 0, 1, Cint, s1, s2, s3);
    i=0;
    backwords=0;
    while i<length(PvaluesC)
        i=i+1;
        P= PvaluesC(i)
        OldSolutions=Solutions;
        Solutions=FindCfCytokines(P, 0, 1, Cint, s1, s2, s3);
        if rem(length(Solutions),2)==0
            Solutions=FindCfCytokines(P, 0, 1, Cint/10, s1, s2, s3);
            if rem(length(Solutions),2)==0
                PvaluesC(i)=[];
                try
                    CvaluesC(i,:)=[];
                catch
                end
                disp('deleted value, even number of solutions')
                i=i-1;
                Solutions=OldSolutions;
                continue
            end
        end
        try
            if length(Solutions)~=length(OldSolutions)
                if backwords==0
                    PvaluesC=[PvaluesC(1:i-2),PvaluesC(i-1):0.001:PvaluesC(i),PvaluesC(i+1:end)];
                    CvaluesC=[CvaluesC;zeros(length(PvaluesC(i-1):0.001:PvaluesC(i))-2,4)];
                    i=i-1;
                    backwords=1;
                    Solutions=OldSolutions;
                else
                    backwords=0;
                end
            end
        catch
        end
        Solutions=sort(Solutions);
        CvaluesC(i,1)= PvaluesC(i);
        for j=1:length(Solutions)
            CvaluesC(i,j+1)=Solutions(j);
        end
    end
    CVALUES=CvaluesC;
    try
        save([pwd '/NullclineData/GMCSF = ' num2str(s1) ' MCSF = ' num2str(s2) ' G-CSF = ' num2str(s3) '.mat'], 'PVALUES', 'CVALUES')
    catch
        save([pwd 'GMCSF = ' num2str(s1) ' MCSF = ' num2str(s2) ' G-CSF = ' num2str(s3) '.mat'], 'PVALUES', 'CVALUES')
    end
else
    CvaluesC=CVALUES;
end

%Compile C/EBPf Solutions
[Cvalues1,Cvalues2]= GraphMatrixCompiler(CvaluesC);

%Plot C/EBP solutions
if isempty(Cvalues2)
    hold(ax1, 'on')
    c1 = plot(ax1, Cvalues1(:,2),Cvalues1(:,1), 'r', 'LineWidth',2);
    C=CalculateCdC(Cvalues1(:,2),Cvalues1(:,1),[s1,s2,s3]); %Calculates the results for C/EBPtotal
    hold(ax2, 'on')
    plot(ax2, C,Cvalues1(:,1), 'r', 'LineWidth',2);
else
    hold(ax1, 'on')
    c1 = plot(ax1, Cvalues1(:,2),Cvalues1(:,1), 'r', 'LineWidth',2); 
    c2= plot(ax1, Cvalues2(:,2),Cvalues2(:,1), 'r', 'LineWidth',2);
    C=CalculateCdC(Cvalues1(:,2),Cvalues1(:,1),[s1,s2,s3]); %Calculates the results for C/EBPtotal
    hold(ax2, 'on')
    plot(ax2, C,Cvalues1(:,1), 'r', 'LineWidth',2);
    C=CalculateCdC(Cvalues2(:,2),Cvalues2(:,1),[s1,s2,s3]); %Calculates the results for C/EBPtotal
    plot(ax2, C,Cvalues2(:,1), 'r', 'LineWidth',2);
end



axes(ax1); hold on
box on
set(gca, 'fontsize',13)
ylabel('[PU.1]','fontsize', 16)
xlabel('[C/EBP]_F','fontsize', 16)
if exist('Title','var')
    title(Title, 'FontWeight', 'normal')
else
    title({['GM-CSF = ', num2str(s1)], ['M-CSF = ', num2str(s2)], ['G-CSF = ', num2str(s3)]})
end
% legend([p1,c1],{'PU.1 nullcline', 'C/EBP\alpha_F nullcline'})
xlim([0 1.1])
ylim([0 1.1])
% fig1=figure(1);
% if exist('Title','var')
%     saveas(fig1, [pwd '/NullclinesCf/GMCSF = ' num2str(s1) ' MCSF = ' num2str(s2) ' G-CSF = ' num2str(s3) ' ' Title '.png'])
% else
%     saveas(fig1, [pwd '/NullclinesCf/GMCSF = ' num2str(s1) ' MCSF = ' num2str(s2) ' G-CSF = ' num2str(s3) '.png'])
% end

axes(ax2); hold on
box on
set(gca, 'fontsize',13)
ylabel('[PU.1]', 'fontsize', 16)
xlabel('[C/EBP]_T','fontsize', 16)
if exist('Title','var')
    title(Title, 'FontWeight', 'normal')
else
    title({['GM-CSF = ', num2str(s1)], ['M-CSF = ', num2str(s2)], ['G-CSF = ', num2str(s3)]})
end
% legend([p1,c1],{'PU.1 nullcline', 'C/EBP\alpha_F nullcline'})
xlim([0 1.1])
ylim([0 1.1])
% fig2=figure(2);
% if exist('Title','var')
%     saveas(fig2, [pwd '/NullclinesCt/GMCSF = ' num2str(s1) ' MCSF = ' num2str(s2) ' G-CSF = ' num2str(s3) ' ' Title '.png'])
% else
%     saveas(fig2, [pwd '/NullclinesCt/GMCSF = ' num2str(s1) ' MCSF = ' num2str(s2) ' G-CSF = ' num2str(s3) '.png'])
% end


%Stochastic Population Generation and Ploting
S1 = 0; S2 = 0; S3 = 0;
Ci=0.142; Pi=0.133; Gi= 0.088; Ii=0.1494; Ei=0.091; GMRi=0.0676; MRi=0.0445; GRi=0.091;
[~,y] = ode45(@SystemODE, [0 50], [Ci, Pi, Gi, Ii, Ei, GMRi, MRi, GRi]);
Ci=y(end,1); Pi=y(end,2); Gi=y(end,3); Ii=y(end,4); Ei=y(end,5); GMRi=y(end,6);  MRi=y(end,7); GRi=y(end,8);

GMHcount = zeros(1, 3);
Data = zeros(NumberOfCells, 7);
GMHcount = zeros(1, 3);
for i=1:length(Data(:,1))
    tCi = normrnd(Ci, StandDev*Ci);
    tPi = normrnd(Pi, StandDev*Pi);
    tGi = normrnd(Gi, StandDev*Gi);
    tIi = normrnd(Ii, StandDev*Ii);
    tEi = normrnd(Ei, StandDev*Ei);
    tGMRi = normrnd(GMRi, StandDev*GMRi);
    tMRi = normrnd(MRi, StandDev*MRi);
    tGRi = normrnd(GRi, StandDev*GRi);
    S1=s1; S2=s2; S3=s3;
    [~,y2] = ode45(@SystemODE, [5 125], [tCi, tPi, tGi, tIi, tEi, tGMRi, tMRi, tGRi]);
    Data(i,:) = [y2(1,1),y2(1,2),tCi, tPi, tGi, tIi, tEi];
    b=y2(:,4)-y2(:,1)+(1/K);
    c=-1*y2(:,1)/K;
    CF=(-1*b+(((b.*b)-4*c).^(1/2)))/2;
    if CF(end,1)/y2(end,2) > 3
        hold(ax1, 'on')
        plot(ax1, CF(:,1),y2(:,2),'--r')
        hold(ax2, 'on')
        plot(ax2, y2(:,1),y2(:,2),'--r')
        GMHcount(1,1) = GMHcount(1,1) + 1;
    elseif CF(end,1)/y2(end,2) < 1/3
        hold(ax1, 'on')
        plot(ax1, CF(:,1),y2(:,2),'--b')
        hold(ax2, 'on')
        plot(ax2, y2(:,1),y2(:,2),'--b')
        GMHcount(1,2) = GMHcount(1,2) + 1;
    else
        hold(ax1, 'on')
        plot(ax1, CF(:,1),y2(:,2),'--g')
        hold(ax2, 'on')
        plot(ax2, y2(:,1),y2(:,2),'--g')
        if y2(end,1) > 0.4
            GMHcount(1,3) = GMHcount(1,3) + 1; %hybrid cells
        end
    end
    hold(ax1, 'on')
    plot(ax1, CF(1,1),y2(1,2),'*c')
    plot(ax1, CF(end,1),y2(end,2),'*k')
    hold(ax2, 'on')
    plot(ax2, y2(1,1),y2(1,2),'*c')
    plot(ax2, y2(end,1),y2(end,2),'*k')
end

% fig1=figure(1);
% if exist('Title','var')
%     saveas(fig1, [pwd '/NullclinesCf/GMCSF = ' num2str(S1) ' MCSF = ' num2str(S2) ' G-CSF = ' num2str(S3) ' ' Title '_projections.png'])
% else
%     saveas(fig1, [pwd '/NullclinesCf/GMCSF = ' num2str(S1) ' MCSF = ' num2str(S2) ' G-CSF = ' num2str(S3) '_projections.png'])
% end
% hold off
% fig2=figure(2);
% if exist('Title','var')
%     saveas(fig2, [pwd '/NullclinesCt/GMCSF = ' num2str(S1) ' MCSF = ' num2str(S2) ' G-CSF = ' num2str(S3) ' ' Title '_projections.png'])
% else
%     saveas(fig2, [pwd '/NullclinesCt/GMCSF = ' num2str(S1) ' MCSF = ' num2str(S2) ' G-CSF = ' num2str(S3) '_projections.png'])
% end
% hold off
S1=0; S2=0; S3=0;
end

