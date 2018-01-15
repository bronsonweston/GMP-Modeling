function [] = Nullclines3D(s1, s2, s3, Load, Title)
global S1 S2 S3 PVALUES CVALUES 
%NullTrajSystem(Title)
%Plots the nullclines of system 2 and plots trajectory lines given the
%inital conditions specified 
NumberOfCells=50;
StandDev=0.2;
K=7.5;
Pint=0.01;
Cint=0.01;


Cvalues = -0.02:0.01:1.11;
Pvalues = zeros(length(Cvalues),4);
Solutions=FindPCytokines(Cvalues(1), 0, 1, Pint, s1, s2, s3 );
S1=s1; S2=s2; S3=s3;
i=0;
backwords=0;

fig1 = figure; ax1=gca;
fig2 = figure; ax2=gca;
fig3 = figure; ax3=gca;
fig4 = figure; ax4=gca;
fig5 = figure; ax5=gca;

% hold on

if strcmp(Load,'yes')
    try 
        load([pwd '/NullclineData/GMCSF = ' num2str(s1) ' MCSF = ' num2str(s2) ' G-CSF = ' num2str(s3) '.mat'], 'PVALUES', 'CVALUES')
    catch
        Load='no';
    end
else
    Load='no';
end

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
                    Pvalues(i,:)=[];
                    disp('deleted value, even number of solutions')
                    i=i-1;
                    Solutions=OldSolutions;
                    continue
                end
            end
%             Solutions
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
[Pvalues1,Pvalues2]= GraphMatrixCompiler(Pvalues);

if isempty(Pvalues2)
    plot(ax1, Pvalues1(:,1),Pvalues1(:,2), 'b', 'LineWidth',2);
%     global CT P G
    [C,G,I]=CalculateCdPGI(Pvalues1(:,1),Pvalues1(:,2),[s1,s2,s3]);
    plot(ax2, C,Pvalues1(:,2), 'b', 'LineWidth',2);
%     figure(3); hold on;
    plot3(ax3, Pvalues1(:,1),Pvalues1(:,2), G, 'b', 'LineWidth',2);
%     figure(4); hold on;
    plot3(ax4, C ,Pvalues1(:,2), G, 'b', 'LineWidth',2);
    plot3(ax5, C ,Pvalues1(:,2), I, 'b', 'LineWidth',2);
    P = Pvalues1(:,2);
    CT=Pvalues1(:,1);
    
else
    stopping='stop... more to do'
    return
    figure(1); hold on;
    plot(Pvalues1(:,1),Pvalues1(:,2), 'b', 'LineWidth',2);
    plot(Pvalues2(:,1),Pvalues2(:,2), 'b', 'LineWidth',2);
    [C, G]=CalculateCG(Pvalues1(:,1),Pvalues1(:,2),[s1,s2,s3]);
    figure(3); hold on;
    plot3(Pvalues1(:,1),Pvalues1(:,2), G, 'b', 'LineWidth',2);
    figure(4); hold on;
    plot3(C,Pvalues1(:,2), G, 'b', 'LineWidth',2);
    figure(2); hold on;
    plot(C,Pvalues1(:,2), 'b', 'LineWidth',2);
    [C, G]=CalculateC(Pvalues2(:,1),Pvalues2(:,2),[s1,s2,s3]);
    plot(C,Pvalues2(:,2), 'b', 'LineWidth',2);
    figure(3); hold on;
    plot3(Pvalues1(:,1),Pvalues2(:,2), G, 'b', 'LineWidth',2);
    figure(4); hold on;
    plot3(C,Pvalues2(:,2), G, 'b', 'LineWidth',2);
end


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
    save([pwd '/NullclineData/GMCSF = ' num2str(s1) ' MCSF = ' num2str(s2) ' G-CSF = ' num2str(s3) '.mat'], 'PVALUES', 'CVALUES')
else
    CvaluesC=CVALUES;
end
[Cvalues1,Cvalues2]= GraphMatrixCompiler(CvaluesC);

if isempty(Cvalues2)
%     figure(1); hold on;
    hold(ax1, 'on')
    c1 = plot(ax1,Cvalues1(:,2),Cvalues1(:,1), 'r', 'LineWidth',2);
    [C,G,I]= CalculateCdCGI(Cvalues1(:,2),Cvalues1(:,1),[s1,s2,s3]);
%     figure(2); hold on;
    hold(ax2, 'on')
    plot(ax2,C,Cvalues1(:,1), 'r', 'LineWidth',2);
%     figure(3); hold on;
    hold(ax3, 'on')
    plot3(ax3,Cvalues1(:,2),Cvalues1(:,1), G, 'r', 'LineWidth',2);
%     figure(4); hold on;
    hold(ax4, 'on')
    plot3(ax4,C,Cvalues1(:,1), G, 'r', 'LineWidth',2);
    hold(ax5, 'on')
    plot3(ax5, C ,Cvalues1(:,1), I, 'r', 'LineWidth',2);
end


% return
axes(ax1); hold on;
box on
set(gca, 'fontsize',13)
ylabel('[PU.1]','fontsize', 16)
xlabel('[C/EBP]_F','fontsize', 16)
if exist('Title','var')
    title(Title, 'FontWeight', 'normal')
else
    title({['GM-CSF = ', num2str(s1)], ['M-CSF = ', num2str(s2)], ['G-CSF = ', num2str(s3)]})
end
xlim([0 1.1])
ylim([0 1.1])

axes(ax2); hold on;
box on
set(gca, 'fontsize',13)
ylabel('[PU.1]', 'fontsize', 16)
xlabel('[C/EBP]_T','fontsize', 16)
if exist('Title','var')
    title(Title, 'FontWeight', 'normal')
else
    title({['GM-CSF = ', num2str(s1)], ['M-CSF = ', num2str(s2)], ['G-CSF = ', num2str(s3)]})
end
xlim([0 1.1])
ylim([0 1.1])


axes(ax3); hold on;
box on
set(gca, 'fontsize',13)
ylabel('[PU.1]', 'fontsize', 16)
xlabel('[C/EBP]_F','fontsize', 16)
zlabel('[Gfi-1]','fontsize', 16) 
xlim([0 1.1])
ylim([0 1.1])
zlim([0 1.1])
hold off

axes(ax4); hold on;
box on
set(gca, 'fontsize',13)
ylabel('[PU.1]', 'fontsize', 16)
xlabel('[C/EBP]_T','fontsize', 16)
zlabel('[Gfi-1]','fontsize', 16) 
xlim([0 1.1])
ylim([0 1.1])
zlim([0 1.1])


axes(ax5); hold on;
box on
set(gca, 'fontsize',16)
ylabel('[PU.1]', 'fontsize', 18)
xlabel('[C/EBP]_T','fontsize', 18)
zlabel('[IRF8]','fontsize', 18) 
xlim([0 1.1])
ylim([0 1.1])
zlim([0 1.1])
S1 = 0; S2 = 0; S3 = 0;
end

