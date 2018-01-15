function [ ] = Bifurcation( Cytokine )
%Generates bifurcation figures published in "Dynamic Analysis of Cytokine-
%Induced Differentiation of Granulocyte-Monocyte Progenitor Cells" by
%Bronson R. Weston et al. 
global S1 S2 S3 Var constantP
Int=0.01;
%--------------------------------------------------------------------------

S1=0; S2=0; S3=0;

%Initial Conditions to start bifurcation solver
naive= [0.08, 0.133];
gran=   [1, 0.05];
mono=   [0.05, 1];
unstableSS1 = [0.5, 0.473];
unstableSS2= [0.05, 0.3];
options= optimset('TolFun', 1e-7, 'TolX', 1e-7);
Save='on';
box on
%% Cytokine Settings
if strcmp(Cytokine,'M-CSF')
    Var='M-CSF';
    MinVar=0; Step=0.01; MaxVar=1.5;
    iSS=[unstableSS1; naive; gran; mono];
elseif strcmp(Cytokine,'G-CSF')
    Var='G-CSF';
    MinVar=0; Step=0.01; MaxVar=1.5;
    iSS=[unstableSS2; unstableSS1; naive; gran; mono];
elseif strcmp(Cytokine,'GM-CSF')
    Var='GM-CSF';
    MinVar=0; Step=0.005; MaxVar=1.5;
    iSS=[unstableSS1; naive; gran; mono];
else
    disp('Cytokine Input Invalid')
    return
end
%%

fig= figure('Position', [25, 50, 1300, 500]);
ax1=subplot(1,3,1);
ax2=subplot(1,3,2);
ax3=subplot(1,3,3);

%Find solution
for i=1:length(iSS(:,1))
    lastParm=iSS(i,:);
    backwards=0;
    Data=[];
    var=MinVar-Step;
    while var < MaxVar && var >= MinVar-Step;
        if backwards== 0
            var= var+Step;
        else
            var= var-Step;
        end
        var
        ChangeVar(Var,var);
        if ~isempty(lastParm)
            parm0=[lastParm(1),lastParm(2)];
            parm=fminsearch(@SSQ,parm0,options);
            [dCdt, dPdt] = dcfdp( parm(1), parm(2), Int);
            if abs(dCdt) <1e-4 && abs(dPdt) <1e-4
                lastParm = parm;
                Data=  [Data; [var,CalculateC(lastParm(1),lastParm(2)),lastParm(1),lastParm(2)]];
            else
                cChange=Data(end,3)-Data(end-1,3); pChange=Data(end,4)-Data(end-1,4);
                parm0=[Data(end,3)+cChange,var];
                constantP=Data(end,4)+pChange;
                parm=fminsearch(@SSQ2,parm0,options);
                var=parm(2);
                [dCdt, dPdt] = dcfdp( parm(1), constantP, Int);
                if abs(dCdt) <1e-4 && abs(dPdt) <1e-4
                    lastParm=[parm(1),constantP];
                    Data=  [Data; [var,CalculateC(lastParm(1),lastParm(2)),lastParm(1),lastParm(2)]];
                    if var<Data(end-1,1)
                        backwards=1;
                    else
                        backwards=0;
                    end
                else
                    break;
                end
            end
        end
    end
    
    currentStability=Stable([Data(1,3),Data(1,4)]);
    pointer=0;
    while pointer<length(Data)
        pointer=pointer+1;
        ChangeVar(Var,Data(pointer,1))
        stability=Stable([Data(pointer,3),Data(pointer,4)]);
        if stability~=currentStability
            if currentStability==1
                linetype='';
            else
                linetype='--';
            end
            hold(ax1, 'on')
            plot(ax1, Data(1:pointer-1,1),Data(1:pointer-1,2),[linetype 'r'],'LineWidth',2.5)
            hold(ax2, 'on')
            plot(ax2, Data(1:pointer-1,1),Data(1:pointer-1,3),[linetype 'k'],'LineWidth',2.5)
            hold(ax3, 'on')
            plot(ax3, Data(1:pointer-1,1),Data(1:pointer-1,4),[linetype 'b'],'LineWidth',2.5)
            Data= Data(pointer:end,:);
            pointer=0;
            currentStability=stability;
            hold on
        end
    end
    if currentStability==1
        linetype='';
    else
        linetype='--';
    end
    
    hold(ax1, 'on')
    plot(ax1, Data(:,1),Data(:,2),[linetype 'r'],'LineWidth',2.5)
    set(gca, 'fontsize',15)
    hold(ax2, 'on')
    plot(ax2, Data(:,1),Data(:,3),[linetype 'k'],'LineWidth',2.5)
    set(gca, 'fontsize',15)
    hold(ax3, 'on')
    plot(ax3, Data(:,1),Data(:,4),[linetype 'b'],'LineWidth',2.5)
    set(gca, 'fontsize',15)
    hold on
end

%--------------- Label Figure ----------------------------------

hold(ax1, 'on')
xlabel(ax1, Var, 'fontsize',20)%,'fontweight','b')
ylabel(ax1, '[C/EBP]_T', 'fontsize',20)%,'fontweight','b')
set(ax1, 'fontsize',18)%,'fontweight','b')
title(ax1, 'A', 'fontsize',20)%,'fontweight','b')
ylim(ax1, [0 1.05])
xlim(ax1, [MinVar MaxVar])
box on

hold(ax2, 'on')
% xlabel(ax2, Var, 'fontsize',18)%,'fontweight','b')
% ylabel(ax2, '[C/EBP]_F', 'fontsize',18)%,'fontweight','b')
set(ax2, 'fontsize',18)%,'fontweight','b')
xlabel(ax2, Var, 'fontsize',20)%,'fontweight','b')
ylabel(ax2, '[C/EBP]_F', 'fontsize',20)%,'fontweight','b')
ylim(ax2, [0 1.05])
title(ax2, 'B', 'fontsize',20)%,'fontweight','b')
xlim(ax2, [MinVar MaxVar])
box on

hold(ax3, 'on')
xlabel(ax3,Var, 'fontsize',20)%,'fontweight','b')
ylabel(ax3,'[PU.1]', 'fontsize',20)%,'fontweight','b')
set(ax3, 'fontsize',18)%,'fontweight','b')
title(ax3, 'C', 'fontsize',20)%,'fontweight','b')
ylim(ax3,[0 1.05])
xlim(ax3,[MinVar MaxVar])
box on
ax2.Box = 'on';
ax1.Box = 'on';
%--------------- Save Figure ----------------------------------
% if strcmp(Save,'on')
%     saveas(fig, [pwd '/Bifurcation/' Var ' Bifurcation.jpg'])
% end

end 
