function [] = HeatMapping( Num_Cells, Rangex, x_axis, Rangey, y_axis, Receptor, Cell_Type, constantThirdStimulant)
%[] = HeatMapping(Num_Cells, Rangex, x_axis, Rangey, y_axis, Receptor, Cell_Type, constantThirdStimulant)
% Generates a heatmap under the specified conditions.
% constantThirdStimulant is the concentration of a third cytokine, given
% that only two are entered into the x and y axis. For example, if 'G-CSF'
% is entered for the x axis, and 'M=CSF' is entered as the y axis, then the
% value entered into the constantThirdStimulant variable will be a constant
% value of GM-CSF throughout all simulations.
% axis options = MCSF, GCSF, and GMCSF
% range specified as [min, step size, max]
% Receptor options = MCSFR, GCSFR, GMCSFR, off or OFF... reports the
% concentration of the receptor within a default cell
% Cell_Type options = Gran, Mono, Gran+Mono, off or OFF... reports the
% fraction of the population that has differentiated into the specified
% phenotype.

% global save
% if isempty(save)
%     save = 'yes';
% end
StandDev=0.2;

if ~exist('constantThirdStimulant','var')
    constantThirdStimulant=0;
end

if strcmp(Receptor,'off') || strcmp(Receptor,'OFF')
    rinst = 0;
elseif strcmp(Receptor,'M-CSFR')
    rinst = 7;
elseif strcmp(Receptor,'GM-CSFR')
    rinst = 6;
elseif strcmp(Receptor,'G-CSFR')
    rinst = 8;
else
    error('Invalid receptor entry')
end

if strcmp(Cell_Type,'off') || strcmp(Cell_Type,'OFF')
    cinst = 0;
elseif strcmp(Cell_Type,'Gran')
    cinst = 1;
elseif strcmp(Cell_Type,'Mono')
    cinst = 2;
elseif strcmp(Cell_Type,'Gran+Mono')
    cinst = 3;
else
    error('Invalid cell type entry')
end

xax = Rangex(1):Rangex(2):Rangex(3);
yax = Rangey(3):-Rangey(2):Rangey(1);
Data1 = zeros(length(yax), length(xax),3);
Data2 = zeros(length(yax), length(xax));
if rinst ~=0
    Data3 = zeros(length(yax), length(xax));
end
if cinst~= 0
    Data4 = zeros(length(yax), length(xax));
    if cinst==3
        Data5 = zeros(length(yax), length(xax));
    end
end

for i=1:length(xax)
    for j=1:length(yax)
        [Values, DefaultConc] = ValueGenerator_HM(Num_Cells,StandDev, x_axis, xax(i), y_axis, yax(j),constantThirdStimulant);
        value = [Values(1)/max(Values), Values(3)/max(Values), Values(2)/max(Values)];
        Data1(j,i,:)= value.^(1/3);
        Data2(j,i)= DefaultConc(1)/DefaultConc(2);
        if rinst ~= 0
            Data3(j,i)= DefaultConc(rinst);
        end
        if cinst ~= 0
            Data4(j,i)= Values(cinst);
            if cinst==3
                Data4(j,i)= Values(1);
                Data5(j,i)= Values(2);
            end
        end
    end
    xax(i)
end
xticklabels={};
yticklabesl={};
for i=1:length(xax) %generate tick marks (so they increase by increments of 0.1
    if rem(round(xax(i),2),0.2)==0
        xticklabels{i} = num2str(xax(i));
    else
        xticklabels{i} = '';
    end
end
for i=1:length(yax) %generate tick marks (so they increase by increments of 0.1
    if rem(round(yax(i),2),0.2)==0 
        yticklabels{i} = num2str(yax(i));
    else
        yticklabels{i} = '';
    end
end

fig1= figure
image(Data1);        % draw image and scale colormap to values range
set(gca, 'fontsize', 15);
xlabel(x_axis, 'fontsize', 16)
ylabel(y_axis, 'fontsize', 16)
set(gca,'XTick',1:length(xax));
set(gca,'XTickLabel',xticklabels );
set(gca,'YTick',1:length(yax));
set(gca,'YTickLabel',yticklabels);
% title('Population Map')
% if strcmp(save, 'yes')
%     saveas(fig1, [pwd '/Heatmap/' x_axis '+' y_axis '.png'])
%     saveas(fig1, [pwd '/Heatmap/' x_axis '+' y_axis '.fig'])
% end

figure
colormap('hot');   % set colormap
imagesc(Data2, [0 2]);        % draw image and scale colormap to values range
colorbar;          % show color scale
set(gca, 'fontsize', 15);
xlabel(x_axis, 'fontsize', 16)
ylabel(y_axis, 'fontsize', 16)
set(gca,'XTick',1:length(xax));
set(gca,'XTickLabel',xticklabels);
set(gca,'YTick',1:length(yax));
set(gca,'YTickLabel',yticklabels);
title('Default Cell C:P ratio')
hold off
if rinst ~= 0
    figure
    colormap('hot');   % set colormap
    imagesc(Data3, [0 1]);        % draw image and scale colormap to values range
    colorbar;          % show color scale
    %im(2,3,:) = [0,1,0];
    set(gca, 'fontsize', 15);
    xlabel(x_axis, 'fontsize', 16)
    ylabel(y_axis, 'fontsize', 16)
    set(gca,'XTick',1:length(xax));
    set(gca,'XTickLabel',xticklabels);
    set(gca,'YTick',1:length(yax));
    set(gca,'YTickLabel',yticklabels);
    title(strcat(Receptor, ' Concentration'))
    %saveas(gca,'testingPixels.png');
end

if cinst == 1 || cinst == 2
    figure
    colormap('hot');   % set colormap
    imagesc(100*Data4/Num_Cells, [0 100]);        % draw image and scale colormap to values range
    colorbar;          % show color scale
    %im(2,3,:) = [0,1,0];
    c.Label.String = '% of Population';
    c.Label.Rotation = -90;
    c.Label.Position = [2.8, 50, 0];
    c.Label.FontSize = 11;    %im(2,3,:) = [0,1,0];
    set(gca, 'fontsize', 15);
    xlabel(x_axis, 'fontsize', 16)
    ylabel(y_axis, 'fontsize', 16)
    set(gca,'XTick',1:length(xax));
    set(gca,'XTickLabel',xticklabels);
    set(gca,'YTick',1:length(yax));
    set(gca,'YTickLabel',yticklabels);
    if cinst == 1
        title('Granulocyte Count')
    elseif cinst == 2
        title('Monocyte Count')
    else
        error('In Cell Count Figure')
    end   
    %saveas(gca,'testingPixels.png');
end
if cinst == 3
    figure
    colormap('hot');   % set colormap
    imagesc(100*Data4/Num_Cells, [0 100]);        % draw image and scale colormap to values range
    c = colorbar;          % show color scale
    c.Label.String = '% of Population';
    c.Label.Rotation = -90;
    c.Label.Position = [2.8, 50, 0];
    c.Label.FontSize = 11;
    %im(2,3,:) = [0,1,0];
    set(gca, 'fontsize', 15);
    xlabel(x_axis, 'fontsize', 16)
    ylabel(y_axis, 'fontsize', 16)
    set(gca,'XTick',1:length(xax));
    set(gca,'XTickLabel',xticklabels);
    set(gca,'YTick',1:length(yax));
    set(gca,'YTickLabel',yticklabels);
    title('Granulocytes')
    
    figure
    colormap('hot');   % set colormap
    imagesc(100*Data5/Num_Cells, [0 100]);       % draw image and scale colormap to values range
    c = colorbar;          % show color scale
    c.Label.String = '% of Population';
    c.Label.Rotation = -90;
    c.Label.Position = [2.8, 50, 0];
    c.Label.FontSize = 11;    %im(2,3,:) = [0,1,0];
    set(gca, 'fontsize', 15);
    xlabel(x_axis, 'fontsize', 16)
    ylabel(y_axis, 'fontsize', 16)
    set(gca,'XTick',1:length(xax));
    set(gca,'XTickLabel',xticklabels);
    set(gca,'YTick',1:length(yax));
    set(gca,'YTickLabel',yticklabels);
    title('Monocytes')
end
end

