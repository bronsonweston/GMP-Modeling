function [ output_args ] = BarChart( GMCSF, MCSF, GCSF )
%BarChart( GMCSF, MCSF, GCSF ) 
%Returns a bar chart of the protein activity levels at the end of a
%simulation of 100 time units. Requires SystemODE.m to run.
%   Detailed explanation goes here
global S1 S2 S3
S1= 0; S2=0; S3=0; K=7.5; 
Ci=0.14242; Pi=0.13313; Gi= 0.088018; Ii=0.1494; Ei=0.091107; GMRi=0.0686; MRi=0.0449; GRi=0.07725;
[t,y] = ode45(@SystemODE, [0 20], [Ci, Pi, Gi, Ii, Ei, GMRi, MRi, GRi]);
Ci=y(end,1); Pi=y(end,2); Gi=y(end,3); Ii=y(end,4); Ei=y(end,5); GMRi=y(end,6);  MRi=y(end,7); GRi=y(end,8);
S1=GMCSF; S2=MCSF; S3=GCSF;
[t,y] = ode45(@SystemODE, [0 100], [y(end,1), y(end,2), y(end,3), y(end,4), y(end,5), y(end,6), y(end,7), y(end,8)]);
b=y(:,4)-y(:,1)+(1/K);
c=-1*y(:,1)/K;
CF=(-1*b+(((b.*b)-4*c).^(1/2)))/2;
barData1= [y(end,1),CF(end), y(end,3),y(end,8), 0, 0, 0, 0, 0];
x1=1:9;
barData2= [y(end,2), y(end,5), y(end,4), y(end,7)];
x2=5:8;
barData3= [y(end,6)]; 
x3=9;
Labels ={'C/EBP_T', 'C/EBP_F', 'Gfi-1', 'G-CSFR', 'PU.1', 'Egr-2', 'IRF8', 'M-CSFR', 'GM-CSFR'};
figure('Position', [150, 250, 750, 450]);
%figure('Position', [150, 250, 1050, 350]);
% figure()
p1 = bar(x1,barData1);
set(gca,'xticklabel',Labels, 'FontSize', 18, 'ylim', [0 1.1], 'XTickLabelRotation', 45)
ylabel('Protein Activity')
hold on;
p2 = bar(x2,barData2);
p3 = bar(x3,barData3);
set(p1,'FaceColor','r');
set(p2,'FaceColor','b');
set(p3,'FaceColor',[0.8 0 0.8]);

S1= 0; S2=0; S3=0;
end

