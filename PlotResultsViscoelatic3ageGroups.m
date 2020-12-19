% This script was developed for the following publication:
% "Multi-scale study of the architecture, topography and mechanics of the
% human ovary from prepuberty to menopause: a blueprint for next-generation bioengineering and diagnosis"
% Ouni et al., currently under review in Nature communication.
% This script permits to plot the final analysed data obtained from viscoelastic response of a tissue
% from AFM measurements.
% Before Running this script, please load 'CombineData3ageGroupsAllsamples.mat' file into
% Matlab workspace. This file contains final analysed data for all age
% groups and all analysed samples.
% Agegroups (Prepub), Menopausal (Menop), and Reproductive (Repro).
% The raw AFM data, as well as all the stages of analysis with Matlab, and the required file can
% be found here:

%% define figure properties
width = 2.5;% Width in inches
height = 3; %height in inches
alw=1.5;% AxesLineWidth
fsz=10; % Fontsize
lw=2;%  % LineWidth
msz=15; % MarkerSize
fszl=8; % Fontsize legend
fontname='Helvetica'; % Fontsize name
tickdir = 'out';
ticksize= [0.02 0.035]; % Fontsize legend
colp = [228,117,0]/255;
colB = [106, 54, 217]/255;
colG = [106, 54, 217]/255;
mc = [33,30,40]/255;
medc = mc;
fcol = [colB;colB;colB];
fcoledge = [colG;colG;colG];
%% reproductive stage all fited data
% the columns in DataForMatLabS/1/2/3 correspond to the following values
% obtained from the fit%
% ViscoElas1 (slow), Viscoelas2 (fast), A/DeltaH1, A/DeltaH2, Residuals (from the fit), Lincomp1, lincomp2

% getting slow viscoelastic constant (Viscoelas2) by selecting third column:
%% reproductive stage all fited data
% viscoelas2
ReproV2 = DataForMatLabS1repro(:,2);
ReproV2(find(isnan(ReproV2)))=[];
%% *************** prepuberty all fited data
% viscoelas2
PrePubV2 = DataForMatLabS2prepub(:,2);
PrePubV2(find(isnan(PrePubV2)))=[];
%% Menopausal data
% viscoelas2
MenoV2 = DataForMatLabS3menop(:,2);
MenoV2(find(isnan(MenoV2)))=[];
%%
% ViscoElas1, Viscoelas2, A/DeltaH1, A/DeltaH2, Res, Lincomp1, lincomp2
outlier = [30 70]; % This step allows to filter out wrongly fitted data 
% obtained with ForceTimeSpetroscopy_AFMviscoelasticModel.m function
%% *******Plot ViscoElas1 viscoelastic constant 1 (slow response)
datviolin2{1} = rmoutliers(PrePubV2,'percentiles',outlier);
datviolin2{2} = rmoutliers(ReproV2,'percentiles',outlier); 
datviolin2{3} = rmoutliers(MenoV2,'percentiles',outlier);
% ylab = '\tau_{viscoelas}^{fast} (s^{-1})';
ylab = '\tau_{viscoelas}^{slow} (s^{-1})';
legendlab = {'Prepub','Repro','Menop'};
xtixklab = legendlab;
bw = [];
grpos = 1:numel(legendlab);
scaleFactor = 0.1;
ylim1 = [-0.1,3];
ngr = numel(legendlab);
normfac = 0.45;

%
[~,~,l1,~,~,~,~,med,avg] = violonkd(datviolin2,bw,normfac,grpos,fcol,fcoledge,mc,medc,0.75,...
width,height,alw,fsz,lw,fszl,fontname,[],ylab,xtixklab);
hold on
legend off
for ii = 1:ngr
nd = numel(datviolin2{ii});

x = grpos(ii) + scaleFactor * randn(nd,1);% * exp(randomPhaseAngles' * i);
plot(x,datviolin2{ii},'.','Color',colp,'MarkerSize',5)%[255 178 33]/255
end
% plot(grpos,avg,':','MarkerSize',msz,'Color',mc)
plot(grpos,med,'.','MarkerSize',msz,'Color',medc)
% ylim(ylim1)
xtickangle(45)
%
ngr = numel(legendlab);
group = [];
datboxplot = [];
for i = 1:ngr
  nigr = numel(datviolin2{i});  
  group = [group;repmat({legendlab{i}}, nigr, 1)];
  datboxplot = [datboxplot;datviolin2{i}]; 
end

[p,tbl,stats] = kruskalwallis(datboxplot,group);
c = multcompare(stats,'CType','tukey-kramer','Alpha',0.05);
%% %% *******Plot ViscoElas2 viscoelastic constant 1 (fast response)
datviolin1{1} = rmoutliers(PrePubV1,'percentiles',outlier);
datviolin1{2} = rmoutliers(ReproV1,'percentiles',outlier); 
datviolin1{3} = rmoutliers(MenoV1,'percentiles',outlier);
ylab = '\tau_{viscoelas}^{fast} (s^{-1})';

legendlab = {'Prepub','Repro','Menop'};
xtixklab = legendlab;
bw = [];
grpos = 1:numel(legendlab);
scaleFactor = 0.1;
ylim1 = [-0.1,3];
ngr = numel(legendlab);
normfac = 0.45;

%
[~,~,l1,~,~,~,~,med,avg] = violonkd(datviolin1,bw,normfac,grpos,fcol,fcoledge,mc,medc,0.75,...
width,height,alw,fsz,lw,fszl,fontname,[],ylab,xtixklab);
hold on
legend off
for ii = 1:ngr
nd = numel(datviolin1{ii});

x = grpos(ii) + scaleFactor * randn(nd,1);% * exp(randomPhaseAngles' * i);
plot(x,datviolin1{ii},'.','Color',colp,'MarkerSize',5)%[255 178 33]/255
end
% plot(grpos,avg,':','MarkerSize',msz,'Color',mc)
plot(grpos,med,'.','MarkerSize',msz,'Color',medc)
% ylim(ylim1)
xtickangle(45)
%
ngr = numel(legendlab);
group = [];
datboxplot = [];
for i = 1:ngr
  nigr = numel(datviolin1{i});  
  group = [group;repmat({legendlab{i}}, nigr, 1)];
  datboxplot = [datboxplot;datviolin1{i}]; 
end

[p,tbl,stats] = kruskalwallis(datboxplot,group);
c = multcompare(stats,'CType','tukey-kramer','Alpha',0.05);
