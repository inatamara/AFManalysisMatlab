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
% viscoelas1
ReproV1 = DataForMatLabS1{:,2};
ReproV1(find(isnan(ReproV1)))=[];
% viscoelas2
ReproV2 = DataForMatLabS1{:,3};
ReproV2(find(isnan(ReproV2)))=[];
% a_on_deltaH1
ReproV3 = DataForMatLabS1{:,4};
ReproV3(find(isnan(ReproV3)))=[];
% a_on_deltaH2
ReproV4 = DataForMatLabS1{:,5};
ReproV4(find(isnan(ReproV4)))=[];
% residuals
ReproV5 = DataForMatLabS1{:,6};
ReproV5(find(isnan(ReproV5)))=[];
% lincomp1
ReproV6 = DataForMatLabS1{:,7};
ReproV6(find(isnan(ReproV6)))=[];
% lincomp2
ReproV7 = DataForMatLabS1{:,8};
ReproV7(find(isnan(ReproV7)))=[];

%% *************** prepuberty all fited data
% viscoelas1
PrePubV1 = DataForMatLabS2{:,2};
PrePubV1(find(isnan(PrePubV1)))=[];
% viscoelas2
PrePubV2 = DataForMatLabS2{:,3};
PrePubV2(find(isnan(PrePubV2)))=[];
% a_on_deltaH1
PrePubV3 = DataForMatLabS2{:,4};
PrePubV3(find(isnan(PrePubV3)))=[];
% a_on_deltaH2
PrePubV4 = DataForMatLabS2{:,5};
PrePubV4(find(isnan(PrePubV4)))=[];
% residuals
PrePubV5 = DataForMatLabS2{:,6};
PrePubV5(find(isnan(PrePubV5)))=[];
% lincomp1
PrePubV6 = DataForMatLabS2{:,7};
PrePubV6(find(isnan(PrePubV6)))=[];
% lincomp2
PrePubV7 = DataForMatLabS2{:,8};
PrePubV7(find(isnan(PrePubV7)))=[];

%% Menopausal data
% viscoelas1
MenoV1 = DataForMatLabS3{:,2};
MenoV1(find(isnan(MenoV1)))=[];
% viscoelas2
MenoV2 = DataForMatLabS3{:,3};
MenoV2(find(isnan(MenoV2)))=[];
% a_on_deltaH1
MenoV3 = DataForMatLabS3{:,4};
MenoV3(find(isnan(MenoV3)))=[];
% a_on_deltaH2
MenoV4 = DataForMatLabS3{:,5};
MenoV4(find(isnan(MenoV4)))=[];
% residuals
MenoV5 = DataForMatLabS3{:,6};
MenoV5(find(isnan(MenoV5)))=[];
% lincomp1
MenoV6 = DataForMatLabS3{:,7};
MenoV6(find(isnan(MenoV6)))=[];
% lincomp2
MenoV7 = DataForMatLabS3{:,8};
MenoV7(find(isnan(MenoV7)))=[];

%%
% ViscoElas1, Viscoelas2, A/DeltaH1, A/DeltaH2, Res, Lincomp1, lincomp2
outlier = [25 75]; % This step allows to filter out wrongly fitted data 
% obtained with ForceTimeSpetroscopy_AFMviscoelasticModel.m function
datviolin1{1} = rmoutliers(PrePubV2,'percentiles',outlier);
datviolin1{2} = rmoutliers(ReproV2,'percentiles',outlier); %1./
datviolin1{3} = rmoutliers(MenoV2,'percentiles',outlier);
ylab = '\tau_{viscoelas}^{slow} (s^{-1})';
legendlab = {'Prepub','Repro','Menop'};
xtixklab = legendlab;
bw = [];
grpos = 1:numel(legendlab);
scaleFactor = 0.1;
ylim1 = [-0.1,3];
ngr = numel(legendlab);
normfac = 0.45;

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
