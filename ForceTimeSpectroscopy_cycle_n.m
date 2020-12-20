% This script was developed for the following publication:
% "Multi-scale study of the architecture, topography and mechanics of the
% human ovary from prepuberty to menopause: a blueprint for next-generation bioengineering and diagnosis"
% Ouni et al., currently under review in Nature communication.
% This script permits to calculate the viscoelastic response of a tissue
% from AFM measurements.
% The raw AFM data use in this study can be accessed here:

% This file works on the '.OUT' AFM file format from JPK.
% This scripts lets you select one file at a time and analyse several force
% cycles.
% The analysis is not automatic but interactive. You can run the analysis
% of this file and perform all the steps along.
% Please, read all the comments along the text before starting.
% Before using this file, it is recomendet to consult Matlab documentation
% for 'getpts' and 'pause' functions, in order to properly use this script.
% For more details, please Contact Kalina Tamara Haas: kalina.haas@inrae.fr
clear all
clc
%% get a folder containing files which you want to process
d = uigetdir(pwd, 'Select a folder for treated samples');
addpath(d)
cd(d)
files = dir('*.out');
fnames = {files.name}';
str = {fnames{:}}';
%% choose proper colums
num_lines1 = 1;
col_t=1; % time column
col_force=4; %  force column from your file
col_hight=3; %  height column from your file
%% display all the availabe files to chose for batch processing
s = listdlg('PromptString','Select files for analysis:',...
                'SelectionMode','multiple',...
                'ListString',str,...
                'Name', 'Chose a files you want to analyse for treated samples',...
                'ListSize',[500,500]);
ns = numel(s); % num of files selected
ChorNamList = cell(ns,1);
Energy = zeros([ns,1]);
%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);%number of columns
% Specify range and delimiter
opts.DataLines = [7, Inf]; % row range
opts.Delimiter = " ";%delimeter
% Specify column names and types
opts.VariableTypes = ["double", "double", "double", "double", "double"]; % data type
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";
%%
    clear T;
    clear Force;
    clear NoseData;
% Import the data
    FileName = str{s(1)};
    data = table2array(readtable(FileName, opts));
    leg = str{s(1)};
    clear T Force Height
    T(:,1) = data(:,col_t);
    Force(:,1) = data(:,col_force); % force from AFM measurements
    Height(:,1) = data(:,col_hight); % cantilever hight
%% peak detection parameters
    minPeakFN = 0.25; % min peak haight in the normalized Force to detect the force decay
    minpeaknumFN = 10; % minimum number of Force decay peaks, estimate value, better to give slightly bigger value
    minPeakFNdecend = 0.06; % a thershold to detect the Tend of the Force decays. Last peak is alway wrong, it detects the last big 
    minPeakDist = numel( Force)/minpeaknumFN;
% ploting force to determine number of cycles in Force-Time plot
figure;
plot(Force/max(Force));
%% number of cycles you wish to analyse
prompt = {'Enter numebr of cycles to analyse:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'5'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
nc = str2num(answer{1}); % number of cyclec to analyse
%% Loop over number of cycles you wish to analyse
for i = 1:nc
% select and get force spikes maxima and and of exponential decay (points B
% and C in Figure 1) for i-th cycle
close all
f1 = figure;
plot(Force/max(Force));
pause; %*****************here matlab code execution is paused.
% The user can zoom in the matlab plot to place data point in a proper
% position. Please zoom in on the nth cycle (in above publication we analysed
% the first cycle only). 
% After zooming in, press any key, and then place data cursor
% with mouse left-click.
% Select two data points:
% (1) first point click on the Force maximum that corresponds to the
% begining of exponential decays (point B on figure 1)
% (2) second point, select the end of exponential decay (Point C in figure
% 1)
% Please, read below tips on using Matlab getpts function. To know more,
% consult matlab documentation.
%Please select the position of a maxima of the first cycle .

[x1,y1] = getpts(f1); % gete selected data coordinates
locsFmax{i} = round(x1(1)); % position of a maxima of the first (nth) cycle
locsFFdecend{i} = round(x1(2)); % position of a end of the first (nth) cycle
close(f1) % close the figure
%% *************** TIPS ON USING GETPTS***************
% getpts Use normal button clicks to add points. A shift-, right-, 
% or double-click adds a final point and ends the selection. Pressing Return 
% or Enter ends the selection without adding a final point. Pressing Backspace 
% or Delete removes the previously selected point.
%% ploting the selected data + few points before force maxima.
% Now we have o select point A shown in Figure 1.
% It corresponds to a x value where force starts to increase (and height
% decrease).
scale = mean(Force(1:50));
stscale = std(Force(1:50));
par = 20; % exlude few last points from the decay  
fig = figure;
hold on;
plot(T(1:locsFFdecend{i}(1)-par),Force(1:locsFFdecend{i}(1)-par))
% axis([ T(1) T(locsFmax{i}(1)) scale-100*stscale scale+100*stscale])
pause %*****************here matlab code execution is paused.
% The user can zoom in the matlab plot to place data point in a proper
% position. Please zoom in on the nth cycle 
% After zooming in, press any key, and then place data cursor
% with mouse left-click.
[xi{i},yi] = getpts(fig);
deltaH(i) = Height(round(xi{i})) - Height(locsFmax{i}(1))
close(fig)
%% fitting to selected data to recove viscoelastic constant
% model fitted: sum of two exponential decays and a linear trend
% the two exponential correspond to fast and slow viscoelastic constant
funlist = {1, @(coef,xdata) exp(-xdata*coef(1)), @(coef,xdata) exp(-xdata*coef(2)), @(coef,xdata) xdata*coef(3)};
% initialization of fit parameters
NLPstart = [1,1,1]; % initial value
clear coef ABC
times = T(locsFmax{i}(1):locsFFdecend{i}(1)-par);
forces = Force(locsFmax{i}(1):locsFFdecend{i}(1)-par);
times = times - times(1);
figure;
hold on;
plot(times,forces)
% try
[INLP,ILP] = fminspleas(funlist,NLPstart,times,forces,[],[],[]);%,options
coef(i,:) = INLP;
ABC(i,:) = ILP;

fun = @(x) ILP(1) + ILP(2)*exp(-INLP(1)*x)+ ILP(3)*exp(-INLP(2)*x)+ ILP(4)*INLP(3)*x;
valfun = fun(T(locsFmax{i}(1):locsFFdecend{i}(1)-par));
residuals(i,1) = sum(valfun - Force(locsFmax{i}(1):locsFFdecend{i}(1)-par)).^2;
residualsall{i,1} = valfun - Force(locsFmax{i}(1):locsFFdecend{i}(1)-par);
ezplot(fun,times)

title(str{s(1)})
viscoelas(i,1) = INLP(1); % viscoelastic constant 1, from exponential decay 1
viscoelas(i,2) = INLP(2); % viscoelastic constant 2, from exponential decay 2
a_on_deltaH(i,1) = ILP(2)/deltaH(i);
a_on_deltaH(i,2) = ILP(3)/deltaH(i);
lincomp(i,1)  = ILP(1);
lincomp(i,2)  = ILP(4)*INLP(3);
fitfun{i} = fun;
title(FileName)
savefig(strcat(FileName,'.fig'))
i
pause %*****************here matlab code execution is paused, press any key to continue
end  %close the loop over all force cycles
%% plotting the results
figure
nr = 3;
ha = tight_subplot(nr,1,[.01 .03],[.1 .01],[.01 .01]);
axes(ha(1)); plot(viscoelas(:,1),'.b');hold on;plot(viscoelas(:,2),'.r')
% subtitle('Viscoelesticity')
axes(ha(2)); plot(a_on_deltaH(:,1),'.b');hold on;plot(a_on_deltaH(:,2),'.r')
% subtitle('a/\DeltaX')
axes(ha(2)); plot(lincomp(:,1),'.b');hold on;plot(lincomp(:,2),'.r')
% subtitle('Linear component')
%% saving
Sample_names = str(s);
% save to mat file
uisave({'residualsall','viscoelas','a_on_deltaH','lincomp','Sample_names',...
    'funlist','deltaH','residuals','fitfun','ABC','coef','Force','Height','locsFFdecend','locsFmax'},strcat(Sample_names{1}(1:end-4),'.mat'))
%
viscoelas1 = viscoelas(:,1);
viscoelas2 = viscoelas(:,2);
a_on_deltaH1 = a_on_deltaH(:,1);
a_on_deltaH2 = a_on_deltaH(:,2);
lincomp1 = lincomp(:,1);
lincomp2 = lincomp(:,2);
%%
% save to excel sheet
T = table(viscoelas1,viscoelas2,a_on_deltaH1,a_on_deltaH2,residuals,lincomp1,lincomp2);
writetable(T,strcat(Sample_names{1}(1:end-4),'.xls'),'Sheet',1)
