% Atlan et al 2024 Code:
% Reproduces figures for Atlan et al 2024, Nat Comms

% Requires the following matlab classes:
% MouseArray (holds all mice together)
% MouseSummary (holds individual mice handles)
% RegModel (opens regression model results)
%% Start here always:
PATH = 'H:\My Drive\Matlab_Code'; %Edit this according to the path of the code
data_path = 'D:\Atlan2021 data\'; %Edit this according to the path where the data is located

%% Run during the first time only - reset path to match your folder locations:
array_path = [PATH, '\full_array_25mice_local.mat']; % Path to MouseArray file.
X = MouseArray(array_path); % load array;
X.reset_path(data_path);
% run to manually reset data path for each object in the array:
for ii = 1:numel(X.MOUSE_ARRAY)
    X.MOUSE_ARRAY(ii).PATH = data_path;
end
X.save_array(); % save the array with the new path - can overwrite original array

% For DREADD array (figure 3) do the same:
cno_array_path = [PATH, '\hm3dq_array_sep_zscore_local.mat']; % Path to MouseArray file.
Y = MouseArray(cno_array_path);
Y.reset_path(data_path);
% run to manually reset data path for each object in the array:
for ii = 1:numel(Y.MOUSE_ARRAY)
    Y.MOUSE_ARRAY(ii).PATH = data_path;
end
Y.save_array(); % save the array with the new path - can overwrite original array
%% For subsequent runs, start here:
% Figure 2:
array_path = [PATH, '\full_array_25mice_local.mat']; % Path to MouseArray file.
X = MouseArray(array_path);

% Figure 2C
X.behavior_curves(1);
% Indivudal trials heatmaps and saliency responses:
% All individual ACC trials (Figure 2E-H)
X.aligned_hms(0, 10, 5, 1); % All ACCp trials (will take a while)
% Trial onset and cue responses
X.saliency_mock(10, 1, 1);

% Trial onset aligned by outcome (Fig 2H-I; S4F-G)
X.baseline_by_outcome_sum();
% dmerge = mergedots(ax1, ax2); % Merge any two scatter plots given axes ax1 and ax2
% Pre-trial activity correlations with outcome (Fig 2J-K; S4H-I)
X.baseline_quant(30, 'trials', 't_info', 0); % additionally, plots a lot of individual data. 

%Figure S2E
summary_fig = free_rec_comp(PATH);

% All individual OFC trials (Figure S4A)
X.aligned_hms(0, 10, 5, 0); % All OFC trials (long)
X.saliency_mock(10, 1, 0);

%% Chemogenetic manipulation (Fig 3)
cno_array_path = [PATH, '\hm3dq_array_sep_zscore_local.mat']; % Path to MouseArray file.
Y = MouseArray(cno_array_path);
Y.cno_graphs(-1); % Figure 3D-E; S5D
Y.split_array([], 0, 1,'outcome_hist', 'hm3dq') %Figure 3F
strategy_analysis(data_path, 'dreadds', 0); % Plot prem vs cloud idx - Figure 3G
Y.cno_graphs(0); % Figure S5B
Y.cno_graphs(1); % Figure S5C
Y.split_array([], 0, 1,'hypovig', 'hm3dq') %For figure S5E - needs pruning.
% control behavior = Figure S5F-G
cno_controls(data_path); %control select the following sessions: B1,B2,SAL,CNO for each mouse when prompted.

%% Behavioral Strategies (Figure 4)
strategy_analysis(data_path, 'main', 1); % Plot grouped prem vs cloud idx - Figure 4A, 4C; S6A
X.split_array([], 0, 1,'psych', 'prem_2_groups') % psychometric curves by group Figure 4B
X.split_array([], 1, 0,'RT', 'prem_2_groups') % Figure 4D
X.split_array([], 0, 1,'baseline_summary', 'prem_2_groups_no_ofc') %Figure 4F-G; J-K
% dmerge = mergedots(ax1, ax2); % Merge any two scatter plots given axes ax1 and ax2
X.split_array([], 1, 1,'baseline', 'prem_2_groups') % Figure 4 H-I; L-M (Pre-trial correlations by strategy)
X.split_array([], 0, 1,'hypovig', 'prem_2_groups_no_ofc')  %Figure S6B RT in hits
X.split_array([], 1, 0,'hit_RT', 'prem_2_groups') %Figure S6C RT in hits


%% In progress: 
%Plotting model results

r = RegModel('label', 0); % make sure right target file is used (CNO or regular mice)
% plotting cue response (with visual aid):
r.UC.Properties.VariableNames{5} = 'AudStimVis';
r.UC = movevars(r.UC, 'AudStimVis', 'Before', 'BBN');
r.SVM.Properties.VariableNames{5} = 'AudStimVis';
r.SVM = movevars(r.SVM, 'AudStimVis', 'Before', 'BBN');
r.plot_label_by_group('AudStim', 'area', 1, 0)

% need to create a minimal version of baseline...

%% Functions
function sum_f = free_rec_comp(path)

ofc_mice = load([path '\OFC_free_rec.mat']);
acc_mice = load([path '\ACC_free_rec.mat']);
sum_f = figure;

F = fieldnames(acc_mice);
for ii = 1:size(F, 1)
    tr = zscore(acc_mice.(F{ii}));
    [pks, locs, w, p] = findpeaks((tr), 1000, 'MinPeakProminence', 1,  'MinPeakDist', 2, 'Annotate','extents', 'WidthReference','halfprom');
    pk_freq(ii) = length(pks) / (length(tr)/1000);
    pk_amp(ii) = mean(pks);
    pk_w(ii) = mean(w);
    medf(ii) = mad(tr, 1);
end
h1 = subplot(1,4,1);
h1.XLim = [0.5, 2.5];
h1.XTick = [1;2];
h1.XTickLabel = ['ACC'; 'OFC'];
hold on
plot(h1, ones(1,size(F, 1)), pk_freq, 'go')

line(h1, [0.9,1.1], [mean(pk_freq), mean(pk_freq)], 'Color', 'k')
title(h1, 'Event Rate')
ylabel(h1, 'Hz')
set(h1, 'YLim',[0, 0.4])

h2 = subplot(1,4,2);
h2.XLim = [0.5, 2.5];
h2.XTick = [1;2];
h2.XTickLabel = ['ACC'; 'OFC'];
hold on
plot(h2, ones(1,size(F, 1)), pk_w, 'go')
line(h2, [0.9,1.1], [mean(pk_w), mean(pk_w)], 'Color', 'k')
title(h2, 'Event Width (Half Max)')
ylabel(h2, 'Time (s)')
set(h2, 'YLim',[0, 4])

h3 = subplot(1,4,3);
h3.XLim = [0.5, 2.5];
h3.XTick = [1;2];
h3.XTickLabel = ['ACC'; 'OFC'];
hold on
plot(h3, ones(1,size(F, 1)), pk_amp, 'go')
line(h3, [0.9,1.1], [mean(pk_amp), mean(pk_amp)],'Color', 'k')
title(h3, 'Event Amplitude')
ylabel(h3, 'STD')
set(h3, 'YLim',[0.4, 2])

h4 = subplot(1,4,4);
h4.XLim = [0.5, 2.5];
h4.XTick = [1;2];
h4.XTickLabel = ['ACC'; 'OFC'];
hold on
plot(h4, ones(1,size(F, 1)), medf, 'go')
plot(h4, [0.9,1.1], [mean(medf), mean(medf)], 'Color', 'k')
title(h4,'\DeltaF/F MAD')
ylabel(h4, 'STD')
set(h4, 'YLim', [0.2, 1.2])


F = fieldnames(ofc_mice);
for ii = 1:size(F, 1)
    tr = zscore(ofc_mice.(F{ii}));
    [pks, locs, w, p] = findpeaks((tr), 1000, 'MinPeakProminence', 1,  'MinPeakDist', 2, 'Annotate','extents', 'WidthReference','halfprom');
    opk_freq(ii) = length(pks) / (length(tr)/1000);
    opk_amp(ii) = mean(pks);
    opk_w(ii) = mean(w);
    omedf(ii) = mad(tr, 1);
end

plot(h1, 2*ones(1,size(F, 1)), opk_freq, 'bo')
line(h1, [1.9,2.1], [mean(opk_freq), mean(opk_freq)],'Color', 'k')

plot(h2, 2*ones(1,size(F, 1)), opk_w, 'bo')
line(h2, [1.9,2.1], [mean(opk_w), mean(opk_w)],'Color', 'k')


plot(h3, 2*ones(1,size(F, 1)), opk_amp, 'bo')
line(h3, [1.9,2.1], [mean(opk_amp), mean(opk_amp)],'Color', 'k')


plot(h4, 2*ones(1,size(F, 1)), omedf, 'bo')
line(h4, [1.9,2.1], [mean(omedf), mean(omedf)],'Color', 'k')
end
function [output] = sem (data) 
output=nanstd(data)/sqrt(size(data,1)); %Standard error of the mean
end
function newfigax = mergedots(ax1, ax2)
% used to merge and replot scatters with lines between them.
prem_dots = [];
miss_dots = [];
lines1 = [];
lines2 = [];
dots1 = findall(ax1, 'Type', 'Line');
dots2 = findall(ax2, 'Type', 'Line');
for ii = 1:numel(dots1)
    if sum(contains(fieldnames(dots1(ii)), 'YData')) > 0 && numel(dots1(ii).YData) < 2 && dots1(ii).XData ~=0
        prem_dots = [prem_dots, double(dots1(ii).YData)];
    else
        lines1 = [lines1, dots1(ii)];
    end
    if sum(contains(fieldnames(dots2(ii)), 'YData')) > 0  && numel(dots2(ii).YData) < 2 && dots2(ii).XData ~=0
        miss_dots = [miss_dots, double(dots2(ii).YData)];
    else
        lines2 = [lines2, dots2(ii)]; %hide lines
    end
end

histprem = findall(ax1, 'Type', 'Histogram');
set(histprem, 'FaceColor', 'r')
histmiss = findall(ax2, 'Type', 'Histogram');
newfigax = axes(figure);
plot(newfigax, [-0.1, -0.05], [prem_dots; miss_dots], 'o-k', 'MarkerFaceColor', 'r')
copyobj(histprem, newfigax)
copyobj(histmiss, newfigax)
copyobj(lines1, newfigax)
copyobj(lines2, newfigax)
set(newfigax, 'XLim', [-0.15, 0.15], 'YLim', [-0.2, 0.2], 'XTick', [-0.1, -0.05], 'XTickLabel', {'Imp.', 'Miss'})
newfigax.Parent.Position(3) = 0.5*newfigax.Parent.Position(2);
newfigax.YLabel.String = 'Zscored \DeltaF/F';
end
function cno_controls (path)
paths{1} = [path, '2_from430'];
paths{2} = [path, '3_from430'];
paths{3} = [path, '4_from430'];
paths{4} = [path, '4_fromI2'];
paths{5} = [path, '6_fromI2'];

prem_rate = [];
t_info = [];
for ii = 1:numel(paths)
path = paths{ii};
   
   files = dir(path);
   files(~[files.isdir]) = [];
   folders = {files.name}; 
   folders(1:2) = []; %delete '.' and '..';
   [is, idx] = ismember('passive', folders);
   if is
    folders(idx) = []; %delete 'passive' folder
   end
   
   [idx, tf] = listdlg('PromptString','Select Sessions:',...
                           'SelectionMode','multiple',...
                           'ListString', folders);
   
   select_folders = folders(idx);
   
  for tt = 1:numel(select_folders)  
      X = matfile(fullfile(path, select_folders{tt}, 'Behavior_Data'));
      try
        tbl = X.TaskAnalysis;
      catch
         X = matfile(fullfile(path, select_folders{tt}, 't_info'));
         tbl.Data = X.t_info;
      end
      if height(tbl.Data) > 360
          tbl.Data(361:end, :) = []; %equalize to 360 trials
      end
      prem_rate(ii, tt) = numel(find(tbl.Data.plot_result == -1)) / numel(tbl.Data.delay); % get premature rate
      miss_rate(ii, tt) = numel(find(abs(tbl.Data.plot_result) ~= 1)) / numel(tbl.Data.delay); % get miss rate
      hit_rate(ii, tt) = numel(find(tbl.Data.plot_result == 1)) / numel(tbl.Data.delay); % get hit rate
  end 
end
% plot cumulative:
miss = sum(miss_rate)./(360*5);
hit = sum(hit_rate)./(360*5);
prem = sum(prem_rate)./(360*5);
ha = axes(figure); 
bar(mean(hit_rate(:, 3:4)) + mean(miss_rate(:, 3:4)) + mean(prem_rate(:, 3:4)), 'FaceColor', 'r')
hold on
bar(mean(hit_rate(:, 3:4)) + mean(miss_rate(:, 3:4)), 'FaceColor', [0.5, 0.5, 0.5])
bar(mean(hit_rate(:, 3:4)), 'FaceColor', 'g')
ylabel('Avg Rate');
xlabel('Session Type');
ha.XTickLabel = {'Sal', 'CNO'};

% plot normalized
norm_rate = prem_rate./mean(prem_rate(:,1:2), 2); % normalize to two baseline days
norm_miss_rate = miss_rate./mean(miss_rate(:,1:2), 2); % normalize to two baseline days
GREY = [0.5 0.5 0.5];
ha = axes(figure);
bar(1:2, [mean(norm_rate(:, end)), mean(norm_miss_rate(:, end))])
hold on
errorbar(1, mean(norm_rate(:,end)), sem(norm_rate(:,end)), 'Color', 'k');
errorbar(2, mean(norm_miss_rate(:,end)), sem(norm_miss_rate(:,end)), 'Color', 'k');

hp = plot(ha, ones(1, size(norm_rate, 1)), norm_rate(:, end), 'o', 'Color', GREY, 'MarkerEdgeColor', GREY);
hp1 = plot(ha, 2*ones(1, size(norm_rate, 1)), norm_miss_rate(:, end), 'o', 'Color', GREY, 'MarkerEdgeColor', GREY);
xlabel('Error Type')
ylabel('norm. rate under CNO')

ha.XLim = [0.5, 2.5];
ha.YLim = [0, 2];
ha.XTick = [1 2];
ha.XTickLabel = {'Imp', 'Miss'};

end
function strategy_analysis(path, type, cluster)
switch type
    case 'dreadds'
        path = [path, 'DREADDS\'];
        mice = {'4_from700', '4_from800', '7_from700', '9_from700', '1_from800',...
            '4_from700C', '4_from800C', '7_from700C', '9_from700C', '1_from800C'};
    otherwise
        mice = {'1_from400', '5_from400', '8_from402'...
            '1_from404','4_from404','2_from404', '3_from404', '5_from404', '6_from404', '7_from404'...
            '1_from406', '2_from406', '6_from406', '3_from406'...
            '3_from410', '4_from410', '4_from410L'...
            '2_from430', '3_from430',...
            '1_from500', '2_from500', '3_from500',...
            '1_from440', '2_from440', '3_from440'} ;
end
P = {'cloudM','light_isM', 'AudM', 'prem', 'nc_prem', 'c_prem', 'srate'};
Ptitle = {'Effect of cloud ', 'Effect of light ', 'Effect of cue'};
% Find parameters for each mouse
for i = 1:numel(P)
    Full.(P{i}) = NaN(1,numel(mice));
end

for m = 1:length(mice)
    file = [path, mice{m}, '\CueInCloud_comb_t_onset.mat'];
    mat = matfile(file);
    t_info = mat.t_info;

    All.(['Mouse',mice{m}]).t_info = t_info;
    
    All.(['Mouse',mice{m}]).strat = rig_strat_clean(All.(['Mouse',mice{m}]).t_info);
    fnames = fieldnames(All);
    
    for ii = 1:length(P)
        for i = 1:size(fields(All),1)
            Params.(P{ii})(i) = All.(fnames{i}).strat.(P{ii});
        end
    end
    
end

rMice = mice'; % transposed mice array
% Plot ungrouped mice scatter
sz = 80 ;
h_fig = figure();
s = scatter(axes(h_fig), Params.cloudM, Params.prem, sz, Params.light_isM, 'filled', 'MarkerEdgeColor', [0.55, 0.7, 0.98], 'LineWidth', 2);
row = dataTipTextRow('Mouse',fnames);
s.DataTipTemplate.DataTipRows(end+1) = row;

map = bone ;
colormap(map(40:end-20,:))
colorbar;
caxis([-1, 1]);
ylabel(colorbar, 'Visual modulation index')
xlabel('Cloud modulation index')
ylabel('Premature rate')
title('Attentional Strategy')

if cluster
    %% Devision to 2 groups
exploiters = {'7_from404', '5_from404',  '2_from406', '6_from406', ...
    '3_from406', '1_from406', '2_from500', '2_from440', '5_from400'};

explorers = {'2_from430', '1_from440','3_from500', '3_from410', '1_from500', '3_from430',...
    '1_from404', '1_from400', '3_from404', '6_from404',...
    '2_from404', '4_from410', '4_from410L', '4_from404', '3_from440','8_from402'};

groups = {exploiters,explorers};
%% -- Plot -- %%
sz = 80 ;

h_fig3 = figure() ;
h_ax = subplot(1, 1, 1, 'Parent', h_fig3) ;
hold(h_ax, 'on') ;
clrs = {[0.55, 0.7, 0.98], [0.248, 0.392, 0.54], [0.516, 0.352, 0.54]};
for ii = 1:numel(groups)
s = scatter(h_ax, Params.cloudM(ismember(mice, groups{ii})), Params.prem(ismember(mice, groups{ii})), sz, Params.light_isM(ismember(mice, groups{ii})), 'filled', 'MarkerEdgeColor', clrs{ii}, 'LineWidth', 2);
end
map = bone ;
colormap(map(40:end-20,:))
colorbar;
caxis([-0.5, 0.5]);
ylabel(colorbar, 'Visual modulation index')
xlabel('Cloud modulation index')
ylabel('Premature rate')
title('Attentional Strategy')

h_ax.XTick = -0.1: 0.1 : 0.45 ;
%% Plot impulsive response vs cloud
fmain = figure;
ax = axes(fmain);
hold(ax, 'on');
clrs = {[0.55, 0.7, 0.98], [0.248, 0.392, 0.54], [0.516, 0.352, 0.54]};
for ii = 1:numel(groups)
    plot([(2*ii -1); (2*ii)], [Params.nc_prem(ismember(mice, groups{ii})); Params.c_prem(ismember(mice, groups{ii}))], 'o-', 'Color', clrs{ii});
    plot([(2*ii -1); (2*ii)], [mean(Params.nc_prem(ismember(mice, groups{ii}))); mean(Params.c_prem(ismember(mice, groups{ii})))], 'o-k');
end
 ax.XLim = [0.8, 4.2];
 
 %% Plot success rate
fmain = figure;
ax = axes(fmain);
hold(ax, 'on');
clrs = {[0.55, 0.7, 0.98], [0.248, 0.392, 0.54], [0.516, 0.352, 0.54]};
for ii = 1:numel(groups)
    plot(ii/2, Params.srate(ismember(mice, groups{ii})), 'o', 'Color', clrs{ii});
    plot(ii/2, mean(Params.srate(ismember(mice, groups{ii}))), 'ok');
end
set(ax, 'XLim', [0.3, 1.2], 'YLim', [0, 0.7], 'XTick', [0.5, 1], 'XTickLabel', {'Exploiters', 'Explorers'});
end
%% -- Functions -- %%

function param = rig_strat_clean(t_info)
% param is a structure that contains parameters of the mouse's strategeis,
% regarding dependence on light, effect of cloud and cue attenuation.

% auditory no cloud no visual, ratio of 2:0.02 ---- 0.5
param.Aud_002 = height(t_info((t_info.cue_int == 0.5 | t_info.cue_int == 0.1) & t_info.is_light == 0 & t_info.cloud == 0 & t_info.plot_result == 1, :))/height(t_info((t_info.cue_int == 0.5 | t_info.cue_int == 0.1) & t_info.is_light == 0 & t_info.cloud == 0 & t_info.plot_result ~= -1, :));
param.Aud_2 = height(t_info(t_info.cue_int == 2 & t_info.is_light == 0 & t_info.cloud == 0 & t_info.plot_result == 1, :))/height(t_info(t_info.cue_int == 2 & t_info.is_light == 0 & t_info.cloud == 0 & t_info.plot_result ~= -1, :));
param.AudM = (param.Aud_2 - param.Aud_002)/(param.Aud_2 + param.Aud_002);

% light vs no light in lowest atten ---- 0.5
param.islight = height(t_info(t_info.cloud == 0 & t_info.is_light == 1 & (t_info.cue_int == 0.5 | t_info.cue_int == 0.1) & t_info.plot_result == 1, :))/height(t_info(t_info.cloud == 0 & t_info.is_light == 1 & (t_info.cue_int == 0.5 | t_info.cue_int == 0.1) & t_info.plot_result ~= -1, :));
param.nolight = height(t_info(t_info.cloud == 0 & t_info.is_light == 0 & (t_info.cue_int == 0.5 | t_info.cue_int == 0.1) & t_info.plot_result == 1, :))/height(t_info(t_info.cloud == 0 & t_info.is_light == 0 & (t_info.cue_int == 0.5 | t_info.cue_int == 0.1) & t_info.plot_result ~= -1, :));
param.light_isM = (param.islight - param.nolight)/(param.islight + param.nolight);

% ratio of success rate in cloud no light trials, through all attenuations
param.iscloud = height(t_info(t_info.cloud == 0.5 & t_info.is_light == 0 & t_info.plot_result == 1, :))/height(t_info(t_info.cloud == 0.5 & t_info.is_light == 0 & t_info.plot_result ~= -1, :));
param.nocloud = height(t_info(t_info.cloud == 0 & t_info.is_light == 0 & t_info.plot_result == 1, :))/height(t_info(t_info.cloud == 0 & t_info.is_light == 0 & t_info.plot_result ~= -1, :));
param.cloudM = (param.nocloud - param.iscloud)/(param.nocloud + param.iscloud);


% general rate of premature licks
param.prem = height(t_info(t_info.plot_result == -1, :))/height(t_info);
%no cloud
param.nc_prem = height(t_info(t_info.cloud == 0 & t_info.plot_result == -1, :))/height(t_info(t_info.cloud == 0, :));
%cloud
param.c_prem = height(t_info(t_info.cloud > 0 & t_info.plot_result == -1, :))/height(t_info(t_info.cloud > 0, :));

param.srate = height(t_info(t_info.plot_result == 1, :))/height(t_info);
end
end
        

% The following helper functions were downloaded from MATLAB FileExchange: 
function varargout=shadedErrorBar(x,y,errBar,lineProps,transparent)
% function H=shadedErrorBar(x,y,errBar,lineProps,transparent)
%
% Purpose 
% Makes a 2-d line plot with a pretty shaded error bar made
% using patch. Error bar color is chosen automatically.
%
% Inputs
% x - vector of x values [optional, can be left empty]
% y - vector of y values or a matrix of n observations by m cases
%     where m has length(x);
% errBar - if a vector we draw symmetric errorbars. If it has a size
%          of [2,length(x)] then we draw asymmetric error bars with
%          row 1 being the upper bar and row 2 being the lower bar
%          (with respect to y). ** alternatively ** errBar can be a
%          cellArray of two function handles. The first defines which
%          statistic the line should be and the second defines the
%          error bar.
% lineProps - [optional,'-k' by default] defines the properties of
%             the data line. e.g.:    
%             'or-', or {'-or','markerfacecolor',[1,0.2,0.2]}
% transparent - [optional, 0 by default] if ==1 the shaded error
%               bar is made transparent, which forces the renderer
%               to be openGl. However, if this is saved as .eps the
%               resulting file will contain a raster not a vector
%               image. 
%
% Outputs
% H - a structure of handles to the generated plot objects.     
%
%
% Examples
% y=randn(30,80); x=1:size(y,2);
% shadedErrorBar(x,mean(y,1),std(y),'g');
% shadedErrorBar(x,y,{@median,@std},{'r-o','markerfacecolor','r'});    
% shadedErrorBar([],y,{@median,@std},{'r-o','markerfacecolor','r'});    
%
% Overlay two transparent lines
% y=randn(30,80)*10; x=(1:size(y,2))-40;
% shadedErrorBar(x,y,{@mean,@std},'-r',1); 
% hold on
% y=ones(30,1)*x; y=y+0.06*y.^2+randn(size(y))*10;
% shadedErrorBar(x,y,{@mean,@std},'-b',1); 
% hold off
%
%
% Rob Campbell - November 2009


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Error checking    
error(nargchk(3,5,nargin))


%Process y using function handles if needed to make the error bar
%dynamically
if iscell(errBar) 
    fun1=errBar{1};
    fun2=errBar{2};
    errBar=fun2(y);
    y=fun1(y);
else
    y=y(:)';
end

if isempty(x)
    x=1:length(y);
else
    x=x(:)';
end


%Make upper and lower error bars if only one was specified
if length(errBar)==length(errBar(:))
    errBar=repmat(errBar(:)',2,1);
else
    s=size(errBar);
    f=find(s==2);
    if isempty(f), error('errBar has the wrong size'), end
    if f==2, errBar=errBar'; end
end

if length(x) ~= length(errBar)
    error('length(x) must equal length(errBar)')
end

%Set default options
defaultProps={'-k'};
if nargin<4, lineProps=defaultProps; end
if isempty(lineProps), lineProps=defaultProps; end
if ~iscell(lineProps), lineProps={lineProps}; end

if nargin<5, transparent=0; end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Plot to get the parameters of the line 
H.mainLine=plot(x,y,lineProps{:});


% Work out the color of the shaded region and associated lines
% Using alpha requires the render to be openGL and so you can't
% save a vector image. On the other hand, you need alpha if you're
% overlaying lines. There we have the option of choosing alpha or a
% de-saturated solid colour for the patch surface .

col=get(H.mainLine,'color');
edgeColor=col+(1-col)*0.55;
patchSaturation=0.15; %How de-saturated or transparent to make patch
if transparent
    faceAlpha=patchSaturation;
    patchColor=col;
    set(gcf,'renderer','openGL')
else
    faceAlpha=1;
    patchColor=col+(1-col)*(1-patchSaturation);
    set(gcf,'renderer','painters')
end

    
%Calculate the error bars
uE=y+errBar(1,:);
lE=y-errBar(2,:);


%Add the patch error bar
holdStatus=ishold;
if ~holdStatus, hold on,  end


%Make the patch
yP=[lE,fliplr(uE)];
xP=[x,fliplr(x)];

%remove nans otherwise patch won't work
xP(isnan(yP))=[];
yP(isnan(yP))=[];


H.patch=patch(xP,yP,1,'facecolor',patchColor,...
              'edgecolor','none',...
              'facealpha',faceAlpha);


%Make pretty edges around the patch. 
H.edge(1)=plot(x,lE,'-','color',edgeColor);
H.edge(2)=plot(x,uE,'-','color',edgeColor);

%Now replace the line (this avoids having to bugger about with z coordinates)
delete(H.mainLine)
H.mainLine=plot(x,y,lineProps{:});


if ~holdStatus, hold off, end


if nargout==1
    varargout{1}=H;
end
end
function [ffit, curve] = FitPsycheCurveWH(xAxis, yData, varargin)
% Mean (u): The mean value of the distribution representing subject bias.
% Standard deviation (v): The variation of the distribution representing the subjects discrimination sensitivity.
% Guess rate (g) and lapse rate (l): Two additional parameters representing the subjects fallibility (ie. potential inability to ever reach 100% performance) at each end of the distribution/stimulus spectrum. 
% The parameters g and l constrain the limits of the cumulative distribution that provides the sigmoid shape for the psychometric curve:
% 
% y = g+(1-g-l)*dist where dist is a cumulitive gaussian distribution
% function





% Start points and limits
if  isempty(varargin) || varargin{:} ~=1
    useLims=0;
elseif numel(varargin{:})>1
    useLims=1;
    UL=varargin{1}(1,:);
    SP=varargin{1}(2,:);
    LM=varargin{1}(3,:);
elseif varargin{:} == 1
    useLims=1;
    UL = [0.9500    0.9000   3   10.0000];
    SP = [0.1000    0.6000    2    1.0000];
    LM = [0         0         0         0];
end

% Transpose if necessary
if size(xAxis,1)<size(xAxis,2)
    xAxis=xAxis';
end
if size(yData,1)<size(yData,2)
    yData=yData';
end

% Check range of data
if min(yData)<0 || max(yData)>1  
     % Attempt to normalise data to range 0 to 1
     yData = yData/(mean(yData)*2);
end
    
% Prepare fitting function
F=@(g,l,u,v,x) g+(1-g-l)*0.5*(1+erf((x-u)/sqrt(2*v^2)));

% Fit using fit function from fit toolbox
if useLims==1
    % SPs and limits specified, use while fitting
    ffit=fit(xAxis,yData,F,'StartPoint',SP,'Upper',UL,'Lower',LM);
else
    % Fits not specified, don't use while fitting
    ffit=fit(xAxis,yData,F);
end

% Create a new xAxis with higher resolution
fineX = linspace(min(xAxis),max(xAxis),numel(xAxis)*50);
% Generate curve from fit
curve = feval(ffit, fineX);
curve = [fineX', curve];
end
function x = beeswarm(x,y,varargin)
%function xbee = beeswarm(x,y)
%
% Input arguments:
%   x               column vector of groups (only tested for integer)
%   y               column vector of data
%
% Optional input arguments:
%   sort_style      ('nosort' - default | 'up' | 'down' | 'fan' | 'rand' | 'square' | 'hex')
%   corral_style    ('none' default | 'gutter' | 'omit' | 'rand')
%   dot_size        relative. default=1
%   overlay_style   (false default | 'box' | 'sd' | 'ci')
%   use_current_axes (false default | true)
%   colormap        (lines default | 'jet' | 'parula' | 'r' | Nx3 matrix of RGB values]
%
% Output arguments:
%   xbee            optimized layout positions
%
% Known Issues:
%       x locations depend on figure aspect ratio. resizing the figure window and rerunning may give different results
%       setting corral to 'none' still has a gutter when the width is large
%
% Usage example:
% 	x = round(rand(150,1)*5);
%   y = randn(150,1);
%   beeswarm(x,y,3,'sort_style','up','overlay_style','ci')
%
% % Ian Stevenson, CC-BY 2019

p = inputParser;
addRequired(p,'x')
addRequired(p,'y')
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addOptional(p,'sort_style','nosort')
addOptional(p,'corral_style','none')
addOptional(p,'dot_size',11/sqrt(length(x)),validScalarPosNum)
addOptional(p,'overlay_style',false)
addOptional(p,'use_current_axes',false)
addOptional(p,'colormap','lines')
addOptional(p,'MarkerFaceColor','')
addOptional(p,'MarkerFaceAlpha',.5)
addOptional(p,'MarkerEdgeColor','none')
parse(p,x,y,varargin{:});

% extra parameters
rwid = .05; % width of overlay box/dash

dcut=8; % spacing factor
nxloc=512; % resolution for optimization
chanwid = .9; % percent width of channel to use
yl = [min(y) max(y)]; % default y-limits
asp_rat = 1;
keep_hold = false;

% get aspect ratio for a figure window
if isfinite(p.Results.dot_size)
    if ~p.Results.use_current_axes
        % make new axes
        s=scatter(x,y);
        xl=[min(x)-.5 max(x)+.5];
    else
        xl=xlim();
    end
    yl=ylim();
    pasp_rat = get(gca,'PlotBoxAspectRatio');
    dasp_rat = get(gca,'DataAspectRatio');
    asp_rat = pasp_rat(1)/pasp_rat(2);
    
    % pix-scale
    pf = get(gcf,'Position');
    pa = get(gca,'Position');
    as = pf(3:4).*pa(3:4); % width and height of panel in pixels
    dcut = dcut*sqrt(p.Results.dot_size)/as(1)*(range(unique(x))+1);
    if ~ishold
        cla
    else
        keep_hold = true;
    end
end

% sort/round y for different plot styles
yorig=y;
switch lower(p.Results.sort_style)
    case 'up'
        [y,sid]=sort(y);
    case 'fan'
        [~,sid]=sort(abs(y-mean(y)));
        sid=[sid(1:2:end); sid(2:2:end)];
        y=y(sid);
    case 'down'
        [y,sid]=sort(y,'descend');
    case 'rand'
        sid=randperm(length(y));
        y=y(sid);
    case 'square'
        nxloc=.9/dcut;
%         [~,e,b]=histcounts(y,ceil((range(x)+1)*chanwid*nxloc/2/asp_rat));
        edges = linspace(min(yl),max(yl),ceil((range(x)+1)*chanwid*nxloc/asp_rat));
        [~,e,b]=histcounts(y,edges);
        y=e(b)'+mean(diff(e))/2;
        [y,sid]=sort(y);
    case 'hex'
        nxloc=.9/dcut;
%         [~,e,b]=histcounts(y,ceil((range(x)+1)*chanwid*nxloc/2/sqrt(1-.5.^2)/asp_rat));
        edges = linspace(min(yl),max(yl),ceil((range(x)+1)*chanwid*nxloc/sqrt(1-.5.^2)/asp_rat));
        [n,e,b]=histcounts(y,edges);
        oddmaj=0;
        if sum(mod(n(1:2:end),2)==1)>sum(mod(n(2:2:end),2)==1),
            oddmaj=1;
        end
        y=e(b)'+mean(diff(e))/2;
        [y,sid]=sort(y);
        b=b(sid);
    otherwise
        sid=1:length(y);
end
x=x(sid);
yorig=yorig(sid);
[ux,~,ic] = unique(x);
% rmult=(range(ux)+1)*2;
rmult=5;

% for each group...
for i=1:length(ux)
    fid = find(ic==i);   
    
    % set of possible x locations
    xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i);

    % rescale y to that things are square visually
    zy=(y(fid)-min(yl))/(max(yl)-min(yl))/asp_rat*(range(ux)+1)*chanwid;
    
    % precalculate y distances so that we only worry about nearby points
    D0=squareform(pdist(zy))<dcut*2;    
    
    if length(fid)>1
        % for each data point in the group sequentially...
        for j=1:length(fid)
            if strcmp(lower(p.Results.sort_style),'hex')
                xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i);
                if mod(b(fid(j)),2)==oddmaj
                    xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i)+mean(diff(xi))/2;
                end
            end
            zid = D0(j,1:j-1);
            e = (xi-ux(i)).^2; % cost function
            if ~strcmp(lower(p.Results.sort_style),'hex') && ~strcmp(lower(p.Results.sort_style),'square')
                if sum(zid)>0
                    D = pdist2([xi ones(length(xi),1)*zy(j)], [x(fid(zid)) zy(zid)]);
                    D(D<=dcut)=Inf;
                    D(D>dcut & isfinite(D))=0;
                    e = e + sum(D,2) + randn(1)*10e-6; % noise to tie-break
                end
            else
                if sum(zid)>0
                    D = pdist2([xi ones(length(xi),1)*zy(j)], [x(fid(zid)) zy(zid)]);
                    D(D==0)=Inf;
                    D(D>dcut & isfinite(D))=0;
                    e = e + sum(D,2) + randn(1)*10e-6; % noise to tie-break
                end
            end

            if strcmp(lower(p.Results.sort_style),'one')
                e(xi<ux(i))=Inf;
            end
            [~,mini] = min(e);
            if mini==1 && rand(1)>.5, mini=length(xi); end
            x(fid(j)) = xi(mini);
        end
    end
%     x(fid)=x(fid)-median(x(fid))+ux(i); % center x locations by median
end

if strcmp(lower(p.Results.sort_style),'randn')
    x=ux(ic)+randn(size(ic))/4;
end

% corral any points outside of the channel
out_of_range = abs(x-ux(ic))>chanwid/2;
switch lower(p.Results.corral_style)
    case 'gutter'
        id = (x-ux(ic))>chanwid/2;
        x(id)=chanwid/2+ux(ic(id));
        id = (x-ux(ic))<-chanwid/2;
        x(id)=-chanwid/2+ux(ic(id));
    case 'omit'
        x(out_of_range)=NaN;
    case 'random'
        x(out_of_range)=ux(ic(out_of_range))+rand(sum(out_of_range),1)*chanwid-chanwid/2;
end

% plot groups and add overlay
if isfinite(p.Results.dot_size)
    if isnumeric(p.Results.colormap)
        cmap=p.Results.colormap;
    else
        cmap = feval(p.Results.colormap,length(ux));
    end
    for i=1:length(ux)
        if isempty(p.Results.MarkerFaceColor')
            scatter(x(ic==i),y(ic==i),p.Results.dot_size*36,'filled','MarkerFaceAlpha',p.Results.MarkerFaceAlpha,'MarkerEdgeColor',p.Results.MarkerEdgeColor,'MarkerFaceColor',cmap(i,:))
        else
            scatter(x(ic==i),y(ic==i),p.Results.dot_size*36,'filled','MarkerFaceAlpha',p.Results.MarkerFaceAlpha,'MarkerEdgeColor',p.Results.MarkerEdgeColor,'MarkerFaceColor',p.Results.MarkerFaceColor)
        end
        hold on
        iqr = prctile(yorig(ic==i),[25 75]);
        switch lower(p.Results.overlay_style)
            case 'box'
                rectangle('Position',[ux(i)-rwid iqr(1) 2*rwid iqr(2)-iqr(1)],'EdgeColor','k','LineWidth',2)
                line([ux(i)-rwid ux(i)+rwid],[1 1]*median(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
            case 'sd'
                line([1 1]*ux(i),mean(yorig(ic==i))+[-1 1]*std(yorig(ic==i)),'Color',cmap(i,:),'LineWidth',2)
                line([ux(i)-2*rwid ux(i)+2*rwid],[1 1]*mean(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
            case 'ci'
                line([1 1]*ux(i),mean(yorig(ic==i))+[-1 1]*std(yorig(ic==i))/sqrt(sum(ic==i))*tinv(0.975,sum(ic==i)-1),'Color',cmap(i,:),'LineWidth',2)
                line([ux(i)-2*rwid ux(i)+2*rwid],[1 1]*mean(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
        end
        
    end
    hold off
    if keep_hold
        hold on
    end
    xlim(xl)
    ylim(yl)
end

% unsort so that output matches the original y data
x(sid)=x;
end