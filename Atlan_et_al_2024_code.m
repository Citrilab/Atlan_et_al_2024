% Atlan et al 2024 Code:
% Reproduces figures for Atlan et al 2024, Nat Comms

% Requires the following matlab classes:
% MouseArray (holds all mice together)
% MouseSummary (holds individual mice handles)
% RegModel (opens regression model results)
%% Start here always:
PATH = 'C:\Users\galat\Documents\MATLAB\atlanetal2024\'; %Edit this according to the path of the code
data_path = 'C:\Users\galat\Desktop\New folder\'; %Edit this according to the path where the data is located
cno_data_path = 'C:\Users\galat\Desktop\New folder\DREADDs\'; %Edit this according to the path where the data is located
sleep_path = [data_path, '\sleep_data'];
%% Run during the first time only - reset path to match your folder locations:
array_path = [PATH, '\full_array_25mice_local.mat']; % Path to MouseArray file.
X = MouseArray(array_path); % load array;
X.reset_path(data_path);
X.get_tags(0.5);
% run to manually reset data path for each object in the array:
for ii = 1:numel(X.MOUSE_ARRAY)
    X.MOUSE_ARRAY(ii).PATH = data_path;
end
X.save_array(); % save the array with the new path - can overwrite original array

% For DREADD array (figure 3) do the same:
cno_array_path = [PATH, '\hm3dq_array_sep_zscore_local.mat']; % Path to MouseArray file.
Y = MouseArray(cno_array_path);
Y.reset_path(cno_data_path);
Y.get_tags(0.5);
% run to manually reset data path for each object in the array:
for ii = 1:numel(Y.MOUSE_ARRAY)
    Y.MOUSE_ARRAY(ii).PATH = data_path;
end
Y.save_array(); % save the array with the new path - can overwrite original array

%% For subsequent runs, start here:
%% Figure 1 - Sleep
claustrum_state_activity(sleep_path); %Figure 1C, S1E
SEA_traces(sleep_path); %Figure 1D, S1F - unzip SEA_zscoredDff_traces in sleep data folder first
bootstrap_sleep(sleep_path, 'SEA_ACCp_zscored'); %Figure 1E
bootstrap_sleep(sleep_path, 'percentNREM_hm3dq'); %Figure 1I
bootstrap_sleep(sleep_path, 'percentAwakening_hm3dq'); %Figure 1L
bootstrap_sleep(sleep_path, 'SEA_OFCp_zscored'); %Figure S1G
bootstrap_sleep(sleep_path, 'percentNREM_ctrl'); %Figure S1I
bootstrap_sleep(sleep_path, 'percentAwakening_ctrl'); %Figure S1J
%% Figure 2:
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

summary_fig = free_rec_comp(PATH); %Figure S2E

% All individual OFC trials (Figure S4A)
X.aligned_hms(0, 10, 5, 0); % All OFC trials (long)
X.saliency_mock(10, 1, 0);

%% Chemogenetic manipulation (Fig 3)
cno_array_path = [PATH, '\hm3dq_array_sep_zscore_local.mat']; % Path to MouseArray file.
Y = MouseArray(cno_array_path);
Y.cno_graphs(-1); % Figure 3D-E; S5D
Y.split_array([], 0, 1,'outcome_hist', 'hm3dq') %Figure 3F
strategy_analysis(cno_data_path, 'dreadds', 0); % Plot prem vs cloud idx - Figure 3G
Y.cno_graphs(0); % Figure S5B
Y.cno_graphs(1); % Figure S5C
Y.split_array([], 0, 1,'hypovig', 'hm3dq') %For figure S5E
%% control behavior = Figure S5F-G
cno_controls(cno_data_path); %control select the following sessions: B1,B2,SAL,CNO for each mouse when prompted.

%% Behavioral Strategies (Figure 4)
strategy_analysis(data_path, 'main', 1); % Plot grouped prem vs cloud idx - Figure 4A, 4C; S6A
X.split_array([], 0, 1,'psych', 'prem_2_groups') % psychometric curves by group Figure 4B
X.split_array([], 1, 0,'RT', 'prem_2_groups') % Figure 4D
model_plot(PATH, 'main', 'BBN', 'prem_2_groups'); %Plotting model results for BBN response (4E)
X.split_array([], 0, 1,'baseline_summary', 'prem_2_groups_no_ofc') %Figure 4F-G; J-K
% dmerge = mergedots(ax1, ax2); % Merge any two scatter plots given axes ax1 and ax2
X.split_array([], 1, 1,'baseline', 'prem_2_groups') % Figure 4 H-I; L-M (Pre-trial correlations by strategy)
X.split_array([], 0, 1,'hypovig', 'prem_2_groups_no_ofc')  %Figure S6B miss streak comparison
X.split_array([], 1, 0,'hit_RT', 'prem_2_groups') %Figure S6C RT in hits

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
path(end-7:end) = []; %remove dreadds from path
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
ylabel('Sucess Rate')
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
function model_plot(path, type, param, group)
% type - CNO or main cohort
% param - what to plot (BBN, Aud, etc...)
% group: area, prem_2_groups, cno
r = RegModel(path, 'label', 0, type); % make sure right target file is used (CNO or regular mice)
% Some prep
r.UC.Properties.VariableNames{5} = 'AudStimVis';
r.UC = movevars(r.UC, 'AudStimVis', 'Before', 'BBN');
r.SVM.Properties.VariableNames{5} = 'AudStimVis';
r.SVM = movevars(r.SVM, 'AudStimVis', 'Before', 'BBN');
switch param
    case 'BBN'
         r.plot_by_group('BBN', group, 1, 0);
    case 'Aud'
        r.plot_label_by_group('AudStim', group, 1, 0);
end
end