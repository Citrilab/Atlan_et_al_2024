classdef MouseArray < matlab.mixin.Copyable % allows the handle to be shallow copied (for taking subsets of it)

    properties
        param_data = struct();
        plot_handles = struct();
        MOUSE_ARRAY
        PATH
    end
    
    methods
        function obj = MouseArray(path)
            obj.PATH = 'C:\Users\galat\Documents\MATLAB\atlanetal2024'; %default path if not specficied
            if nargin > 0
                tmp = load(path);
                obj.PATH = fileparts(path);
                obj.MOUSE_ARRAY = tmp.mse_array;
                obj.param_data.zscore = 0;
            else
                load_flg = questdlg('Load Existing Array?', 'Loard Array or Create New', 'Yes', 'No', 'Yes');
                
                switch load_flg
                    case 'Yes'
                        [filename, foldername] = uigetfile(obj.PATH, ...
                            'Select StimArray file', ...
                            'MultiSelect', 'off');
                        if isequal(filename, 0) || isequal(foldername, 0)
                            disp('No file selected');                                          % Case where no file selected
                            return
                        end
                        fullFileName = fullfile(foldername, filename);
                        tmp = load(fullFileName);
                        obj.MOUSE_ARRAY = tmp.mse_array;
                        
                    case 'No'
                        obj.MOUSE_ARRAY = obj.load_array();
                        load_passive = questdlg('Do you wish to load passive responses?', 'Array Loaded', 'Yes', 'No', 'No');
                        switch load_passive
                            case 'Yes'
                                obj.load_passive_stim();
                        end
                        
                        now_what = questdlg('Do you wish to save?', 'Array Loaded', 'Yes', 'Calculate Tags', 'Yes');
                        switch now_what
                            case 'Yes'
                                obj.save_array();
                            case 'No'
                            case 'Calculate Tags'
                                sd = 0.5;
                                obj.get_tags(sd);
                                obj.save_array();
                        end
                end
                sz = questdlg('Zscore when plotting?', 'Choose yes if data is not already zscored', 'Yes', 'No', 'Yes'); % add parameter for zscoring
                switch sz
                    case 'Yes'
                        obj.param_data.zscore = 1;
                    case 'No'
                        obj.param_data.zscore = 0;
                end
            end
        end
        
        function get_tags(obj, sd_thresh)
            if nargin < 2
                sd_thresh = 0.5; % default sets 30/30/40% split
            end
            for ii = 1:size(obj.MOUSE_ARRAY, 2)
                obj.MOUSE_ARRAY(ii).get_state_tag_by_baseline(sd_thresh);
                
                %
                info  = obj.MOUSE_ARRAY(ii).INFO.t_info;                                      % set any field called attenuation (old mice) to opto (new name)
                [flg, loc]  = ismember('attenuation', info.Properties.VariableNames);
                if flg
                    obj.MOUSE_ARRAY(ii).INFO.t_info.Properties.VariableNames(loc) = {'opto'};
                    obj.MOUSE_ARRAY(ii).INFO.l_info.Properties.VariableNames(loc) = {'opto'};
                    obj.MOUSE_ARRAY(ii).INFO.cue_info.Properties.VariableNames(loc) = {'opto'};
                end
            end
        end
        
        function reset_path(obj, path)
            %used to change the path for all mice in the array (to their data) 
            for ii = 1:length(obj.MOUSE_ARRAY)
                obj.MOUSE_ARRAY(ii).PATH = path;
                obj.MOUSE_ARRAY(ii).get_data;
            end
        end
               
        function load_passive_stim(obj)
            % function LOAD_PASSIVE_STIM loads the passive stimuli matrix
            % if it exists
            for ii = 1:length(obj.MOUSE_ARRAY)
                array = obj.MOUSE_ARRAY(ii);
                array.get_passive_data();
            end
        end
        
        function info = combine_info(obj, nme)
            % combines an info matrix from all mice
            info = [];
            for ii = 1:length(obj.MOUSE_ARRAY)
                tmp = obj.MOUSE_ARRAY(ii).INFO.(nme);
                if strcmp(tmp.Properties.VariableNames(end), 'activity_tag_ch1')
                    tmp(:, end + 1) = {nan};
                    tmp.Properties.VariableNames(end) = {'activity_tag_ch2'};
                end
                tmp(:, end+1) = {obj.MOUSE_ARRAY(ii).MOUSE_NM};
                tmp.Properties.VariableNames(end) = {'mse_nm'};
                info = [info; tmp];
            end
        end
        
        function [p, stats] = run_analysis(obj, epoch, param, plot_individual, cond)
            if nargin < 5
                h1 = obj.compare_param(epoch, plot_individual);
            else
                h1 = obj.compare_param(epoch, plot_individual, cond);
            end
            sz = size(h1);
            root_find = get(h1, 'Type');
            h1(strcmp(root_find, 'root')) = [];                         % remove empty rows (from single channel)
            h1 = reshape(h1, [], sz(2));                                % organize back
            h2 = obj.merge_ind_plots(h1);                           % combine plots from the same mouse;
            obj.plot_handles.(epoch).individual_plots.(param) = h2;
            [obj.plot_handles.(epoch).summary_plots.(param), a_sum, o_sum] = obj.plot_summary(epoch, param, 1, 'ACC', 'OFC');
            if plot_individual
                obj.set_plot_scales(obj.plot_handles.(epoch).individual_plots.(param));     % scales all individual plots
            end
            [p, tbl, stats] = obj.do_statistics(epoch, a_sum, o_sum);
            
        end
        
        function [h_ind] = compare_param(obj, epoch, plot_flg, cond)
            array = obj.MOUSE_ARRAY;
            param_struct = obj.param_data;
            
            if ~plot_flg
                set(0, 'DefaultFigureVisible','off');                           % block plot visibilty
            end
            for ii = 1:length(array)
                data = array(ii).normalize_data(epoch, 1);
                switch epoch
                    case 'passive'
                        t_vec = linspace(-1,4, length(data));
                    otherwise
                        t_vec = linspace(-5,15, length(data));
                end
                [nme, num_run, data_param, t_win] = obj.get_plot_settings(epoch);
                if nargin < 4
                    [info, data_idx, params] = obj.get_trial_subset(array(ii).INFO, epoch);
                else
                    [info, data_idx, params] = obj.get_trial_subset(array(ii).INFO, epoch, cond);
                end
                data = data(data_idx, :); % can use this to create sorted heatmap
                for jj = 1:num_run
                    sub_data = data(info.(data_param) == params(jj), :);
                    [ha_ind(ii, jj), q] = array(ii).analyze_data_segment(t_vec, sub_data, param_struct.zscore, t_win, 100);  % change zscore to 0 to stop zscoring
                    sgtitle([strrep(array(ii).MOUSE_NM, '_', ' '), ' ', array(ii).PROJ_TO{1}, ' ', strrep(nme, '_', ' '), ' ', num2str(jj)]);
                    if ~strcmp(epoch, 'passive')
                         param_struct.(array(ii).PROJ_TO{1}).(['m_', array(ii).MOUSE_NM]).(nme).([data_param, num2str(jj)]) = q;
                        param_struct.(array(ii).PROJ_TO{1}).(['m_', array(ii).MOUSE_NM]).(nme).([data_param, num2str(jj)]).auc = ...
                                                                 obj.compare_response_auc(downsample(sub_data', 100)', [1:50], [51:65]);
                    else
                        param_struct.(array(ii).PROJ_TO{1}).(['m_', array(ii).MOUSE_NM]).(nme).([data_param, num2str(jj)]) = q;
                        param_struct.(array(ii).PROJ_TO{1}).(['m_', array(ii).MOUSE_NM]).(nme).([data_param, num2str(jj)]).auc = ...
                                                                 obj.compare_response_auc(downsample(sub_data', 100)', [1:10], [11:30]);
                    end
                end
                
                if ~isempty(array(ii).PROJ_TO{2})
                    data = array(ii).normalize_data(epoch, 2);
                    data = data(data_idx, :);
                    for jj = 1:num_run
                        sub_data = data(info.(data_param) == params(jj), :);
                        [ho_ind(ii, jj), q] = array(ii).analyze_data_segment(t_vec, sub_data, param_struct.zscore, t_win, 100); % change zscore to 0 to stop zscoring
                        title([strrep(array(ii).MOUSE_NM, '_', ' '), ' ', array(ii).PROJ_TO{2}, ' ', strrep(nme, '_', ' '), ' ', num2str(jj)]);
                        if ~strcmp(epoch, 'passive')
                             param_struct.(array(ii).PROJ_TO{2}).(['m_', array(ii).MOUSE_NM]).(nme).([data_param, num2str(jj)]) = q;
                            param_struct.(array(ii).PROJ_TO{2}).(['m_', array(ii).MOUSE_NM]).(nme).([data_param, num2str(jj)]).auc = ...
                                                                 obj.compare_response_auc(downsample(sub_data', 100)', [1:50], [51:65]);
                        else
                            param_struct.(array(ii).PROJ_TO{2}).(['m_', array(ii).MOUSE_NM]).(nme).([data_param, num2str(jj)]) = q;
                            param_struct.(array(ii).PROJ_TO{2}).(['m_', array(ii).MOUSE_NM]).(nme).([data_param, num2str(jj)]).auc = ...
                                obj.compare_response_auc(downsample(sub_data', 100)', [1:10], [11:30]);
                        end
                    end
                else
                      ho_ind = [];
                end
                h_ind = [ha_ind; ho_ind];
            end
            obj.param_data = param_struct;
            set(0, 'DefaultFigureVisible','on');                                                        % back to enable plotting visibility
        end
        
        function [h_sum, a_pks, o_pks] = plot_summary(obj, epoch, param, plot_flg, target1, target2)
            param_struct = obj.param_data;
            [nme, num_run, data_param, ~, plot_type] = obj.get_plot_settings(epoch);
            h = axes(figure);
            h_sum = h.Parent;
            hold (h, 'on')
            for jj = 1:num_run
                fnames = fieldnames(param_struct.(target1));
                for ii = 1:length(fnames)
                    if ~isfield(param_struct.(target1).(fnames{ii}).(nme).([data_param, num2str(jj)]),(param)) || isempty(param_struct.(target1).(fnames{ii}).(nme).([data_param, num2str(jj)]).(param)) %|| param_struct.(target1).(fnames{ii}).(nme).([data_param, num2str(jj)]).(param) == 0
                        a_pks(ii, jj) = 0;
                    else
                        a_pks(ii, jj) = param_struct.(target1).(fnames{ii}).(nme).([data_param, num2str(jj)]).(param);
                    end
                end
                a_table = table(a_pks, fnames);
                fnames = fieldnames(param_struct.(target2));
                for ii = 1:length(fnames)
                    if ~isfield(param_struct.(target2).(fnames{ii}).(nme).([data_param, num2str(jj)]),(param)) || isempty(param_struct.(target2).(fnames{ii}).(nme).([data_param, num2str(jj)]).(param)) %|| param_struct.(target2).(fnames{ii}).(nme).([data_param, num2str(jj)]).(param) == 0
                        o_pks(ii, jj) = 0;
                    else
                        o_pks(ii, jj) = param_struct.(target2).(fnames{ii}).(nme).([data_param, num2str(jj)]).(param);
                    end
                end
                o_table = table(o_pks, fnames);
            end
            if plot_flg
                switch plot_type
                    case 'stand alone'
                        for tt = 1:numel(o_pks)
                            tmp_h = plot(h, 1, o_pks(tt), 'ob');
                            tmp_h.Tag = o_table.fnames{tt};
                        end
                        for tt = 1:numel(a_pks)
                            tmp_h = plot(h, 2, a_pks(tt), 'ok');
                            tmp_h.Tag = a_table.fnames{tt};
                        end
                        set(h, 'XLim', [0.5, 2.5], 'XTick', ([1,2]), 'XTickLabel', {target2, target1});
                    case 'vs'
                        for tt = 1:size(o_pks, 1)
                            tmp_h = plot(h, [1:num_run(end)], o_pks(tt, :), 'ob');
                            tmp_h.Tag = o_table.fnames{tt};
                            if tt ~= 1
                                tmp_h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                            end
                        end
                        for tt = 1:size(a_pks, 1)
                            tmp_h = plot(h, [1:num_run(end)], a_pks(tt, :), 'ok');
                            tmp_h.Tag = a_table.fnames{tt};
                            if tt ~= 1
                                tmp_h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                            end
                        end
                        legend(target2, target1)
                        set(h, 'XLim', [0.5, num_run(end) + 0.5], 'XTick', ([1:num_run(end)]));
                        if strcmp(epoch, 'lick')
                            h.XTickLabel = {'Unrewarded', 'Rewarded'};
                            set(h.Children, 'LineStyle' ,'-');
                        elseif strcmp(epoch, 'vis')
                            h.XTickLabel = {'Aud', 'Vis + Aud'};
                            set(h.Children, 'LineStyle' ,'-');
                        elseif strcmp(epoch, 'passive')
                            h.XTickLabel = {'Pre', 'Post'};
                            set(h.Children, 'LineStyle' ,'-');
                        end
                    case 'cross-mice'
                        h1 = subplot(2,1,1);
                        hold(h1, 'on')
                        clr = {'r', 'g', 'k'};
                        for tt = 1:size(a_pks, 1)
                            for kk = 1:3
                                tmp_h = scatter(h1, kk, a_pks(tt, kk), clr{kk});
                                tmp_h.Tag = a_table.fnames{tt};
                            end
                        end
                        bar(h1, 1, mean(a_pks(:, 1)), 'FaceColor', [1, 0, 0], 'FaceAlpha', 0.5)
                        bar(h1, 2, mean(a_pks(:, 2)), 'FaceColor', [0, 1, 0], 'FaceAlpha', 0.5)
                        bar(h1, 3, mean(a_pks(:, 3)), 'FaceColor', [1, 1, 1], 'FaceAlpha', 0.5)
                        h1.XTick = 1:3;
                        h1.XTickLabel = {'Premature', 'Hit', 'Miss'};
                        set(h1.Title, 'String', target1)
                        
                        h2 = subplot(2,1,2);
                        hold(h2, 'on')
                        for tt = 1:size(o_pks, 1)
                            for kk = 1:3
                                tmp_h = scatter(h2, kk, o_pks(tt, kk), clr{kk});
                                tmp_h.Tag = o_table.fnames{tt};
                            end
                        end
                        bar(h2, 1, mean(o_pks(:, 1)), 'FaceColor', [1, 0, 0], 'FaceAlpha', 0.5)
                        bar(h2, 2, mean(o_pks(:, 2)), 'FaceColor', [0, 1, 0], 'FaceAlpha', 0.5)
                        bar(h2, 3, mean(o_pks(:, 3)), 'FaceColor', [1, 1, 1], 'FaceAlpha', 0.5)
                        h2.XTick = 1:3;
                        h2.XTickLabel = {'Premature', 'Hit', 'Miss'};
                        
                        set([h1,h2], 'YLim', [min(h1.YLim(1), h2.YLim(1)), max(h1.YLim(2), h2.YLim(2))])
                        set(h2.XLabel, 'String', param)
                        set(h2.Title, 'String', target2)
                end
                sgtitle( [epoch, ' quantification']);
                ax_h = findall(h_sum, 'Type', 'Axes');
                
                if strcmp(epoch, 'baseline')
                    set([ax_h.YLabel], 'String', '\DeltaRate');
                    set(h_sum, 'Position',  [680   382   344   716])
                else
                    set([ax_h.YLabel], 'String', param);
                end
                
                dcm_obj = datacursormode(h_sum);
                set(dcm_obj,'UpdateFcn',@myupdatefcn)
            else
                close(h.Parent); % if not plotting, close figure
            end
            function txt = myupdatefcn(~,event_obj)
                % Customizes text of data tips
                try
                    pos = get(event_obj,'Position');
                    plot_obj = get(event_obj, 'Target');
                    obj_name = strrep(plot_obj.Tag(3:end), '_', ' ');
                    txt = {num2str(pos(2)), obj_name};              % add to the highlighted tag the name of the mouse
                catch
                    keyboard();                                 % in case something goes wrong
                end
            end
        end
        
        function [r, p, del_ch, sumfig_h] = baseline_quant(obj, binnum, align, inf_name, cond)
            % binnum  - number of bins
            % epoc to align to (trials etc)
            % inf_name (corresponding info name)
            % cond - optional intersectioning     
            array = obj.MOUSE_ARRAY;
            if isempty(binnum)
                binnum = 20;    % default if not defined
            end
            for ii = 1:length(array)
                info = array(ii).INFO.(inf_name);
                if nargin > 4 && cond ~= 0
                    idx = eval(cond);
                else
                    idx = true(1, size(array(ii).INFO.(inf_name), 1));
                end
                session_indices = array(ii).get_session_subset(340, 't_info');  %get 300 indices per session
                
                sm_win = 10; % how much to smooth baseline
                shuff = 1000; % how many times to shuffle
                cut_win = 20; % cut out first 20 trials (with drift)
                
                for kk = 1:array(ii).DATA_CH
                    if kk > 1
                        baseline1 = baseline;
                    end
                    baseline = array(ii).data.(['ch', num2str(kk),'_baseline']);
                    data = array(ii).normalize_data(align, kk);
                    [r(kk, ii), p(kk, ii), top{kk, ii}, bot{kk, ii}, del_ch(ii).(['ch', num2str(kk)]), probs] = obj.base_outcome_corr(baseline(idx), info(idx,:), binnum, data(idx, :));

                    obj.param_data.(array(ii).PROJ_TO{kk}).(['m_', array(ii).MOUSE_NM]).baseline.outcome1.low = del_ch(ii).(['ch', num2str(kk)])(1, 1); % save paramters (low, high; prem, cor, om)
                    obj.param_data.(array(ii).PROJ_TO{kk}).(['m_', array(ii).MOUSE_NM]).baseline.outcome2.low = del_ch(ii).(['ch', num2str(kk)])(2, 1);
                    obj.param_data.(array(ii).PROJ_TO{kk}).(['m_', array(ii).MOUSE_NM]).baseline.outcome3.low = del_ch(ii).(['ch', num2str(kk)])(3, 1);
                    obj.param_data.(array(ii).PROJ_TO{kk}).(['m_', array(ii).MOUSE_NM]).baseline_probs = probs;

                    obj.param_data.(array(ii).PROJ_TO{kk}).(['m_', array(ii).MOUSE_NM]).baseline.outcome1.high = del_ch(ii).(['ch', num2str(kk)])(1, 2);
                    obj.param_data.(array(ii).PROJ_TO{kk}).(['m_', array(ii).MOUSE_NM]).baseline.outcome2.high = del_ch(ii).(['ch', num2str(kk)])(2, 2);
                    obj.param_data.(array(ii).PROJ_TO{kk}).(['m_', array(ii).MOUSE_NM]).baseline.outcome3.high = del_ch(ii).(['ch', num2str(kk)])(3, 2);
                    
                    top_idx = find(top{kk, ii});
                    bot_idx = find(bot{kk, ii});
                    top_idx(top_idx == 1) = []; % remove trial #1 because no previous trial
                    bot_idx(bot_idx == 1) = [];
                    
                    info.plot_result(info.plot_result == 2) = 0; % lates are misses
                    out_before_bot = info.plot_result(bot_idx - 1);
                    out_before_top = info.plot_result(top_idx - 1);
                    obj.param_data.(array(ii).PROJ_TO{kk}).(['m_', array(ii).MOUSE_NM]).pre_bot = hist(out_before_bot,unique(out_before_bot))/length(out_before_bot);
                    obj.param_data.(array(ii).PROJ_TO{kk}).(['m_', array(ii).MOUSE_NM]).pre_top = hist(out_before_top,unique(out_before_top))/length(out_before_top);
                end
                
            end
            %plot summary graph for all mice
            obj.plot_baseline_summary('avg');
        end
        

        function hf = analyze_rt(obj, comb, cond1, cond2)
            % compares RT across 2 conditions
            % (e.g. 'info.cue_int > 0.5 & info.is_light == 1', 'info.cue_int < 1 & info.is_light == 1')
            % comb flag indicates whether to analyze within mice or across
            if comb == 1
                for ii = 1:length(obj.MOUSE_ARRAY)
                    array = obj.MOUSE_ARRAY(ii);
                    info = array.INFO.t_info;
                    if strcmp(cond1, 'info.plot_result == -1')
                        hf(ii) = figure;
                        sgtitle(array.MOUSE_NM, 'Interpreter', 'none')
                        subplot(1, 2, 1)
                        hold on
                        
                        if nargin < 4
                            RT1 = info.first_lick(info.plot_result == -1);
                            plot(RT1, info.delay(info.plot_result == -1), 'r.')
                            xlabel('RT')
                            ylabel('Delay')
                            subplot(1, 2, 2)
                            hold on
                            histogram(RT1, 20, 'FaceColor', 'r');
                            histogram(info.delay(info.plot_result == -1), 20, 'FaceColor', 'm');
                            legend('Prem. RT', 'Delay')
                            xlabel('RT')
                        else
                            RT1 = info.first_lick(info.plot_result == -1 & eval(cond2));
                            plot(RT1, info.delay(info.plot_result == -1 & eval(cond2)), 'r.')
                            
                            RT2 = info.first_lick(info.plot_result == -1 & ~eval(cond2));
                            plot(RT2, info.delay(info.plot_result == -1 & ~eval(cond2)), 'm.')
                            xlabel('RT')
                            ylabel('Delay')
                            subplot(1, 2, 2)
                            hold on
                            histogram(RT1, 20, 'FaceColor', 'r');
                            histogram(RT2, 20, 'FaceColor', 'm');
                            xlabel('RT')
                        end
                        
                        
                    else
                        hf = figure;
                        title(array.MOUSE_NM, 'Interpreter', 'none')
                        subplot(1, 2, 1)
                        hold on
                        RT1 = info.first_lick(info.plot_result == 1 & eval(cond1)) - info.delay(info.plot_result == 1 & eval(cond1));
                        plot(RT1, info.delay(info.plot_result == 1 & eval(cond1)), '.')
                        
                        RT2 = info.first_lick(info.plot_result == 1 & eval(cond2)) - info.delay(info.plot_result == 1 & eval(cond2));
                        plot(RT2, info.delay(info.plot_result == 1 & eval(cond2)), '.')
                        
                        xlabel('RT')
                        ylabel('Delay')

                        subplot(1, 2, 2)
                        hold on
                        histogram(RT1, 20);
                        histogram(RT2, 20);
                        legend(cond1, cond2, 'Interpreter', 'none')
                        xlabel('RT')
                    end
                end
            elseif comb == 0
                info = obj.combine_info('t_info');
                if strcmp(cond1, 'info.plot_result == -1')
                    hf = figure;
                    sgtitle('All Mice', 'Interpreter', 'none')
                    subplot(1, 2, 1)
                    hold on
                    
                    if nargin < 4
                        RT1 = info.first_lick(info.plot_result == -1);
                        scatter(RT1, info.delay(info.plot_result == -1), 'r', 'filled')
                        xlabel('RT')
                        ylabel('Delay')
                        subplot(1, 2, 2)
                        hold on
                        histogram(RT1, 60, 'FaceColor', 'r');
                        histogram(info.delay(info.plot_result == -1), 60, 'FaceColor', 'm');
                        legend('Prem. RT', 'Delay')
                        xlabel('RT')
                    else
                        RT1 = info.first_lick(info.plot_result == -1 & eval(cond2));
                        plot(RT1, info.delay(info.plot_result == -1 & eval(cond2)), 'r.')
                        
                        RT2 = info.first_lick(info.plot_result == -1 & ~eval(cond2));
                        plot(RT2, info.delay(info.plot_result == -1 & ~eval(cond2)), 'm.')
                        xlabel('RT')
                        ylabel('Delay')
                        subplot(1, 2, 2)
                        hold on
                        histogram(RT1, 60, 'FaceColor', 'r');
                        histogram(RT2, 60, 'FaceColor', 'm');
                        xlabel('RT')
                    end
                    
                    
                else
                    hf = figure;
                    title('All Mice', 'Interpreter', 'none')
                    subplot(1, 2, 1)
                    hold on
                    RT1 = info.first_lick(info.plot_result == 1 & eval(cond1)) - info.delay(info.plot_result == 1 & eval(cond1));
                    plot(RT1, info.delay(info.plot_result == 1 & eval(cond1)), '.')
                    
                    RT2 = info.first_lick(info.plot_result == 1 & eval(cond2)) - info.delay(info.plot_result == 1 & eval(cond2));
                    plot(RT2, info.delay(info.plot_result == 1 & eval(cond2)), '.')
                    
                    xlabel('RT')
                    ylabel('Delay')
                    
                    
                    ha = subplot(1, 2, 2);
                    hold on
                    h1 = histogram(RT1, 60);
                    h2 = histogram(RT2, 60);
                    legend(cond1, cond2, 'Interpreter', 'none')
                    xlabel('RT')
                end
            elseif comb == 2
                info = obj.combine_info('t_info');
                for kk = 1:numel(unique(info.day))
                    hf(kk) = figure;
                    sgtitle('All Mice', 'Interpreter', 'none')
                    subplot(1, 2, 1)
                    hold on
                    
                    RT1 = info.first_lick(info.plot_result == -1 & info.day == kk);
                    scatter(RT1, info.delay(info.plot_result == -1 & info.day == kk), 'r', 'filled')
                    xlabel('RT')
                    ylabel('Delay')
                    subplot(1, 2, 2)
                    hold on
                    histogram(RT1, 60, 'FaceColor', 'r');
                    histogram(info.delay(info.plot_result == -1 & info.day == kk), 60, 'FaceColor', 'm');
                    legend('Prem. RT', 'Delay')
                    xlabel('RT')
                end
            end
        end
        
        
        function h = saliency_mock(obj, ds_fac, comb_flg, target)
            if ~comb_flg
                for ii = 1:length(obj.MOUSE_ARRAY)
                h(ii) = figure;
                array = obj.MOUSE_ARRAY(ii);
                sgtitle(array.MOUSE_NM, 'Interpreter', 'none')
                t_info = array.INFO.t_info;
                cue_info = array.INFO.cue_info;
                
                bbn_data = downsample(array.normalize_data('trials', 1)', ds_fac)';
                cue_data = downsample(array.normalize_data('cue', 1)', ds_fac)';
                
                
                bbn_data = bbn_data - mean(bbn_data(:, floor(1:(5086/ds_fac))), 2);
                cue_data = cue_data - mean(cue_data(:, floor(1:(5086/ds_fac))), 2);
                
                t_vec = linspace(-5, 15, size(bbn_data, 2));
                subplot(array.DATA_CH, 2, 1)
                hold on
                clr = [41, 41, 116] ./255;
                shadedErrorBar(t_vec, bbn_data(t_info.is_light == 1 & t_info.first_lick > 1 & t_info.cue_int == 2, :), {@mean, @sem}, {'Color', 'c'});
                shadedErrorBar(t_vec, bbn_data(t_info.is_light == 0 & t_info.first_lick > 1 & t_info.cue_int < 0.5, :), {@mean, @sem}, {'Color', [clr, 0.4]});
                shadedErrorBar(t_vec, bbn_data(t_info.is_light == 0 & t_info.first_lick > 1 & t_info.cue_int == 0.5, :), {@mean, @sem}, {'Color', [clr, 0.6]});
                shadedErrorBar(t_vec, bbn_data(t_info.is_light == 0 & t_info.first_lick > 1 & t_info.cue_int == 1, :), {@mean, @sem}, {'Color', [clr, 0.8]});
                shadedErrorBar(t_vec, bbn_data(t_info.is_light == 0 & t_info.first_lick > 1 & t_info.cue_int == 2, :), {@mean, @sem}, {'Color', [clr, 1]});
                ylabel('DeltaF/F')
                xlim([-0.2, 1.5])
                subplot(array.DATA_CH, 2, 2)
                hold on
                shadedErrorBar(t_vec, cue_data(cue_info.is_light == 1 & cue_info.plot_result == 0 & cue_info.cue_int == 2, :), {@mean, @sem}, {'Color', 'c'});
                shadedErrorBar(t_vec, cue_data(cue_info.is_light == 0 & cue_info.plot_result == 0 & cue_info.cue_int < 0.5, :), {@mean, @sem}, {'Color', [clr, 0.4]});
                shadedErrorBar(t_vec, cue_data(cue_info.is_light == 0 & cue_info.plot_result == 0 & cue_info.cue_int == 0.5, :), {@mean, @sem}, {'Color', [clr, 0.6]});
                shadedErrorBar(t_vec, cue_data(cue_info.is_light == 0 & cue_info.plot_result == 0 & cue_info.cue_int == 1, :), {@mean, @sem}, {'Color', [clr, 0.8]});
                shadedErrorBar(t_vec, cue_data(cue_info.is_light == 0 & cue_info.plot_result == 0 & cue_info.cue_int == 2, :), {@mean, @sem}, {'Color', [clr, 1]});
                xlim([-0.2, 2])
                
                if ~isempty(array.PROJ_TO{2})
                    bbn_data = downsample(array.MATFILE.trials.af_trials', ds_fac)';
                    cue_data = downsample(array.MATFILE.cue.af_trials', ds_fac)';
                    bbn_data = bbn_data - mean(bbn_data(:, floor(1:(5086/ds_fac))), 2);
                    cue_data = cue_data - mean(cue_data(:, floor(1:(5086/ds_fac))), 2);
                    subplot(array.DATA_CH, 2, 3)
                    ylabel('DeltaF/F')
                    xlim([-0.2, 1.5])
                    hold on
                    shadedErrorBar(t_vec, bbn_data(t_info.first_lick > 1 & t_info.cue_int < 0.5, :), {@mean, @sem}, {'Color', [0.7, 0.7, 0.7, 0.7]});
                    shadedErrorBar(t_vec, bbn_data(t_info.first_lick > 1 & t_info.cue_int == 0.5, :), {@mean, @sem}, {'Color', [0.5, 0.5, 0.5, 0.7]});
                    shadedErrorBar(t_vec, bbn_data(t_info.first_lick > 1 & t_info.cue_int == 1, :), {@mean, @sem}, {'Color', [0.25, 0.25, 0.25, 0.7]});
                    shadedErrorBar(t_vec, bbn_data(t_info.first_lick > 1 & t_info.cue_int == 2, :), {@mean, @sem}, {'Color', [0, 0, 0, 0.7]});
                    subplot(array.DATA_CH, 2, 4)
                    hold on
                    shadedErrorBar(t_vec, cue_data(cue_info.is_light == 1 & cue_info.plot_result == 0 & cue_info.cue_int == 2, :), {@mean, @sem}, {'Color', 'c'});
                    shadedErrorBar(t_vec, cue_data(cue_info.is_light == 0 & cue_info.plot_result == 0 & cue_info.cue_int < 0.5, :), {@mean, @sem}, {'Color', [0.7, 0.7, 0.7, 0.7]});
                    shadedErrorBar(t_vec, cue_data(cue_info.is_light == 0 & cue_info.plot_result == 0 & cue_info.cue_int == 0.5, :), {@mean, @sem}, {'Color', [0.5, 0.5, 0.5, 0.7]});
                    shadedErrorBar(t_vec, cue_data(cue_info.is_light == 0 & cue_info.plot_result == 0 & cue_info.cue_int == 1, :), {@mean, @sem}, {'Color', [0.25, 0.25, 0.25, 0.7]});
                    shadedErrorBar(t_vec, cue_data(cue_info.is_light == 0 & cue_info.plot_result == 0 & cue_info.cue_int == 2, :), {@mean, @sem}, {'Color', [0, 0, 0, 0.7]});
                    xlim([-0.2, 4])
                end
                end
            else
                [acc_tch, ofc_tch, acc_tinfo, ofc_tinfo] = obj.combine_by_target('trials', 0);
                [acc_cch, ofc_cch, acc_cinfo, ofc_cinfo] = obj.combine_by_target('cue', 0);
                if target
                    bbn_data = downsample(acc_tch', ds_fac)';
                    cue_data = downsample(acc_cch', ds_fac)';
                    t_info = acc_tinfo;
                    cue_info = acc_cinfo;
                else
                    bbn_data = downsample(ofc_tch', ds_fac)';
                    cue_data = downsample(ofc_cch', ds_fac)';
                    t_info = ofc_tinfo;
                    cue_info = ofc_cinfo;
                end
                
                bbn_data = bbn_data - mean(bbn_data(:, floor(1:(5086/ds_fac))), 2);
                cue_data = cue_data - mean(cue_data(:, floor(1:(5086/ds_fac))), 2);
                
                t_vec = linspace(-5, 15, size(bbn_data, 2));
                subplot(1, 2, 1)
                hold on
                clr = [41, 41, 116] ./255;
                shadedErrorBar(t_vec, bbn_data(t_info.is_light == 1 & t_info.first_lick > 1 & t_info.cue_int == 2, :), {@mean, @sem}, {'Color', 'c'});
                shadedErrorBar(t_vec, bbn_data(t_info.is_light == 0 & t_info.first_lick > 1 & t_info.cue_int < 0.5, :), {@mean, @sem}, {'Color', [clr, 0.4]});
                shadedErrorBar(t_vec, bbn_data(t_info.is_light == 0 & t_info.first_lick > 1 & t_info.cue_int == 0.5, :), {@mean, @sem}, {'Color', [clr, 0.6]});
                shadedErrorBar(t_vec, bbn_data(t_info.is_light == 0 & t_info.first_lick > 1 & t_info.cue_int == 1, :), {@mean, @sem}, {'Color', [clr, 0.8]});
                shadedErrorBar(t_vec, bbn_data(t_info.is_light == 0 & t_info.first_lick > 1 & t_info.cue_int == 2, :), {@mean, @sem}, {'Color', [clr, 1]});
                ylabel('DeltaF/F')
                xlim([-0.2, 1.5])
                subplot(1, 2, 2)
                hold on
                shadedErrorBar(t_vec, cue_data(cue_info.is_light == 1 & cue_info.plot_result == 0 & cue_info.cue_int == 2, :), {@mean, @sem}, {'Color', 'c'});
                shadedErrorBar(t_vec, cue_data(cue_info.is_light == 0 & cue_info.plot_result == 0 & cue_info.cue_int < 0.5, :), {@mean, @sem}, {'Color', [clr, 0.4]});
                shadedErrorBar(t_vec, cue_data(cue_info.is_light == 0 & cue_info.plot_result == 0 & cue_info.cue_int == 0.5, :), {@mean, @sem}, {'Color', [clr, 0.6]});
                shadedErrorBar(t_vec, cue_data(cue_info.is_light == 0 & cue_info.plot_result == 0 & cue_info.cue_int == 1, :), {@mean, @sem}, {'Color', [clr, 0.8]});
                shadedErrorBar(t_vec, cue_data(cue_info.is_light == 0 & cue_info.plot_result == 0 & cue_info.cue_int == 2, :), {@mean, @sem}, {'Color', [clr, 1]});
                xlim([-0.2, 2])
            end
            
        end
        
        function hypovig(obj, plot_flg)
            for ii = 1:length(obj.MOUSE_ARRAY)
                array = obj.MOUSE_ARRAY(ii); 
                    
                if contains(obj.MOUSE_ARRAY(1).MATFILE.trials.Properties.Source, 'DREADD')
                    %if running for DREADDs, eliminate first 3 saline days from comp.
                    if numel(unique(array.INFO.t_info.day)) > 3
                        array.INFO.t_info(array.INFO.t_info.day < 4, :) = [];
                    end
                end
                    
                [env_trials, ~, seq_trials] = array.id_hypovig(5);
                info = array.INFO.t_info;
                                
                if plot_flg == 3
                    % plot baseline streak vs rest                    
                    acc_ch = find(strcmp(array.PROJ_TO, 'ACC'));
                    acc_ch = num2str(acc_ch);
                    if ~isempty(acc_ch)
                    end
                    info.plot_result(info.plot_result == 2) = 0;    %late is miss
                    D = diff([false,diff(info.plot_result)'==0,false]); % to get all misses by size of streak
                    misses_by_size = arrayfun(@(a,b)info.plot_result(a:b)',find(D>0),find(D<0),'uni',0);   % gets consective outcomes
                    misses_by_size(cellfun(@sum, misses_by_size)~=0) = [];                      %keeps only misses (sum = 0)
                    miss_size{ii} = cellfun(@numel, misses_by_size);                                % how many consecutive misses
                    miss_size{ii}(miss_size{ii} > 50) = [];
                    
                end
            end
       
            if plot_flg == 3
         
                figure;
                cdfplot(cell2mat(miss_size));
                figure;
                histogram(cell2mat(miss_size), 'BinEdges', linspace(1,50, 50), 'Normalization', 'count')
            end
        end
        

        function plot_sig(obj, cond)
            array = obj.MOUSE_ARRAY;
            for ii = 1:length(array)
                
                ch1 = array(ii).normalize_data('trials', 1);
                t_vec = linspace(-5, 15, size(ch1, 2));
                figure
                title([array(ii).MOUSE_NM, array(ii).PROJ_TO{1}], 'Interpreter', 'none')
                hold on
                t_info = array(ii).INFO.t_info;
                info = t_info(eval(cond), :);
                data = ch1(eval(cond), :);
                obj.MOUSE_ARRAY(ii).plot_by_outcome(t_vec, data, info, 1, 100);
                
                if ~isempty(array(ii).PROJ_TO{2})
                    ch2 = array(ii).normalize_data('trials', 2);
                    figure
                    title([array(ii).MOUSE_NM, array(ii).PROJ_TO{2}], 'Interpreter', 'none')
                    hold on
                    t_info = array(ii).INFO.t_info;
                    info = t_info(eval(cond), :);
                    data = ch2(eval(cond), :);
                    obj.MOUSE_ARRAY(ii).plot_by_outcome(t_vec, data, info, 1, 100);
                end
                
            end
        end
        
        function [data1, data2, info, info2] = combine_array_data(obj, data_to_combine)
            data1 = [];
            data2 = [];
            info = table;
            info2 = table;
            for ii = 1:length(obj.MOUSE_ARRAY)
                array = obj.MOUSE_ARRAY(ii);
                switch data_to_combine
                    case 'baseline'
                        data1{ii} = array.data.ch1_baseline;
                        if array.DATA_CH == 2
                            data2{ii} = array.data.ch2_baseline;
                        end
                    case 'passive'
                        data1 = [data1; array.MATFILE.passive.all_trials];
                        tmp_info = array.INFO.passive_info;
                        tmp_info.mouse(:) = ii;
                        info = [info; tmp_info];

                        if array.DATA_CH == 2
                            data2 = [data2; array.MATFILE.passive.af_trials];
                            tmp_info2 = array.INFO.passive_info;
                            tmp_info2.mouse(:) = ii;
                            info2 = [info2; tmp_info2];
                        end
                    case 'all'
                        data1 = [data1; array.MATFILE.trials.all_trials];
                        tmp_info = array.INFO.t_info;
                        tmp_info.mouse(:) = ii;
                        info = [info; tmp_info];
                        if array.DATA_CH == 2
                            data2 = [data2; array.MATFILE.trials.af_trials];
                            tmp_info2 = array.INFO.t_info;
                            tmp_info2.mouse(:) = ii;
                            info2 = [info2; tmp_info2];
                        end
                end
            end
        end
        
        function [r, p] = ch_correlation(obj, type, ds_fac)
            for ii = 1:length(obj.MOUSE_ARRAY)
                array = obj.MOUSE_ARRAY(ii);
                switch type
                    case 'passive'
                        info = array.INFO.passive_info;
                        ch1 = downsample(array.normalize_data('passive', 1)', ds_fac)';
                        ch2 = downsample(array.normalize_data('passive', 2)', ds_fac)';
                        
                        pre1 = reshape(ch1(strcmp(info.cond, 'pre_awake'), :)', 1, []);
                        pre2 = reshape(ch2(strcmp(info.cond, 'pre_awake'), :)', 1, []);
                        [r(1, ii), p(1, ii)] = corr(pre1(:), pre2(:));
                        
                        post1 = reshape(ch1(strcmp(info.cond, 'post_awake'), :)', 1, []);
                        post2 = reshape(ch2(strcmp(info.cond, 'post_awake'), :)', 1, []);
                        [r(2, ii), p(2, ii)] = corr(post1(:), post2(:));
                        
                    case 'task'
                        info = array.INFO.t_info;
                        ch1 = downsample(array.normalize_data('trials', 1)', ds_fac)';
                        ch2 = downsample(array.normalize_data('trials', 2)', ds_fac)';
                        
                        ch1 = reshape(ch1', 1, []);
                        ch2 = reshape(ch2', 1, []);
                        [r(ii), p(ii)] = corr(ch1(:), ch2(:));
                end
                
            end
        end
        
        function [ax, avg_data] = aligned_hms (obj, channel, ds_fac, zero_point, target)
            % If channel == 0 combine all data, otherwise plot indivudally
            % target 1 - ACC, 0 - OFC
            if channel~= 0
            for ii = 1:length(obj.MOUSE_ARRAY)
                array = obj.MOUSE_ARRAY(ii);
                if channel == 2 & array.DATA_CH == 1
                    % do nothing - skip mouse
                else
                ax(ii) = figure;
                sgtitle(array.MOUSE_NM, 'Interpreter', 'none')
                curr_ax = subplot(1, 4, 3);           
                ch = array.normalize_data('trials', channel);
                fs = size(ch, 2)/ 20;
                zero_sample = fs * zero_point;
                win = [-0.2, 1.5];
                s_win = floor(zero_sample + (win * fs)); % sample to start from, end on
                cut_ch = ch(:, s_win(1):s_win(2));
                cut_ch = cut_ch - mean(cut_ch(:, 1:floor(abs(win(1))*fs)), 2);
                
                ds_ch = obj.smooth_sig(cut_ch, ds_fac);
                t_vec = linspace(win(1), win(2), size(ds_ch, 2));
                
                sample_norm = ds_ch;
                avg_data{ii, 3} = mean(sample_norm) ./ max(mean(sample_norm));
                
                % sorting:
                [srted_idx, srted_info] = array.sort_info_by('trials', 'cue_int');
                srted_idx(srted_info.first_lick < 1) = []; % remove trials with early lick
                srted_data = sample_norm(srted_idx, :);
                array.plot_heatmap(srted_data, t_vec, curr_ax);

                line([0,0], [0, size(srted_data, 1)], 'Color', 'c', 'LineWidth', 3)
                title('Trial Onset')

                colormap([linspace(1, 41/255, 100)', linspace(1, 41/255, 100)', linspace(1, 116/255, 100)'])
                xlabel('Time (s)')
                ylabel('Trial #')
                curr_ax = subplot(1, 4, 4);
                ch = array.normalize_data('cue', channel);
                win = [-0.2, 1.5];
                s_win = floor(zero_sample + (win * fs)); % sample to start from, end on
                cut_ch = ch(:, s_win(1):s_win(2));
                cut_ch = cut_ch - mean(cut_ch(:, 1:floor(abs(win(1))*fs)), 2);
                ds_ch = obj.smooth_sig(cut_ch, ds_fac);
                ds_ch = ds_ch - mean(ds_ch(1:floor(abs(win(1)) * fs/ds_fac)));
                t_vec = linspace(win(1), win(2), size(ds_ch, 2));

                sample_norm = ds_ch;
                avg_data{ii, 4} = mean(sample_norm) ./ max(mean(sample_norm));
                % sorting:
                [srted_idx, srted_info] = array.sort_info_by('cue', 'cue_int');
                srted_data = sample_norm(srted_idx, :);
                srted_info = [srted_info(srted_info.is_light == 0, :); srted_info(srted_info.is_light > 0, :)];
                srted_data = [srted_data(srted_info.is_light == 0, :); srted_data(srted_info.is_light > 0, :)];
                srted_data(srted_info.plot_result~=0,:) = []; % keep only miss trials
                
                
                array.plot_heatmap(flipud(srted_data), t_vec, curr_ax);         % start with high intensities first
                line([0,0], [0, size(srted_data, 1)], 'Color', 'c', 'LineWidth', 3)
                title('Cue Onset')
                xlabel('Time (s)')
                ylabel('Trial #')
                
                curr_ax = subplot(1, 4, 1);
                ch = array.normalize_data('trials', channel);
                win = [-0.2, 10];
                s_win = floor(zero_sample + (win * fs)); % sample to start from, end on
                cut_ch = ch(:, s_win(1):s_win(2));
                cut_ch = cut_ch - mean(cut_ch(:, 1:floor(abs(win(1))*fs)), 2);
                ds_ch = obj.smooth_sig(cut_ch, ds_fac);
                t_vec = linspace(win(1), win(2), size(ds_ch, 2));
                sample_norm = ds_ch;
                avg_data{ii, 1} = mean(sample_norm)./ max(mean(sample_norm));
                [srted_idx, srted_info] = array.sort_info_by('trials', 'delay');
                srted_data = sample_norm(srted_idx, :);
                srted_data = srted_data(srted_info.plot_result == 1, :);
              
                array.plot_heatmap(srted_data, t_vec, curr_ax);
                title('Hit Trials')
                xlabel('Time (s)')
                ylabel('Trial #')
                lick_times = srted_info.first_lick(srted_info.plot_result == 1);
                line(repmat([lick_times'; lick_times'], [2, 1]), repmat([[1:size(lick_times)]-0.8; [1:size(lick_times)]+0.8], [2, 1]), 'Color', [62,98,83] ./ 255, 'LineWidth', 2)
                 
                 curr_ax = subplot(1, 4, 2);
                 win = [-0.2, 10];
                 s_win = floor(zero_sample + (win * fs)); % sample to start from, end on
                 cut_ch = ch(:, s_win(1):s_win(2));
                 cut_ch = cut_ch - mean(cut_ch(:, 1:floor(abs(win(1))*fs)), 2);
                 ds_ch = obj.smooth_sig(cut_ch, ds_fac);
                 t_vec = linspace(win(1), win(2), size(ds_ch, 2));
                 sample_norm = ds_ch;
                 avg_data{ii, 2} = mean(sample_norm)./ max(mean(sample_norm)); 
                 [srted_idx, srted_info] = array.sort_info_by('trials', 'first_lick');
                 srted_data = sample_norm(srted_idx, :);
                 srted_data = srted_data(srted_info.plot_result == -1, :);
                 array.plot_heatmap(srted_data, t_vec, curr_ax);
                 title('Premature Trials')
                 xlabel('Time (s)')
                 ylabel('Trial #')
                 lick_times = srted_info.first_lick(srted_info.plot_result == -1);
                 line(repmat([lick_times'; lick_times'], [2, 1]), repmat([[1:size(lick_times)]-0.8; [1:size(lick_times)]+0.8], [2, 1]), 'Color', [129,88,135]./ 255, 'LineWidth', 2)
                 
                 hms = findall(ax(ii), 'Type', 'axes');
                 set(hms, 'CLim', [-0.5, 1.5])
                 set(ax(ii), 'Position', [ 6, 328, 1672, 770]);
                end
            end
            else
                ax = figure;
                sgtitle('All Mice')
                curr_ax = subplot(1, 3, 3);
                
                [acc_ch, ofc_ch, acc_info, ofc_info] = obj.combine_by_target('cue', 0);
                if target
                    ch = acc_ch;
                    info = acc_info;
                else
                    ch = ofc_ch;
                    info = ofc_info;
                end
                
                fs = size(ch, 2)/ 20;
                zero_sample = fs * zero_point;
                
                win = [-0.2, 2];
                s_win = floor(zero_sample + (win * fs)); % sample to start from, end on
                cut_ch = ch(:, s_win(1):s_win(2));
                cut_ch = cut_ch - mean(cut_ch(:, 1:floor(abs(win(1))*fs)), 2);
                ds_ch = obj.smooth_sig(cut_ch, ds_fac);
                ds_ch = ds_ch - mean(ds_ch(1:floor(abs(win(1)) * fs/ds_fac)));
                t_vec = linspace(win(1), win(2), size(ds_ch, 2));
                
                sample_norm = ds_ch;
                avg_data{4} = mean(sample_norm) ./ max(mean(sample_norm));
                
                % sorting:
                [~, d_idx] = sort(info.('cue_int'));
                srted_idx = d_idx;
                srted_info = info(d_idx, :);
                
                % plot cue breaks
                srted_info.cue_int(srted_info.cue_int == 0.1) = 0.5;
                
                srted_data = sample_norm(srted_idx, :);
                srted_info = [srted_info(srted_info.is_light == 0, :); srted_info(srted_info.is_light > 0, :)];
                srted_data = [srted_data(srted_info.is_light == 0, :); srted_data(srted_info.is_light > 0, :)];
                srted_data(srted_info.plot_result~=0,:) = []; % keep only miss trials
                
               
                
                obj.MOUSE_ARRAY(1).plot_heatmap(flipud(srted_data), t_vec, curr_ax);         % start with high intensities first
                colorbar('off')
                line([0,0], [0, size(srted_data, 1)], 'Color', 'k', 'LineWidth', 2)
                
                srted_info(srted_info.plot_result ~= 0, :) = [];
                c_breaks = find(diff(srted_info.cue_int) ~= 0);
                first_light = find(srted_info.is_light == 1, 1, 'first');
                c_breaks(c_breaks > first_light) = [];
                c_breaks = abs(c_breaks - size(srted_data, 1));
                line(curr_ax, repmat(win, [4,1])', repmat(c_breaks, [1,2])', 'Color', 'k')

                title('Cue Onset')
                xlabel('Time (s)')

                curr_ax = subplot(1, 3, 1);
               [acc_ch, ofc_ch, acc_info, ofc_info] = obj.combine_by_target('trials', 0);
                if target
                    ch = acc_ch;
                    info = acc_info;
                else
                    ch = ofc_ch;
                    info = ofc_info;
                end
                if target 
                    win = [-0.2, 6];
                else
                    win = [-0.2, 10];
                end
                
                s_win = floor(zero_sample + (win * fs)); % sample to start from, end on
                cut_ch = ch(:, s_win(1):s_win(2));
                cut_ch = cut_ch - mean(cut_ch(:, 1:floor(abs(win(1))*fs)), 2);
                ds_ch = obj.smooth_sig(cut_ch, ds_fac);
                t_vec = linspace(win(1), win(2), size(ds_ch, 2));
sample_norm = ds_ch;
avg_data{1} = mean(sample_norm)./ max(mean(sample_norm));                
                 [~, d_idx] = sort(info.('delay'));
                srted_idx = d_idx;
                srted_info = info(d_idx, :);

                srted_data = sample_norm(srted_idx, :);
                srted_data = srted_data(srted_info.plot_result == 1, :);
              
                obj.MOUSE_ARRAY(1).plot_heatmap(srted_data, t_vec, curr_ax);
                colorbar('off')
                line([0,0], [0, size(srted_data, 1)], 'Color', 'k', 'LineWidth', 2)
                title('Hit Trials')
                xlabel('Time (s)')
                ylabel('Trial #')
                lick_times = srted_info.first_lick(srted_info.plot_result == 1);
                line(curr_ax, repmat([lick_times'; lick_times'], [2, 1]), repmat([[1:size(lick_times)]-0.8; [1:size(lick_times)]+0.8], [2, 1]), 'Color', [62,98,83] ./ 255)%, 'LineWidth', 2)

                 
                 curr_ax = subplot(1, 3, 2);
                 if target
                     win = [-0.2, 6];
                 else
                     win = [-0.2, 10];
                 end
                 s_win = floor(zero_sample + (win * fs)); % sample to start from, end on
                 cut_ch = ch(:, s_win(1):s_win(2));
                 cut_ch = cut_ch - mean(cut_ch(:, 1:floor(abs(win(1))*fs)), 2);
                 ds_ch = obj.smooth_sig(cut_ch, ds_fac);
                 t_vec = linspace(win(1), win(2), size(ds_ch, 2));
                 sample_norm = ds_ch;
                 avg_data{2} = mean(sample_norm)./ max(mean(sample_norm)); 
                  [~, d_idx] = sort(info.('first_lick'));
                srted_idx = d_idx;
                srted_info = info(d_idx, :);

                 srted_data = sample_norm(srted_idx, :);
                 srted_data = srted_data(srted_info.plot_result == -1, :);
                 
                 obj.MOUSE_ARRAY(1).plot_heatmap(srted_data, t_vec, curr_ax);
                 colorbar('off')
                 line([0,0], [0, size(srted_data, 1)], 'Color', 'k', 'LineWidth', 2)
                 title('Premature Trials')
                 xlabel('Time (s)')
                 lick_times = srted_info.first_lick(srted_info.plot_result == -1);
                 line(curr_ax, repmat([lick_times'; lick_times'], [2, 1]), repmat([[1:size(lick_times)]-0.8; [1:size(lick_times)]+0.8], [2, 1]), 'Color', [129,88,135]./ 255) %, 'LineWidth', 2)
                 
                 hms = findall(ax, 'Type', 'axes');       
                 ax.Position(3) = 1.5*ax.Position(3);
                 ax.Position(4) = 2 * ax.Position(4);
                 set(hms, 'CLim', [-0.5, 1.5])
                 colormap([linspace(1, 41/255, 100)', linspace(1, 41/255, 100)', linspace(1, 116/255, 100)']) % set color to blue
            end
        end
        
        function plot_baseline_summary(obj, what_to_plot)
            % plots a summary plot and fit for all mice in array
            proj_targets = fieldnames(obj.param_data);
            proj_targets(strcmp(proj_targets, 'zscore')) = [];
            hf = figure;
            clr = {'.r', '.g', '.k'};
            nme = {'prem', 'corr', 'miss'};
            
            for tt = 1:numel(proj_targets)
                prem = [];
                corr = [];
                miss = [];
                data = struct;
                bseline = [];
                
                fnames = fieldnames(obj.param_data.(proj_targets{tt}));
                for ii = 1:numel(fnames)
                    params = obj.param_data.(proj_targets{tt}).(fnames{ii}).baseline_probs;
                    prem(ii, :) = params(1, :);
                    corr(ii, :) = params(2, :);
                    miss(ii, :) = params(3, :);
                    
                    mse_idx = find(strcmp({obj.MOUSE_ARRAY.MOUSE_NM}, fnames{ii}(3:end))); % find where this mouse is in the object
                    rel_ch = find(cell2mat((cellfun(@strcmp, {obj.MOUSE_ARRAY(mse_idx).PROJ_TO}, proj_targets(tt), 'UniformOutput', false)))); % which channel
                    switch what_to_plot
                        case 'avg'
                            bseline(ii) = mean(obj.MOUSE_ARRAY(mse_idx).data.(['ch', num2str(rel_ch),'_baseline'])); %get avg baseline
                        case 'var'
                            bseline(ii) = var(obj.MOUSE_ARRAY(mse_idx).data.(['ch', num2str(rel_ch),'_baseline'])); %get baseline variance
                    end
                end
                data.prem = prem;
                data.corr = corr;
                data.miss = miss;
               
                
                for kk = 1:3
                    figure(hf)
                    subplot(numel(proj_targets), 3, tt + 2*(tt-1) + (kk-1))
                    hold on
                    obj.plot_corr(data.(nme{kk}), clr{kk}, nme{kk});
                    if kk == 2
                        title(proj_targets{tt})
                    end
                end
            end
            
            ha = findall(hf, 'Type', 'Axes');
            set(ha, 'YLim', [0, max([ha.YLim])], 'XTick', [0,1], 'XTickLabel', {'1', num2str(length(prem))}); % sprintfc('%g',1:length(pre_bot))
            hlab = [ha.XLabel];
            set(hlab, 'String', 'Quantile')
        end
        
        function [acc_ch, ofc_ch, acc_info, ofc_info] = combine_by_target(obj, data_to_combine, axon_data)
            switch data_to_combine
                case 'trials'
                    info_nm = 't_info';
                case 'lick'
                    info_nm = 'l_info';
                case 'cue'
                    info_nm = 'cue_info';
            end
            acc_ch = [];
            acc_info = [];
            ofc_ch = [];
            ofc_info = [];
            array = obj.MOUSE_ARRAY;
            for ii = 1:length(array)
                if axon_data
                    %for inputs, use:
                    acc = find(strcmp(array(ii).PROJ_TO, 'ACC_input'));
                    ofc = find(strcmp(array(ii).PROJ_TO, 'AUD_input'));
                else
                    acc = find(strcmp(array(ii).PROJ_TO, 'ACC'));
                    ofc = find(strcmp(array(ii).PROJ_TO, 'OFC'));
                end
                if ~isempty(acc)
                    acc_ch = [acc_ch; array(ii).normalize_data(data_to_combine, acc)];
                    info = array(ii).INFO.(info_nm);
                    if strcmp(info.Properties.VariableNames(end), 'activity_tag_ch1')
                        info(:, end + 1) = {nan};
                        info.Properties.VariableNames(end) = {'activity_tag_ch2'};
                    end
                    info(:, end+1) = table({array(ii).MOUSE_NM}, 'VariableName' , {'mouse'});
                    acc_info = [acc_info; info];
                end
                if ~isempty(ofc)
                    ofc_ch = [ofc_ch; array(ii).normalize_data(data_to_combine,ofc)];
                    info = array(ii).INFO.(info_nm);
                    if strcmp(info.Properties.VariableNames(end), 'activity_tag_ch1')
                        info(:, end + 1) = {nan};
                        info.Properties.VariableNames(end) = {'activity_tag_ch2'};
                    end
                    info(:, end+1) = table({array(ii).MOUSE_NM}, 'VariableName' , {'mouse'});
                    ofc_info = [ofc_info; info];
                    
                end
            end
        end
        
        
        function  arg_out = split_array(obj, param, calc_flg, plot_flg, epoch, group_by)
            % compare groups by behavior, calc flg gets numbers, plot flg
            % plots data, can be independent -CURRENTLY WORKS ONLY TOGETHER
            % 
            groups = obj.group_members(group_by);
             
            for ii = 1:numel(groups)
                grpnmes = groups{ii};
                tmp_obj = copy(obj);
                tmp_obj.MOUSE_ARRAY(~ismember({tmp_obj.MOUSE_ARRAY.MOUSE_NM}, grpnmes)) = [];                
                f_names = fieldnames(tmp_obj.param_data);
                f_names(strcmp(f_names, 'zscore')) = [];
                if numel(f_names) > 0
                    for kk = 1:numel(f_names)
                        mse_nmes = fieldnames(tmp_obj.param_data.(f_names{kk}));
                        tmp_obj.param_data.(f_names{kk}) = rmfield(tmp_obj.param_data.(f_names{kk}), (mse_nmes(~contains(mse_nmes, grpnmes))));
                    end
                end
               
                if calc_flg
                    switch epoch
                        case 'baseline'
                             [r, p, null] = tmp_obj.baseline_quant(30, 'trials', 't_info', 0);
                        case 'RT'
                            if ~strcmp(group_by, 'hm3dq')
                                hf(ii) = tmp_obj.analyze_rt(0, 'info.plot_result == -1');
                            else
                                tmp_obj.analyze_rt(2, 'info.plot_result == -1');
                            end
                        case 'hit_RT'
                            info = tmp_obj.combine_info('t_info');
                            info = info(info.plot_result == 1, :);
                            num_mice = unique(info.mse_nm);
                            RT_all = info.first_lick - info.delay;
                            for jj = 1:numel(num_mice)
                                RT(jj) = mean(RT_all(find(strcmp(num_mice{jj}, info.mse_nm))));
                            end
                            figure; plot(ii, RT, 'ob')
                            hold on
                            plot(ii, mean(RT), 'ok')
                            title(['group ', num2str(ii)])
                        case 'outcome'
                            tmp_obj.grouped_outcomes(0, 1)
                    end
                    if isfield('ACC', obj.param_data)
                        obj.param_data = [obj.param_data, tmp_obj.param_data];
                    else
                        obj.param_data = tmp_obj.param_data;
                    end
                    
                    if ~(strcmp(epoch, 'baseline') | strcmp(epoch, 'cloud') | strcmp(epoch, 'hit_RT') | strcmp(epoch, 'RT')  | strcmp(epoch, 'outcome'))
                        sz = size(h1);
                        root_find = get(h1, 'Type');
                        h1(strcmp(root_find, 'root')) = [];                         % remove empty rows (from single channel)
                        h1 = reshape(h1, [], sz(2));                                % organize back
                        h2 = tmp_obj.merge_ind_plots(h1);                           % combine plots from the same mouse;
                        set(h2, 'UserData' ,ii);                                    % label as group number
                        obj.plot_handles.(epoch).individual_plots.(param) = h2;
                        
                    end
                end
                if plot_flg
                    switch epoch
                        case 'outcome_hist'
                            
                            for tt = 1:numel(tmp_obj.MOUSE_ARRAY)
                                if contains(group_by, 'hm')
                                    if ii == 1
                                        bkp_info = tmp_obj.MOUSE_ARRAY(tt).INFO.t_info;
                                        tmp_info = (tmp_obj.MOUSE_ARRAY(tt).INFO.t_info(tmp_obj.MOUSE_ARRAY(tt).INFO.t_info.day > 3, :));
                                        tmp_obj.MOUSE_ARRAY(tt).INFO.t_info = tmp_info;
                                        tmp_obj.MOUSE_ARRAY(tt).plot_outcome_dist;
                                        tmp_obj.MOUSE_ARRAY(tt).INFO.t_info = bkp_info;
                                    else
                                        tmp_obj.MOUSE_ARRAY(tt).plot_outcome_dist;
                                    end
                                    title(tmp_obj.MOUSE_ARRAY(tt).MOUSE_NM, 'Interpreter', 'none')
                                end
                            end
                            info = tmp_obj.combine_info('t_info');
                            if contains(group_by, 'hm')
                                if ii == 1
                                    info = info(info.day > 3, :);       % take only day 4-5-6 for SAL data
                                end
                            else
                            end
                            pr = info.trial_number(info.plot_result == -1);
                            om = info.trial_number(info.plot_result == 0 | info.plot_result == 2);
                            cr = info.trial_number(info.plot_result == 1);
                            h_hist = axes(figure);
                            hold (h_hist, 'on')
                            histfit(cr, 50, 'kernel')
                            histfit(om, 50, 'kernel')
                            histfit(pr, 50, 'kernel')
                            line_h = flipud(findall(h_hist, 'Type', 'Line'));          % order is reversed (last one drawn is first)
                            bar_h = flipud(findall(h_hist, 'Type', 'Bar'));
                            legend(h_hist, line_h, {'correct', 'miss', 'premature'})
                            clrs = [204, 204, 204; 62, 98, 83; 129, 88, 135] ./255;
                            line_h(1).Color = clrs(2,:); line_h(2).Color = clrs(1,:); line_h(3).Color = clrs(3,:);
                            bar_h(1).FaceColor = clrs(2,:); bar_h(1).FaceAlpha = 0.9;
                            bar_h(2).FaceColor = clrs(1,:); bar_h(2).FaceAlpha = 0.9;
                            bar_h(3).FaceColor = clrs(3,:); bar_h(3).FaceAlpha = 0.9;
                            title(['group ', num2str(ii), ' summary'])

                        case 'baseline'
                            tmp_obj.plot_summary('baseline', 'low', 1, 'ACC', 'OFC');
                            tmp_obj.plot_summary('baseline', 'high', 1, 'ACC', 'OFC');
                        
                        case 'hypovig'
                            tmp_obj.hypovig(3);
                        case 'baseline_summary'
                            [acc1(ii), acc2(ii)] = tmp_obj.baseline_by_outcome_sum();
                            arg_out = [acc1; acc2];
                        case 'psych'
                            tmp_obj.behavior_curves(1);
                        case 'bbn_response'
                            tmp_obj.trial_onset_summary;
                            sgtitle(['Group ', num2str(ii)])
                        otherwise
                            [obj.plot_handles.trials.summary_plots.(param), a_sum, o_sum] = tmp_obj.plot_summary(epoch, param, 1, 'ACC', 'OFC');
                            fighandles = get(groot, 'Children');
                            sum_f = findobj(fighandles, 'Type', 'Figure', 'Visible', 'on');
                            close(fighandles(~ismember(fighandles, sum_f)))
                            
                            sum_f = sum_f(cellfun(@isnumeric, {sum_f.UserData}) & ~cellfun(@isempty, {sum_f.UserData}));    % handles for individual figures
                            sum_f = sum_f([sum_f.UserData] == ii);                                                          % handles of current group figure
                            
                            tmp_obj.set_plot_scales(sum_f);
                            [p, tbl, stats] = tmp_obj.do_statistics((epoch), a_sum, o_sum);
                            disp(p)
                    end
                end
            end
            % group summary code:
            if strcmp(epoch, 'RT') & ~strcmp(group_by, 'hm3dq')
                arg_out = obj.create_RT_summary(hf);
            end
        end
        
        function ax = hit_rt_vs_baseline(obj, split_flg)
            if ~split_flg % ALL mice
            all_info = obj.combine_info('t_info');
            [data1, data2] = obj.combine_array_data('baseline');
            data1 = cellfun(@transpose, data1, 'UniformOutput', false);
            data1 = [data1{:}];
            data2 = cellfun(@transpose, data2, 'UniformOutput', false);
            data2 = [data2{:}];
            crrct = all_info(all_info.plot_result ==1, :);
            hitRT = crrct.first_lick - crrct.delay;
            ax = axes(figure); 
            plot(ax, data1(all_info.plot_result == 1), hitRT, '.');
            mice = unique(all_info.mse_nm, 'Stable');
            single_ch_mice = [obj.MOUSE_ARRAY.DATA_CH] ~= 2;
            dual_ch_info = all_info(~contains(all_info.mse_nm, mice(single_ch_mice)), :);   % get just the info for dual channels
            hold on
            plot(ax, data2(dual_ch_info.plot_result == 1), hitRT(~contains(crrct.mse_nm, mice(single_ch_mice))), '.');
            ylabel('RT (sec)')
            xlabel('Pretrial Activity')
            else
                [acc_ch, ofc_ch, acc_info, ofc_info] = X.combine_by_target('trials', 0);
                bseline = mean(acc_ch(acc_info.plot_result == 1, 1000:5000));
                bseline = mean(acc_ch(acc_info.plot_result == 1, 1000:5000), 2);
                RTs = acc_info.first_lick(acc_info.plot_result == 1) - acc_info.delay(acc_info.plot_result == 1);
                plot(bseline, RTs, '.')
                ylabel('RT (sec)')
                xlabel('Pretrial Activity')
            end
        end
        
        function save_array(obj)
            mse_array = obj.MOUSE_ARRAY;
            [baseFileName, folder] = uiputfile('mse_array.mat', 'Save array to file');
            if baseFileName == 0
                % User clicked the Cancel button.
            end
            fullFileName = fullfile(folder, baseFileName);
            save(fullFileName, 'mse_array', '-v7.3');
        end
                % ----------------------------- Panel graphs ----------- %
                
        function behavior_curves(obj, merge)
            % plots psychometric curves for all mice in
            % array
            for ii = 1:length(obj.MOUSE_ARRAY)
                array = obj.MOUSE_ARRAY(ii);
                obj.MOUSE_ARRAY(ii).INFO.t_info.cue_int(obj.MOUSE_ARRAY(ii).INFO.t_info.cue_int == 0.1) = 0.5;% Fix wrongly labeld attenuations
                info = obj.MOUSE_ARRAY(ii).INFO.t_info;
              
                if merge
                    sum_outcme(ii, :, :) = array.analyze_behavior(0); % colums are [light+cloud, light-cloud, cloud, -cloud];               
                    prm(ii, :, :) = [sum(info.plot_result == -1 & info.cloud == 0) / sum(info.cloud == 0), ...
                                    sum(info.plot_result == -1 & info.cloud > 0) / sum(info.cloud > 0)];
                else
                    array.analyze_behavior(1);
                end
            end
            
            if merge
                avg_psych = squeeze(mean(sum_outcme, 1));
                sem_psych = squeeze(sem(sum_outcme));
                h_ax = axes(figure); hold on
                
                cloud_clr = [0.85,0.33,0.10];
                no_cloud_clr = [0.00,0.35,0.74];
                
                errorbar(avg_psych(:, [1,3]), sem_psych(:, [1,3]), 'o', 'Color', cloud_clr)
                errorbar(avg_psych(:, [2,4]), sem_psych(:, [2,4]), 'o', 'Color', no_cloud_clr)
                for ii = 1:4
                    [~, curve] = FitPsycheCurveWH(1:4, avg_psych(:, ii), 1);
                    tmp_p = plot(curve(:, 1), curve(:, 2));
                    if mod(ii, 2) == 0
                        set(tmp_p, 'Color', no_cloud_clr);
                    else
                        set(tmp_p, 'Color', cloud_clr);
                    end
                    if ii > 2
                        set(tmp_p, 'LineStyle', '--')
                    end
                end
                set(h_ax, 'YLim', [0,1], 'XTick', [1:4], 'XTickLabel', {'0.02', '0.5', '1', '2'}, 'XLim', [0.9, 4.1])
                ylabel(h_ax, 'Hit Rate %')
                xlabel(h_ax, 'Cue Intensity')
                
                ax_2 = axes(h_ax.Parent, 'Position', [0.6, 0.16, 0.27, 0.27], 'XTick', [1:2], 'XTickLabel', {'no cloud', 'cloud'});
                hold(ax_2, 'on')
                ylabel(ax_2, 'Premature Rate')
                bar(ax_2, 1, nanmean(squeeze(prm(:, :, 1))), 'FaceColor', no_cloud_clr)
                bar(ax_2, 2, nanmean(squeeze(prm(:, :, 2))), 'FaceColor', cloud_clr)
                errorbar(nanmean(squeeze(prm(:, :, 1))), sem(squeeze(prm(:, :, 1))), 'Color', no_cloud_clr)
                errorbar(2, nanmean(squeeze(prm(:, :, 2))), sem(squeeze(prm(:, :, 2))), 'Color', cloud_clr)
                
            end
        end
        
        function grouped_outcomes(obj, show_ind, groupby)
            info = obj.combine_info('t_info');
            mice = unique(info.mse_nm);
            for ii = 1:numel(mice)
                sub_info = info(strcmp(info.mse_nm, mice{ii}), :);
                res(1, ii) = sum(sub_info.plot_result == 1) / height(sub_info);
                res(2, ii) = sum(sub_info.plot_result == -1) / height(sub_info);
                res(3, ii) = sum(abs(sub_info.plot_result) ~= 1) / height(sub_info); % misses
            end
            clr = [63, 98, 83; 129, 89, 135; 113, 114, 114] ./ 255;
            
            if ~show_ind
                f_main = figure;
                hold on
                for ii = 1:3
                    hb = bar(ii, mean(res(ii, :)), 'FaceColor', clr(ii, :), 'LineStyle', 'none');
                    if ii == 3
                        hb.FaceAlpha = 0.1;
                    end
                    errorbar(ii, mean(res(ii, :)), sem(res(ii, :)'), 'k', 'LineWidth', 2)
                end
            else
                
                f_main = figure;
                hold on
                hp = plot(1:3, res, '.', 'Color', [0.8, 0.8, 0.8]);
                for ii = 1:3
                    errorbar(ii, mean(res(ii, :)), sem(res(ii, :)'), 'ko', 'MarkerFaceColor', clr(ii, :), 'LineWidth', 2)
                end
                
                if nargin > 2
                    group = obj.group_members(groupby);
                    for ii = 1:numel(group)
                        idx = contains(mice, group{ii});
                        set(hp(idx), 'Color', clr(ii, :));
                        f(ii) = figure;
                        copyobj(hp(idx), axes(f(ii)))
                        set(f(ii).Children.Title, 'String', ['Group ', num2str(ii)])
                    end
                end
            end
            figure(f_main);
            xlim([0.25, 3.75])
            xticks(1:3);
            xticklabels({'Hit', 'Premature', 'Miss'})
            ylabel('Proportion of Trials %')
            title('Task Outcomes')
           
        end

        function [ax_acc, ax_acc2, ax_ofc, ax_ofc2] = baseline_by_outcome_sum(obj)
            [acc_ch, ofc_ch, acc_info, ofc_info] = obj.combine_by_target('trials', 0);           
            t_vec = linspace(-5, 15, size(acc_ch, 2));
            acc_info.plot_result(acc_info.plot_result == 2) = 0;
            acc_ax = axes(figure);
            obj.MOUSE_ARRAY(1).plot_by_outcome(t_vec, acc_ch, acc_info, 1, 100);
            a = acc_ax.Children;
            acc_subax = axes(acc_ax.Parent, 'Position', [0.6, 0.53, 0.29, 0.35]);
            copyobj(a, acc_subax)
            acc_subax.XLim = [-4.5, -0.5];
            set(acc_ax.Legend, 'Position', [0.1517, 0.8047, 0.1911, 0.1048])
            rectangle(acc_ax, 'Position', [-4.5, -0.1, 4, 0.2])
            title(acc_ax, 'ACCp')           
            ylabel(acc_ax, 'Zscored \DeltaF/F');
            xlabel(acc_ax, 'Time (s)');
%             
            acc_bseline = mean(acc_ch(:, 1000:5000), 2);
            
            
            acc_tbl = table(double(acc_bseline), categorical(acc_info.plot_result), ...
                            categorical(acc_info.day), categorical(acc_info{:,end}), ...
                            'VariableNames', {'baseline', 'outcome', 'day', 'mouse'}); 
                        
            ofc_ax = axes(figure);
            if ~isempty (ofc_info) && numel(unique(ofc_info.Var15)) > 1
                
                ofc_info.plot_result(ofc_info.plot_result == 2) = 0;
                obj.MOUSE_ARRAY(1).plot_by_outcome(t_vec, ofc_ch, ofc_info, 1, 100);
                o = ofc_ax.Children;
                ofc_subax = axes(ofc_ax.Parent, 'Position', [0.6, 0.53, 0.29, 0.35]);
                copyobj(o, ofc_subax)
                ofc_subax.XLim = [-4.5, -0.5];
                set(ofc_ax.Legend, 'Position', [0.1517, 0.8047, 0.1911, 0.1048])
                rectangle(ofc_ax, 'Position', [-4.5, -0.1, 4, 0.2])
                xlabel(ofc_ax, 'Time (s)');
                ylabel(ofc_ax, 'Zscored \DeltaF/F');
                title(ofc_ax, 'OFCp')
                ofc_bseline = mean(ofc_ch(:, 1000:5000), 2);            
                ofc_tbl = table(double(ofc_bseline), categorical(ofc_info.plot_result), ...
                            categorical(ofc_info.day), categorical(ofc_info{:,end}), ...
                            'VariableNames', {'baseline', 'outcome', 'day', 'mouse'});
                                      
              ofc_info.trial_result(ofc_info.plot_result ~= -1) = {'not-prem'};
               ofc_info.trial_result(ofc_info.plot_result == 0) = {'miss'};
               ofc_info.trial_result(ofc_info.plot_result == 1) = {'hit'};
              
              tmp_tbl_ofc_not_prem = table(ofc_bseline(ofc_info.plot_result ~= -1), ofc_info.(ofc_info.Properties.VariableNames{end})(ofc_info.plot_result ~= -1),  'VariableNames', {'Values', 'mouse'});
              tmp_tbl_ofc_miss = table(ofc_bseline(ofc_info.plot_result == 0), ofc_info.(ofc_info.Properties.VariableNames{end})(ofc_info.plot_result == 0),  'VariableNames', {'Values', 'mouse'});
              tmp_tbl_ofc_corr = table(ofc_bseline(ofc_info.plot_result == 1), ofc_info.(ofc_info.Properties.VariableNames{end})(ofc_info.plot_result == 1),  'VariableNames', {'Values', 'mouse'});
              tmp_tbl_ofc_prem = table(ofc_bseline(ofc_info.first_lick < 0.5), ofc_info.(ofc_info.Properties.VariableNames{end})(ofc_info.first_lick < 0.5),  'VariableNames', {'Values', 'mouse'});
              statarray = grpstats(tmp_tbl_ofc_miss, 'mouse');
              statarray1 = grpstats(tmp_tbl_ofc_corr, 'mouse');
              statarray2 = grpstats(tmp_tbl_ofc_prem, 'mouse');
              statarray3 = grpstats(tmp_tbl_ofc_not_prem, 'mouse');
                            tmp_tbl_ofc = table([statarray1{:,3}; statarray{:,3}; statarray2{:,3}], ...
                  [repmat({'hit'}, [height(statarray1), 1]); ...
                  repmat({'miss'}, [height(statarray1), 1]); ...
                  repmat({'prem'}, [height(statarray1), 1])], ...
                  'VariableNames', {'Values', 'Identifiers'});

              ax_ofc = obj.bootstrap_plot(tmp_tbl_ofc, 'prem', 5000);
              ax_ofc2 = obj.bootstrap_plot(tmp_tbl_ofc, 'miss', 5000);
            else
                ax_ofc = axes(figure);
                ax_ofc2 = axes(figure);
            end


                        
              acc_info.trial_result(acc_info.plot_result ~= -1) = {'not-prem'};
              acc_info.trial_result(acc_info.plot_result == 0) = {'miss'};
              acc_info.trial_result(acc_info.plot_result == 1) = {'hit'};
              
              tmp_tbl_acc_not_prem = table(acc_bseline(acc_info.plot_result ~= -1), acc_info.(acc_info.Properties.VariableNames{end})(acc_info.plot_result ~= -1),  'VariableNames', {'Values', 'mouse'});
              tmp_tbl_acc_miss = table(acc_bseline(acc_info.plot_result == 0), acc_info.(acc_info.Properties.VariableNames{end})(acc_info.plot_result == 0),  'VariableNames', {'Values', 'mouse'});
              tmp_tbl_acc_corr = table(acc_bseline(acc_info.plot_result == 1), acc_info.(acc_info.Properties.VariableNames{end})(acc_info.plot_result == 1),  'VariableNames', {'Values', 'mouse'});
              tmp_tbl_acc_prem = table(acc_bseline(acc_info.first_lick < 0.5), acc_info.(acc_info.Properties.VariableNames{end})(acc_info.first_lick < 0.5),  'VariableNames', {'Values', 'mouse'});
              statarray = grpstats(tmp_tbl_acc_miss, 'mouse');
              statarray1 = grpstats(tmp_tbl_acc_corr, 'mouse');
              statarray2 = grpstats(tmp_tbl_acc_prem, 'mouse');
              statarray3 = grpstats(tmp_tbl_acc_not_prem, 'mouse');
              

              tmp_tbl_acc = table([statarray1{:,3}; statarray{:,3}; statarray2{:,3}], ...
                                  [repmat({'hit'}, [height(statarray1), 1]); ...
                                   repmat({'miss'}, [height(statarray1), 1]); ...
                                   repmat({'prem'}, [height(statarray1), 1])], ...
                                    'VariableNames', {'Values', 'Identifiers'});
                                      
              ax_acc = obj.bootstrap_plot(tmp_tbl_acc, 'prem', 5000);
              set([ax_ofc, ax_acc], 'YLim', [min(ax_acc.YLim(1), ax_ofc.YLim(1)), ...
                  max(ax_acc.YLim(2), ax_ofc.YLim(2))]);
              title(ax_ofc, 'OFCp prem')
              title(ax_acc, 'ACCp prem')
              
              
              ax_acc2 = obj.bootstrap_plot(tmp_tbl_acc, 'miss', 5000);
              set([ax_ofc2, ax_acc2], 'YLim', [min(ax_acc2.YLim(1), ax_ofc2.YLim(1)), ...
                  max(ax_acc2.YLim(2), ax_ofc2.YLim(2))]);
              title(ax_ofc2, 'OFCp miss')
              title(ax_acc2, 'ACCp miss')
        end
        
        
        function trial_onset_summary(obj)
            [acc_ch, ofc_ch, acc_info, ofc_info] = obj.combine_by_target('trials', 0);
            acc_ch(acc_info.first_lick < 1, :) = [];
            ofc_ch(ofc_info.first_lick < 1, :) = [];
            
            t_vec = linspace(-5, 15, size(acc_ch, 2));
            f = figure; 
            ax1 = subplot(1,2,1);
            shadedErrorBar(downsample(t_vec', 10)', downsample(acc_ch', 10)', {@mean, @sem}, 'k');
            ylabel('zscored \DeltaF/F')
            xlabel('Time relative to trial onset (s)')
            title('ACCp')
            ax2 = subplot(1,2,2);
            shadedErrorBar(downsample(t_vec', 10)', downsample(ofc_ch', 10)', {@mean, @sem}, 'k');
            ylabel('zscored \DeltaF/F')
            xlabel('Time relative to trial onset (s)')
            title('OFCp')
            set([ax1, ax2], 'XLim', [-0.1, 1], 'YLim', [-0.07, 0.2]);
        end
        
        
        function cno_graphs(obj, outcome)
            if nargin < 2
                outcome = -1; % default is impulsive errors
            end
            % For each mouse, combine info for saline / cno in right order
            for ii = 1:length(obj.MOUSE_ARRAY)/2
                array = obj.MOUSE_ARRAY(ii);
                nme = array.MOUSE_NM;
                array_cno = obj.MOUSE_ARRAY(strcmp([nme, 'C'], {obj.MOUSE_ARRAY.MOUSE_NM}));
                info_sal = array.INFO.t_info;
                info_sal.day(info_sal.day == 6) = 8;  
                info_sal.day(info_sal.day == 5) = 6;
                info_cno = array_cno.INFO.t_info;
                info_cno.day(info_cno.day == 3) = 9;
                info_cno.day(info_cno.day == 2) = 7;
                info_cno.day(info_cno.day == 1) = 5;
                
                info = [info_sal; info_cno];
                for tt = 1:9
                prem_rate(ii, tt) = numel(find(info.plot_result(info.day == tt) == -1)) / numel(find((info.day == tt)));
                miss_rate(ii, tt) = numel(find(abs(info.plot_result(info.day == tt)) ~= 1)) / numel(find((info.day == tt)));
                hit_rate(ii, tt) = numel(find((info.plot_result(info.day == tt)) == 1)) / numel(find((info.day == tt)));
                RT(ii, tt) =  mean(info.first_lick(info.plot_result == 1 & info.day == tt) - info.delay(info.plot_result == 1 & info.day == tt));
                end
            end
            switch outcome
                case -1
                    rate = prem_rate;
                case 1
                    rate = hit_rate;
                case 0
                    rate = miss_rate;
            end
            avg_rate = mean(rate(:, 1:3), 2);
            norm_rate = rate(:, 4:end) ./ avg_rate;
            stat_tbl = table(norm_rate(:), reshape(repmat([1:6], [5, 1]), [], 1), ...
                             repmat([1:5]', [6 ,1]), reshape(repmat({'sal', 'cno'}, [5, 3]), [], 1), ... 
                             'VariableNames', {'FA', 'day', 'mouse', 'cond'});
            a = axes(figure);
            bar(a, 1:2:6, mean(norm_rate(:, 1:2:6)), 'FaceColor', 'c', 'BarWidth', 0.25)
            hold(a, 'on')
            bar(a, 2:2:6, mean(norm_rate(:, 2:2:6)), 'FaceColor', [0.8, 0.8, 0.8], 'BarWidth', 0.25)
            plot(a, repmat(1:6, [5,1])', norm_rate', 'ko-') 
            a.XTick = 1:6;
            a.XTickLabel = {'Sal', 'CNO','Sal', 'CNO','Sal', 'CNO'};
            ylabel('Normalized rate')
            xlabel('Session Type')
            a.XLim = [0.2, 6.8];
            a.YLim(2) = 1.6;
            
            % plot cumulative graph without 3 Saline baseline days
            miss = sum(miss_rate(:, 4:end))./(360*5);
            hit = sum(hit_rate(:, 4:end))./(360*5);
            prem = sum(prem_rate(:, 4:end))./(360*5);
            a1 = axes(figure);
            bar(a1, (hit) + (miss) + (prem), 'FaceColor', 'r')
            hold(a1, 'on');
            bar(a1, (hit) + (miss), 'FaceColor', [0.5, 0.5, 0.5])
            bar(a1, (hit) ,'FaceColor', 'g')
            ylabel('Avg Rate');
            xlabel('Session Type');
            set(a1, 'XTickLabel', {'Sal', 'CNO'})
            
            % plot RT:
            a = axes(figure);
            boxplot(a, RT(:, 4:end))
            hold(a, 'on')
            plot(a, 1:6, RT(:, 4:end), 'ok')
            a.XTickLabel = {'Sal', 'CNO'};
            ylabel('RT(sec)')
            xlabel('Session Type')
            a.XLim = [0.2, 6.8];
            a.YLim(2) = 0.8;
        end

    end
    methods (Static = true)
        
        function mse_array = load_array()
            end_flg = false;
            mse_array = [];
            while ~end_flg
                X = MouseSummary;       % loads data, goes through it, get day coefficients based on GLM and movement data
                
                mse_array = [mse_array, X];
                answer = questdlg('Load another mouse?', '', 'Yes','No' ,'No');
                switch answer
                    case 'Yes'
                        end_flg = false;
                    case 'No'
                        end_flg = true;
                end
            end
        end
        
        function [nme, num_run, data_param, t_win, plot_type] = get_plot_settings(epoch)
            switch epoch
                case 'trials'
                    nme = 'bbn_res';
                    num_run = 1;
                    data_param = 'BBN';
                    t_win = [4000, 7500];
                    plot_type = 'stand alone';
                case 'cue'
                    nme = 'cue_res';
                    data_param = 'cue_int';
                    num_run = 4;
                    t_win = [4000, 9000];
                    plot_type = 'vs';
                case 'vis'
                    nme = 'cue_res';
                    data_param = 'is_light';
                    num_run = 2;
                    t_win = [4000, 9000];
                    plot_type = 'vs';
                case 'lick'
                    nme = 'lick_res';
                    num_run = 2;
                    data_param = 'plot_result';
                    t_win = [4000, 20000];
                    plot_type = 'vs'; 
                case 'mov'
                    nme = 'movement_res';
                    data_param = 'move';
                    num_run = 1;
                    t_win = [1000, 20000];
                    plot_type = 'stand alone';
                case 'passive'
                    nme = 'passive_res';
                    num_run = 2;
                    data_param = 'cond_num';
                    t_win = [1, 5086];
                    plot_type = 'vs';
                case 'baseline'
                    nme = 'baseline';
                    num_run = 3;
                    data_param = 'outcome';
                    t_win = [1, 5000];
                    plot_type = 'cross-mice';
                case 'cloud'                        % FOR PLOTTING PRE-CUE BASELINE
                    nme = 'cloud_res';
                    data_param = 'cloud';
                    num_run = 2;
                    t_win = [1, 5000];
                    plot_type = 'vs';
            end
            
        end
        
        function [info, data_idx, params] = get_trial_subset(info_struct, epoch, cond)
            if nargin < 3               % no cond parameter
                switch epoch
                    case 'trials'
                        info = info_struct.t_info(info_struct.t_info.plot_result ~= -1, :);
                        data_idx = info_struct.t_info.plot_result ~= -1;
                        info.BBN = ones(height(info), 1);
                        data_param = 'BBN';
                        params = unique(info.(data_param));
                    case 'cue'
                        info = info_struct.cue_info; 
                        data_idx = 1:height(info); % only no lick trials (also rule out misses with lick)
                        data_param = 'cue_int';
                        params = unique(info.(data_param));
                    case 'lick'
                        info = info_struct.l_info;
                        info.plot_result(info.plot_result ~= 1) = 0;        % treat all non-correct as non-rewearded (0)
                        data_idx = ~isnan(info_struct.l_info.first_lick);
                        data_param = 'plot_result';
                        params = unique(info.(data_param));
                    case 'vis'
                        info = info_struct.cue_info; 
                        data_idx = 1:height(info); % only no lick trials (also rule out misses with lick)
                        data_param = 'is_light';
                        params = unique(info.(data_param));
                    case 'mov'
                        info = info_struct.mov_info;
                        info.move = ones(height(info), 1);                     % dummy paramater to use
                        data_idx = ~isnan(info_struct.mov_info.closest_trial); % dummy condition, takes all move events
                        params = 1;
                    case 'cloud'
                        info = info_struct.cue_info;
                        data_idx = 1:height(info);
                        data_param = 'cloud';
                        params = unique(info.(data_param));
                end
            else                    % cond parameter exists
                switch epoch
                    case 'trials'
                        t_info = info_struct.t_info;
                        info = info_struct.t_info(info_struct.t_info.plot_result ~= -5 & eval(cond), :); 
                        info.BBN = ones(height(info), 1);                                           % dummy paramater to use
                        data_idx = find(info_struct.t_info.plot_result ~= -5 & eval(cond));
                        data_param = 'BBN';
                        params = unique(info.(data_param));
                    case 'cue'
                        cue_info = info_struct.cue_info;
                        info = info_struct.cue_info(eval(cond), :); 
                        data_idx = find(eval(cond));
                        data_param = 'cue_int';
                        params = unique(info.(data_param));
                    case 'lick'
                        l_info = info_struct.l_info;
                        info = info_struct.l_info(eval(cond), :);
                        info.plot_result(info.plot_result ~= 1) = 0;                            % treat all non-correct as non-rewearded (0)
                        data_idx = find(~isnan(info_struct.l_info.first_lick) & eval(cond));
                        data_param = 'plot_result';
                        params = unique(info.(data_param));
                    case 'vis'
                        cue_info = info_struct.cue_info;
                        info = info_struct.cue_info(eval(cond), :);
                        data_idx = find(eval(cond));
                        data_param = 'is_light';
                        params = unique(info.(data_param));
                    case 'mov'
                        mov_info = info_struct.mov_info;
                        info = info_struct.mov_info(eval(cond), :);
                        info.move = ones(height(info), 1);      % dummy paramater to use
                        data_idx = find(~isnan(info_struct.mov_info.closest_trial) & eval(cond)); % dummy condition, takes all move events
                        params = 1;
                    case 'cloud'
                        cue_info = info_struct.cue_info;
                        info = info_struct.cue_info(eval(cond), :);
                        data_idx = find(eval(cond));
                        data_param = 'cloud';
                        params = unique(info.(data_param));
                end
            end
        end
        function [p, tbl, stats] = do_statistics(epoch, acc_data, ofc_data, aud_data)
            % use [] instead of acc/ofc array to check only ofc/aud
            switch epoch
                case 'cue'
                    param_num = 4;
                case 'trials'
                    param_num = 1;
                case 'mov'
                    param_num = 1;
                otherwise
                    param_num = 2;
            end
            aT = table;
            acc_mice = size(acc_data, 1);
            aT.peaks = acc_data(:);
            param = repmat(1:param_num, [acc_mice, 1]);
            aT.param = param(:);
            
            if nargin > 2
                ofc_mice = size(ofc_data, 1);
                param = repmat(1:param_num, [ofc_mice, 1]);
                oT = table;
                oT.peaks = ofc_data(:);
                oT.param = param(:);
            end
            if nargin  > 3
                aud_mice = size(aud_data, 1);
                param = repmat(1:param_num, [aud_mice, 1]);
                dT = table;
                dT.peaks = aud_data(:);
                dT.param = param(:);
            else
                dT = table();
            end
            switch epoch
                case 'trials'
                    %[p,tbl,stats] = anova1(aT.peaks, oT.peaks, aT.peaks);
                    [~, p, tbl, stats] = ttest2(aT.peaks, oT.peaks); % two sample ttest
                case 'mov'
                    [~, p, tbl, stats] = ttest2(aT.peaks, oT.peaks); % two sample ttest
                case 'cue'
                    T = [aT; oT; dT];
                    [p,tbl,stats] = anova1(T.peaks, T.param);
                    %[p,tbl,stats] = kruskalwallis(T.peaks, T.param);
                    xlabel(epoch)
                otherwise
                    T = [aT; oT; dT];
                    [~, p, tbl, stats] = ttest(T.peaks(T.param == 1), T.peaks(T.param == 2)); % paired ttest %(use peaks for licks)
            end
        end
        function h_merge = merge_ind_plots(h_ind)
            % trim empty rows (from single channel mice)
            h_merge = h_ind(:, 1);
            if size(h_ind, 2) > 1
                for ii = 1:size(h_merge, 1)
                    for jj = 2:size(h_ind, 2)
                        clr = [0.8, 0.8, 0.8] - jj * [0.1, 0.1, 0.1];
                        lnes = findobj(h_ind(ii, jj), 'type', 'Line');
                        ptches = findobj(h_ind(ii, jj), 'type', 'Patch');
                        set(ptches, 'FaceColor', clr);
                        gos = [lnes; ptches];
                        copyobj(gos, findobj(h_merge(ii), 'type', 'axes'));
                        close(h_ind(ii, jj)); % close uncessecary plots
                    end
                end
            end
        end
        function [r, p, top_ten_idx, bottom_ten_idx, subs, probs] = base_outcome_corr(norm_r_base, t_info, binnum, data)
            [~, edges, bin] = histcounts(norm_r_base, prctile(norm_r_base, [1:100/binnum:100]));
            prem_prob = zeros(1, binnum);
            corr_prob = zeros(1, binnum);
            miss_prob = zeros(1, binnum);
            
            t_info.plot_result(t_info.plot_result == 2) = 0; % turn late into miss
            
            for ii = 1: binnum
                prem_prob(ii) = sum(t_info.plot_result(bin == ii - 1, :) == -1) / size(t_info.plot_result(bin == ii - 1, :) == -1, 1); % % to do premature for all trials in a bin.
                
                corr_prob(ii) = sum(t_info.plot_result(bin == ii - 1, :) == 1) / size(t_info.plot_result(bin == ii - 1, :) == 1, 1); % % to do premature for all trials in a bin.
                
                miss_prob(ii) = sum(t_info.plot_result(bin == ii - 1, :) == 0) / size(t_info.plot_result(bin == ii - 1, :) == 0, 1); % % to do premature for all trials in a bin.
                  
            
            end
            prem_prob(isnan(prem_prob)) = 0; %If there are NaNs, means there were no premature with this baseline, change prob to zero.
            corr_prob(isnan(corr_prob)) = 0;
            miss_prob(isnan(miss_prob)) = 0;

            max_prem = mean(prem_prob);
            max_corr = mean(corr_prob);
            max_miss = mean(miss_prob);
            
            prem_prob = prem_prob ./ max_prem;    % normalize rates
            corr_prob = corr_prob ./ max_corr;
            miss_prob = miss_prob ./ max_miss;
           
            margin_to_take = 25; %in percent
            perbin = margin_to_take / (100 / binnum);
            top_ten_idx = ismember(bin, [round((binnum - perbin)):binnum]);
            bottom_ten_idx = ismember(bin, [0: (perbin - 1)]);

            sub_prem_top = (sum(t_info.plot_result(top_ten_idx) == -1) / sum(top_ten_idx)) / (mean(prem_prob) * max_prem);       % prem prob for top baseline, normalized by mean overall prem_prob
            sub_prem_bot = (sum(t_info.plot_result(bottom_ten_idx) == -1) / sum(bottom_ten_idx)) / (mean(prem_prob) * max_prem); % prem prob for bottom baseline, normalized by mean overall prem_prob
            
            sub_corr_top = (sum(t_info.plot_result(top_ten_idx) == 1) / sum(top_ten_idx)) /(mean(corr_prob) * max_corr);
            sub_corr_bot = (sum(t_info.plot_result(bottom_ten_idx) == 1) / sum(bottom_ten_idx)) /(mean(corr_prob) * max_corr);
            
            sub_miss_top = (sum(t_info.plot_result(top_ten_idx) == 0) / sum(top_ten_idx)) / (mean(miss_prob) * max_miss);
            sub_miss_bot = (sum(t_info.plot_result(bottom_ten_idx) == 0) / sum(bottom_ten_idx)) / (mean(miss_prob) * max_miss);
            
            sub_prem = [sub_prem_bot - 1, sub_prem_top - 1];
            sub_cor = [sub_corr_bot - 1, sub_corr_top - 1];
            sub_miss = [sub_miss_bot - 1, sub_miss_top - 1];
            subs = [sub_prem; sub_cor; sub_miss]; % column 1 is low, 2 is high (prem, cor, miss)
            probs = [prem_prob; corr_prob; miss_prob];
            [r, p] = corrcoef(edges, prem_prob);
            r = r(2);
            p = p(2);
        end
        function plot_corr(data, clr, txt)
            if ~isempty(data)
                x = linspace(0, 1, size(data, 2));
                y = mean(data, 1);
                plot(x, y,  clr)
                if size(data, 1) > 1
                    errorbar(x, y, sem(data),clr(1), 'Color', clr(2))
                end
                ylabel(['norm. ', txt, ' Rate'])
                [bls,bint,r,rint,stats] = regress(y',[x',ones(length(x),1)]);
                
                mdl = fitlm(x,y);
                Xnew = linspace(min(x), max(x), 100)';
                [ytag,yci] = predict(mdl, Xnew);

                plot(Xnew, ytag, 'b-')
                plot(Xnew, yci, '--k')
                rsq = sprintf('%1.2f', stats(1));
                pval = sprintf('%1.5f', stats(3));
                eq = ['y = ', sprintf('%1.2f', bls(1)), ' x+ ', sprintf('%1.2f', bls(2))];
                corr_text = sprintf('%s\n R^2 = %s\n p = %s', eq, rsq, pval);
                text(0.9, 0.5, corr_text,'HorizontalAlignment','right')
            end
        end
        function set_plot_scales(handles, min_y, max_y)
            % set scales by min max of all plots, or by user defined
            % numbers
            handles = findall(handles, 'Type', 'Figure'); % remove dummy figure handles.
            for ii = 1:length(handles)
                tool_bar = findall(handles(ii), 'Type', 'AxesToolbar');
                tmp = findall(handles(ii), 'Type', 'Axes');
                h(ii) = tmp(~ismember(tmp, tool_bar));
                lim(ii, :) = h(ii).YLim;
            end
            
            for ii = 1:length(handles)
                figure(handles(ii))
                if nargin == 1
                    h(ii).YLim = [min(lim(:, 1)), max(lim(:, 2))];
                elseif nargin == 2
                    h(ii).YLim = [min_y, max(lim(:, 2))];
                elseif nargin == 3
                    if ~isempty(min_y)
                        h(ii).YLim = [min_y, max_y];
                    else
                        h(ii).YLim = [min(lim(:, 1)), max_y];
                    end
                end
            end
            addToolbarExplorationButtons(handles)       % bring back the buttons
        end
        function fig_struct = copy_plot(h, dims, pos)
            fig_struct = struct();
            fig_struct.fig = figure;
            fig_struct.sub_h(pos) = subplot(dims(1), dims(2), pos);
            hold(fig_struct.sub_h(pos), 'on');
            labels = findall(h, 'Type', 'Text');
            
            if isa(h, 'matlab.ui.Figure')
                tmp = h.Children.Children;
            else
                tmp = h.Children;
            end
            copyobj(tmp, fig_struct.sub_h(pos));
            copyobj(labels, fig_struct.sub_h(pos));
        end
        function sm = smooth_sig(sig, ds_fac)
            sm = arrayfun(@(i) mean(sig(:, i:i+ds_fac-1)'),1:ds_fac:size(sig, 2)-ds_fac+1, 'UniformOutput', false);
            sm = cell2mat(sm);
            sm = reshape(sm, size(sig, 1), []);
        end

        function h_new = create_RT_summary(hf)
            % organizes grouped RT figure
            clrs = [140, 179, 250; 62,98,83; 129, 88, 135] ./255;
            nmes = {'chill', 'middle', 'stressed'};
            lick_hists = findall(hf, 'Type', 'Histogram', 'FaceColor', [1, 0, 0]);
            delay_hists = findall(hf, 'Type', 'Histogram', 'FaceColor', [1, 0, 1]);
            dots = findall(hf, 'Type', 'Scatter');
            for ii = 1:numel(dots)
                set(dots(ii), 'MarkerFaceColor', clrs(ii, :), 'MarkerFaceAlpha', 0.3, 'SizeData', 1)
                set(lick_hists(ii), 'FaceColor', clrs(ii, :), 'Normalization', 'probability') % can use cdf to look at what percent of licks occur before time t
                set(delay_hists(ii), 'FaceColor', clrs(ii, :))
            end
            h_new = figure;
            % plot lick histograms
            ax1 = axes(h_new);
            hold(ax1, 'on')
            
            copyobj(lick_hists, ax1)
            hl = legend(ax1, ax1.Children, nmes, 'AutoUpdate', 'off', 'Orientation', 'horizontal');
            xlabel(ax1, 'RT (sec)'); ylabel(ax1, 'Trial Count'); title(ax1, 'Premature RT by Group');
            for ii = 1:numel(dots)
                lh(ii) = line(ax1, [median(lick_hists(ii).Data), median(lick_hists(ii).Data)], [ax1.YLim(1), ax1.YLim(2)], ...
                    'Color', clrs(ii, :), 'LineWidth', 1, 'LineStyle', '--'); %plot mean line
            end
            uistack(lh, 'bottom')
            uistack(hl, 'top')
            
            % plot delay subhistograms
            ax2 = axes('Position', [0.5, 0.2095, 0.3250, 0.2096]);
            hold(ax2, 'on')
            copyobj(delay_hists, ax2)
            bin_cts = {delay_hists.BinCounts};
            bin_edges = delay_hists(1).BinEdges(2:end);
            dum1 = {};
            dum2 = {};
            for ii = 1:numel(dots)
                dum1{ii} = bin_edges;
                dum2{ii} = 1;
            end
            ps = cellfun(@polyfit, dum1, bin_cts, dum2, 'UniformOutput', false); % get linear fits
            yhats =  cellfun(@polyval, ps, dum1, 'UniformOutput', false);  % get y estimates
            for ii = 1:numel(dots)
                plot(ax2, bin_edges, yhats{ii}, '-', 'Color', clrs(ii, :), 'LineWidth', 2);        % plot slopes
                text(ax2, ax2.XLim(2), yhats{ii}(end), sprintf('%1.2f', ps{ii}(1)), 'Color', clrs(ii, :))
            end
            xlabel(ax2, 'Delay (sec)'); ylabel(ax2, 'Trial Count'); title(ax2, 'Premature Delay dist.');
            text(ax2, ax2.XLim(2), ax2.Title.Position(2), 'Slope:', 'Color', 'k', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'FontSize', ax2.Title.FontSize)
            
            % plot lick scatters
            for ii = 1:numel(dots)
                ax3(ii) = axes('Position', [ax2.Position(1) + ax2.Position(3)/2.65 * (ii-1) , ax2.Position(2)*2.65, ax2.Position(3)/4, ax2.Position(4)]);
                copyobj(dots(ii), ax3(ii));
                xlabel('RT')
                ylabel('Delay')
                if ii ==2
                    title('Premature Licks')
                end
            end
            hl.Position = [0.49, 0.86, 0.35, 0.036];
            h_new.Position = [474,   503,   766,   595];
            close(hf)
        end
        
        function groups = group_members(groupby)
            switch groupby                
                case 'prem_2_groups'
                    exploiters = {'7_from404', '5_from404',  '2_from406', '6_from406', ...
                        '3_from406', '1_from406', '2_from500', '2_from440', '5_from400'};
                    
                    explorers = {'2_from430', '1_from440','3_from500', '3_from410', '1_from500', '3_from430',...
                        '1_from404', '1_from400', '3_from404', '6_from404',...
                        '2_from404', '4_from410', '4_from410L', '4_from404', '3_from440','8_from402'};
                    
                    groups = {exploiters,explorers};
                case 'prem_2_groups_no_ofc'
                                        exploiters = {'2_from406', '6_from406', ...
                        '3_from406', '1_from406', '2_from500', '2_from440'};
                    
                    explorers = {'2_from430', '1_from440','3_from500', '3_from410', '1_from500', '3_from430',...
                        '1_from404', '1_from400', '3_from404',...
                        '2_from404', '4_from410', '4_from410L', '4_from404', '3_from440'};
                    
                    groups = {exploiters,explorers};
                    
                case 'hm3dq'
                    %CNO
                    hm3dq_sal = {'4_from700', '7_from700', '9_from700', '1_from800', '4_from800'};
                    hm3dq_cno = {'4_from700C', '7_from700C', '9_from700C', '1_from800C', '4_from800C'};
                    groups = {hm3dq_sal, hm3dq_cno};
            end
        end
        
        function [per, prob, period_fft] = shuffle_baseline(shuff_num, sub_sig, sm_win, f_ax)
            sub_ax1_a = subplot(2, 3, 4, axes(f_ax));
            hold(sub_ax1_a, 'on')
            for ii = 1:size(sub_sig, 2)
                sess_sig = sub_sig(:, ii);
                sm_sess_sig = smooth(sess_sig, sm_win) - mean(smooth(sess_sig, sm_win));     % smooth and zero mean per session
                [period(ii)] = fit_sine(1:size(sess_sig, 1), sm_sess_sig'); % fit sine
                [f, fft_data(ii, :)] = quick_fft([], sess_sig, 1, 0);                        % fft session
                sm_shuff_sig = [];
                for tt = 1:shuff_num
                    shuff_sig = sess_sig(randperm(length(sess_sig)));       % randomized session
                    sm_shuff_sig(:, tt) = smooth(shuff_sig, sm_win) - mean(smooth(shuff_sig, sm_win)); % smooth and zero mean per session
                    [period_shuff(ii, tt)] = fit_sine(1:size(sess_sig, 1), sm_shuff_sig(:, tt)');   % fit sine
                    [~, fft_data_shuff(tt,:)] = quick_fft([], shuff_sig, 1, 0);              % fft shuff
                end
                [perm_c(ii, :), ~] = xcorr(shuff_sig, size(sub_sig, 1)/2, 'coeff');   % cross-corr non smoothed
                [c(ii, :), lags]= xcorr(sess_sig, size(sub_sig, 1)/2, 'coeff');   % cross-corr
                fft_shuff_sess(ii, :) = max(fft_data_shuff);                
            end
            per = mean(period);
            histogram(sub_ax1_a, period_shuff(:), shuff_num, 'Normalization', 'pdf')
            line(sub_ax1_a, [per, per], [0, sub_ax1_a.YLim(2)], 'Color', 'k')
            
            pd = fitdist(period_shuff(:), 'Lognormal');
            prob = 1 - cdf(pd, per);
            plot(sub_ax1_a, pdf(pd, 0:max(period_shuff(:))))
            text(sub_ax1_a, max(period_shuff(:))/2, sub_ax1_a.YLim(2)/2, sprintf('period: %1.2f\nprob: %1.4f', per, prob))
            
            sub_ax1_b = subplot(2, 3, 5, axes(f_ax)); 
            hold(sub_ax1_b, 'on')
            plot(sub_ax1_b, lags, mean(c))
            plot(sub_ax1_b, lags, mean(perm_c))
            legend('smoothed baseline', 'shuffled baseline')
            
            sub_ax1_c = subplot(2, 3, 6, axes(f_ax)); 
            hold(sub_ax1_c, 'on')
            avg_fft = mean(fft_data);
            avg_shuff_fft = mean(fft_shuff_sess);
            plot(f, avg_fft)
            plot(f, avg_shuff_fft)
            
            tresh = median(avg_shuff_fft);
            if sum(avg_fft>tresh) == 0
                period_fft = 0;
            else
                period_fft = f(avg_fft == max(avg_fft(avg_fft>tresh))); % peak of significant power 0 if none (but gives most mice same number)
            end
            text(sub_ax1_c, 0.01, sub_ax1_c.YLim(2)/2, sprintf('period: %1.5f', period_fft))
            
            legend('baseline', 'shuffled')
            xlim([0,0.01]);
            
            function [f, pspec] = quick_fft(ax, sig, Fs, plot_flg)
                if nargin < 3
                    ax = axes(figure);
                    hold(ax, 'on')
                end
                %L = length(sig);% Length of signal
                Y = fft(sig, 1000);
                L = length(Y);
                P2 = abs(Y/L);
                pspec = P2(1:L/2+1);
                pspec(2:end-1) = 2*pspec(2:end-1);
                f = Fs*(0:(L/2))/L;
                f = f./20;          % trial into seconds
                if plot_flg
                    plot(ax, f, pspec)
                    title('Single-Sided Amplitude Spectrum of X(t)')
                    xlabel('f (Hz)')
                    ylabel('|P1(f)|')
                end
            end
            function [per] = fit_sine(x, sig)
                yu = max(sig);
                yl = min(sig);
                yr = (yu-yl);                               % Range of �y�
                yz = sig-yu+(yr/2);
                zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
                per = 2*mean(diff(zx));                     % Estimate period

            end
        end
        function ax = bootstrap_plot(data, nme, bstrp_num)
            prem_data = data.Values(strcmp(data.Identifiers, nme));
            ax = axes(figure); 
            hold(ax, 'on');
            prem_btstrp = bootstrp(bstrp_num, @mean, prem_data);
            plot(ax, -0.2, prem_data, 'ok')
            line(ax, [-0.2, 0.2], [0, 0], 'Color', 'k', 'LineStyle', '--')
            line(ax, [-0.2, 0.2], [mean(prem_data), mean(prem_data)], 'Color', 'k', 'LineStyle', '--')
            histogram(ax, prem_btstrp, 'Orientation', 'horizontal', 'FaceColor', 'k', 'Normalization', 'probability', 'LineStyle', 'none')
            plot(ax, 0, median(prem_btstrp), 'ok', 'MarkerFaceColor', 'k')
            line(ax, [0, 0], [prctile(prem_btstrp, 5), prctile(prem_btstrp, 95)], 'Color', 'k', 'LineWidth', 2)
            ax.XLim = [-0.25, 0.25];
            
        end
    end
end

