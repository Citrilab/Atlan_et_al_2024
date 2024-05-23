classdef MouseSummary < handle
    
    properties
        data = struct();
        MOUSE_NM
        PROJ_TO
        INFO
        PATH
        MATFILE
        DATA_CH
    end
    
    methods
        function obj = MouseSummary ()
            % This path needs to match the location of your data:
            obj.PATH = 'D:\Atlan2021 data\';
            
            obj.MOUSE_NM = char(inputdlg('Enter mouse name:'));              % mouse name;
            obj.PROJ_TO = inputdlg({'Enter ID of Primary Channel', ...
                'Enter ID of Secondary Channel'}, ...
                'Set Mouse Type: ACC/OFC/AUD', [1,35]);
            obj.DATA_CH = 1 + 1 * ~sum(cellfun(@isempty, obj.PROJ_TO));        % 1 if once channel left empty (assumed secondary) 2 otherwise.
            obj.get_data();
        end
        
        function get_data(obj)
            % Loads info files into the INFO property
            path = obj.PATH;                                                 % tmp variable names
            mse_nm = obj.MOUSE_NM;                                           % tmp variable names
            Files = dir(fullfile([path, mse_nm], 'Cue*.mat'));               % get file names
            
            for ii = 1:length({Files.name})
                a_type = Files(ii).name;
                us = strfind(a_type, '_');
                X = matfile([path, mse_nm, '\', Files(ii).name]);
                
                switch(a_type(us(2) + 1:end-4))
                    case 't_onset'
                        obj.INFO.t_info = X.t_info;
                        if ismember('lick_trialid', fieldnames(X))
                            obj.INFO.lick_trialid = X.lick_trialid;
                        end
                        obj.MATFILE.trials = X;
                    case 'lick'
                        obj.MATFILE.lick = X;
                    case 'cloud'
                        obj.MATFILE.cloud = X;
                    case 'cue'
                        obj.INFO.cue_info = X.t_info;
                        obj.INFO.cue_info(obj.INFO.cue_info.plot_result==-1, :) =[]; % info excluding premature trials
                        obj.MATFILE.cue = X;
                    case 'movement'
                        obj.INFO.mov_info = X.mov_info;
                        obj.INFO.AF = X.AF_var;
                        obj.MATFILE.mov = X;
                end
            end
            obj.data.estimates(: ,1) = obj.get_day_difference(1);
            if obj.DATA_CH > 1
                obj.data.estimates(: ,2) = obj.get_day_difference(0);
            end
            obj.get_movement_data();
            obj.get_state_tag_by_outcome();
        end

        function estimates = get_day_difference(obj, flg, ext)
            % function GET_DAY_DIFFERENCE
            %
            % Inputs
            %   flg - 1 will normalize by all_trials, otherwise by af_trials
            %   ext - used to normalized to external data set (struct:
            %   ext.data & ext.info
            if nargin > 2 % normalize to external data
                obj.INFO.t_info = [ext.info; obj.INFO.t_info];
                 if flg == 1
                    trials = [ext.data; obj.MATFILE.trials.all_trials];
                else
                    trials = [ext.data; obj.MATFILE.trials.af_trials];
                end
            else
                if flg == 1
                    trials = obj.MATFILE.trials.all_trials;
                else
                    trials = obj.MATFILE.trials.af_trials;
                end
            end
            obj.INFO.l_info = obj.INFO.t_info(~isnan(obj.INFO.t_info.first_lick), :); %lick trials (including omissions)
            obj.INFO.cue_info = obj.INFO.t_info(obj.INFO.t_info.plot_result >= 0, :);
            
            t_breaks = find(obj.INFO.t_info.trial_number == 1);
            t_breaks = [t_breaks; size(obj.INFO.t_info, 1) + 1];
            
            l_breaks = find(obj.INFO.l_info.trial_number(2:end) < obj.INFO.l_info.trial_number (1:end -1));
            l_breaks = [1; l_breaks + 1; size(obj.INFO.l_info, 1) + 1];
            
            r_day = [];
            l_day = [];
            
            for ii = 1:(length(t_breaks) - 1)
                r_day(t_breaks(ii):t_breaks(ii + 1) - 1) = ii;      % tag each recording day
            end
            
            for ii = 1:(length(l_breaks) - 1)
                l_day(l_breaks(ii):l_breaks(ii + 1) -1) = ii;
            end
            
            r_day = categorical(r_day');
            r_outcome = categorical(obj.INFO.t_info.trial_result);           % tag outcome
            
            r_base = double(mean(trials(:, 1000:5000), 2));                  % 4 seconds preceeding T_onset
            
            r_set = table(r_base, r_day, r_outcome, 'VariableNames', {'baseline', 'day', 'outcome'});
            G = fitlme(r_set, 'baseline ~ outcome + day');
            estimates = G.Coefficients.Estimate;
            obj.INFO.t_info.day = double(r_day);
            obj.INFO.cue_info.day = double(r_day(obj.INFO.t_info.plot_result >= 0));
            obj.INFO.l_info.day = double(l_day)';
            obj.data.t_breaks = t_breaks;
            for ii = 1:double(r_day(end))
                obj.INFO.l_info.trial_number(l_day == ii) =  obj.INFO.l_info.trial_number(l_day == ii) + t_breaks(ii) - 1;         %get corresponding trial numbers
            end
            obj.data.anova = anova(G);
            if nargin > 2
                estimates(1) = []; %remove correction for dummy day since it will be removed from data
                obj.INFO.t_info(obj.INFO.t_info.day == 1, :) = []; %delete first day from data
                obj.INFO.cue_info(obj.INFO.t_info.day == 1, :) = []; %delete first day from data
                obj.INFO.l_info(obj.INFO.t_info.day == 1, :) = []; %delete first day from data
                obj.data.t_breaks = []; %delete first day from data
                obj.INFO.t_info.day = obj.INFO.t_info.day - 1; %reset days
                obj.INFO.l_info.day = obj.INFO.l_info.day - 1; %reset days
                obj.INFO.cue_info.day = obj.INFO.cue_info.day - 1; %reset days
            end
        end
        function get_movement_data(obj)
            % get movement data
            m_day = [];
            m_breaks = find(obj.INFO.mov_info.closest_trial(2:end) < obj.INFO.mov_info.closest_trial(1:end -1));
            m_breaks = [1; m_breaks + 1; size(obj.INFO.mov_info, 1) + 1];
            for ii = 1:(length(m_breaks) - 1)
                if strcmp(obj.MOUSE_NM, '1_from800') && ii > 1
                    m_day(m_breaks(ii):m_breaks(ii + 1) -1) = ii + 1;
                else
                    m_day(m_breaks(ii):m_breaks(ii + 1) -1) = ii;
                end
            end
            movement_thresh = 5;                % how much around BBN should a movement occur to be considered a movement trial
            mov_trials = obj.INFO.mov_info.closest_trial;
            for jj = 2:length(m_breaks) - 1
                if strcmp(obj.MOUSE_NM, '1_from800') && jj > 1
                    mov_trials(m_breaks(jj): m_breaks(jj + 1) - 1) = mov_trials(m_breaks(jj): m_breaks(jj + 1) - 1) + obj.data.t_breaks(jj+1) - 1;     % 1_800 second session has no movements - this takes care of it
                else
                    mov_trials(m_breaks(jj): m_breaks(jj + 1) - 1) = mov_trials(m_breaks(jj): m_breaks(jj + 1) - 1) + obj.data.t_breaks(jj) - 1;     % each closest trial becomes relative to all trials numbers
                end
            end
            mov_trials(obj.INFO.mov_info.distance_from_trial < -movement_thresh ...
                | obj.INFO.mov_info.distance_from_trial > movement_thresh) = [];
            obj.data.mov_trials = mov_trials;       % keep move trial idxes (with duplicates)
            
            mov_trials = unique(mov_trials);        % get rid of trials labeled twice
            
            obj.INFO.t_info.move = zeros(height(obj.INFO.t_info), 1); % add column of move or not
            obj.INFO.t_info.move(mov_trials) = 1;
            
            obj.INFO.cue_info.move = obj.INFO.t_info.move(obj.INFO.t_info.plot_result ~= -1,:);
            obj.INFO.l_info.move = obj.INFO.t_info.move(~isnan(obj.INFO.t_info.first_lick));
            obj.INFO.mov_info.day = m_day';
        end
        function get_state_tag_by_outcome(obj)
            outcome = obj.INFO.t_info.plot_result;
            outcome(outcome == -1) = -5;
            outcome(outcome == 0) = 5;
            outcome(outcome == 2) = 2;
            avg_outcome = movmean(outcome, 10);
            T = zeros(height(obj.INFO.t_info), 1);
            T(avg_outcome < -1) = -1;
            T(avg_outcome > 1) = 1;
            obj.INFO.t_info.Tag = T;
            obj.INFO.l_info.Tag = obj.INFO.t_info.Tag(~isnan(obj.INFO.t_info.first_lick));
            obj.INFO.cue_info.Tag = obj.INFO.t_info.Tag(obj.INFO.t_info.plot_result ~= -1);      
        end
        function get_state_tag_by_baseline(obj, stds)
            % get state tag from baseline of data (or of both channels if there are two)
            data = obj.normalize_data('trials', 1);
            norm_r_base = double(mean(data(:, 1000:5000), 2));
            z_base = norm_r_base - mean(norm_r_base);           % normalize around zero
            T = zeros(size(z_base, 1), 1);
            T(z_base < -std(z_base)* stds) = -1;
            T(z_base > std(z_base) * stds) = 1;
            nme = 'activity_tag_ch1';
            obj.INFO.t_info.(nme) = T;
            obj.INFO.l_info.(nme) = T(~isnan(obj.INFO.t_info.first_lick));
            obj.INFO.cue_info.(nme) = T(obj.INFO.t_info.plot_result ~= -1);
            obj.data.ch1_baseline = norm_r_base;        % store mean subtracted baseline
            
            if obj.DATA_CH > 1
                data = obj.normalize_data('trials', 2);
                norm_r_base = double(mean(data(:, 1000:5000), 2));
                z_base = norm_r_base - mean(norm_r_base);           % normalize around zero
                T = zeros(size(z_base, 1), 1);
                T(z_base < -std(z_base)* stds) = -1;
                T(z_base > std(z_base) * stds) = 1;
                nme = 'activity_tag_ch2';
                obj.INFO.t_info.(nme) = T;
                obj.INFO.l_info.(nme) = T(~isnan(obj.INFO.t_info.first_lick));
                obj.INFO.cue_info.(nme) = T(obj.INFO.t_info.plot_result ~= -1);
                obj.data.ch2_baseline = norm_r_base; % store mean subtracted baseline
            end
        end
        function plot_outcome_dist(obj)
            pr = obj.INFO.t_info.trial_number(obj.INFO.t_info.plot_result == -1);
            om = obj.INFO.t_info.trial_number(obj.INFO.t_info.plot_result == 0 | obj.INFO.t_info.plot_result == 2);
            cr = obj.INFO.t_info.trial_number(obj.INFO.t_info.plot_result == 1);
            h_hist = axes(figure);
            hold (h_hist, 'on')
            histfit(cr, 200, 'kernel')
            histfit(om, 200, 'kernel')
            histfit(pr, 200, 'kernel')
            line_h = flipud(findall(h_hist, 'Type', 'Line'));          % order is reversed (last one drawn is first)
            bar_h = flipud(findall(h_hist, 'Type', 'Bar'));
            legend(h_hist, line_h, {'correct', 'miss', 'premature'})
            line_h(1).Color = 'g'; line_h(2).Color = [0.4, 0.4, 0.4]; line_h(3).Color = 'r';
            bar_h(1).FaceColor = [0.5, 0.8, 0.5]; bar_h(1).FaceAlpha = 0.9;
            bar_h(2).FaceColor = [0.8, 0.8, 0.8]; bar_h(2).FaceAlpha = 0.9;
            bar_h(3).FaceColor = [0.64, 0.08, 0.1]; bar_h(3).FaceAlpha = 0.9;
        end
        
        function norm_data = normalize_data(obj, data_name, ch_num)
            % function NORMALIZE_DATA uses the day_estimates to normalize
            % data by day
            if ch_num == 1
                ch = 'all_trials';
            else
                ch = 'af_trials';
            end
            
            switch data_name
                case 'cloud'
                    info = obj.INFO.t_info;
                case 'cue'
                    info = obj.INFO.cue_info;
                case 'lick'
                    info = obj.INFO.l_info;
                case 'mov'
                    info = obj.INFO.mov_info;
                case 'trials'
                    info = obj.INFO.t_info;
                case 'vis'
                    info = obj.INFO.cue_info;
                    data_name = 'cue';
                case 'passive'
                    info = obj.INFO.passive_info;
            end
            
            trials = obj.MATFILE.(data_name).(ch);
            if ~strcmp(data_name, 'passive')
                norm_data = trials - obj.data.estimates(1, ch_num); % subtract intercept (1st day / correct)
                %norm_data = trials;        %if normalize by day only 
                for ii = 2:length(unique(info.day)) 
                    day_es = obj.data.estimates(ii, ch_num);
                    norm_data(info.day == ii, :) = norm_data(info.day == ii, :) - day_es; % remove each days' intercept: Use - mean(mean(norm_data(info.day == ii, 1000:5000), 2)); if by day only%
                end
            else
                norm_data = trials - mean(trials(:, 1:1000), 2);   % subtract first second
            end
        end
        
        function [srted_idx, srted_info] = sort_info_by(obj, data_name, sortby)            
            switch data_name
                case 'cloud'
                    info = obj.INFO.t_info;
                case 'cue'
                    info = obj.INFO.cue_info;
                case 'lick'
                    info = obj.INFO.l_info;
                case 'mov'
                    info = obj.INFO.mov_info;
                case 'trials'
                    info = obj.INFO.t_info;
            end
                
            [~, d_idx] = sort(info.(sortby));
            srted_idx = d_idx;
            srted_info = info(d_idx, :);
        end
        
        function new_vec = fs_transfer (obj, vec, fs)
            if nargin < 3
                fs = obj.TRACE.streams.(obj.STREAM_STORE2).fs;
            end
            new_vec = round(vec * fs + 1);
        end
        
        function [envelope_trials, other_trials, seq_trials] = id_hypovig(obj, seq_len)
            info = obj.INFO.t_info;
            info.plot_result(info.plot_result == 2) = 0; % count late as miss
            len = zeros(1, seq_len);    % sequence of omissions
            seq_idx = strfind(info.plot_result', len); % idxs of trials in sequence
            seq_idx(seq_idx < 2 * seq_len | seq_idx > height(info) - seq_len * 2) = []; % remove idxes which are too early or too late
            start_idx = seq_idx - 1;
            start_idx(ismember(start_idx, seq_idx)) = [];
            end_idx = seq_idx + 1;
            end_idx(ismember(end_idx, seq_idx)) = [];
            end_idx = end_idx + (seq_len - 1);
            
            envelope_trials = [];
            seq_trials = [];
            for ii = 1:numel(start_idx)
                envelope_trials = [envelope_trials, (start_idx(ii) - seq_len):start_idx(ii)];
                seq_trials = [seq_trials, (start_idx(ii) + 1):(end_idx(ii) -1)]; % there can be some overlap between envelope and seq trials                
            end
            other_trials = 1:height(info);
            other_trials(union(envelope_trials, seq_trials)) = [];
        end
        
        function sum_outcme = analyze_behavior(obj, full)
            t_info = obj.INFO.t_info;
            no_p_res = unique(t_info.plot_result);
            no_p_res(no_p_res == -1) = [];
            cases = {'light on', 'light off'};
            cue_ints = unique(t_info.cue_int);
            if full
                for ii = 1:length(no_p_res)
                    switch num2str(no_p_res(ii))
                        case '0'
                            cond_name = 'omitted';
                        case '1'
                            cond_name = 'correct';
                        case '2'
                            cond_name = 'late';
                    end
                    
                    tot_behav_res_ax(ii) = axes(figure('Name', ['Overall %', cond_name]));
                    tot_behav_res_ax(ii).YLim = [0, 1];
                    hold(tot_behav_res_ax(ii), 'on');
                    h_tot_leg = legend(tot_behav_res_ax(ii), '-DynamicLegend', 'Location', 'best');
                    for tt = 1:2
                        trial_ns = find(t_info.plot_result == no_p_res(ii));
                        for jj = 1:length(cue_ints)
                            cue_locs{jj} = find(t_info.cue_int == cue_ints(jj));     % should be same number
                        end
                        h1 = obj.per_graphs(tot_behav_res_ax(ii), intersect(trial_ns, find(t_info.is_light == (tt == 1))), intersect(find(t_info.cloud > 0), intersect(find(t_info.plot_result ~= -1), find(t_info.is_light == (tt == 1)))), cue_locs);
                        h2 = obj.per_graphs(tot_behav_res_ax(ii), intersect(trial_ns, find(t_info.is_light == (tt == 1))), intersect(find(t_info.cloud == 0), intersect(find(t_info.plot_result ~= -1), find(t_info.is_light == (tt == 1)))), cue_locs);
                        h1.HandleVisibility = 'off';
                        h2.HandleVisibility = 'off';
                        h1.Color = [0.85,0.33,0.10];
                        h2.Color = [0.00,0.35,0.74];
                        l_style = {'-', '--'};
                        
                        if no_p_res(ii) == 1 || no_p_res(ii) == 0
                            tmp1c = h1.Color;
                            tmp2c = h2.Color;
                            h1.LineStyle = 'none';
                            h2.LineStyle = 'none';
                            
                            ydata1 = h1.YData;
                            ydata2 = h2.YData;
                            ydata1(isnan(ydata1)) = 0;
                            ydata2(isnan(ydata2)) = 0; % replace NaNs with 0
                            
                            [ffit, curve] = FitPsycheCurveWH(h1.XData, ydata1, 1);
                            plot(tot_behav_res_ax(ii), curve(:, 1), curve(:, 2), 'Color', tmp1c, 'LineStyle', l_style{tt}, 'DisplayName', ['cloud - ', cases{tt}])
                            if no_p_res == 1
                                cl_fit{tt, 1} = ffit; %correct cloud (light and no light)
                            end
                            [ffit, curve] = FitPsycheCurveWH(h1.XData, ydata2, 1);
                            plot(tot_behav_res_ax(ii), curve(:, 1), curve(:, 2), 'Color', tmp2c, 'LineStyle', l_style{tt}, 'DisplayName', ['no cloud - ', cases{tt}])
                            if no_p_res(ii) == 1
                                nocl_fit{tt, 1} = ffit; %correct no cloud (light and no light)
                            end
                        end
                    end
                    tot_behav_res_ax(ii).YLim = [0, 1];
                    tot_behav_res_ax(ii).XTick = 1:numel(cue_ints);
                    tot_behav_res_ax(ii).XTickLabel = cue_ints;
                    tot_behav_res_ax(ii).Title.String = ['% ', cond_name, obj.MOUSE_NM];
                    xlabel(tot_behav_res_ax(ii), 'cue intensity');
                    ylabel(tot_behav_res_ax(ii), '% of trials');
                end
            else
                cond_name = 'correct';
                no_p_res = 1;
                sum_outcme = [];
                tot_behav_res_ax = axes(figure('Name', ['Overall %', cond_name]));
                tot_behav_res_ax.YLim = [0, 1];
                hold(tot_behav_res_ax, 'on');
                h_tot_leg = legend(tot_behav_res_ax, '-DynamicLegend', 'Location', 'best');
                for tt = 1:2
                    trial_ns = find(t_info.plot_result == no_p_res);
                    for jj = 1:length(cue_ints)
                        cue_locs{jj} = find(t_info.cue_int == cue_ints(jj));     % should be same number
                    end
                    h1 = obj.per_graphs(tot_behav_res_ax, intersect(trial_ns, find(t_info.is_light == (tt == 1))), intersect(find(t_info.cloud > 0), intersect(find(t_info.plot_result ~= -1), find(t_info.is_light == (tt == 1)))), cue_locs);
                    h2 = obj.per_graphs(tot_behav_res_ax, intersect(trial_ns, find(t_info.is_light == (tt == 1))), intersect(find(t_info.cloud == 0), intersect(find(t_info.plot_result ~= -1), find(t_info.is_light == (tt == 1)))), cue_locs);
                    
                    sum_outcme = [sum_outcme, h1.YData', h2.YData'];
                end
                close(tot_behav_res_ax.Parent)
                
            end
        end
        
        function idx = get_session_subset(obj, trial_num, nme)
            % gets indices for trials by session where trials num is the
            % number of trials per session, sessions with less trials are
            % excluded
            
            info = obj.INFO.(nme);
            idx = [];
            for ii = 1:numel(unique(info.day))
                idx_to_add = find(info.day == ii);
                if length(idx_to_add) < trial_num
                    idx_to_add =[];
                elseif length(idx_to_add) > trial_num
                    idx_to_add(trial_num+1:end) = [];
                end
                idx = [idx, idx_to_add];
            end
        end
        
        function auto_corr = create_autocorr(obj, proj_target, smooth_fac, max_lag, ds_fac)
            auto_corr = nan; % if no ACC channel keep nan
            for ii = 1:obj.DATA_CH
                if strcmp(obj.PROJ_TO{ii}, proj_target)
                    sig = obj.normalize_data('trials', ii);
                    sig = reshape(sig', 1, []);
                    sig = smooth(sig', smooth_fac)';
                    sig = downsample(sig', ds_fac)';
                    auto_corr = xcorr(sig, max_lag*(1000/ds_fac), 'coeff');
                end
            end
        end
        
        function [lick_rise, mov_rise] = mov_lick_comp(obj, plot_flg)
            %compares 75% rise time between lick and movement, plots if flg
            lick_rise = nan(2,1);
            mov_rise = nan(2,1);
            pure_lick_idx = find(obj.INFO.l_info.move == 0);
            info = obj.INFO.l_info(pure_lick_idx, :);
            pure_mov_idx = find(cellfun(@isempty, (cellfun(@(x) find(-5<x & x<5), obj.INFO.mov_info.lick_times, 'UniformOutput', false)))); % locomotion events with no licks 5 s around them
            mov_info = obj.INFO.mov_info(pure_mov_idx, :);
            
            
            for ii = 1:obj.DATA_CH
                lick_data = double(obj.normalize_data('lick', ii));
                lick_data = downsample(lick_data(pure_lick_idx, :)', 100)';
                mov_data = double(obj.normalize_data('mov', ii));
                mov_data = downsample(mov_data(pure_mov_idx, :)', 100)';
                t_vec = linspace(-5,15,size(lick_data, 2));
                
                [low_lick, high_lick, lick_rise(ii)] = obj.get_perctiles(mean(lick_data(:, 1:70)), [0.2, 0.95]);
                [low_mov, high_mov, mov_rise(ii)] = obj.get_perctiles(mean(mov_data(:, 1:70)), [0.2, 0.95]);
                if (isempty(lick_rise(ii)) || isempty(mov_rise(ii)))
                    keyborard();
                end
               
                if plot_flg
                    figure; hold on
                    shadedErrorBar(t_vec, lick_data, {@mean, @sem}, {'Color', [0.4, 0, 0.3]})
                    shadedErrorBar(t_vec, mov_data, {@mean, @sem}, {'Color', [0.3, 0.3, 0.3]})
                    xlabel('Time (s)')
                    ylabel('Zscored \DeltaF/F')
                    title([obj.MOUSE_NM, ' ', obj.PROJ_TO(ii)])
                    mean_lick = mean(lick_data);
                    mean_mov = mean(mov_data);
                    line(repmat(t_vec([low_lick, high_lick]), [2, 1]), [0, 0; mean_lick(low_lick), mean_lick(high_lick)],'LineStyle', '--', 'Color', [0.4, 0, 0.3]);
                    line(repmat(t_vec([low_mov, high_mov]), [2, 1]), [-0.05, -0.05; mean_mov(low_mov), mean_mov(high_mov)], 'LineStyle', '--', 'Color', [0.3, 0.3, 0.3]);
                    
                    line(t_vec([low_lick, high_lick]), [0, 0],'LineWidth', 4, 'Color', [0.4, 0, 0.3]);
                    line(t_vec([low_mov, high_mov]), [-0.05, -0.05],'LineWidth', 4, 'Color', [0.3, 0.3, 0.3]);
                end
            end 
        end

    end

    methods (Static = true)
        function [lo, hi, rise_time] =  get_perctiles(data, range)
            data = data - min(data);    % set minimum to zero
            [~, maxi] = max(data);      % get idx of maximum
            if maxi <= 20
                [~, new_maxi] = max(data(21:end));
                maxi = new_maxi + 21;
            end
            values = range*data(maxi);
            lo = find(data(1:maxi)<values(1), 1, 'last');   
            if isempty(lo)
                [~, lo] = min(data(maxi-20:maxi)); % take minimum of 2 sec back from max activity
                lo = lo + maxi-20;
            end
            hi = find(data(21:end)>values(2), 1, 'first');
            hi = hi + 21;
            if isempty(hi) || hi < lo
                rise_time = nan;
            end
            rise_time = (hi-lo) / 10;         % data is downsampled to 10 hz
        end
        function plot_by_outcome (t_vec, trials, info, late_flag, ds_fac)
            % function PLOT_BY_OUTCOME plots event-aligned data according
            % to the outcome of the trial.
            % The function does not label axes or adds title!
            %
            % Inputs:
            %   t_vec - time vector (length = 1 trial)
            %   trials - data to plot (each row: trial)
            %   info - table with trial data (must match size of trials)
            %   late_flag - 0: ignore late trials, 1: merge with omission 2: plot late
            %   ds_fac - downsample factor, 1 or no value for no downsampling.
            
            if nargin < 5
                ds_fac = 1;         % set downsample to 1 if not specified (no downsample)
            end
            
            if late_flag == 2
                leg_tags = {'Omitted', 'Correct', 'Premature', 'Late'};
            elseif late_flag == 1
                leg_tags = {'Miss', 'Hit', 'Premature'};
                info.plot_result(info.plot_result == 2) = 0;                % change late to omission
            else
                leg_tags = {'Omitted', 'Correct', 'Premature'};
                trials(info.plot_result == 2, :) = [];                      % delete late trials
                info(info.plot_result == 2, :) = [];                        % delete late trial info
            end
            
            hold(gca, 'on')
            clr = {'k', [62,98,83]./255, [129,88,135]./255};
            outcomes = [0, 1, -1];
            for ii = 1:3
                if sum(info.plot_result == outcomes(ii)) > 1
                    shadedErrorBar(downsample(t_vec', ds_fac)', downsample(trials(info.plot_result == outcomes(ii), :)', ds_fac)', {@mean, @sem}, {'Color', clr{ii}});
                end
            end
            if late_flag == 2 && sum(info.plot_result == 2) > 1
                shadedErrorBar(downsample(t_vec', ds_fac)', downsample(trials(info.plot_result == 2, :)', ds_fac)', {@mean, @sem}, 'b');
            end
            
            patch_h = flipud(findall(gca, 'Type', 'patch'));          % order is reversed (last one drawn is first)
            legend(gca, patch_h, leg_tags(ismember(outcomes,unique(info.plot_result))))
        end
        
        function hrates = plot_heatmap(data, x_axis, h_axes)
            if nargin < 3 
                hrates = axes(figure);
            else
                hrates = h_axes;
            end
            
            hm = imagesc(hrates, data);
            hm.XData = x_axis;
            hrates.XLim = ([x_axis(1), x_axis(end)]);
            caxis([prctile(data(:), 1), prctile(data(:), 99)])
            colorbar
        end
        
        function add_marks(h, mark)
            mark_x = repmat(mark, [1, 2]);
            mark_y = [(0:(length(mark) - 1))' + 0.5, (1:length(mark))' + 0.5];
            for jj = 1:size(mark_x, 1)
                line(h, mark_x(jj,:), mark_y(jj,:), 'Color','k', 'LineWidth',4)
            end
        
        end
        
        function [h_f, quant] = analyze_data_segment(t_vec, trials, norm_flg, window, ds_fac)
            if nargin < 5
                ds_fac = 1;                         % optional to downsample
            end
            if norm_flg 
                sig = zscore(trials(:, window(1):window(2))')';
                y_l = 'zscored \DeltaF/F';
            else
                sig = trials(:, window(1):window(2));
                y_l = '\DeltaF/F';
            end
            
            sig = sig - mean(sig(:, 1:(round(size(trials, 2) / 4.5) - window(1))), 2); % remove each trial baseline (up to ~5 s) before downsampling
            
            ds_tvec = downsample(t_vec(window(1):window(2))', ds_fac)';       % downsample t_vec
            ds_sig =  downsample(sig', ds_fac)';                                % downsample sig

            h_f = figure; shadedErrorBar(ds_tvec, ds_sig, {@mean, @sem}, 'k');
            xlabel('Time relative to epoch (s)')
            ylabel(y_l)
            
            [pks, locs, w, p] = findpeaks(mean(ds_sig), ds_tvec', 'MinPeakProminence',0.05,  'MinPeakDist', 2, 'Annotate','extents', 'WidthReference','halfprom', 'NPeaks', 1);
            
            ind = 2;            % peak index
            while locs < 0
                % fix problem where initial peak is found before time zero
                [pks, locs, w, p] = findpeaks(mean(ds_sig), ds_tvec', 'MinPeakProminence',0.05, 'Annotate','extents', 'WidthReference','halfprom');
                if numel(locs) == 1
                    pks = []; locs = []; w = []; p = [];                    % no peaks after time zero
                else
                    if ind <= numel(locs)
                        pks = pks(ind); locs = locs(ind); w = w(ind); p = p(ind);   % move through additional peaks until find one after zero
                    else
                        pks = []; locs = []; w = []; p = [];                    % no peaks after time zero 
                    end
                end
                ind = ind + 1;
            end
            if ~isempty(sig)
                quant.peaks = pks;
                quant.location = locs;
                quant.width = w;
                quant.prominence = p;
            else
                quant.peaks = nan;
                quant.location = nan;
                quant.width = nan;
                quant.prominence = nan;
            end
        end
        function h = per_graphs(h_ax, trial_ns, cond, int_fact)
            % function PLOT_PER_GRAPHS() plots psychometric curves based on the
            % behavioral data and the outcome chosen
            %
            % Inputs:
            %   h_ax       - axes handle to plot on
            %   trial_ns   - vector of indices of trials with chosen outcome (correct, premature etc...)
            %   cond       - vector of indices for separating condition (cloud, light, etc)
            %   int_factor - vector of indices for gradient facor trials (e.g. cue int)
            %   trialidx   - all trials and their results
            %
            % Output:
            %   h - plot handle
            for ii = 1:length(int_fact)
                select_trials = intersect(int_fact{ii}, trial_ns);                      % trial outcome (crossed with condition) cross with factor
                per_trials(ii) = size(intersect(select_trials, cond), 1) / ...
                    size(intersect(cond, int_fact{ii}), 1); % ratio between outcome and condition mix
            end
            h = plot(h_ax, per_trials, 'o--');                                          % draw
        end
    end
end