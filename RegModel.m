% Model Analysis
classdef RegModel < handle
    
    properties
        data = struct();
        SVM
        UC
        GSVM
        GUC
        GROUPS
        PATH
    end
    methods
    function obj = RegModel(path, scope, merge_aud, source)
        % scope = 'full' : Takes CVR2 data from entire signal
        %         'label': Takes CVR2 data limited to the scope of that label
                 
            obj.PATH = path;

            switch source
                case 'cno'
                    load(strcat(obj.PATH, '\model_params_session_avg_080421_full_model_noruns_DREADDs_twoT_shuffled'));
                case 'main'
                    load(strcat(obj.PATH, '\model_params_session_avg_130622_full_model_window_noruns_allmice_twoT_fixed_SHUFFLED_PreBL_FINAL'));
            end
            
           
            switch scope
                case 'full'
                    obj.SVM = SVM_full;
                    obj.UC = UC_full;
                    obj.GSVM = GSVM;
                    obj.GUC = GUC;
% %                     % Normalize by window size
                    switch width(SVM_full)-1
                        case 17
                            obj.data.win_length = [5.5, 5.5, 5.5, 8, 11, 3, 3, 3, 3, 3, 1.5, 10, 10, 10, 5, 5, 5];
                        case 16
                            obj.data.win_length = [5.5, 5.5, 5.5, 5.5, 3, 3, 3, 3, 3, 1.5, 10, 10, 10, 5, 5, 5];
                        case 15
                            obj.data.win_length = [5.5, 5.5, 5.5, 5.5, 3, 3, 3, 3, 3, 1.5, 10, 10, 5, 5, 5];
                    end
                    obj.SVM{:, 1:width(SVM_full)-1} = obj.SVM{:, 1:width(SVM_full)-1} ./ obj.data.win_length;
                    obj.UC{:, 1:width(UC_full)-1} = obj.UC{:, 1:width(UC_full)-1} ./ obj.data.win_length;
                    % End normalize by window size
                    
                case 'label'
                    obj.SVM = SVM_eventwindow;
                    obj.UC = UC_eventwindow;
                    obj.GSVM = GSVM_eventwindow;
                     obj.GUC = GUC_eventwindow;
                    switch width(SVM_full)-1
                        case 17
                            obj.data.win_length = [5.5, 5.5, 5.5, 8, 11, 3, 3, 3, 3, 3, 1.5, 10, 10, 10, 5, 5, 5];
                        case 16
                            obj.data.win_length = [5.5, 5.5, 5.5, 5.5, 3, 3, 3, 3, 3, 1.5, 10, 10, 10, 5, 5, 5];
                        case 15
                            obj.data.win_length = [5.5, 5.5, 5.5, 5.5, 3, 3, 3, 3, 3, 1.5, 10, 10, 5, 5, 5];
                    end
                    obj.SVM{:, 1:width(SVM_full)-1} = obj.SVM{:, 1:width(SVM_full)-1} ./ obj.data.win_length;
                    obj.UC{:, 1:width(UC_full)-1} = obj.UC{:, 1:width(UC_full)-1} ./ obj.data.win_length;
                    % End normalize by window size
            end
            %normalize by total SVM
            obj.SVM{:, 1:end-1} =  obj.SVM{:, 1:end-1} ./  obj.SVM{:, end};
            obj.UC{:, 1:end-1} =  obj.UC{:, 1:end-1} ./  obj.UC{:, end};
            obj.GSVM{:, 1:end-1} =  obj.GSVM{:, 1:end-1} ./  obj.GSVM{:, end};
            obj.GUC{:, 1:end-1} =  obj.GUC{:, 1:end-1} ./  obj.GUC{:, end};
            
            
            
            
            obj.GROUPS = grouping;
            nmes = fieldnames(grouping);
            for ii = 1:numel(nmes)
                obj.data.group_idx.(nmes{ii}) = obj.tag_by_group(grouping.(nmes{ii}));
            end
            obj.data.cmap = [140, 179, 250; 62, 98, 83; 129, 88, 135] ./255;
            
            if merge_aud
                obj.merge_aud;
            end
    end
    
    function merge_aud(obj)
        data_nms = {'SVM', 'UC'};
        for ii = 1:numel(data_nms)
            tmp = obj.(data_nms{ii});
            colidx = find(contains(tmp.Properties.VariableNames, 'Aud'));
            visidx = find(contains(tmp.Properties.VariableNames, 'Vis'));
            tmp_tbl = table(mean(tmp{:, colidx(1:2)},2), ...
                            mean(tmp{:, colidx(3:4)},2), ...
                            tmp{:, visidx},...
                            'VariableNames', ...
                            {'AudLow', 'AudHigh', 'AudVis'});
            tmp(:, [colidx, visidx]) = [];
            tmp = [tmp(:, 1:(min([colidx, visidx])-1)), tmp_tbl, tmp(:, min([colidx, visidx]):end)];
            obj.(data_nms{ii}) = tmp; 
        end
    end
    function [ha, svm_plot, uc_plot] = plot_by_group(obj, nme, divide_by, exclude, log_scale)
        % returns handles for the plots by paramter and by group
        % divide by: "cloud", "orig"
        % if exclude = 1: show only ACC and OFC; if exclude = 2: show only ACC
        % if log scale - log scale (not really beneficial)
        if strcmp(divide_by, 'area')
            groups = [{obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'ACC'))}, ...
                {obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'OFC'))}];
            grp_nms = {'ACC', 'OFC'};
        elseif strcmp(divide_by, 'cno')
            groups = [{obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'C_'))}, ...
                {obj.SVM.Properties.RowNames(~contains(obj.SVM.Properties.RowNames, 'C_'))}];
            grp_nms = {'Saline', 'CNO'};
        elseif strcmp(divide_by, 'input')
            groups = [{obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'ACCin'))}, ...
                {obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'AUDin'))}];
            grp_nms = {'ACCin', 'AUDin'};
        else
            groups = obj.GROUPS.([divide_by, '_sep']);
            grp_nms = {'exploiters', 'explorers', ' '};
        end
              
        f = figure;
        f.Position = [253,  0,   987,   684];
        
        for ii = 1:numel(groups)
            subgroup = groups{ii};
            [sub_svm, sub_uc, sub_gsvm, sub_guc] = obj.get_data(nme, subgroup);
            switch exclude
                case 1
                    exclude_idx = find(contains(sub_svm.Properties.RowNames, {'in', 'AUD'}));
                    sub_svm(exclude_idx, :) = [];
                    sub_uc(exclude_idx, :) = [];
                    sub_gsvm(exclude_idx, :) = [];
                    sub_guc(exclude_idx, :) = [];
                case 2
                    exclude_idx = find(contains(sub_svm.Properties.RowNames, {'in', 'AUD', 'OFC'}));
                    sub_svm(exclude_idx, :) = [];
                    sub_uc(exclude_idx, :) = [];
                    sub_gsvm(exclude_idx, :) = [];
                    sub_guc(exclude_idx, :) = [];
                case 'input'
                    exclude_idx = find(~contains(sub_svm.Properties.RowNames, {'in'}));
                    sub_svm(exclude_idx, :) = [];
                    sub_uc(exclude_idx, :) = [];
                    sub_gsvm(exclude_idx, :) = [];
                    sub_guc(exclude_idx, :) = [];
            end
            svm_plot{ii} = table2array(sub_svm);
            uc_plot{ii} = table2array(sub_uc);
            xax = repmat([1:size(svm_plot{ii}, 2)], [size(sub_svm, 1), 1]);
            
            ha(1, ii) = subplot(2, numel(groups), ii);
            hold on
            title(grp_nms{ii})
            ylabel('CVR^2')
            obj.bootstrap_plot(svm_plot{ii}, '', 5000, ha(1, ii));
            for tt = 1:size(sub_svm, 1)
                tmp_h = plot(xax(tt, :), svm_plot{ii}(tt, :), '--ob');
                tmp_h.Tag = sub_svm.Properties.RowNames{tt};
            end
            %boxplot(svm_plot{ii})
            
            ha(2, ii) = subplot(2, numel(groups), numel(groups) + ii);
            hold on
            ylabel('\DeltaR^2')
            obj.bootstrap_plot(uc_plot{ii}, '', 5000, ha(2, ii));
            for tt = 1:size(sub_uc, 1)
                tmp_h = plot(xax(tt, :), uc_plot{ii}(tt, :), '--ob');
                tmp_h.Tag = sub_uc.Properties.RowNames{tt};
            end
            
        end
        for ii = 1:2
        ylims = [ha(ii, :).YLim];
        set(ha(ii,:), 'YLim', [min(ylims(:)), max(ylims(:))])
        set(ha(ii,:), 'XTick', [1:size(sub_svm, 2)], ...
            'XTickLabel', sub_svm.Properties.VariableNames(:), ...
            'TickLabelInterpreter', 'none')
        if log_scale
            set(ha(ii, :), 'YScale', 'log')
        end
        end
        dcm_obj = datacursormode(f);
        set(dcm_obj,'UpdateFcn',@myupdatefcn)
        
        function txt = myupdatefcn(~,event_obj)
            % Customizes text of data tips
            try
                pos = get(event_obj,'Position');
                plot_obj = get(event_obj, 'Target');
                obj_name = strrep(plot_obj.Tag(1:end), '_', ' ');
                txt = {num2str(pos(2)), obj_name};              % add to the highlighted tag the name of the mouse
            catch
                keyboard();                                 % in case something goes wrong
            end
        end
    end

    function h_bar = model_group_summary(obj, data_type, divide_by, avg_flg, exclude)
        if strcmp(divide_by, 'area')
            groups = [{obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'ACC'))}, ...
                {obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'OFC'))}];
            grp_nms = {'ACC', 'OFC'};
        elseif strcmp(divide_by, 'input')
            groups = [{obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'ACCin'))}, ...
                {obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'AUDin'))}];
            grp_nms = {'ACCin', 'AUDin'};
        elseif strcmp(divide_by, 'cno')
            groups = [{obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'C_'))}, ...
                {obj.SVM.Properties.RowNames(~contains(obj.SVM.Properties.RowNames, 'C_'))}];
            grp_nms = {'Saline', 'CNO'};
        else
            groups = obj.GROUPS.([divide_by, '_sep']);
            grp_nms = {'exploiters', 'explorers', ' '};
        end
        h_bar = figure;
        for ii = 1:numel(groups)
            subgroup = groups{ii};
            [sub_svm, sub_uc, sub_gsvm, sub_guc] = obj.get_data([], subgroup);
            switch exclude
                case 1
                    exclude_idx = find(contains(sub_svm.Properties.RowNames, {'in', 'AUD'}));
                    sub_svm(exclude_idx, :) = [];
                    sub_uc(exclude_idx, :) = [];
                    sub_gsvm(exclude_idx, :) = [];
                    sub_guc(exclude_idx, :) = [];
                case 2
                    exclude_idx = find(contains(sub_svm.Properties.RowNames, {'in', 'AUD', 'OFC'}));
                    sub_svm(exclude_idx, :) = [];
                    sub_uc(exclude_idx, :) = [];
                    sub_gsvm(exclude_idx, :) = [];
                    sub_guc(exclude_idx, :) = [];
            end
            switch data_type
                case 'SVM'
                    values = sub_svm(:, 1:end-1);
                case 'UC'
                    values = sub_uc(:, 1:end-1);
                case 'GSVM'
                    values = sub_gsvm(:, 1:end-1);
                case 'GUC'
                    values = sub_guc(:, 1:end-1);
            end
            subplot(numel(groups), 1, ii)
            hold on
            if avg_flg
                bar(mean(table2array(values)', 2))
                errorbar(mean(table2array(values)', 2), sem(table2array(values)), '.')
            else
                 bar(table2array(values)')
            end
            title(grp_nms{ii})
            xticks(1:numel(values.Properties.VariableNames));
            xticklabels(values.Properties.VariableNames);
            ylabel(data_type);
            
        end
        set(h_bar.Children, 'YLim', [min([h_bar.Children.YLim]), max([h_bar.Children.YLim])])
    end
    
    function comp_tbl = get_stat_tbl(obj, param, nme, divide_by, exclude)
        % get table for statistics on parameters
        % param - SVM or UC
        % nme - columns to take
        % divide by 'orig' 'cloud'
        % exclude AUD and inputs
        comp_tbl = obj.(param)(:, find(contains(obj.(param).Properties.VariableNames, nme)));
        comp_tbl.group = obj.data.group_idx.([divide_by, '_sep']);
        switch exclude
            case 1
                exclude_idx = find(contains(comp_tbl.Properties.RowNames, {'in', 'AUD'}));
                comp_tbl(exclude_idx, :) = [];
            case 2
                exclude_idx = find(contains(comp_tbl.Properties.RowNames, {'in', 'AUD', 'OFC'}));
                comp_tbl(exclude_idx, :) = [];
            case 'input'
                    exclude_idx = find(~contains(sub_svm.Properties.RowNames, {'in'}));
                    comp_tbl(exclude_idx, :) = [];
        end
    end
        
    function [sub_svm, sub_uc, sub_gsvm, sub_guc] = get_data(obj, nme, group)
        rowidx = find(contains(obj.SVM.Properties.RowNames, group));
        if ~isempty(nme)
            colidx = find(contains(obj.SVM.Properties.VariableNames, nme));
        else 
            colidx = 1:width(obj.SVM);
        end
        
        sub_svm = obj.SVM(rowidx, colidx);
        sub_uc = obj.UC(rowidx, colidx);
        sub_gsvm = obj.GSVM(rowidx, :);
        sub_guc = obj.GUC(rowidx, :);
    end
     
    function tags = tag_by_group(obj, group)
        tags = zeros(height(obj.SVM), 1);
        for ii = 1:numel(group)
        tags(contains(obj.SVM.Properties.RowNames, group{ii})) = ii;
        end
    end

    function h_bar = vert_summary(obj, data_type)
        groups = [{obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'ACC'))}, ...
            {obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'OFC'))}];
        grp_nms = {'ACC', 'OFC'};
        clrs = [139, 197, 184; 189, 191, 224] ./ 255;
        spcing = 0.2;
        h_bar = figure;
        avgs = [];
        for ii = 1:numel(groups)
            subgroup = groups{ii};
            [sub_svm, sub_uc, sub_gsvm, sub_guc] = obj.get_data([], subgroup);
            exclude_idx = find(contains(sub_svm.Properties.RowNames, {'in', 'AUD'}));
            sub_svm(exclude_idx, :) = [];
            sub_uc(exclude_idx, :) = [];
            sub_gsvm(exclude_idx, :) = [];
            sub_guc(exclude_idx, :) = [];
            
            hold on            
            switch data_type
                case 'SVM'
                    values = sub_svm(:, 1:end-1);
                case 'UC'
                    values = sub_uc(:, 1:end-1);
                case 'GSVM'
                    values = sub_gsvm(:, 1:end-1);
                case 'GUC'
                    values = sub_guc(:, 1:end-1);
            end
            
            nums = fliplr(table2array(values)); % flip order
            plot(nums, repmat(1:size(nums, 2)', [size(nums, 1), 1])+(-1)^ii*spcing, '*', 'Color', [clrs(ii, :), 0.2])
            avgs(:, ii) = median(nums)';
            hb(ii) = barh((1:size(nums, 2))'+(-1)^ii*spcing, avgs(:, ii), 'FaceColor', clrs(ii, :), 'BarWidth', spcing, 'FaceAlpha', 0.5);
            errs(:, ii) = sem(nums);
            errorbar(hb(ii).YData', (1:size(nums, 2))'+(-1)^ii*spcing, prctile(nums, 25), prctile(nums, 75),'horizontal', 'k.')
        end
        labels = fliplr(values.Properties.VariableNames);         
        yticks(1:numel(values.Properties.VariableNames));
        yticklabels(labels);
        xlabel(data_type);
        title('Model Summary')
    end    
    
    function [ha] = plot_label_by_group(obj, nme, divide_by, exclude, log_scale)
        % returns handles for the plots by paramter and by group
        % divide by: "cloud", "orig"
        % if exclude = 1: show only ACC and OFC; if exclude = 2: show only ACC
        % if log scale - log scale (not really beneficial)
        if strcmp(divide_by, 'area')
            groups = [{obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'ACC'))}, ...
                {obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'OFC'))}];
            grp_nms = {'ACC', 'OFC'};
            cmap = [139, 197, 184; 189, 191, 224];
        elseif strcmp(divide_by, 'cno')
            groups = [{obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'C_'))}, ...
                {obj.SVM.Properties.RowNames(~contains(obj.SVM.Properties.RowNames, 'C_'))}];
            grp_nms = {'Saline', 'CNO'};
            cmap = [12, 114, 186; 148, 29, 94];
        elseif strcmp(divide_by, 'input')
            groups = [{obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'ACCin'))}, ...
                {obj.SVM.Properties.RowNames(contains(obj.SVM.Properties.RowNames, 'AUDin'))}];
            grp_nms = {'ACCin', 'AUDin'};
            cmap = [139, 197, 184; 189, 191, 224];
        else
            groups = obj.GROUPS.([divide_by, '_sep']);
            grp_nms = {'exploiters', 'explorers', ' '};
            cmap = obj.data.cmap .* 255;
        end
        
        f = figure;
        f.Position = [73   256   900   684];
        
        svm_plot = table;
        uc_plot = table;
        
        for ii = 1:numel(groups)
            subgroup = groups{ii};
            [sub_svm, sub_uc, sub_gsvm, sub_guc] = obj.get_data(nme, subgroup);
            switch exclude
                case 1
                    exclude_idx = find(contains(sub_svm.Properties.RowNames, {'in', 'AUD'}));
                    sub_svm(exclude_idx, :) = [];
                    sub_uc(exclude_idx, :) = [];
                    sub_gsvm(exclude_idx, :) = [];
                    sub_guc(exclude_idx, :) = [];
                case 2
                    exclude_idx = find(contains(sub_svm.Properties.RowNames, {'in', 'AUD', 'OFC'}));
                    sub_svm(exclude_idx, :) = [];
                    sub_uc(exclude_idx, :) = [];
                    sub_gsvm(exclude_idx, :) = [];
                    sub_guc(exclude_idx, :) = [];
                case 'input'
                    exclude_idx = find(~contains(sub_svm.Properties.RowNames, {'in'}));
                    sub_svm(exclude_idx, :) = [];
                    sub_uc(exclude_idx, :) = [];
                    sub_gsvm(exclude_idx, :) = [];
                    sub_guc(exclude_idx, :) = [];
            end
            
            sub_svm.group(:) = ii;
            svm_plot = [svm_plot; sub_svm];
            sub_uc.group(:) = ii;
            uc_plot = [uc_plot; sub_uc];
        end
        
        num_labels = width(sub_svm) - 1;
        num_pts = height(svm_plot);
        plot_tbl = table(svm_plot{:,1:num_labels}(:), ...
            uc_plot{:,1:num_labels}(:), ...
            repmat(svm_plot.group, [num_labels, 1]), ...
            reshape((1:num_labels) .* ones(num_pts, 1), [], 1), ...
            repmat((1:num_pts)', [num_labels, 1]),...
            'VariableNames', {'svm_data', 'uc_data', 'group', 'label', 'mouse'});        
 
        sgtitle(f, [nme, ' Comparison'])
        ha(1) = subplot(2, 1, 1);
        hold(ha(1), 'on')
        ylabel(ha(1), 'SVM')
        
        ha(2) = subplot(2, 1, 2);
        hold(ha(2), 'on')
        ylabel(ha(2), 'UC')
        
        drift_size = 0.5*num_labels/4;
        drift = linspace(-drift_size, drift_size, num_labels);
        spacing = 1 + num_labels/4;
        
        if num_labels > 1
            for tt=1:num_labels
                htmp = axes(figure);
                clrs = cmap;
                beeswarm(spacing * plot_tbl.group(plot_tbl.label == tt) + drift(tt), plot_tbl.svm_data(plot_tbl.label == tt), 'sort_style', 'up', 'corral_style', 'gutter', 'colormap', clrs./255, 'MarkerFaceAlpha', 0.7, 'overlay_style', 'sd');
                dots = htmp.Children;
                copyobj(dots, ha(1))
                close(htmp.Parent);
            end
            
            
            for tt=1:num_labels
                htmp = axes(figure);
                clrs = cmap;
                beeswarm(spacing * plot_tbl.group(plot_tbl.label == tt) + drift(tt), plot_tbl.uc_data(plot_tbl.label == tt), 'sort_style', 'down', 'corral_style', 'gutter', 'colormap', clrs./255, 'MarkerFaceAlpha', 0.7, 'overlay_style', 'sd');
                dots = htmp.Children;
                copyobj(dots, ha(2))
                close(htmp.Parent);
            end
            xs = sort([drift' * [ones(1,numel(grp_nms))] + spacing * [1:numel(grp_nms)]]);
            for ii = 1:2
                ylims = [ha(ii).YLim];
                set(ha(ii), 'YLim', [min(ylims(:)), max(ylims(:))], 'XLim', [xs(1)-xs(1)/2, xs(end)+xs(1)/2])
                set(ha(ii), 'XTick', xs(:), ...
                    'XTickLabel', sub_svm.Properties.VariableNames(1:num_labels), ...
                    'XTickLabelRotation', 45,   ...
                    'TickLabelInterpreter', 'none')
                if log_scale
                    set(ha(ii), 'YScale', 'log')
                end
            end
        else
            htmp = axes(figure);
            clrs = cmap;
            beeswarm(spacing * plot_tbl.group(plot_tbl.label == 1), plot_tbl.svm_data(plot_tbl.label == 1), 'sort_style', 'rand', 'colormap', clrs./255, 'MarkerFaceAlpha', 0.7, 'overlay_style', 'sd');
            dots = htmp.Children;
            copyobj(dots, ha(1))
            close(htmp.Parent);
            
            htmp = axes(figure);
            clrs = cmap;
            beeswarm(spacing * plot_tbl.group(plot_tbl.label == 1), plot_tbl.uc_data(plot_tbl.label == 1), 'sort_style', 'rand', 'colormap', clrs./255, 'MarkerFaceAlpha', 0.7, 'overlay_style', 'sd');
            dots = htmp.Children;
            copyobj(dots, ha(2))
            close(htmp.Parent);
            
            xs = sort([drift' * [ones(1,numel(grp_nms))] + spacing * [1:numel(grp_nms)]]);
            for ii = 1:2
                ylims = [ha(ii).YLim];
                set(ha(ii), 'YLim', [min(ylims(:)), max(ylims(:))])
                set(ha(ii), 'XTick', xs(:), ...
                    'XTickLabel', grp_nms, ...
                    'TickLabelInterpreter', 'none')
                if log_scale
                    set(ha(ii), 'YScale', 'log')
                end
            end
        end
        
    end
    
    % ----------------------------------------------------------------------%
    end
    methods (Static = true)
        function ax = bootstrap_plot(data, nme, bstrp_num, hax)
            if istable(data)
                prem_data = data(data.group == nme, 1);
            elseif size(data, 2) > 1
                prem_data = data(strcmp(data(:, 2), nme), 1);
            else
                prem_data = data;
            end
            
            if nargin < 4
                ax = axes(figure);
            else
                ax = hax;
            end
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