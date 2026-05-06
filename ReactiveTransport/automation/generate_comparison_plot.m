function generate_comparison_plot(inversion_results, global_time, global_porosity, sample_indices, output_dir, folder_name, show_nmr)
    % generate_comparison_plot - Generate comparison plot between NMR inversion water content and original porosity
    %
    % Inputs:
    %   inversion_results - Struct containing inversion results
    %       .timesteps        - Array of timestep indices
    %       .water_contents   - Array of inverted water contents
    %       .timestep_strings - Cell array of timestep strings
    %   global_time       - Time data from global_evolution.xlsx
    %   global_porosity   - Porosity data from global_evolution.xlsx
    %   sample_indices    - Indices of sampled timesteps
    %   output_dir        - Output directory for the plot
    %   folder_name       - Name of the dissolution folder
    %   show_nmr          - (Optional) Whether to show NMR inversion results (default: true)
    
    % Handle optional parameter
    if nargin < 7 || isempty(show_nmr)
        show_nmr = true;  % Default: show NMR results
    end
    
    % === 可视化参数设置 ===
    fontSize = 20;           % 坐标轴标签和刻度字体大小
    titleFontSize = 22;      % 标题字体大小
    legendFontSize = 14;     % 图例字体大小
    fontName = 'Helvetica';  % 字体类型
    
    % Create figure - 只保留一个图
    fig = figure('Visible', 'off', 'Position', [100 100 1200 600]);
    
    % ========== Single Plot: Porosity Evolution ==========
    hold on;
    
    % Plot original porosity from global_evolution.xlsx
    if isnumeric(global_time)
        h1 = plot(global_time, global_porosity, 'b-', 'LineWidth', 2);
    else
        h1 = plot(1:length(global_porosity), global_porosity, 'b-', 'LineWidth', 2);
    end
    
    legendHandles = h1;
    legendEntries = {'Simulation'};
    
    % Plot inverted water content at sampled timesteps (only if show_nmr is true)
    if show_nmr && ~isempty(inversion_results.timesteps)
        % Match timesteps to global_time
        inv_times = zeros(size(inversion_results.timesteps));
        for k = 1:length(inversion_results.timesteps)
            ts_idx = inversion_results.timesteps(k);
            if ts_idx <= length(global_time)
                if isnumeric(global_time)
                    inv_times(k) = global_time(ts_idx);
                else
                    inv_times(k) = ts_idx;
                end
            else
                inv_times(k) = ts_idx;
            end
        end
        
        h2 = plot(inv_times, inversion_results.water_contents, 'ro-', ...
            'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        legendHandles = [legendHandles, h2];
        legendEntries = [legendEntries, {'NMR Inversion'}];
    end
    
    hold off;
    xlabel('Time (s)', 'FontSize', fontSize, 'FontName', fontName);
    ylabel('Porosity / Water Content', 'FontSize', fontSize, 'FontName', fontName);
    
    % Adjust title based on whether NMR is shown
    if show_nmr
        title('Porosity Evolution: Simulation vs NMR Inversion', 'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
    else
        title('Porosity Evolution: Simulation', 'FontSize', titleFontSize, 'FontName', fontName, 'FontWeight', 'normal');
    end
    
    legend(legendHandles, legendEntries, 'Location', 'best', 'FontSize', legendFontSize, 'FontName', fontName);
    set(gca, 'FontSize', fontSize-2, 'FontName', fontName);
    grid on;
    box on;
    
    % ========== Save figure ==========
    output_png = fullfile(output_dir, 'porosity_comparison.png');
    print(fig, output_png, '-dpng', '-r300');
    close(fig);
    
    fprintf('        [OK] Comparison plot saved: %s\n', output_png);
    
    % ========== Save comparison data to Excel ==========
    if show_nmr && ~isempty(inversion_results.timesteps)
        try
            valid_indices = inversion_results.timesteps(inversion_results.timesteps <= length(global_porosity));
            n_valid = length(valid_indices);
            
            if n_valid > 0
                original_vals = global_porosity(valid_indices(1:n_valid));
                inverted_vals = inversion_results.water_contents(1:n_valid);
                rel_error = (inverted_vals - original_vals) ./ original_vals * 100;
                
                % Create table
                comparison_table = table(...
                    valid_indices(:), ...
                    original_vals(:), ...
                    inverted_vals(:), ...
                    rel_error(:), ...
                    'VariableNames', {'TimeStep', 'Original_Porosity', 'NMR_WaterContent', 'RelativeError_Percent'});
                
                % Save to Excel
                output_xlsx = fullfile(output_dir, 'porosity_comparison.xlsx');
                writetable(comparison_table, output_xlsx);
                fprintf('        [OK] Comparison data saved: %s\n', output_xlsx);
            end
        catch ME
            fprintf('        [WARN] Failed to save comparison data: %s\n', ME.message);
        end
    end
end