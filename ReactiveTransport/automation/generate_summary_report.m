function generate_summary_report(data_root, folders, folder_params)
    % generate_summary_report - 生成批处理总结报告
    %
    % 输入:
    %   data_root     - 数据根目录
    %   folders       - 文件夹列表
    %   folder_params - 解析的参数列表
    
    % 创建报告文件
    report_file = fullfile(data_root, 'batch_logs', ...
        sprintf('summary_report_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
    
    fid = fopen(report_file, 'w');
    
    fprintf(fid, '═══════════════════════════════════════════════════════════════════════\n');
    fprintf(fid, '              NMR T2 反演批处理总结报告\n');
    fprintf(fid, '═══════════════════════════════════════════════════════════════════════\n\n');
    fprintf(fid, '生成时间: %s\n\n', datestr(now));
    
    fprintf(fid, '处理的文件夹:\n');
    fprintf(fid, '───────────────────────────────────────────────────────────────────────\n');
    
    for i = 1:length(folders)
        params = folder_params{i};
        folder_path = fullfile(data_root, folders(i).name);
        
        fprintf(fid, '\n[%d] %s\n', i, folders(i).name);
        fprintf(fid, '    参数: Da=%.4f, Pe=%.4f, L=%.4f\n', params.Da, params.Pe, params.L);
        fprintf(fid, '    尺寸: X=%.4f cm, Y=%.4f cm\n', params.lengthXAxis, params.lengthYAxis);
        fprintf(fid, '    类型: %s\n', params.layoutType);
        
        % 统计DXF文件
        dxf_pore_dir = fullfile(folder_path, 'dxf_pore');
        if exist(dxf_pore_dir, 'dir')
            pore_files = dir(fullfile(dxf_pore_dir, '*.dxf'));
            fprintf(fid, '    DXF文件数: %d\n', length(pore_files));
        end
        
        % 统计COMSOL结果
        comsol_dir = fullfile(folder_path, 'comsol_results');
        if exist(comsol_dir, 'dir')
            xlsx_files = dir(fullfile(comsol_dir, '*.xlsx'));
            fprintf(fid, '    COMSOL结果: %d 个\n', length(xlsx_files));
        else
            fprintf(fid, '    COMSOL结果: 无\n');
        end
        
        % 统计反演结果
        inversion_dir = fullfile(folder_path, 'inversion_results');
        if exist(inversion_dir, 'dir')
            png_files = dir(fullfile(inversion_dir, '*.png'));
            fprintf(fid, '    反演图像: %d 个\n', length(png_files));
        else
            fprintf(fid, '    反演图像: 无\n');
        end
    end
    
    fprintf(fid, '\n───────────────────────────────────────────────────────────────────────\n');
    fprintf(fid, '报告结束\n');
    
    fclose(fid);
    
    fprintf('总结报告已保存: %s\n', report_file);
end
