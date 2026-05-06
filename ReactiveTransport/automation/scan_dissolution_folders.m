function folders = scan_dissolution_folders(data_root)
    % scan_dissolution_folders - 扫描所有RTM结果文件夹
    %
    % 输入:
    %   data_root - 数据根目录
    %
    % 输出:
    %   folders - 包含所有dissolution文件夹信息的结构体数组
    
    % 获取所有子目录
    all_items = dir(data_root);
    
    % 筛选出RTM结果文件夹
    folders = [];
    
    for i = 1:length(all_items)
        item = all_items(i);
        
        % 跳过非目录
        if ~item.isdir
            continue;
        end
        
        % 跳过 . 和 ..
        if strcmp(item.name, '.') || strcmp(item.name, '..')
            continue;
        end
        
        % 跳过batch_logs目录
        if strcmp(item.name, 'batch_logs')
            continue;
        end
        
        folder_path = fullfile(data_root, item.name);
        is_legacy_name = startsWith(item.name, 'dissolution_results-');
        has_metadata = exist(fullfile(folder_path, 'run_metadata.json'), 'file') == 2;
        has_dxf = exist(fullfile(folder_path, 'dxf_pore'), 'dir') && ...
                  exist(fullfile(folder_path, 'dxf_solid'), 'dir');

        if has_dxf && (is_legacy_name || has_metadata)
            folders = [folders; item];
        end
    end
    
    % 按名称排序
    if ~isempty(folders)
        [~, idx] = sort({folders.name});
        folders = folders(idx);
    end
end
