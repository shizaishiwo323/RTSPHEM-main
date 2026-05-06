function [pore_files, solid_files] = get_dxf_files(folder_path)
    % get_dxf_files - 获取配对的孔隙和固体DXF文件
    %
    % 输入:
    %   folder_path - dissolution结果文件夹路径
    %
    % 输出:
    %   pore_files  - 孔隙DXF文件列表
    %   solid_files - 固体DXF文件列表 (与pore_files一一对应)
    
    pore_dir = fullfile(folder_path, 'dxf_pore');
    solid_dir = fullfile(folder_path, 'dxf_solid');
    
    pore_files = [];
    solid_files = [];
    
    % 检查目录存在性
    if ~exist(pore_dir, 'dir') || ~exist(solid_dir, 'dir')
        return;
    end
    
    % 获取所有孔隙DXF文件
    pore_all = dir(fullfile(pore_dir, 'pore_t*.dxf'));
    
    if isempty(pore_all)
        return;
    end
    
    % 按名称排序
    [~, idx] = sort({pore_all.name});
    pore_all = pore_all(idx);
    
    % 查找配对的固体文件
    for i = 1:length(pore_all)
        pore_name = pore_all(i).name;
        
        % 从pore_t0001.dxf提取时间步编号
        timestep = regexp(pore_name, 'pore_t(\d+)\.dxf', 'tokens');
        if isempty(timestep)
            continue;
        end
        
        % 构造对应的solid文件名
        solid_name = sprintf('solid_t%s.dxf', timestep{1}{1});
        solid_path = fullfile(solid_dir, solid_name);
        
        % 检查solid文件是否存在
        if exist(solid_path, 'file')
            pore_files = [pore_files; pore_all(i)];
            solid_files = [solid_files; dir(solid_path)];
        end
    end
end
