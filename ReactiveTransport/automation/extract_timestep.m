function timestep = extract_timestep(filename)
    % extract_timestep - 从文件名提取时间步编号
    %
    % 输入:
    %   filename - 文件名 (如 pore_t0001.dxf)
    %
    % 输出:
    %   timestep - 时间步字符串 (如 '0001')
    
    match = regexp(filename, 't(\d+)', 'tokens');
    if ~isempty(match)
        timestep = match{1}{1};
    else
        timestep = '0000';
    end
end
