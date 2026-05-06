function success = run_comsol_processing(mph_file, pore_dxf, solid_dxf, lengthXAxis, lengthYAxis, excel_output, config)
    % run_comsol_processing - 运行COMSOL几何处理和求解
    %
    % 输入:
    %   mph_file     - COMSOL模型文件路径
    %   pore_dxf     - 孔隙DXF文件路径
    %   solid_dxf    - 固体DXF文件路径
    %   lengthXAxis  - X轴长度 (用于几何缩放验证)
    %   lengthYAxis  - Y轴长度 (用于边界计算)
    %   excel_output - Excel输出路径
    %   config       - 配置对象
    %
    % 输出:
    %   success - 处理是否成功
    
    % 记录几何参数供日志使用
    fprintf('      几何尺寸: X=%.6f cm, Y=%.6f cm\n', lengthXAxis, lengthYAxis);
    
    % 检查 lengthYAxis 是否有效
    if isnan(lengthYAxis) || lengthYAxis <= 0
        fprintf('      ⚠ 警告: lengthYAxis 无效 (%.6f)，将使用默认估算\n', lengthYAxis);
        % 尝试从 lengthXAxis 估算 (假设 Y ≈ 2/3 X)
        if ~isnan(lengthXAxis) && lengthXAxis > 0
            lengthYAxis = lengthXAxis * 0.667;
            fprintf('      → 从 lengthXAxis 估算 lengthYAxis = %.6f\n', lengthYAxis);
        end
    end
    
    success = false;
    
    % 检查是否启用COMSOL
    if ~config.enable_comsol
        fprintf('      ⚠ COMSOL处理已禁用\n');
        success = true;
        return;
    end
    
    % 检查输入文件
    if ~exist(mph_file, 'file')
        fprintf('      × 模型文件不存在: %s\n', mph_file);
        return;
    end
    
    if ~exist(pore_dxf, 'file')
        fprintf('      × 孔隙DXF不存在: %s\n', pore_dxf);
        return;
    end
    
    if ~exist(solid_dxf, 'file')
        fprintf('      × 固体DXF不存在: %s\n', solid_dxf);
        return;
    end
    
    % 检查是否已存在结果且不需要覆盖
    % 同时检查短文件名和长文件名格式
    [output_dir, output_name, output_ext] = fileparts(excel_output);
    
    % 从文件名中提取时间步
    timestep_match = regexp(output_name, 't(\d{4})', 'tokens');
    if ~isempty(timestep_match)
        timestep_str = timestep_match{1}{1};
        % 检查短文件名格式
        short_filename = sprintf('T2_t%s%s', timestep_str, output_ext);
        short_filepath = fullfile(output_dir, short_filename);
        
        if exist(excel_output, 'file') && ~config.overwrite_existing
            fprintf('      ⚠ 结果已存在,跳过COMSOL处理\n');
            success = true;
            return;
        elseif exist(short_filepath, 'file') && ~config.overwrite_existing
            fprintf('      ⚠ 结果已存在(短文件名),跳过COMSOL处理\n');
            success = true;
            return;
        end
    else
        % 没有找到时间步，只检查原始文件名
        if exist(excel_output, 'file') && ~config.overwrite_existing
            fprintf('      ⚠ 结果已存在,跳过COMSOL处理\n');
            success = true;
            return;
        end
    end
    
    try
        %% 连接COMSOL
        fprintf('      连接COMSOL...\n');
        
        % 尝试连接到本地服务器
        fprintf('=== 步骤1: 初始化COMSOL ===\n');
        model = comsol_init();
        %% 加载模型
        fprintf('      加载模型...\n');
        model = mphload(mph_file);
        
        %% 更新几何路径
        fprintf('      更新几何路径...\n');
        comp_tag = 'comp1';
        geom_tag = 'geom1';
        
        model.component(comp_tag).geom(geom_tag).feature('imp1').set('filename', pore_dxf);
        model.component(comp_tag).geom(geom_tag).feature('imp2').set('filename', solid_dxf);
        
        %% 分析DXF获取周长
        fprintf('      分析DXF周长...\n');
        pore_total_perimeter = analyze_dxf_perimeter(pore_dxf);
        
        %% 构建几何
        fprintf('      构建几何...\n');
        
        % 获取几何对象
        geom = model.component(comp_tag).geom(geom_tag);
        
        % 设置缩放系数
        scale_factor = config.scale_factor;
        fprintf('        缩放系数: %d\n', scale_factor);
        
        % 更新缩放特征的缩放系数
        scale_updated = false;
        
        % 方法1: 直接尝试访问 sca1
        try
            geom.feature('sca1').set('factor', num2str(scale_factor));
            fprintf('        ✓ 缩放系数已更新为 %d (通过sca1)\n', scale_factor);
            scale_updated = true;
        catch
            % 方法2: 遍历查找缩放特征
            try
                feature_tags = geom.feature.tags;
                for i = 1:length(feature_tags)
                    ftag = char(feature_tags(i));
                    if contains(lower(ftag), 'sca')
                        try
                            geom.feature(ftag).set('factor', num2str(scale_factor));
                            fprintf('        ✓ 缩放系数已更新为 %d (通过%s)\n', scale_factor, ftag);
                            scale_updated = true;
                            break;
                        catch
                        end
                    end
                end
            catch
            end
        end
        
        if ~scale_updated
            fprintf('        ⚠ 警告: 无法更新缩放系数,将使用模型中原有的设置\n');
        end
        
        % 直接运行完整几何序列 - 这是最可靠的方法
        % COMSOL会自动按顺序执行所有几何操作
        geom.run;
        fprintf('        ✓ 几何序列已构建\n');
        
        %% 配置边界
        fprintf('      配置边界...\n');
        configure_boundaries(model, comp_tag, lengthYAxis, pore_total_perimeter);
        
        %% 求解
        fprintf('      求解模型...\n');
        study_tags = model.study.tags;
        if length(study_tags) > 0
            study_tag = char(study_tags(1));
            model.study(study_tag).run;
        end
        
        %% 导出结果
        fprintf('      导出结果...\n');
        
        % 确保输出目录存在
        [output_dir, ~, ~] = fileparts(excel_output);
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
        
        % 查找导出节点并设置路径
        export_tags = model.result.export.tags;
        if length(export_tags) > 0
            xlsx_tag = char(export_tags(1));
            model.result.export(xlsx_tag).set('filename', excel_output);
            model.result.export(xlsx_tag).run;
            
            % 验证文件生成
            if exist(excel_output, 'file')
                fprintf('      ✓ 结果已导出: %s\n', excel_output);
                success = true;
            else
                fprintf('      × 导出文件未生成\n');
            end
        else
            fprintf('      × 未找到导出节点\n');
        end
        
        %% 导出mph文件 (如果启用)
        % 安全获取 export_mph 配置值
        should_export_mph = false;
        try
            export_mph_val = config.export_mph;
            % 转换为标量逻辑值
            if ~isempty(export_mph_val)
                should_export_mph = all(logical(export_mph_val(:)));
            end
        catch
            % 如果获取失败，默认不导出
            should_export_mph = false;
        end
        
        if should_export_mph
            fprintf('      导出mph文件...\n');
            
            try
                % 获取输入文件夹路径 (从pore_dxf获取)
                % 确保路径是字符数组
                pore_dxf_char = char(pore_dxf);
                excel_output_char = char(excel_output);
                
                [dxf_folder, ~, ~] = fileparts(pore_dxf_char);
                [parent_folder, ~, ~] = fileparts(dxf_folder);  % 返回上一级目录
                
                % 创建mph输出目录
                mph_output_dir = fullfile(parent_folder, 'mph_models');
                if exist(mph_output_dir, 'dir') == 0
                    mkdir(mph_output_dir);
                    fprintf('        创建mph输出目录: %s\n', mph_output_dir);
                end
                
                % 生成mph文件名 (使用与excel相同的命名规则)
                [~, excel_name, ~] = fileparts(excel_output_char);
                mph_filename = [excel_name, '.mph'];
                mph_output = fullfile(mph_output_dir, mph_filename);
                
                mphsave(model, mph_output);
                if exist(mph_output, 'file') == 2
                    file_info = dir(mph_output);
                    fprintf('        ✓ mph文件已导出: %s (%.2f MB)\n', mph_filename, file_info.bytes/1024/1024);
                else
                    fprintf('        × mph文件未生成\n');
                end
            catch ME_mph
                fprintf('        × mph导出失败: %s\n', ME_mph.message);
            end
        end
        
    catch ME
        fprintf('      × COMSOL处理错误: %s\n', ME.message);
        success = false;
    end
end

function perimeter = analyze_dxf_perimeter(dxf_file)
    % 分析DXF文件中的POLYLINE总周长
    
    perimeter = 0;
    
    try
        fid = fopen(dxf_file, 'r');
        if fid == -1
            return;
        end
        
        lines = {};
        while ~feof(fid)
            lines{end+1} = fgetl(fid);
        end
        fclose(fid);
        
        i = 1;
        while i <= length(lines)
            if strcmp(strtrim(char(lines{i})), 'POLYLINE')
                vertices = [];
                is_closed = 0;
                
                % 检查是否闭合
                for j = i+1:min(i+20, length(lines))
                    try
                        code = str2double(strtrim(char(lines{j})));
                        if code == 70
                            flag = str2double(strtrim(char(lines{j+1})));
                            is_closed = mod(flag, 2);
                            break;
                        elseif code == 0
                            break;
                        end
                    catch
                    end
                end
                
                % 提取顶点
                j = i + 1;
                while j <= min(length(lines), i + 5000)
                    try
                        if strcmp(strtrim(char(lines{j})), 'VERTEX')
                            for k = j+1:min(j+30, length(lines))
                                code = str2double(strtrim(char(lines{k})));
                                if code == 10
                                    x = str2double(strtrim(char(lines{k+1})));
                                    for m = k+2:min(k+6, length(lines))
                                        y_code = str2double(strtrim(char(lines{m})));
                                        if y_code == 20
                                            y = str2double(strtrim(char(lines{m+1})));
                                            vertices(end+1, :) = [x, y];
                                            break;
                                        end
                                    end
                                    break;
                                elseif code == 0
                                    break;
                                end
                            end
                        elseif strcmp(strtrim(char(lines{j})), 'SEQEND')
                            break;
                        end
                    catch
                    end
                    j = j + 1;
                end
                
                % 计算周长
                if size(vertices, 1) >= 2
                    poly_perimeter = 0;
                    for k = 1:size(vertices, 1)-1
                        seg_len = sqrt((vertices(k+1,1)-vertices(k,1))^2 + ...
                                      (vertices(k+1,2)-vertices(k,2))^2);
                        poly_perimeter = poly_perimeter + seg_len;
                    end
                    
                    if is_closed
                        seg_len = sqrt((vertices(1,1)-vertices(end,1))^2 + ...
                                      (vertices(1,2)-vertices(end,2))^2);
                        poly_perimeter = poly_perimeter + seg_len;
                    end
                    
                    perimeter = perimeter + poly_perimeter;
                end
            end
            i = i + 1;
        end
        
    catch
        perimeter = 0;
    end
end

function configure_boundaries(model, comp_tag, lengthYAxis, pore_total_perimeter)
    % 配置通量/源边界选择
    % 完全借鉴 process_geometry.m 中的逻辑
    
    fprintf('        配置通量/源边界选择...\n');
    
    try
        % === 步骤1: 找到正确的物理场标签 ===
        physics_tags = model.component(comp_tag).physics.tags;
        
        % 尝试多个可能的物理场标签
        possible_physics_tags = {'c', 'C', 'coefficient', 'pde', 'cf', 'cfeq'};
        physics_tag = '';
        
        for i = 1:length(physics_tags)
            ptag = char(physics_tags(i));
            if any(strcmp(ptag, possible_physics_tags))
                physics_tag = ptag;
                break;
            end
        end
        
        % 如果没找到,使用第一个物理场
        if isempty(physics_tag) && length(physics_tags) > 0
            physics_tag = char(physics_tags(1));
        end
        
        if isempty(physics_tag)
            fprintf('        ⚠ 未找到任何物理场\n');
            return;
        end
        
        fprintf('        物理场: %s\n', physics_tag);
        
        % === 步骤2: 找到通量/源特征 ===
        feature_tags = model.component(comp_tag).physics(physics_tag).feature.tags;
        
        % 尝试多个可能的通量/源特征标签
        possible_flux_tags = {'src1', 'flux1', 'fs1', 'source1', 'src', 'flux'};
        flux_tag = '';
        
        for i = 1:length(feature_tags)
            ftag = char(feature_tags(i));
            if contains(lower(ftag), 'src') || contains(lower(ftag), 'flux') || contains(lower(ftag), 'source')
                flux_tag = ftag;
                break;
            end
        end
        
        % 如果没找到,尝试默认标签
        if isempty(flux_tag)
            for i = 1:length(possible_flux_tags)
                try
                    temp = model.component(comp_tag).physics(physics_tag).feature(possible_flux_tags{i});
                    flux_tag = possible_flux_tags{i};
                    break;
                catch
                end
            end
        end
        
        if isempty(flux_tag)
            fprintf('        ⚠ 未找到通量/源特征\n');
            return;
        end
        
        fprintf('        通量/源特征: %s\n', flux_tag);
        
        % === 步骤3: 获取所有边界编号 ===
        % 先设置为all
        model.component(comp_tag).physics(physics_tag).feature(flux_tag).selection.all;
        
        % 获取选择的边界
        sel_obj = model.component(comp_tag).physics(physics_tag).feature(flux_tag).selection;
        
        try
            entities = sel_obj.entities(1); % 1D边界
            all_boundary_ids = double(entities);
            % 确保是行向量，便于后续处理
            all_boundary_ids = all_boundary_ids(:)';
        catch
            fprintf('        ⚠ 无法获取边界实体\n');
            return;
        end
        
        if isempty(all_boundary_ids)
            fprintf('        ⚠ 未识别到任何边界\n');
            return;
        end
        
        fprintf('        边界总数: %d\n', length(all_boundary_ids));
        fprintf('        边界范围: %d - %d\n', min(all_boundary_ids), max(all_boundary_ids));
        
        % === 步骤4: 计算需要排除的左右边界数量 ===
        % 公式: lengthYAxis / pore_total_perimeter = exclude_count / total_boundary_count
        fprintf('        === 边界排除计算 ===\n');
        fprintf('        lengthYAxis: %.6f cm\n', lengthYAxis);
        fprintf('        pore_total_perimeter: %.6f\n', pore_total_perimeter);
        fprintf('        边界总数: %d\n', length(all_boundary_ids));
        
        if pore_total_perimeter > 0 && ~isnan(lengthYAxis) && lengthYAxis > 0
            exclude_count_float = (lengthYAxis / pore_total_perimeter) * length(all_boundary_ids);
            exclude_count = round(exclude_count_float);
            fprintf('        计算比例: %.6f\n', lengthYAxis / pore_total_perimeter);
            fprintf('        计算排除数量: %.2f → %d\n', exclude_count_float, exclude_count);
        else
            fprintf('        ⚠ 无法计算排除比例,使用默认值\n');
            % 使用默认排除数量 (约占边界总数的5%)
            exclude_count = round(length(all_boundary_ids) * 0.05);
            fprintf('        默认排除数量: %d (5%%)\n', exclude_count);
        end
        
        % 确保排除数量在合理范围内 (不超过边界总数的1/3)
        exclude_count = max(1, min(exclude_count, floor(length(all_boundary_ids)/3)));
        fprintf('        最终排除数量: %d (每侧)\n', exclude_count);
        
        % === 步骤5: 确定左右边界范围 ===
        % 注意：这里是基于数组索引，而不是边界ID值
        % 因为边界ID可能不连续 (如 1-170, 172-181, ...)
        
        % 左边界: 数组的前 exclude_count 个元素
        left_end_idx = min(exclude_count, length(all_boundary_ids));
        left_boundaries = all_boundary_ids(1:left_end_idx);
        
        % 右边界: 数组的后 exclude_count 个元素
        right_start_idx = max(1, length(all_boundary_ids) - exclude_count + 1);
        right_boundaries = all_boundary_ids(right_start_idx:end);
        
        fprintf('        左边界: 索引 1-%d, ID范围 %d-%d\n', left_end_idx, min(left_boundaries), max(left_boundaries));
        fprintf('        右边界: 索引 %d-%d, ID范围 %d-%d\n', right_start_idx, length(all_boundary_ids), min(right_boundaries), max(right_boundaries));
        
        % === 步骤6: 应用边界排除 ===
        % 合并左右边界ID
        combined_exclude = [left_boundaries, right_boundaries];
        
        % 计算保留的边界
        boundaries_to_keep = setdiff(all_boundary_ids, combined_exclude);
        
        fprintf('        排除左边界数: %d\n', length(left_boundaries));
        fprintf('        排除右边界数: %d\n', length(right_boundaries));
        fprintf('        保留边界数: %d\n', length(boundaries_to_keep));
        
        % 尝试设置边界选择
        success = false;
        
        % 方法1: 先选择所有，再移除左右边界
        try
            model.component(comp_tag).physics(physics_tag).feature(flux_tag).selection.all;
            % 注意：remove 需要数值数组
            model.component(comp_tag).physics(physics_tag).feature(flux_tag).selection.remove(combined_exclude);
            fprintf('        ✓ 边界选择已应用 (使用selection.remove)\n');
            success = true;
        catch ME_remove
            fprintf('        × selection.remove失败: %s\n', ME_remove.message);
            
            % 方法2: 直接设置保留的边界
            try
                model.component(comp_tag).physics(physics_tag).feature(flux_tag).selection.set(boundaries_to_keep);
                fprintf('        ✓ 边界选择已应用 (使用selection.set)\n');
                success = true;
            catch ME_set
                fprintf('        × selection.set也失败: %s\n', ME_set.message);
            end
        end
        
        if success
            fprintf('        ✓ 通量/源边界配置完成\n');
        else
            fprintf('        ⚠ 边界选择设置失败，请手动检查\n');
        end
        
    catch ME
        fprintf('      ⚠ 边界配置警告: %s\n', ME.message);
    end
end
