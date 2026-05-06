%% extract_comsol_info.m
% 提取 COMSOL .mph 文件中的全部模拟信息并导出为 Markdown 文件
%
% 核心策略:
%   策略A - mphsave 导出 .m 脚本（最可靠，包含完整模型定义）
%   策略B - 直接 API 提取（补充信息：网格统计等运行时数据）
%
% 用法:
%   extract_comsol_info()
%   extract_comsol_info('path/to/model.mph')
%   extract_comsol_info('path/to/model.mph', 'output.md')

function extract_comsol_info(mph_file, output_md)

    %% 默认路径
    if nargin < 1 || isempty(mph_file)
        config = AutomationConfig();
        mph_file = config.mph_file;
    end
    if nargin < 2 || isempty(output_md)
        [fdir, fname, ~] = fileparts(mph_file);
        output_md = fullfile(fdir, [fname, '_info.md']);
    end

    fprintf('====================================================\n');
    fprintf('  COMSOL 模型信息提取工具\n');
    fprintf('====================================================\n');
    fprintf('模型路径: %s\n', mph_file);
    fprintf('输出路径: %s\n\n', output_md);

    if ~exist(mph_file, 'file')
        error('模型文件不存在: %s', mph_file);
    end

    %% 初始化 COMSOL 连接
    fprintf('[步骤 1/4] 连接 COMSOL...\n');
    comsol_init_local();

    %% 加载模型
    fprintf('[步骤 2/4] 加载模型...\n');
    import com.comsol.model.*
    import com.comsol.model.util.*
    model = mphload(mph_file);
    model_tag = char(model.tag);
    fprintf('  ✓ 模型已加载 (tag = %s)\n\n', model_tag);

    %% 策略A: 导出为 .m 文件（最完整的方式）
    fprintf('[步骤 3/4] 导出模型为 .m 脚本...\n');
    [fdir, fname, ~] = fileparts(mph_file);
    m_file = fullfile(fdir, [fname, '_export.m']);
    try
        mphsave(model, m_file);
        fprintf('  ✓ 已导出: %s\n', m_file);
        m_content = fileread(m_file);
    catch ME
        fprintf('  ✗ mphsave失败: %s\n', ME.message);
        fprintf('  → 尝试用 ModelUtil.showProgress 导出...\n');
        m_content = '';
    end

    %% 策略B: 直接 API 提取（带详细调试输出）
    fprintf('[步骤 4/4] 提取模型信息...\n');

    lines = {};
    lines = add(lines, '# COMSOL 模型信息报告');
    lines = add(lines, '');
    lines = add(lines, sprintf('**模型文件**: `%s`', mph_file));
    lines = add(lines, sprintf('**模型标签**: `%s`', model_tag));
    lines = add(lines, sprintf('**提取时间**: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS')));
    lines = add(lines, '');
    lines = add(lines, '---');
    lines = add(lines, '');

    % ── 0. 模型树概览 ──────────────────────────────────────────
    lines = add(lines, '## 0. 模型树概览');
    lines = add(lines, '');
    lines = dump_model_tree(lines, model);

    % ── 1. 全局参数 ─────────────────────────────────────────────
    lines = add(lines, '## 1. 全局参数 (Parameters)');
    lines = add(lines, '');
    lines = section_parameters(lines, model);

    % ── 2. 组件与几何 ───────────────────────────────────────────
    lines = add(lines, '## 2. 几何 (Geometry)');
    lines = add(lines, '');
    lines = section_geometry(lines, model);

    % ── 3. 定义（变量/探针/选择等）──────────────────────────────
    lines = add(lines, '## 3. 定义 (Definitions)');
    lines = add(lines, '');
    lines = section_definitions(lines, model);

    % ── 4. 材料 ─────────────────────────────────────────────────
    lines = add(lines, '## 4. 材料 (Materials)');
    lines = add(lines, '');
    lines = section_materials(lines, model);

    % ── 5. 物理场 ───────────────────────────────────────────────
    lines = add(lines, '## 5. 物理场与控制方程 (Physics)');
    lines = add(lines, '');
    lines = section_physics(lines, model);

    % ── 6. 网格 ─────────────────────────────────────────────────
    lines = add(lines, '## 6. 网格划分 (Mesh)');
    lines = add(lines, '');
    lines = section_mesh(lines, model);

    % ── 7. 算例/求解器 ──────────────────────────────────────────
    lines = add(lines, '## 7. 算例与求解器 (Study / Solver)');
    lines = add(lines, '');
    lines = section_study(lines, model);

    % ── 8. 结果 ─────────────────────────────────────────────────
    lines = add(lines, '## 8. 结果与后处理 (Results)');
    lines = add(lines, '');
    lines = section_results(lines, model);

    % ── 9. 导出 ─────────────────────────────────────────────────
    lines = add(lines, '## 9. 数据导出 (Export)');
    lines = add(lines, '');
    lines = section_export(lines, model);

    % ── 10. 完整模型 .m 脚本 ────────────────────────────────────
    if ~isempty(m_content)
        lines = add(lines, '## 10. 完整模型定义脚本 (mphsave 输出)');
        lines = add(lines, '');
        lines = add(lines, '以下是 `mphsave` 导出的完整 MATLAB 脚本，包含模型的所有设置细节:');
        lines = add(lines, '');
        lines = add(lines, '```matlab');
        % 按行写入
        m_lines = strsplit(m_content, '\n');
        for i = 1:length(m_lines)
            lines = add(lines, m_lines{i});
        end
        lines = add(lines, '```');
        lines = add(lines, '');
    end

    %% 写入文件
    fid = fopen(output_md, 'w', 'n', 'UTF-8');
    if fid == -1, error('无法创建输出文件: %s', output_md); end
    for i = 1:length(lines)
        fprintf(fid, '%s\n', lines{i});
    end
    fclose(fid);

    fprintf('\n====================================================\n');
    fprintf('  ✓ 完成！信息已保存到:\n  %s\n', output_md);
    fprintf('====================================================\n');
end

%% ================================================================
%  工具函数
%% ================================================================

function lines = add(lines, s)
    lines{end+1} = s;
end

function c = jtags(java_tags)
    % 安全地将 Java String[] 转为 MATLAB cell of char
    c = {};
    if isempty(java_tags), return; end
    try
        n = java_tags.length;  % Java array
        for i = 1:n
            c{end+1} = char(java_tags(i-1));  % 0-based
        end
    catch
        try
            c = cell(java_tags);
            for i = 1:length(c)
                c{i} = char(c{i});
            end
        catch
            try
                c = cellstr(java_tags);
            catch
                try
                    % 单元素
                    c = {char(java_tags)};
                catch
                    c = {};
                end
            end
        end
    end
end

% ──────────────────────────────────────────────────────────────────
%  0. 模型树概览
% ──────────────────────────────────────────────────────────────────
function lines = dump_model_tree(lines, model)
    % 列出顶层节点
    top_nodes = {'param', 'func', 'component', 'material', ...
                 'physics', 'mesh', 'study', 'sol', 'result'};
    lines = add(lines, '模型顶层节点探测:');
    lines = add(lines, '');
    for i = 1:length(top_nodes)
        nd = top_nodes{i};
        try
            obj = eval(['model.', nd]);
            tags = jtags(obj.tags);
            if ~isempty(tags)
                lines = add(lines, sprintf('- **%s**: %s', nd, strjoin(tags, ', ')));
                fprintf('  [树] %s: %s\n', nd, strjoin(tags, ', '));
            else
                lines = add(lines, sprintf('- **%s**: _(存在但无子节点)_', nd));
                fprintf('  [树] %s: (空)\n', nd);
            end
        catch ME_tree
            lines = add(lines, sprintf('- **%s**: _(不可访问: %s)_', nd, ME_tree.message));
            fprintf('  [树] %s: 不可访问 - %s\n', nd, ME_tree.message);
        end
    end

    % 列出 component 下的子节点
    comp_tags = {};
    try comp_tags = jtags(model.component.tags); catch; end
    for ci = 1:length(comp_tags)
        ctag = comp_tags{ci};
        lines = add(lines, sprintf('- **component/%s** 子节点:', ctag));
        sub_nodes = {'geom', 'mesh', 'physics', 'material', ...
                     'variable', 'probe', 'cpl', 'coordSystem'};
        for si = 1:length(sub_nodes)
            snd = sub_nodes{si};
            try
                sub_obj = model.component(ctag).(snd);
                stags = jtags(sub_obj.tags);
                if ~isempty(stags)
                    lines = add(lines, sprintf('  - %s: %s', snd, strjoin(stags, ', ')));
                end
            catch
            end
        end
    end

    % 同时检测 modelNode（某些 COMSOL 版本使用）
    try
        mn_tags = jtags(model.modelNode.tags);
        if ~isempty(mn_tags)
            lines = add(lines, sprintf('- **modelNode**: %s', strjoin(mn_tags, ', ')));
        end
    catch; end

    lines = add(lines, '');
end

% ──────────────────────────────────────────────────────────────────
%  1. 全局参数
% ──────────────────────────────────────────────────────────────────
function lines = section_parameters(lines, model)
    found_any = false;

    % 方法1: model.param.varnames (默认参数组)
    names = {};
    try
        names = jtags(model.param.varnames);
    catch ME1
        fprintf('  [参数] model.param.varnames 失败: %s\n', ME1.message);
    end

    if ~isempty(names)
        lines = add(lines, '| 参数名 | 表达式 | 描述 |');
        lines = add(lines, '|--------|--------|------|');
        for k = 1:length(names)
            nm = names{k};
            expr = ''; descr = '';
            try expr = char(model.param.get(nm)); catch; end
            try descr = char(model.param.descr(nm)); catch; end
            lines = add(lines, sprintf('| `%s` | `%s` | %s |', nm, expr, descr));
            found_any = true;
        end
        lines = add(lines, '');
    end

    % 方法2: 检查参数组
    try
        param_tags = jtags(model.param.tags);
        for gi = 1:length(param_tags)
            gtag = param_tags{gi};
            gnames = {};
            try gnames = jtags(model.param(gtag).varnames); catch; end
            if isempty(gnames), continue; end
            lines = add(lines, sprintf('### 参数组: `%s`', gtag));
            lines = add(lines, '| 参数名 | 表达式 | 描述 |');
            lines = add(lines, '|--------|--------|------|');
            for k = 1:length(gnames)
                nm = gnames{k};
                expr = ''; descr = '';
                try expr = char(model.param(gtag).get(nm)); catch; end
                try descr = char(model.param(gtag).descr(nm)); catch; end
                lines = add(lines, sprintf('| `%s` | `%s` | %s |', nm, expr, descr));
                found_any = true;
            end
            lines = add(lines, '');
        end
    catch; end

    if ~found_any
        lines = add(lines, '_（未找到全局参数, 可能在 .m 脚本中）_');
        lines = add(lines, '');
    end
end

% ──────────────────────────────────────────────────────────────────
%  2. 几何
% ──────────────────────────────────────────────────────────────────
function lines = section_geometry(lines, model)
    comp_tags = {};
    try comp_tags = jtags(model.component.tags); catch; end
    if isempty(comp_tags)
        lines = add(lines, '_（无组件信息, 见 .m 脚本）_');
        lines = add(lines, '');
        return;
    end

    for ci = 1:length(comp_tags)
        ctag = comp_tags{ci};
        geom_tags = {};
        try geom_tags = jtags(model.component(ctag).geom.tags); catch; end

        for gi = 1:length(geom_tags)
            gtag = geom_tags{gi};
            lines = add(lines, sprintf('### 组件 `%s` → 几何 `%s`', ctag, gtag));

            % 空间维度
            try
                geom_obj = model.component(ctag).geom(gtag);
                ndim = geom_obj.geomRep.numDomains;
                lines = add(lines, sprintf('- 域数量: %d', ndim));
            catch; end

            % 特征序列
            feat_tags = {};
            try feat_tags = jtags(model.component(ctag).geom(gtag).feature.tags); catch; end
            if ~isempty(feat_tags)
                lines = add(lines, '');
                lines = add(lines, '**几何特征序列**:');
                lines = add(lines, '| # | 标签 | 名称 |');
                lines = add(lines, '|---|------|------|');
                for fi = 1:length(feat_tags)
                    ftag = feat_tags{fi};
                    flabel = safe_label(model.component(ctag).geom(gtag).feature(ftag));
                    lines = add(lines, sprintf('| %d | `%s` | %s |', fi, ftag, flabel));
                end
            end
            lines = add(lines, '');
        end
    end
end

% ──────────────────────────────────────────────────────────────────
%  3. 定义（变量、探针、选择等）
% ──────────────────────────────────────────────────────────────────
function lines = section_definitions(lines, model)
    comp_tags = {};
    try comp_tags = jtags(model.component.tags); catch; end
    found_any = false;

    for ci = 1:length(comp_tags)
        ctag = comp_tags{ci};

        % 变量
        try
            var_tags = jtags(model.component(ctag).variable.tags);
            for vi = 1:length(var_tags)
                vtag = var_tags{vi};
                vlabel = safe_label(model.component(ctag).variable(vtag));
                lines = add(lines, sprintf('### 变量组 `%s` (%s)', vtag, vlabel));
                found_any = true;

                vnames = {};
                try vnames = jtags(model.component(ctag).variable(vtag).varnames); catch; end
                if ~isempty(vnames)
                    lines = add(lines, '| 变量名 | 表达式 | 描述 |');
                    lines = add(lines, '|--------|--------|------|');
                    for ni = 1:length(vnames)
                        nm = vnames{ni};
                        expr = ''; descr = '';
                        try expr = char(model.component(ctag).variable(vtag).get(nm)); catch; end
                        try descr = char(model.component(ctag).variable(vtag).descr(nm)); catch; end
                        lines = add(lines, sprintf('| `%s` | `%s` | %s |', nm, expr, descr));
                    end
                end
                lines = add(lines, '');
            end
        catch; end

        % 探针
        try
            prb_tags = jtags(model.component(ctag).probe.tags);
            for pi = 1:length(prb_tags)
                ptag = prb_tags{pi};
                plabel = safe_label(model.component(ctag).probe(ptag));
                lines = add(lines, sprintf('- 探针 `%s`: %s', ptag, plabel));
                found_any = true;
            end
            lines = add(lines, '');
        catch; end
    end

    % 全局函数
    try
        func_tags = jtags(model.func.tags);
        if ~isempty(func_tags)
            lines = add(lines, '### 全局函数');
            for fi = 1:length(func_tags)
                ftag = func_tags{fi};
                flabel = '';
                try flabel = safe_label(model.func(ftag)); catch; end
                lines = add(lines, sprintf('- `%s`: %s', ftag, flabel));
                found_any = true;
            end
            lines = add(lines, '');
        end
    catch; end

    if ~found_any
        lines = add(lines, '_（无额外定义, 见 .m 脚本）_');
        lines = add(lines, '');
    end
end

% ──────────────────────────────────────────────────────────────────
%  4. 材料
% ──────────────────────────────────────────────────────────────────
function lines = section_materials(lines, model)
    comp_tags = {};
    try comp_tags = jtags(model.component.tags); catch; end
    found_any = false;

    for ci = 1:length(comp_tags)
        ctag = comp_tags{ci};
        mat_tags = {};
        try mat_tags = jtags(model.component(ctag).material.tags); catch; end

        for mi = 1:length(mat_tags)
            mtag = mat_tags{mi};
            found_any = true;
            mlabel = safe_label(model.component(ctag).material(mtag));
            lines = add(lines, sprintf('### 材料 `%s` (标签: %s)', mlabel, mtag));

            % 属性组
            prop_tags = {};
            try prop_tags = jtags(model.component(ctag).material(mtag).propertyGroup.tags); catch; end
            for pi = 1:length(prop_tags)
                ptag = prop_tags{pi};
                plabel = safe_label(model.component(ctag).material(mtag).propertyGroup(ptag));
                lines = add(lines, sprintf('**属性组 `%s`** (%s):', ptag, plabel));

                % 尝试用 set/get 提取
                try
                    props = model.component(ctag).material(mtag).propertyGroup(ptag);
                    % 遍历可能的属性名
                    prop_names = {};
                    try prop_names = jtags(props.properties); catch; end
                    for ni = 1:length(prop_names)
                        pn = prop_names{ni};
                        pv = '';
                        try pv = char(props.getString(pn)); catch; end
                        if ~isempty(pv)
                            lines = add(lines, sprintf('- `%s` = `%s`', pn, pv));
                        end
                    end
                catch; end
            end
            lines = add(lines, '');
        end
    end

    if ~found_any
        lines = add(lines, '_（未定义材料, 见 .m 脚本）_');
        lines = add(lines, '');
    end
end

% ──────────────────────────────────────────────────────────────────
%  5. 物理场
% ──────────────────────────────────────────────────────────────────
function lines = section_physics(lines, model)
    comp_tags = {};
    try comp_tags = jtags(model.component.tags); catch; end
    found_any = false;

    for ci = 1:length(comp_tags)
        ctag = comp_tags{ci};
        phys_tags = {};
        try phys_tags = jtags(model.component(ctag).physics.tags); catch; end

        for pi = 1:length(phys_tags)
            ptag = phys_tags{pi};
            found_any = true;
            plabel = safe_label(model.component(ctag).physics(ptag));

            lines = add(lines, sprintf('### 物理场 `%s` (标签: %s)', plabel, ptag));

            % 因变量
            try
                fieldname = char(model.component(ctag).physics(ptag).fieldName);
                lines = add(lines, sprintf('- **因变量**: `%s`', fieldname));
            catch
                try
                    fields = jtags(model.component(ctag).physics(ptag).field.tags);
                    for ffi = 1:length(fields)
                        fn = char(model.component(ctag).physics(ptag).field(fields{ffi}).toString);
                        lines = add(lines, sprintf('- **场**: `%s` → %s', fields{ffi}, fn));
                    end
                catch; end
            end

            % 方程形式
            try
                eq_form = char(model.component(ctag).physics(ptag).prop('EquationForm').getString('form'));
                lines = add(lines, sprintf('- **方程形式**: %s', eq_form));
            catch; end

            % 物理场特征
            feat_tags = {};
            try feat_tags = jtags(model.component(ctag).physics(ptag).feature.tags); catch; end
            if ~isempty(feat_tags)
                lines = add(lines, '');
                lines = add(lines, '**物理场特征（方程 / 边界条件）**:');
                lines = add(lines, '| # | 标签 | 名称 | 关键设置 |');
                lines = add(lines, '|---|------|------|----------|');
                for fi = 1:length(feat_tags)
                    ftag = feat_tags{fi};
                    flabel = safe_label(model.component(ctag).physics(ptag).feature(ftag));
                    ks = extract_feature_settings(model, ctag, ptag, ftag);
                    lines = add(lines, sprintf('| %d | `%s` | %s | %s |', fi, ftag, flabel, ks));
                end
            end
            lines = add(lines, '');
        end
    end

    if ~found_any
        lines = add(lines, '_（无物理场信息, 见 .m 脚本）_');
        lines = add(lines, '');
    end
end

function s = extract_feature_settings(model, ctag, ptag, ftag)
    s = '';
    common_props = {'D', 'u', 'f', 'g', 'q', 'h', 'r', 'N', 'J', ...
                    'c', 'a', 'da', 'f0', 'flux', 'conc', ...
                    'T2', 'T1', 'rho', 'sigma', 'mu', 'D_eff'};
    parts = {};
    feat = model.component(ctag).physics(ptag).feature(ftag);
    for i = 1:length(common_props)
        try
            val = char(feat.getString(common_props{i}));
            if ~isempty(val) && ~strcmp(val, '0')
                parts{end+1} = sprintf('%s=%s', common_props{i}, val);
            end
        catch; end
    end
    if ~isempty(parts)
        s = strjoin(parts, '; ');
        if length(s) > 100, s = [s(1:97), '...']; end
    end
end

% ──────────────────────────────────────────────────────────────────
%  6. 网格
% ──────────────────────────────────────────────────────────────────
function lines = section_mesh(lines, model)
    comp_tags = {};
    try comp_tags = jtags(model.component.tags); catch; end
    found_any = false;

    for ci = 1:length(comp_tags)
        ctag = comp_tags{ci};
        mesh_tags = {};
        try mesh_tags = jtags(model.component(ctag).mesh.tags); catch; end

        for mi = 1:length(mesh_tags)
            mtag = mesh_tags{mi};
            found_any = true;
            mlabel = safe_label(model.component(ctag).mesh(mtag));
            lines = add(lines, sprintf('### 网格 `%s` (标签: %s)', mlabel, mtag));

            % 特征序列
            feat_tags = {};
            try feat_tags = jtags(model.component(ctag).mesh(mtag).feature.tags); catch; end
            if ~isempty(feat_tags)
                lines = add(lines, '| # | 标签 | 名称 |');
                lines = add(lines, '|---|------|------|');
                for fi = 1:length(feat_tags)
                    ftag = feat_tags{fi};
                    flabel = safe_label(model.component(ctag).mesh(mtag).feature(ftag));
                    % 尝试补充尺寸信息
                    size_str = '';
                    try
                        hmax = char(model.component(ctag).mesh(mtag).feature(ftag).getString('hmax'));
                        if ~isempty(hmax), size_str = [size_str, ' hmax=', hmax]; end
                    catch; end
                    try
                        hmin = char(model.component(ctag).mesh(mtag).feature(ftag).getString('hmin'));
                        if ~isempty(hmin), size_str = [size_str, ' hmin=', hmin]; end
                    catch; end
                    if ~isempty(size_str)
                        flabel = [flabel, ' [', strtrim(size_str), ']'];
                    end
                    lines = add(lines, sprintf('| %d | `%s` | %s |', fi, ftag, flabel));
                end
                lines = add(lines, '');
            end

            % 网格统计
            try
                stats = mphmeshstats(model, mtag);
                lines = add(lines, '**网格统计**:');
                if isfield(stats, 'numelem')
                    lines = add(lines, sprintf('- 单元总数: %d', sum(stats.numelem)));
                end
                if isfield(stats, 'numvtx')
                    lines = add(lines, sprintf('- 节点总数: %d', stats.numvtx));
                end
                if isfield(stats, 'types')
                    lines = add(lines, sprintf('- 单元类型: %s', strjoin(stats.types, ', ')));
                end
                if isfield(stats, 'minqual')
                    lines = add(lines, sprintf('- 最小质量: %.4f', stats.minqual));
                end
                if isfield(stats, 'meanqual')
                    lines = add(lines, sprintf('- 平均质量: %.4f', stats.meanqual));
                end
                lines = add(lines, '');
            catch ME_mesh
                fprintf('  [网格统计] %s: %s\n', mtag, ME_mesh.message);
            end
        end
    end

    if ~found_any
        lines = add(lines, '_（无网格信息, 见 .m 脚本）_');
        lines = add(lines, '');
    end
end

% ──────────────────────────────────────────────────────────────────
%  7. 算例与求解器
% ──────────────────────────────────────────────────────────────────
function lines = section_study(lines, model)
    study_tags = {};
    try study_tags = jtags(model.study.tags); catch; end
    if isempty(study_tags)
        lines = add(lines, '_（无算例信息, 见 .m 脚本）_');
        lines = add(lines, '');
        return;
    end

    for si = 1:length(study_tags)
        stag = study_tags{si};
        slabel = safe_label(model.study(stag));
        lines = add(lines, sprintf('### 算例 `%s` (标签: %s)', slabel, stag));

        % 算例步骤
        step_tags = {};
        try step_tags = jtags(model.study(stag).feature.tags); catch; end
        if ~isempty(step_tags)
            lines = add(lines, '**算例步骤**:');
            lines = add(lines, '| # | 标签 | 名称 | 关键参数 |');
            lines = add(lines, '|---|------|------|----------|');
            for fi = 1:length(step_tags)
                ftag = step_tags{fi};
                flabel = safe_label(model.study(stag).feature(ftag));
                ki = get_study_step_info(model, stag, ftag);
                lines = add(lines, sprintf('| %d | `%s` | %s | %s |', fi, ftag, flabel, ki));
            end
            lines = add(lines, '');
        end
    end

    % 求解器配置
    sol_tags = {};
    try sol_tags = jtags(model.sol.tags); catch; end
    if ~isempty(sol_tags)
        lines = add(lines, '### 求解器配置');
        for si = 1:length(sol_tags)
            stag_sol = sol_tags{si};
            sol_label = safe_label(model.sol(stag_sol));
            lines = add(lines, sprintf('**求解器 `%s` (%s)**:', stag_sol, sol_label));

            sfeat_tags = {};
            try sfeat_tags = jtags(model.sol(stag_sol).feature.tags); catch; end
            if ~isempty(sfeat_tags)
                lines = add(lines, '| # | 标签 | 名称 | 额外信息 |');
                lines = add(lines, '|---|------|------|----------|');
                for fi = 1:length(sfeat_tags)
                    ftag = sfeat_tags{fi};
                    flabel = safe_label(model.sol(stag_sol).feature(ftag));
                    extra = get_solver_feature_info(model, stag_sol, ftag);
                    lines = add(lines, sprintf('| %d | `%s` | %s | %s |', fi, ftag, flabel, extra));
                end
                lines = add(lines, '');
            end
        end
    end
end

function s = get_study_step_info(model, stag, ftag)
    s = '';
    props = {'tlist', 'pname', 'plistarr', 'punit', 'freq', 'freqlist', ...
             'rtol', 'timestep', 'endtime'};
    parts = {};
    feat = model.study(stag).feature(ftag);
    for i = 1:length(props)
        try
            val = char(feat.getString(props{i}));
            if ~isempty(val)
                parts{end+1} = sprintf('%s=%s', props{i}, val);
            end
        catch; end
    end
    s = strjoin(parts, '; ');
    if length(s) > 120, s = [s(1:117), '...']; end
end

function s = get_solver_feature_info(model, stag_sol, ftag)
    s = '';
    props = {'tol', 'reltol', 'abstol', 'maxiter', 'tstepsBDF', ...
             'linsolver', 'initstep', 'maxstep', 'timestep'};
    parts = {};
    feat = model.sol(stag_sol).feature(ftag);
    for i = 1:length(props)
        try
            val = char(feat.getString(props{i}));
            if ~isempty(val) && ~any(strcmp(val, {'0', 'none', ''}))
                parts{end+1} = sprintf('%s=%s', props{i}, val);
            end
        catch; end
    end
    if ~isempty(parts)
        s = strjoin(parts, '; ');
        if length(s) > 100, s = [s(1:97), '...']; end
    end
end

% ──────────────────────────────────────────────────────────────────
%  8. 结果
% ──────────────────────────────────────────────────────────────────
function lines = section_results(lines, model)
    res_tags = {};
    try res_tags = jtags(model.result.tags); catch; end

    if isempty(res_tags)
        lines = add(lines, '_（无结果节点, 见 .m 脚本）_');
        lines = add(lines, '');
        return;
    end

    lines = add(lines, '| # | 标签 | 名称 |');
    lines = add(lines, '|---|------|------|');
    for ri = 1:length(res_tags)
        rtag = res_tags{ri};
        rlabel = safe_label(model.result(rtag));
        lines = add(lines, sprintf('| %d | `%s` | %s |', ri, rtag, rlabel));
    end
    lines = add(lines, '');

    % 数据表
    tbl_tags = {};
    try tbl_tags = jtags(model.result.table.tags); catch; end
    if ~isempty(tbl_tags)
        lines = add(lines, '**数据表**:');
        for ti = 1:length(tbl_tags)
            ttag = tbl_tags{ti};
            tlabel = safe_label(model.result.table(ttag));
            lines = add(lines, sprintf('- `%s`: %s', ttag, tlabel));
        end
        lines = add(lines, '');
    end
end

% ──────────────────────────────────────────────────────────────────
%  9. 导出
% ──────────────────────────────────────────────────────────────────
function lines = section_export(lines, model)
    exp_tags = {};
    try exp_tags = jtags(model.result.export.tags); catch; end

    if isempty(exp_tags)
        lines = add(lines, '_（无导出节点, 见 .m 脚本）_');
        lines = add(lines, '');
        return;
    end

    lines = add(lines, '| # | 标签 | 名称 | 文件名 |');
    lines = add(lines, '|---|------|------|--------|');
    for ei = 1:length(exp_tags)
        etag = exp_tags{ei};
        elabel = safe_label(model.result.export(etag));
        efile = '';
        try efile = char(model.result.export(etag).getString('filename')); catch; end
        lines = add(lines, sprintf('| %d | `%s` | %s | `%s` |', ei, etag, elabel, efile));
    end
    lines = add(lines, '');
end

% ──────────────────────────────────────────────────────────────────
%  安全获取 label
% ──────────────────────────────────────────────────────────────────
function lbl = safe_label(obj)
    lbl = '';
    try lbl = char(obj.label); catch; end
    if isempty(lbl)
        try lbl = char(obj.tag); catch; end
    end
    if isempty(lbl)
        try lbl = char(obj.toString); catch; end
    end
end

% ──────────────────────────────────────────────────────────────────
%  本地 COMSOL 初始化（不依赖外部 comsol_init.m 返回值）
% ──────────────────────────────────────────────────────────────────
function comsol_init_local()
    comsol_path = 'C:\Program Files\COMSOL\COMSOL63\Multiphysics';
    mli_path = fullfile(comsol_path, 'mli');
    if exist(mli_path, 'dir')
        addpath(mli_path);
    end

    import com.comsol.model.*
    import com.comsol.model.util.*

    connected = false;

    % 方法1: 连接到本地已运行的服务器（默认端口 2036）
    try
        mphstart('localhost', 2036);
        fprintf('  ✓ 已连接 COMSOL 服务器（端口 2036）\n');
        connected = true;
    catch ME1
        % "Already connected" 视为成功
        if contains(ME1.message, 'Already connected') || ...
           contains(ME1.message, 'already connected')
            fprintf('  ✓ COMSOL 已处于连接状态\n');
            connected = true;
        else
            fprintf('  ✗ 端口 2036 连接失败: %s\n', ME1.message);
        end
    end

    % 方法2: 尝试其他常用端口
    if ~connected
        ports = [2037, 2038, 3036];
        for i = 1:length(ports)
            try
                mphstart('localhost', ports(i));
                fprintf('  ✓ 已连接 COMSOL 服务器（端口 %d）\n', ports(i));
                connected = true;
                break;
            catch ME2
                if contains(ME2.message, 'Already connected') || ...
                   contains(ME2.message, 'already connected')
                    fprintf('  ✓ COMSOL 已处于连接状态\n');
                    connected = true;
                    break;
                else
                    fprintf('  ✗ 端口 %d 连接失败\n', ports(i));
                end
            end
        end
    end

    % 方法3: 尝试直接创建模型（内置模式）
    if ~connected
        try
            ModelUtil.create('TestConn');
            ModelUtil.remove('TestConn');
            fprintf('  ✓ COMSOL 内置模式可用\n');
            connected = true;
        catch
        end
    end

    % 方法4: 启动新的服务器实例
    if ~connected
        try
            mphstart();
            fprintf('  ✓ COMSOL 服务器已启动（mphstart）\n');
            connected = true;
        catch ME4
            if contains(ME4.message, 'Already connected') || ...
               contains(ME4.message, 'already connected')
                fprintf('  ✓ COMSOL 已处于连接状态\n');
                connected = true;
            else
                error('无法连接 COMSOL：%s\n请先启动 COMSOL with MATLAB 服务器', ME4.message);
            end
        end
    end
end
