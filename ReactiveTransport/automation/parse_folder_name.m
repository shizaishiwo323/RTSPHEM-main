function params = parse_folder_name(folder_input, metadata_filename)
    % parse_folder_name - Read RTM run parameters for NMR automation.
    %
    % New runs are driven by run_metadata.json. Legacy folders named like
    % dissolution_results-Da_... are still supported only when JSON metadata
    % is missing or cannot be decoded.

    if nargin < 2 || strlength(string(metadata_filename)) == 0
        metadata_filename = 'run_metadata.json';
    end

    params = default_params();

    folder_path = char(folder_input);
    if exist(folder_path, 'dir')
        [~, folder_name] = fileparts(folder_path);
        metadata_file = fullfile(folder_path, char(metadata_filename));
        if exist(metadata_file, 'file')
            [params, metadata_ok] = parse_metadata_file(metadata_file, params);
            if metadata_ok
                return;
            end
            warning('run_metadata.json 不可用，将回退到旧版文件夹名解析: %s', metadata_file);
        end
    else
        folder_name = char(folder_input);
    end

    params = parse_legacy_folder_name(folder_name, params);
    params.parameterSource = 'folder_name';
end

function params = parse_legacy_folder_name(folder_name, params)
    name = strrep(folder_name, 'dissolution_results-', '');
    num_pattern = '([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)';

    da_match = regexp(name, ['Da_', num_pattern], 'tokens');
    if ~isempty(da_match)
        params.Da = str2double(da_match{1}{1});
    end

    pe_match = regexp(name, ['Pe_', num_pattern], 'tokens');
    if ~isempty(pe_match)
        params.Pe = str2double(pe_match{1}{1});
    end

    l_match = regexp(name, ['_L_', num_pattern], 'tokens');
    if ~isempty(l_match)
        params.L = str2double(l_match{1}{1});
    end

    x_match = regexp(name, ['lengthXAxis_', num_pattern], 'tokens');
    if ~isempty(x_match)
        params.lengthXAxis = str2double(x_match{1}{1});
    end

    y_match = regexp(name, ['lengthYAxis_', num_pattern], 'tokens');
    if ~isempty(y_match)
        params.lengthYAxis = str2double(y_match{1}{1});
    end

    if contains(name, 'square')
        params.layoutType = 'square';
    elseif contains(name, 'hex')
        params.layoutType = 'hex';
    elseif contains(name, 'random')
        params.layoutType = 'random';
    end
end

function params = default_params()
    params = struct();
    params.Da = NaN;
    params.Pe = NaN;
    params.L = NaN;
    params.lengthXAxis = NaN;
    params.lengthYAxis = NaN;
    params.layoutType = 'unknown';
    params.parameterSource = 'unknown';
    params.metadataFile = '';
    params.schemaVersion = '';
    params.runId = '';
end

function [params, ok] = parse_metadata_file(metadata_file, params)
    ok = false;

    fid = fopen(metadata_file, 'r');
    if fid == -1
        error('无法读取元数据文件: %s', metadata_file);
    end
    cleaner = onCleanup(@() fclose(fid));
    raw = fread(fid, '*char')';
    clear cleaner;

    if isempty(strtrim(raw))
        warning('run_metadata.json 为空: %s', metadata_file);
        return;
    end

    try
        metadata = jsondecode(raw);
    catch ME
        warning('run_metadata.json 解析失败: %s (%s)', metadata_file, ME.message);
        return;
    end

    candidates = collect_parameter_structs(metadata);

    params.Da = get_json_number(candidates, { ...
        'Da', 'damkohlerNumber', 'DamkohlerNumber', 'Damkohler', ...
        'damkohler_number', 'damkohlerNumber_value'}, params.Da);
    params.Pe = get_json_number(candidates, { ...
        'Pe', 'pecletNumber', 'PecletNumber', 'Peclet', ...
        'peclet_number', 'pecletNumber_value'}, params.Pe);
    params.L = get_json_number(candidates, { ...
        'characteristicLength_cm', 'characteristicLength', ...
        'characteristic_length_cm', 'characteristicLengthCm', ...
        'L_cm', 'L'}, params.L);
    params.lengthXAxis = get_json_number(candidates, { ...
        'lengthXAxis_cm', 'lengthXAxis', 'length_x_axis_cm', ...
        'lengthXAxisCm', 'x_length_cm', 'domainLengthX_cm'}, params.lengthXAxis);
    params.lengthYAxis = get_json_number(candidates, { ...
        'lengthYAxis_cm', 'lengthYAxis', 'length_y_axis_cm', ...
        'lengthYAxisCm', 'y_length_cm', 'domainLengthY_cm'}, params.lengthYAxis);

    layout_value = get_json_string(candidates, { ...
        'layoutType', 'layout_type', 'geometryLayout', 'geometry', 'layout'}, '');
    if strlength(string(layout_value)) > 0
        params.layoutType = char(layout_value);
    end

    params.parameterSource = 'run_metadata.json';
    params.metadataFile = metadata_file;
    params.schemaVersion = get_json_string({metadata}, {'schema_version', 'schemaVersion'}, '');
    params.runId = get_json_string({metadata}, {'run_id', 'runId'}, '');

    warn_missing_metadata_fields(params, metadata_file);
    ok = true;
end

function candidates = collect_parameter_structs(metadata)
    candidates = {metadata};

    nested_names = { ...
        'parameters', 'params', 'rtm_parameters', 'rtmParameters', ...
        'simulation_parameters', 'simulationParameters', 'config', ...
        'domain', 'geometry_parameters', 'geometryParameters'};

    for i = 1:length(nested_names)
        name = nested_names{i};
        if isfield(metadata, name) && isstruct(metadata.(name))
            candidates{end+1} = metadata.(name); %#ok<AGROW>
        end
    end
end

function value = get_json_number(candidates, field_names, default_value)
    value = default_value;

    for i = 1:length(candidates)
        data = candidates{i};
        for j = 1:length(field_names)
            field_name = field_names{j};
            if isfield(data, field_name) && ~isempty(data.(field_name))
                converted = to_number(data.(field_name));
                if ~isnan(converted)
                    value = converted;
                    return;
                end
            end
        end
    end
end

function value = get_json_string(candidates, field_names, default_value)
    value = default_value;

    for i = 1:length(candidates)
        data = candidates{i};
        for j = 1:length(field_names)
            field_name = field_names{j};
            if isfield(data, field_name) && ~isempty(data.(field_name))
                raw_value = data.(field_name);
                if isstring(raw_value) || ischar(raw_value)
                    value = char(raw_value);
                    return;
                elseif isnumeric(raw_value) || islogical(raw_value)
                    value = num2str(raw_value);
                    return;
                end
            end
        end
    end
end

function value = to_number(raw_value)
    value = NaN;

    if isnumeric(raw_value) && isscalar(raw_value)
        value = double(raw_value);
    elseif islogical(raw_value) && isscalar(raw_value)
        value = double(raw_value);
    elseif isstring(raw_value) || ischar(raw_value)
        value = str2double(char(raw_value));
    end
end

function warn_missing_metadata_fields(params, metadata_file)
    missing = {};

    if isnan(params.Da); missing{end+1} = 'Da'; end
    if isnan(params.Pe); missing{end+1} = 'Pe'; end
    if isnan(params.L); missing{end+1} = 'L/characteristicLength_cm'; end
    if isnan(params.lengthXAxis); missing{end+1} = 'lengthXAxis_cm'; end
    if isnan(params.lengthYAxis); missing{end+1} = 'lengthYAxis_cm'; end
    if strcmp(params.layoutType, 'unknown'); missing{end+1} = 'layoutType'; end

    if ~isempty(missing)
        warning('run_metadata.json 缺少部分参数 (%s): %s', strjoin(missing, ', '), metadata_file);
    end
end
