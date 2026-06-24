function precip_InitializePrecipitationAreaTimeseries(csvFile)
% precip_InitializePrecipitationAreaTimeseries - Create precipitation area CSV.

parentDir = fileparts(csvFile);
if ~isempty(parentDir) && ~exist(parentDir, 'dir')
    mkdir(parentDir);
end

fid = fopen(csvFile, 'w');
if fid == -1
    error('RTSPHEM:Precipitate:AreaTimeseriesOpenFailed', ...
        'Cannot open precipitation area time-series CSV: %s', csvFile);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, ['timestep,time_s,total_net_solid_area_cm2,', ...
    'first_pore_net_solid_area_cm2,first_three_pores_net_solid_area_cm2,', ...
    'total_solid_area_cm2,first_pore_solid_area_cm2,', ...
    'first_three_pores_solid_area_cm2\n']);
clear cleanupObj;
end
