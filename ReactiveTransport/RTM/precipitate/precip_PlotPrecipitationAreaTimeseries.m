function precip_PlotPrecipitationAreaTimeseries(csvFile, pngFile, markerTimesSeconds)
% precip_PlotPrecipitationAreaTimeseries - Plot net precipitation area curves.

if nargin < 3 || isempty(markerTimesSeconds)
    markerTimesSeconds = [];
end
if ~isfile(csvFile)
    error('RTSPHEM:Precipitate:MissingAreaTimeseriesCsv', ...
        'Precipitation area time-series CSV does not exist: %s', csvFile);
end

data = readtable(csvFile);
requiredColumns = {'time_s', 'total_net_solid_area_cm2', ...
    'first_pore_net_solid_area_cm2', 'first_three_pores_net_solid_area_cm2'};
for iColumn = 1:numel(requiredColumns)
    if ~ismember(requiredColumns{iColumn}, data.Properties.VariableNames)
        error('RTSPHEM:Precipitate:InvalidAreaTimeseriesCsv', ...
            'Missing required precipitation area column: %s', requiredColumns{iColumn});
    end
end

parentDir = fileparts(pngFile);
if ~isempty(parentDir) && ~exist(parentDir, 'dir')
    mkdir(parentDir);
end

fig = figure('Visible', 'off', 'Position', [100, 100, 900, 600]);
cleanupObj = onCleanup(@() close(fig));
timeMinutes = data.time_s ./ 60;
plot(timeMinutes, data.total_net_solid_area_cm2, '-o', 'LineWidth', 1.5);
hold on;
plot(timeMinutes, data.first_pore_net_solid_area_cm2, '-s', 'LineWidth', 1.5);
plot(timeMinutes, data.first_three_pores_net_solid_area_cm2, '-^', 'LineWidth', 1.5);
markerTimesMinutes = markerTimesSeconds(:).' ./ 60;
for iMarker = 1:numel(markerTimesMinutes)
    xline(markerTimesMinutes(iMarker), '--', sprintf('%.0f min', markerTimesMinutes(iMarker)), ...
        'LabelVerticalAlignment', 'bottom');
end
grid on;
xlabel('time [min]');
ylabel('net solid area / precipitated area [cm^2]');
legend({'total', 'first pore', 'first three pores'}, 'Location', 'best');
title('Precipitation area time series');
try
    exportgraphics(fig, pngFile, 'Resolution', 200);
catch
    saveas(fig, pngFile);
end
clear cleanupObj;
end
