function diagnostic = precip_TestMixingWidth(fractionInletA, spec)
% precip_TestMixingWidth - Basic split-inlet symmetry and width diagnostic.
%
% Inputs:
%   fractionInletA - numY-by-numX field of CaCl2 inlet fraction.
%   spec           - benchmark spec.
%
% Output:
%   diagnostic     - mixing width and centerline symmetry metrics.

if size(fractionInletA, 1) ~= spec.numY || size(fractionInletA, 2) ~= spec.numX
    error('RTSPHEM:Precipitate:InvalidMixingField', ...
        'fractionInletA must match [numY, numX].');
end

symmetryResidual = fractionInletA + flipud(fractionInletA) - 1;
diagnostic = struct();
diagnostic.symmetryError = max(abs(symmetryResidual(:)));
diagnostic.width_cm = zeros(1, spec.numX);
yCenters = ((1:spec.numY)' - 0.5) .* spec.dy_cm;
for iX = 1:spec.numX
    profile = fractionInletA(:, iX);
    mixed = profile >= 0.1 & profile <= 0.9;
    if any(mixed)
        diagnostic.width_cm(iX) = max(yCenters(mixed)) - min(yCenters(mixed)) + spec.dy_cm;
    end
end
diagnostic.meanWidth_cm = mean(diagnostic.width_cm);
end
