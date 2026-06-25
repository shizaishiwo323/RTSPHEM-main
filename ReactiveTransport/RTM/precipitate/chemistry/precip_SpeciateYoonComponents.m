function [chem, metadata] = precip_SpeciateYoonComponents(samples, spec)
% precip_SpeciateYoonComponents - Dispatch Yoon aqueous speciation backend.
%
% The Yoon benchmark owns the Vaterite kinetic law and Vm update. Chemistry
% backends provide aqueous activities and saturation state only.

if nargin < 2 || isempty(spec)
    spec = precip_ZhangYoonBenchmarkSpec();
end

backend = getFieldOrDefault(spec, 'chemistryBackend', 'yoon_equilibrium');
backend = char(string(backend));
switch lower(backend)
    case 'yoon_equilibrium'
        chem = precip_YoonCarbonateEquilibrium(samples, spec);
        validateChemistryOutput(chem, samples);
        metadata = makeMetadata('yoon_equilibrium', ...
            'aqueous_speciation', 'local_yoon_carbonate_equilibrium');
    case 'iphreeqc_speciation'
        speciationFcn = getFieldOrDefault(spec, 'iphreeqcSpeciationFcn', []);
        if isempty(speciationFcn) || ~isa(speciationFcn, 'function_handle')
            error('RTSPHEM:Precipitate:ChemistryBackendUnavailable', ...
                ['iphreeqc_speciation requires spec.iphreeqcSpeciationFcn ' ...
                'in this explicit interface path. No IPhreeqcCOM server ' ...
                'is created automatically.']);
        end
        [chem, userMetadata] = callSpeciationFunction(speciationFcn, samples, spec);
        validateChemistryOutput(chem, samples);
        source = getFieldOrDefault(userMetadata, 'source', ...
            ['function_handle:' func2str(speciationFcn)]);
        metadata = makeMetadata('iphreeqc_speciation', ...
            'aqueous_speciation', source);
    otherwise
        error('RTSPHEM:Precipitate:InvalidChemistryBackend', ...
            'Unsupported Yoon chemistry backend: %s.', backend);
end
end

function [chem, metadata] = callSpeciationFunction(speciationFcn, samples, spec)
metadata = struct();
try
    outputCount = nargout(speciationFcn);
catch
    outputCount = -1;
end

if outputCount == 1
    chem = speciationFcn(samples, spec);
elseif outputCount == 0
    error('RTSPHEM:Precipitate:InvalidChemistryBackendOutput', ...
        'Injected speciation function must return a chemistry struct.');
else
    [chem, metadata] = speciationFcn(samples, spec);
end

if isempty(metadata)
    metadata = struct();
elseif ~isstruct(metadata)
    error('RTSPHEM:Precipitate:InvalidChemistryBackendOutput', ...
        'Injected speciation metadata must be a struct when provided.');
end
end

function validateChemistryOutput(chem, samples)
if ~isstruct(chem)
    error('RTSPHEM:Precipitate:InvalidChemistryBackendOutput', ...
        'Chemistry backend must return a struct.');
end

requiredFields = {'pH', 'aH', 'aCa', 'aCO3', 'aH2CO3', ...
    'omegaVaterite', 'saturationIndexVaterite'};
expectedSize = size(samples.Ca_total);
for iField = 1:numel(requiredFields)
    fieldName = requiredFields{iField};
    if ~isfield(chem, fieldName)
        error('RTSPHEM:Precipitate:InvalidChemistryBackendOutput', ...
            'Chemistry backend output is missing required field: %s.', ...
            fieldName);
    end
    value = chem.(fieldName);
    if ~isnumeric(value) || ~isequal(size(value), expectedSize) || ...
            any(~isfinite(value(:)))
        error('RTSPHEM:Precipitate:InvalidChemistryBackendOutput', ...
            ['Chemistry backend field %s must be numeric, finite, ' ...
            'and match the component array size.'], fieldName);
    end
end
end

function metadata = makeMetadata(backend, role, source)
metadata = struct();
metadata.backend = backend;
metadata.role = role;
metadata.source = source;
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
