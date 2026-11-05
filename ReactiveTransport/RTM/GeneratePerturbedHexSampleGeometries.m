function summary = GeneratePerturbedHexSampleGeometries(options)
% GeneratePerturbedHexSampleGeometries
% Reproducible 600 x 400 um circular-grain geometries for M14 samples.
% Uses a perturbed hexagonal lattice and Voronoi local-porosity index Iv.
% Kxx is the upstream RTSPHEM periodic Stokes cell-problem permeability.

if nargin < 1, options = struct(); end
rtmDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(rtmDir));
addpath(genpath(fullfile(projectRoot, 'ReactiveTransport')));

outDir = getOpt(options, 'outputDir', fullfile(projectRoot, 'outputs', ...
    'se_geometry_calibration_600x400um'));
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Values transcribed from Wang_SE_experiment/碳酸岩样品/照片/孔隙率.jpg.
targetPorosity = getOpt(options, 'targetPorosity', [0.382; 0.374; 0.342; 0.356]);
validateattributes(targetPorosity, {'numeric'}, {'vector','numel',4,'>',0,'<',1});
targets = table(["M14_24"; "M14_25"; "M14_26"; "M14_34"], ...
    targetPorosity(:), [124.3; 107.6; 59.7; 108.7], ...
    'VariableNames', {'sampleName','targetPorosity','targetKxx_mD'});

Lx = 0.060; Ly = 0.040;                 % 600 x 400 um, cm
radius = getOpt(options, 'circleRadius_cm', 0.003); % 30 um, no scale factor
nx = getOpt(options, 'permeabilityNx', 192);
ny = getOpt(options, 'permeabilityNy', 128);
nSeeds = getOpt(options, 'candidateCount', 60);
phiTol = getOpt(options, 'porosityTolerance', 0.008);
kTol = getOpt(options, 'permeabilityRelativeTolerance', 0.12);
jitterMaxFactor = getOpt(options, 'jitterMaxFactor', 0.50);
minimumThroatFactor = getOpt(options, 'minimumThroatFactor', 0.30);
randomWalkSweeps = getOpt(options, 'randomWalkSweeps', 1);
polydispersityCvMax = getOpt(options, 'polydispersityCvMax', 0.28);
calibrationMode = string(getOpt(options, 'calibrationMode', 'porosity_only_2d_k'));
if ~ismember(calibrationMode, ["porosity_only_2d_k", "joint_phi_kxx"])
    error('RTM:GeometryCalibration:UnknownMode', ...
        'calibrationMode must be porosity_only_2d_k or joint_phi_kxx.');
end
validateattributes(phiTol, {'numeric'}, {'scalar','positive','<=',0.02});
validateattributes(kTol, {'numeric'}, {'scalar','positive','<=',0.25});
validateattributes(radius, {'numeric'}, {'scalar','positive','<',Ly/8});
rng(getOpt(options, 'rngSeed', 20260730), 'twister');

summary = targets;
summary.actualPorosity = nan(height(targets),1);
summary.Kxx_mD = nan(height(targets),1);
summary.Iv = nan(height(targets),1);
summary.seed = nan(height(targets),1);
summary.geometryFile = strings(height(targets),1);
summary.calibrationMode = repmat(calibrationMode, height(targets), 1);
summary.accepted = false(height(targets),1);

for i = 1:height(targets)
    target = targets(i,:);
    best = struct('score', inf);
    validCandidates = {};
    % Porosity is controlled physically by lattice spacing; disorder controls
    % local porosity and preferential-flow contrast. No post-hoc K scaling.
    spacing0 = spacingForPorosity(target.targetPorosity, radius);
    % Finite 600 x 400 um windows quantise the number of grains, so explore
    % a broad physical spacing range instead of trusting infinite-lattice phi.
    spacingGrid = linspace(max(2e-5, 0.15*spacing0), 1.35*spacing0, 13);
    for candidate = 1:nSeeds
        spacing = spacingGrid(1 + mod(candidate-1, numel(spacingGrid)));
        base = hexCenters(Lx, Ly, radius, spacing);
        radii = matchedRadii(base, radius, target.targetPorosity, Lx, Ly, ...
            polydispersityCvMax * rand);
        if isempty(radii), continue; end
        jitter = radius * (0.02 + (jitterMaxFactor-0.02) * rand);
        seed = randi(2^31-1);
        rng(seed, 'twister');
        centers = perturbCenters(base, radii, Lx, Ly, minimumThroatFactor*spacing, jitter, randomWalkSweeps);
        [~, Iv] = voronoiPorosity(centers, radii, Lx, Ly);
        phi = rasterPorosity(centers, radii, Lx, Ly, 600, 400);
        if ~isfinite(phi) || ~isfinite(Iv), continue; end
        % Screen φ before the costly Stokes cell problem.
        if abs(phi-target.targetPorosity) > 2.5*phiTol, continue; end
        K = upstreamKxx(centers, radii, Lx, Ly, nx, ny);
        if ~isfinite(K) || K <= 0, continue; end
        candidateRecord = struct('centers',centers,'spacing',spacing, ...
            'radii',radii,'jitter',jitter,'seed',seed,'phi',phi,'Iv',Iv,'Kxx',K);
        validCandidates{end+1} = candidateRecord; %#ok<AGROW>
        score = abs(phi-target.targetPorosity)/phiTol;
        if calibrationMode == "joint_phi_kxx"
            score = score + abs(log(K/target.targetKxx_mD))/log(1+kTol);
        end
        if score < best.score
            best = candidateRecord;
            best.score = score;
        end
    end
    if calibrationMode == "porosity_only_2d_k" && ~isempty(validCandidates)
        phiError = cellfun(@(x) abs(x.phi-target.targetPorosity), validCandidates);
        eligible = validCandidates(phiError <= phiTol);
        if isempty(eligible), eligible = validCandidates; end
        logK = cellfun(@(x) x.Kxx, eligible);
        logK = log(logK);
        [~, representativeIndex] = min(abs(logK - median(logK)));
        best = eligible{representativeIndex};
    end
    if ~isfield(best, 'centers')
        error('RTM:GeometryCalibration:NoCandidate', ...
            'No valid candidate for %s; increase candidateCount.', target.sampleName);
    end

    circleCenters = best.centers; %#ok<NASGU>
    circleRadii = best.radii; %#ok<NASGU>
    circleRadius = mean(circleRadii); %#ok<NASGU>
    lengthXAxis = Lx; lengthYAxis = Ly; %#ok<NASGU>
    targetInitialPorosity = target.targetPorosity; %#ok<NASGU>
    finalPorosity = best.phi; %#ok<NASGU>
    targetAvgSpacing = best.spacing; %#ok<NASGU>
    finalAvgSpacing = best.spacing; %#ok<NASGU>
    circleSpacing = best.spacing; %#ok<NASGU>
    minThroatRandom = 0.30 * best.spacing; %#ok<NASGU>
    sampleName = target.sampleName; %#ok<NASGU>
    experimentalPorosity = target.targetPorosity; %#ok<NASGU>
    experimentalPermeability_mD = target.targetKxx_mD; %#ok<NASGU>
    Kxx_mD = best.Kxx; %#ok<NASGU>
    Iv = best.Iv; %#ok<NASGU>
    lvRandomnessIndex = Iv; %#ok<NASGU>
    jitterAmplitude_cm = best.jitter; %#ok<NASGU>
    seed = best.seed; %#ok<NASGU>
    calibrationMethod = "perturbed_hex_voronoi_Iv_upstream_RTSPHEM_Kxx"; %#ok<NASGU>
    calibrationModeSaved = calibrationMode; %#ok<NASGU>
    geometryFile = fullfile(outDir, sprintf('%s_perturbed_hex_600x400um.mat', target.sampleName));
    save(geometryFile, 'circleCenters','circleRadii','circleRadius', ...
        'lengthXAxis','lengthYAxis','targetInitialPorosity','finalPorosity', ...
        'targetAvgSpacing','finalAvgSpacing','circleSpacing','minThroatRandom', ...
        'sampleName','experimentalPorosity','experimentalPermeability_mD', ...
        'Kxx_mD','Iv','lvRandomnessIndex','jitterAmplitude_cm','seed', ...
        'calibrationMethod','calibrationModeSaved','-v7');
    writeGeometryFigure(fullfile(outDir, sprintf('%s_perturbed_hex_600x400um.png', target.sampleName)), ...
        circleCenters, circleRadii, Lx, Ly, target, best);
    writeGeometryFigure(fullfile(outDir, sprintf('%s_perturbed_hex_600x400um.svg', target.sampleName)), ...
        circleCenters, circleRadii, Lx, Ly, target, best);

    summary.actualPorosity(i) = best.phi;
    summary.Kxx_mD(i) = best.Kxx;
    summary.Iv(i) = best.Iv;
    summary.seed(i) = best.seed;
    summary.geometryFile(i) = string(geometryFile);
    summary.accepted(i) = abs(best.phi-target.targetPorosity) <= phiTol && ...
        (calibrationMode == "porosity_only_2d_k" || ...
        abs(best.Kxx-target.targetKxx_mD)/target.targetKxx_mD <= kTol);
end
writetable(summary, fullfile(outDir, 'geometry_calibration_summary.csv'));
writetable(summary, fullfile(outDir, 'geometry_calibration_summary.xlsx'));
if ~all(summary.accepted)
    error('RTM:GeometryCalibration:TargetMiss', ...
        'At least one geometry misses the stated phi/Kxx tolerances.');
end
end

function value = getOpt(s, name, fallback)
if isfield(s,name) && ~isempty(s.(name)), value = s.(name); else, value = fallback; end
end

function spacing = spacingForPorosity(phi, r)
% Circle area / triangular lattice cell area: 1-phi = pi*r^2/(sqrt(3)/2*d^2).
d = r * sqrt(2*pi/(sqrt(3)*(1-phi)));
spacing = d - 2*r;
end

function centers = hexCenters(Lx, Ly, r, spacing)
d = 2*r + spacing; dy = sqrt(3)/2*d;
% Keep all physical circles inside the cell. Repetition then has a pore seam
% rather than unmatched half-grains, which is compatible with a periodic cell.
margin = 1.4*r;
xs = margin:d:(Lx-margin); ys = margin:dy:(Ly-margin);
centers = zeros(0,2);
for j = 1:numel(ys)
    row = xs + mod(j-1,2)*d/2;
    row = row(row <= Lx);
    centers = [centers; row(:), repmat(ys(j),numel(row),1)]; %#ok<AGROW>
end
end

function centers = perturbCenters(base, radii, Lx, Ly, minThroat, jitter, sweeps)
centers = base;
for sweep = 1:sweeps
order = randperm(size(base,1));
for i = order
    original = centers(i,:);
    r = radii(i);
    if original(1) <= r || original(1) >= Lx-r || original(2) <= r || original(2) >= Ly-r
        continue
    end
    for trial = 1:100
        candidate = original + jitter*randn(1,2);
        if candidate(1)<r || candidate(1)>Lx-r || candidate(2)<r || candidate(2)>Ly-r, continue; end
        other = centers; other(i,:)=[];
        otherR = radii; otherR(i)=[];
        if all(vecnorm(other-candidate,2,2) >= r + otherR + minThroat)
            centers(i,:) = candidate; break
        end
    end
end
end
end

function [phi, Iv] = voronoiPorosity(centers, radii, Lx, Ly)
% Luisa method: Voronoi local porosity and an internal 4r boundary cut-off.
mirror = [centers; [-centers(:,1), centers(:,2)]; [2*Lx-centers(:,1),centers(:,2)]; ...
    [centers(:,1),-centers(:,2)]; [centers(:,1),2*Ly-centers(:,2)]];
mirror = uniquetol(mirror, 1e-12, 'ByRows', true);
[V,C] = voronoin(mirror);
domain = polyshape([0 Lx Lx 0],[0 0 Ly Ly]);
local = nan(size(centers,1),1); areaCell = nan(size(local));
for i=1:size(centers,1)
    pointId=find(all(abs(mirror-centers(i,:))<1e-12,2),1);
    id=C{pointId}; if any(id==1) || numel(id)<3, continue; end
    p=intersect(polyshape(V(id,1),V(id,2)),domain); if isempty(p.Vertices), continue; end
    a=area(p); if a<=0, continue; end
    local(i)=1-pi*radii(i)^2/a; areaCell(i)=a;
end
valid=isfinite(local); phi=1-sum(pi*radii(valid).^2)/sum(areaCell(valid));
r=max(radii); inner=centers(:,1)>4*r & centers(:,1)<Lx-4*r & centers(:,2)>4*r & centers(:,2)<Ly-4*r;
Iv=std(local(valid & inner),0);
end

function phi = rasterPorosity(centers, radii, Lx, Ly, nx, ny)
x = ((0:nx-1)+0.5)*Lx/nx; y = ((0:ny-1)+0.5)*Ly/ny;
[X,Y] = meshgrid(x,y); solid = false(ny,nx);
for i=1:size(centers,1)
    solid = solid | hypot(X-centers(i,1),Y-centers(i,2)) <= radii(i);
end
phi = 1-mean(solid,'all');
end

function Kxx = upstreamKxx(centers, radii, Lx, Ly, nx, ny)
base=FoldedCartesianGrid(2,[0 Lx 0 Ly],[nx ny]);
coord=base.coordinates; level=-inf(size(coord,1),1);
for i=1:size(centers,1)
    level=max(level,radii(i)-hypot(coord(:,1)-centers(i,1),coord(:,2)-centers(i,2)));
end
try
    K=computePhysicalPermeabilityTensor(coord, base.triangles, level);
    if any(diag(K) <= 0) || norm(K-K.','fro') > 1e-6*max(1,norm(K,'fro'))
        Kxx = NaN;
    else
        Kxx=K(1,1)/9.869233e-12;
    end
catch
    Kxx = NaN;
end
end

function writeGeometryFigure(path, centers, radii, Lx, Ly, target, best)
f=figure('Visible','off','Color','w'); hold on; axis equal; xlim([0 Lx]); ylim([0 Ly]);
t=linspace(0,2*pi,80); for i=1:size(centers,1), plot(centers(i,1)+radii(i)*cos(t),centers(i,2)+radii(i)*sin(t),'k-'); end
xlabel('x (cm)'); ylabel('y (cm)'); title(sprintf('%s: phi %.3f (design %.3f), 2-D Kxx %.1f mD, Iv %.3f', ...
    target.sampleName,best.phi,target.targetPorosity,best.Kxx,best.Iv));
exportgraphics(f,path,'Resolution',200); close(f);
end

function radii = matchedRadii(centers, meanRadius, targetPhi, Lx, Ly, cv)
% Bounded log-normal radii; bisection fixes physical solid area, not K.
z = exp(cv*randn(size(centers,1),1) - 0.5*cv^2);
z = min(1.35, max(0.65, z));
lo = 0.5*meanRadius; hi = 1.5*meanRadius;
for iter=1:18
    scale=(lo+hi)/2; candidate=scale*z;
    phi=rasterPorosity(centers,candidate,Lx,Ly,360,240);
    if phi > targetPhi, lo=scale; else, hi=scale; end
end
radii=((lo+hi)/2)*z;
if any(centers(:,1) < radii | centers(:,1) > Lx-radii | ...
       centers(:,2) < radii | centers(:,2) > Ly-radii), radii=[]; return; end
dx=centers(:,1)-centers(:,1).'; dy=centers(:,2)-centers(:,2).';
dist=hypot(dx,dy); dist(1:size(dist,1)+1:end)=inf;
if any(dist < radii + radii.' + 1e-5, 'all'), radii=[]; end
end
