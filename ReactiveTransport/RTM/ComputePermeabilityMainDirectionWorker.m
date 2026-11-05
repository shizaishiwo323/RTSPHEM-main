function result = ComputePermeabilityMainDirectionWorker(meshPath, levelSet)
% Evaluate the upstream RTSPHEM cell problem in the main flow direction.

mesh = load(meshPath, 'workerCoord', 'workerTriangles');
K_cm2 = computePhysicalPermeabilityTensor( ...
    mesh.workerCoord, mesh.workerTriangles, levelSet, 'MainDirectionOnly');
result = struct('K_cm2', K_cm2, 'K_mD', K_cm2 / 9.869233e-12);
end
