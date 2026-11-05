function solution = SolveStokesLevelWorker(meshPath, levels, flowDirection, ...
        lengthXAxis, lengthYAxis, nxParts, nyParts, inletVelocity, mu)
% Rebuild and solve Stokes exactly from compact numeric inputs.

mesh = load(meshPath, 'workerCoord', 'workerTriangles');
grid = Grid(mesh.workerCoord, mesh.workerTriangles);
edgeX = reshape(grid.coordV(grid.V0E(:), 1), size(grid.V0E));
edgeY = reshape(grid.coordV(grid.V0E(:), 2), size(grid.V0E));
grid.idE(all(edgeY < eps, 2)) = 1;
grid.idE(all(edgeX > lengthXAxis - eps, 2)) = 2;
grid.idE(all(edgeY > lengthYAxis - eps, 2)) = 3;
grid.idE(all(edgeX < eps, 2)) = 4;
stepper = Stepper(0:1);
pressure = Variable(grid, stepper, 'pressure head', 'P1');
velocity = Variable(grid, stepper, 'Stokes Velocity', 'P2P2');
pressure.setdata(0, @(t, x) 0);
velocity.setdata(0, @(t, x) 0);
flowConfig = GetFlowBoundaryConfig( ...
    flowDirection, lengthXAxis, lengthYAxis, nxParts, nyParts);
problem = StokesLEVEL(grid, stepper, 'Stokes problem');
problem.L.setdata(levels);
problem.id2D = num2cell(flowConfig.stokesDirichletIds);
problem.uD.setdata(@(t, x) double(flowConfig.inletPredicate(x)) .* ...
    flowConfig.velocityVector(inletVelocity));
problem.F.setdata(@(t, x) 0);
problem.U = velocity;
problem.P = pressure;
problem.N.setdata(mu);
stepper.next;
problem.computeLevel('s');
solution = struct( ...
    'velocity', problem.U.getdata(1), ...
    'pressure', problem.P.getdata(1));
end
