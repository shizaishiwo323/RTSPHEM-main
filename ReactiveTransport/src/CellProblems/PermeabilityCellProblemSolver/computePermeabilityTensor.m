function [permeability] = computePermeabilityTensor(GridHyPHM, ...
    levelSet, varargin)

mainDirectionOnly = any(ismember(varargin, 'MainDirectionOnly'));
if mainDirectionOnly
    permeability = nan(2, 2);
else
    permeability = zeros(2, 2);
end
st = Stepper(0:1);

%noRefine does not apply local mesh refinement and folding
if ~any(ismember(varargin, 'noRefine'))
    [GridHyPHM, levelSet] = localMeshRefinement(GridHyPHM, ...
        levelSet);
else
    assert(isa(GridHyPHM, 'FoldedGrid'), 'Input must be Folded');
end

%Use P1-bubble/P1 discretization if specified; P2/P1 default
if any(ismember(varargin, 'Bubble'))
    Discretization = 'Bubble';
else
    Discretization = '';
end
problem = Perm(GridHyPHM, st, 'StokesPermeability', Discretization);
problem.solveSecondDirection = ~mainDirectionOnly;
problem.L.setdata(0, levelSet);

st.next; %bU contains integrals of basis functions
bU = problem.computeLevel('s');
velocities1 = problem.U.getdata(1);
% problem.P.visualize();
% problem.U.visualize();
permeability(1, 1) = bU' * velocities1(:, 1); %integration
permeability(1, 2) = bU' * velocities1(:, 2);
if ~mainDirectionOnly
    velocities2 = problem.U2.getdata(1);
    permeability(2, 1) = bU' * velocities2(:, 1);
    permeability(2, 2) = bU' * velocities2(:, 2);
end

if ~mainDirectionOnly && norm(permeability) < eps
    permeability = eye(2) * eps;
elseif mainDirectionOnly && permeability(1, 1) < eps
    permeability(1, 1) = eps;
end

end
