function [SOL] = solveSystemFE(grid, A, rhs, isDoF)
%Solve system of equations given by system matrix A and rhs for the
%unknowns given by isDoF

RHS = rhs(isDoF, :);
%     if( isa( grid, 'FoldedCartesianGrid' ) )
%         RHS( end + 1, : ) = 0;
%     end
if (isa(grid, 'FoldedCartesianGrid'))
    ANEW = A(1:end-1, :);
else ANEW = A;
end
%Remark: last line contained normalization condition --> solution now
%defined uniquely up to constant. PCG will chosen one that is normalized
%later

% sol = A\[RHS; 0 0];
r = symrcm(ANEW);
r = 1:size(ANEW, 1);
ANEW = ANEW(r, r);
RHS = RHS(r, :);
opts.type = 'ict';
opts.droptol = 0.0005;
temp(r) = 1:size(ANEW, 1);

try
    L = ichol(ANEW, opts);
    [sol2, flag2] = pcg(ANEW, RHS(:, 1), 10^(-10), 1000, L, L');
    [sol3, flag3] = pcg(ANEW, RHS(:, 2), 10^(-10), 1000, L, L');
    if flag2 + flag3 ~= 0
        error('solveSystemFE:PCGFailed', ...
            'Diffusion PCG solver failed with flags [%d, %d].', flag2, flag3);
    end
catch ME
    warning('solveSystemFE:DirectFallback', ...
        'Diffusion PCG/ichol solver failed (%s). Falling back to direct solve.', ME.message);
    solDirect = ANEW \ RHS;
    if any(~isfinite(solDirect(:)))
        regularization = max(1, norm(ANEW, 1)) * 1e-10;
        warning('solveSystemFE:RegularizedFallback', ...
            'Direct diffusion solve returned non-finite values. Retrying with diagonal regularization %.3e.', ...
            regularization);
        solDirect = (ANEW + regularization * speye(size(ANEW))) \ RHS;
    end
    sol2 = solDirect(:, 1);
    sol3 = solDirect(:, 2);
end

if any(~isfinite(sol2)) || any(~isfinite(sol3))
    warning('solveSystemFE:NonFiniteSolutionClamped', ...
        'Diffusion solve still returned non-finite values; replacing them with zero.');
    sol2(~isfinite(sol2)) = 0;
    sol3(~isfinite(sol3)) = 0;
end

sol2 = sol2 - sum(sol2) / size(sol2, 1);
sol2 = sol2(temp);
sol3 = sol3 - sum(sol3) / size(sol3, 1);
sol3 = sol3(temp);

%    norm(sol2-sol(:,1))/ norm(sol(:,1)) + norm(sol3-sol(:,2))/ norm(sol(:,2))
sol = [sol2, sol3];
SOL = zeros(grid.nodes, 2);
SOL(isDoF, :) = sol;

if (isa(grid, 'FoldedCartesianGrid'))
    SOL(:, 1) = grid.synchronizeValues(SOL(:, 1));
    SOL(:, 2) = grid.synchronizeValues(SOL(:, 2));
end

end
