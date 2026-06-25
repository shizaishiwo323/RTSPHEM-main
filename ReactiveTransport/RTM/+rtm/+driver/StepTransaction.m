classdef StepTransaction < handle
    %STEPTRANSACTION Snapshot/commit/rollback helper for RTM driver steps.

    properties (Access = private)
        OriginalState
        OriginalGeometry
        OriginalSolverState
        CurrentStatus = 'open'
    end

    methods
        function obj = StepTransaction(state, geometry, solverState)
            if nargin < 3
                solverState = struct();
            end
            obj.OriginalState = state;
            obj.OriginalGeometry = geometry;
            obj.OriginalSolverState = solverState;
        end

        function [state, geometry, solverState] = commit(obj, state, geometry, solverState)
            obj.assertOpen();
            if nargin < 4
                solverState = struct();
            end
            obj.CurrentStatus = 'committed';
        end

        function [state, geometry, solverState] = rollback(obj)
            obj.assertOpen();
            state = obj.OriginalState;
            geometry = obj.OriginalGeometry;
            solverState = obj.OriginalSolverState;
            obj.CurrentStatus = 'rolled_back';
        end

        function value = status(obj)
            value = obj.CurrentStatus;
        end
    end

    methods (Access = private)
        function assertOpen(obj)
            if ~strcmp(obj.CurrentStatus, 'open')
                error('RTSPHEM:Driver:TransactionClosed', ...
                    'StepTransaction is already %s.', obj.CurrentStatus);
            end
        end
    end
end
