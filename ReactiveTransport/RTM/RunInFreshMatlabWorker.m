function result = RunInFreshMatlabWorker(functionHandle, varargin)
% Run one memory-heavy exact solve in a disposable MATLAB process.

inputPath = [tempname, '_rtsphem_input.mat'];
outputPath = [tempname, '_rtsphem_output.mat'];
cleanup = onCleanup(@() cleanupFiles(inputPath, outputPath));
payload = struct('functionName', func2str(functionHandle), ...
    'arguments', {varargin});
save(inputPath, 'payload', '-v7.3');

rtmDir = fileparts(mfilename('fullpath'));
batchCode = sprintf("addpath('%s'); RunMatlabWorkerFile('%s','%s')", ...
    escapeMatlabString(rtmDir), escapeMatlabString(inputPath), ...
    escapeMatlabString(outputPath));
maxAttempts = 3;
for attempt = 1:maxAttempts
    [status, commandOutput] = system(sprintf('matlab -batch "%s"', batchCode));
    if status == 0 && exist(outputPath, 'file') == 2
        break;
    end
    retryableMemoryFailure = contains(commandOutput, ...
        'chol(A(ordering, ordering))') || contains(commandOutput, ...
        'A = sparse(Arow, Acol, Aval') || contains(commandOutput, ...
        'Out of memory', 'IgnoreCase', true) || contains(commandOutput, '内存');
    if ~retryableMemoryFailure || attempt == maxAttempts
        error('RTSPHEM:WorkerIsolation:ProcessFailed', ...
            'Isolated MATLAB solve failed after %d attempt(s) (status %d):\n%s', ...
            attempt, status, commandOutput);
    end
    fprintf(['Isolated solve hit a transient sparse-memory failure; ', ...
        'retrying (%d/%d).\n'], attempt + 1, maxAttempts);
    pause(5);
end
loaded = load(outputPath, 'result');
result = loaded.result;
clear cleanup;
cleanupFiles(inputPath, outputPath);
end

function value = escapeMatlabString(value)
value = strrep(char(value), '''', '''''');
end

function cleanupFiles(varargin)
for i = 1:nargin
    if exist(varargin{i}, 'file') == 2
        delete(varargin{i});
    end
end
end
