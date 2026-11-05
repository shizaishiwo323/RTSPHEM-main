function RunMatlabWorkerFile(inputPath, outputPath)
% File-based entry point for one isolated exact solve.

rtmDir = fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(rtmDir)));
loaded = load(inputPath, 'payload');
functionHandle = str2func(loaded.payload.functionName);
result = functionHandle(loaded.payload.arguments{:});
save(outputPath, 'result', '-v7.3');
end
