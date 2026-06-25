classdef MockIPhreeqcEngine < handle
    properties
        RawOutput
        LoadDatabaseCount = 0
        RunStringCount = 0
        RunFileCount = 0
        LoadedDatabasePaths = strings(0, 1)
        LastInputText = ""
        LastInputPath = ""
        ErrorString = ""
        LoadStatus = 0
        RunStatus = 0
        Deleted = false
    end

    methods
        function obj = MockIPhreeqcEngine(rawOutput)
            if nargin < 1
                rawOutput = {};
            end
            obj.RawOutput = rawOutput;
        end

        function status = LoadDatabase(obj, databasePath)
            obj.LoadDatabaseCount = obj.LoadDatabaseCount + 1;
            obj.LoadedDatabasePaths(end + 1, 1) = string(databasePath);
            status = obj.LoadStatus;
        end

        function status = RunString(obj, inputText)
            obj.RunStringCount = obj.RunStringCount + 1;
            obj.LastInputText = string(inputText);
            status = obj.RunStatus;
        end

        function status = RunFile(obj, inputPath)
            obj.RunFileCount = obj.RunFileCount + 1;
            obj.LastInputPath = string(inputPath);
            status = obj.RunStatus;
        end

        function rawOutput = GetSelectedOutputArray(obj)
            if isa(obj.RawOutput, 'function_handle')
                rawOutput = obj.RawOutput(obj);
            else
                rawOutput = obj.RawOutput;
            end
        end

        function errorString = GetErrorString(obj)
            errorString = char(obj.ErrorString);
        end

        function delete(obj)
            obj.Deleted = true;
        end
    end
end
