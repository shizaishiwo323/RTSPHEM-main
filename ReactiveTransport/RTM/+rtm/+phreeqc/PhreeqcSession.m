classdef PhreeqcSession < handle
    % PhreeqcSession - Persistent IPhreeqc session wrapper for RTM coupling.
    %
    % Keeps one PHREEQC engine alive, loads an exact database once, and runs
    % chemistry text through RunString so callers do not need per-step .phr
    % files. A mock engineFactory can be injected for unit tests.

    properties (Access = private)
        Engine = []
        EngineFactory = []
        OwnsEngine = true
        LoadedDatabasePath = ""
        Closed = false
    end

    properties (SetAccess = private)
        RuntimeConfig = struct()
        DatabasePath = ""
        DatabaseSha256 = ""
        LastLoadStatus = NaN
        LastRunStatus = NaN
        LastErrorString = ""
    end

    methods
        function obj = PhreeqcSession(runtimeConfig)
            if nargin < 1 || isempty(runtimeConfig)
                runtimeConfig = struct();
            end
            obj.RuntimeConfig = runtimeConfig;
            if isfield(runtimeConfig, 'engineFactory') && ~isempty(runtimeConfig.engineFactory)
                obj.EngineFactory = runtimeConfig.engineFactory;
            end
            if isfield(runtimeConfig, 'engine') && ~isempty(runtimeConfig.engine)
                obj.Engine = runtimeConfig.engine;
                obj.OwnsEngine = false;
            end
        end

        function loadDatabaseExact(obj, databasePath)
            obj.assertOpen();
            databasePath = char(databasePath);
            if exist(databasePath, 'file') ~= 2
                error('RTSPHEM:Phreeqc:MissingDatabase', ...
                    'PHREEQC exact database does not exist: %s', databasePath);
            end
            resolvedPath = string(char(java.io.File(databasePath).getCanonicalPath()));
            if obj.LoadedDatabasePath == resolvedPath
                return;
            end

            engine = obj.ensureEngine();
            status = engine.LoadDatabase(char(resolvedPath));
            if isempty(status)
                status = 0;
            end
            obj.LastLoadStatus = double(status);
            obj.LastErrorString = obj.readErrorString();
            if obj.LastLoadStatus ~= 0
                error('RTSPHEM:Phreeqc:LoadDatabaseFailed', ...
                    'PHREEQC LoadDatabase failed for %s: %s', ...
                    char(resolvedPath), char(obj.LastErrorString));
            end

            obj.LoadedDatabasePath = resolvedPath;
            obj.DatabasePath = resolvedPath;
            obj.DatabaseSha256 = obj.computeSha256(char(resolvedPath));
        end

        function rawOutput = runString(obj, inputText)
            obj.assertOpen();
            if strlength(obj.LoadedDatabasePath) == 0
                error('RTSPHEM:Phreeqc:DatabaseNotLoaded', ...
                    'Call loadDatabaseExact before runString.');
            end

            engine = obj.ensureEngine();
            status = engine.RunString(char(inputText));
            if isempty(status)
                status = 0;
            end
            obj.LastRunStatus = double(status);
            obj.LastErrorString = obj.readErrorString();
            if obj.LastRunStatus ~= 0
                error('RTSPHEM:Phreeqc:RunStringFailed', ...
                    'PHREEQC RunString failed: %s', char(obj.LastErrorString));
            end
            rawOutput = engine.GetSelectedOutputArray;
        end

        function status = getLastStatus(obj)
            status = struct( ...
                'loadStatus', obj.LastLoadStatus, ...
                'runStatus', obj.LastRunStatus, ...
                'errorString', obj.LastErrorString);
        end

        function errorString = getLastErrorString(obj)
            errorString = obj.LastErrorString;
        end

        function manifest = getDatabaseManifest(obj)
            manifest = struct( ...
                'databasePath', obj.DatabasePath, ...
                'databaseSha256', obj.DatabaseSha256, ...
                'engineType', string(obj.getRuntimeOption('engineType', 'iphreeqc_com')), ...
                'comProgId', string(obj.getRuntimeOption('comProgId', 'IPhreeqcCOM.Object')));
        end

        function close(obj)
            if obj.Closed
                return;
            end
            if ~isempty(obj.Engine) && obj.OwnsEngine
                try
                    delete(obj.Engine);
                catch
                end
            end
            obj.Closed = true;
        end

        function delete(obj)
            obj.close();
        end
    end

    methods (Access = private)
        function assertOpen(obj)
            if obj.Closed
                error('RTSPHEM:Phreeqc:SessionClosed', ...
                    'PHREEQC session is already closed.');
            end
        end

        function engine = ensureEngine(obj)
            if isempty(obj.Engine)
                if ~isempty(obj.EngineFactory)
                    obj.Engine = obj.EngineFactory();
                else
                    progId = char(obj.getRuntimeOption('comProgId', 'IPhreeqcCOM.Object'));
                    obj.Engine = actxserver(progId);
                end
            end
            engine = obj.Engine;
        end

        function value = getRuntimeOption(obj, fieldName, defaultValue)
            if isfield(obj.RuntimeConfig, fieldName) && ~isempty(obj.RuntimeConfig.(fieldName))
                value = obj.RuntimeConfig.(fieldName);
            else
                value = defaultValue;
            end
        end

        function errorString = readErrorString(obj)
            errorString = "";
            try
                if ~isempty(obj.Engine)
                    errorString = string(obj.Engine.GetErrorString());
                end
            catch
                errorString = "";
            end
        end
    end

    methods (Static, Access = private)
        function digest = computeSha256(pathValue)
            fid = fopen(pathValue, 'r');
            if fid == -1
                error('RTSPHEM:Phreeqc:HashOpenFailed', ...
                    'Cannot hash PHREEQC database: %s', pathValue);
            end
            cleanupObj = onCleanup(@() fclose(fid));
            bytes = fread(fid, Inf, '*uint8');
            clear cleanupObj;

            messageDigest = java.security.MessageDigest.getInstance('SHA-256');
            hashBytes = typecast(messageDigest.digest(bytes), 'uint8');
            digest = string(lower(sprintf('%02x', hashBytes)));
        end
    end
end
