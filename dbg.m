% dbg.m
% Light-weight conditional logger.
% Usage: dbg(params, 'format %d', val);
function dbg(params, varargin)
    if isfield(params,'verbose') && params.verbose
        % ts = datestr(now,'HH:MM:SS.FFF');
        ts = char(datetime('now','Format','HH:mm:ss.SSS'));
        try
            fprintf('[%s] ', ts);
            fprintf(varargin{:});
            if isempty(regexp(varargin{1}, '\n$','once'))
                fprintf('\n');
            end
        catch
            % Fail-silent if bad format
            fprintf('[%s] <dbg format error>\n', ts);
        end
    end
end