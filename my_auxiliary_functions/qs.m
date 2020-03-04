classdef qs < handle
    % creates a script in the tmp directory that can be called also by the qs() call.
    % qs(1) creates the script in pwd instead of tmp. qs(0) or the deletion of the 
    % 'quickscript' handle deletes the script file.
    properties  (SetAccess = 'public')
        script_exists
        script_path
    end
    properties (Access=private)
        in_var
    end
    methods (Access = 'public')
        function obj = qs(varargin)
            obj.script_path = [tempdir 'qs.m'];
            obj.script_exists = exist(obj.script_path,'file');
            if size(varargin)>0
                if varargin{:} == 1 && obj.script_exists
                    copyfile(obj.script_path,pwd);
                    delete(obj.script_path);
                    disp('script moved to pwd');
                    edit([pwd '/qs.m']);
                    return
                elseif varargin{:} == 1 && ~obj.script_exists
                    obj.script_path = [pwd '/qs.m'];
                elseif varargin{:} ~= 1
                    evalin('caller', 'clear(''quickscript'')');
                    return
                end
            end
            obj.in_var = 0;
            if ~evalin( 'base', 'exist(''quickscript'',''var'') == 1' )
                f = fopen(obj.script_path, 'w');
                fclose(f);
                k='script created';
                disp(k);
                edit(obj.script_path);
                obj1 = obj;
                obj1.in_var = 1;
                assignin('base','quickscript',obj1);
            else
                if obj.script_exists
                    text = fileread(obj.script_path);
                    evalin('caller', text);
                else
                    disp('script file doesn''t exist.');
                end
            end
        end
        function delete(obj)
            if obj.in_var
                delete(obj.script_path);
                disp('script deleted');
            end
        end
    end
end