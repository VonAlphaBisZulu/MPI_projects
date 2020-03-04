classdef Progressbar
    % generates a progress bar window that can be updated
    % 
    % To create the window
    % pb = Progressbar(init val , Text , max val) 
    %
    % To update
    % pb.update(val)
    %
    % To close
    % pb.delete
    
    properties
        h
        val
        total
        starttime
        timetable
    end
    
    methods
        function obj = Progressbar(val,text,total)
            obj.val = val;
            switch nargin
                case 1
                    obj.h = waitbar(val);
                    obj.total = 1;
                case 2
                    obj.h = waitbar(val,text);
                    obj.total = 1;
                    obj.h.Name = text;
                case 3
                    obj.h = waitbar(val,text);
                    obj.total = total;
                    obj.h.Name = text;
                otherwise
                    obj.h = waitbar(0);
                    obj.total = total;
                    obj.h.Name = text;
            end
            obj.timetable = repmat({0 0},0,1);
        end
        function obj = update(obj,x)
            obj.timetable(end+1,:) = [{clock} {x/obj.total}];
            % calculate remaining time
            text = [num2str(x) ' / ' num2str(obj.total)];
            if size(obj.timetable,1)>20
                timePassed = cellfun(@etime,obj.timetable((end-19):end,1),repmat({clock},20,1));
                obj.timetable((end-19):end,:);
                P = polyfit(timePassed,cell2mat(obj.timetable((end-19):end,2)),1);
                eta = (1-P(2))/P(1);
                text = [text '  ' obj.secToStrTime(eta) ' remaining'];
            end
            % update waitbar
            waitbar(obj.timetable{end,2},obj.h,text);
        end
        function time = secToStrTime(obj,t)
            t = round(t);
            dys    = floor(t/24/60/60);
            t = t - dys*24*60*60;
            hrs   = floor(t/60/60);
            t = t - hrs*60*60;
            min = floor(t/60);
            sec = t - min*60;
            time = '';
            if dys
                time = [time num2str(dys) 'd '];
            end
            if hrs
                time = [time num2str(hrs) 'h '];
            end
            if min
                time = [time num2str(min) 'm '];
            end
            if sec
                time = [time num2str(sec) 's '];
            end
        end
        function delete(obj)
            obj.h.delete();
        end
    end
    
end

